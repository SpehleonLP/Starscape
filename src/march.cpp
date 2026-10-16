#include "march.h"
#include <algorithm>
#include <cmath>
#include <span>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

// The calibration and the star re-march need a step's extent, its transmittance and its stars; the light of its
// gas is the GPU march's alone.
struct Step
{
	double                        a;
	double                        b;
	glm::dvec3                    transmittance;
	std::span<Contribution const> parts;
};

// Steps every interval at its midpoints, handing each step's contributions and the
// transmittance up to it to `visit`. Returns the total visual extinction in magnitudes.
template<class Visit>
double Walk(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dir, double maxDistance, Visit&& visit)
{
	thread_local std::vector<Interval> intervals;
	thread_local std::vector<Contribution> parts;
	RayIntervals(scene, dir, maxDistance, intervals);
	double depth = 0.0;
	for (Interval const& interval : intervals)
	{
		double s = interval.begin;
		while (s < interval.end)
		{
			double const remaining = interval.end - s;
			double const ds = std::min(std::clamp(params.alpha * s, params.dsMin, params.dsMax), remaining);
			double const mid = s + 0.5 * ds;
			glm::dvec3 const position = scene.origin + mid * dir;
			parts.clear();
			double extinction = 0.0;
			for (PreparedGalaxy const& galaxy : scene.galaxies)
			{
				glm::dvec3 const offset = position - galaxy.position;
				if (glm::dot(offset, offset) > galaxy.boundRadius * galaxy.boundRadius)
				{
					continue;
				}
				GalaxySample sample;
				Evaluate(galaxy, position, sample);
				parts.insert(parts.end(), sample.parts.begin(), sample.parts.begin() + sample.count);
				extinction += sample.extinction;
			}
			visit(Step{s, s + ds, TransmittanceForDepth(depth + 0.5 * extinction * ds), parts});
			depth += extinction * ds;
			s = ds == remaining ? interval.end : s + ds;
		}
	}
	return depth;
}

}

MarchParams ParamsFor(Quality quality)
{
	if (quality == Quality::Quick)
	{
		return {.samplesPerAxis = 1, .alpha = 0.1, .dsMin = 0.01, .dsMax = 50.0};
	}
	return {.samplesPerAxis = 4, .alpha = 0.02, .dsMin = 0.01, .dsMax = 25.0};
}

PreparedScene PrepareScene(Scene const& scene)
{
	PreparedScene prepared;
	for (Population const& population : scene.populations)
	{
		prepared.populations.push_back(Tabulate(population, scene.style.starSaturation));
	}
	for (Galaxy const& galaxy : scene.galaxies)
	{
		PreparedGalaxy ready = PrepareGalaxy(galaxy);
		Disk& disk = ready.galaxy.disk;
		if (disk.armPopulation >= 0)
		{
			// Gas glows at a multiple of its young stars' light, and only their table knows how bright they are.
			disk.gas.youngLight = Luminance(prepared.populations[size_t(disk.armPopulation)].lightBelow[kBins]);
		}
		prepared.galaxies.push_back(std::move(ready));
	}
	prepared.origin = scene.observer.position;
	prepared.skyToScene = SkyToScene(scene.observer);
	prepared.sMin = scene.observer.sMin;
	return prepared;
}

void RayIntervals(PreparedScene const& scene, glm::dvec3 dir, double maxDistance, std::vector<Interval>& out)
{
	out.clear();
	for (PreparedGalaxy const& galaxy : scene.galaxies)
	{
		glm::dvec3 const offset = scene.origin - galaxy.position;
		double const b = glm::dot(offset, dir);
		double const c = glm::dot(offset, offset) - galaxy.boundRadius * galaxy.boundRadius;
		double const discriminant = b * b - c;
		if (discriminant <= 0.0)
		{
			continue;
		}
		double const root = std::sqrt(discriminant);
		double const begin = std::max(-b - root, scene.sMin);
		double const end = std::min(-b + root, maxDistance);
		if (end > begin)
		{
			out.push_back({begin, end});
		}
	}
	std::sort(out.begin(), out.end(), [](Interval const& a, Interval const& b) { return a.begin < b.begin; });
	// Overlapping spheres share their stretch of ray; merging marches each point once.
	size_t merged = 0;
	for (size_t i = 0; i < out.size(); ++i)
	{
		if (merged > 0 && out[i].begin <= out[merged - 1].end)
		{
			out[merged - 1].end = std::max(out[merged - 1].end, out[i].end);
		}
		else
		{
			out[merged++] = out[i];
		}
	}
	out.resize(merged);
}

glm::dvec3 TransmittanceTo(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double distance)
{
	return TransmittanceForDepth(Walk(scene, params, scene.skyToScene * dirSky, distance, [](Step const&) {}));
}

void FluxHistogram::Add(double flux, double stars)
{
	if (!(flux > 0.0) || !(stars > 0.0))
	{
		return;
	}
	int const bucket = std::clamp(int(std::floor((std::log10(flux) - kLogMin) / kBucketWidth)), 0, kBuckets - 1);
	count[size_t(bucket)] += stars;
}

void FluxHistogram::Merge(FluxHistogram const& other)
{
	for (size_t i = 0; i < count.size(); ++i)
	{
		count[i] += other.count[i];
	}
}

// Within a bucket, counts are taken as uniform in log flux.
double FluxHistogram::CountAbove(double flux) const
{
	double const position = (std::log10(flux) - kLogMin) / kBucketWidth;
	if (position >= kBuckets)
	{
		return 0.0;
	}
	int const bucket = std::max(0, int(std::floor(position)));
	double sum = count[size_t(bucket)] * std::clamp(1.0 - (position - bucket), 0.0, 1.0);
	for (int b = bucket + 1; b < kBuckets; ++b)
	{
		sum += count[size_t(b)];
	}
	return sum;
}

std::optional<double> FluxHistogram::FluxForCount(double target) const
{
	double cumulative = 0.0;
	for (int b = kBuckets - 1; b >= 0; --b)
	{
		double const here = count[size_t(b)];
		if (here > 0.0 && cumulative + here >= target)
		{
			double const fractionAbove = (target - cumulative) / here;
			return std::pow(10.0, kLogMin + (b + 1 - fractionAbove) * kBucketWidth);
		}
		cumulative += here;
	}
	return std::nullopt;
}

void HistogramRay(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double sampleSolidAngle,
	FluxHistogram& out)
{
	Walk(scene, params, scene.skyToScene * dirSky, HUGE_VAL, [&](Step const& step)
	{
		// Flux falls as s^-2, so a step spans 2 log10(b/a) decades. CountAbove interpolates linearly inside a bucket, so a
		// piece's mass must be a small fraction of one: quarter-bucket pieces keep that error under the 3% count tolerance
		// (half-bucket pieces still rippled +-4.6%). Each piece is deposited at its own midpoint.
		int const pieces = std::max(1, int(std::ceil(8.0 * std::log10(step.b / step.a) / FluxHistogram::kBucketWidth)));
		for (int piece = 0; piece < pieces; ++piece)
		{
			double const x = step.a + (step.b - step.a) * piece / pieces;
			double const y = step.a + (step.b - step.a) * (piece + 1) / pieces;
			double const middle = 0.5 * (x + y);
			double const perLuminosity = step.transmittance.y * 100.0 * kFluxAt10pc / (middle * middle);
			double const volume = sampleSolidAngle * CubeSpan(x, y);
			for (Contribution const& part : step.parts)
			{
				PopulationTable const& population = scene.populations[size_t(part.population)];
				for (int bin = 0; bin < kBins; ++bin)
				{
					out.Add(population.luminosity[size_t(bin)] * perLuminosity, part.density * population.fraction[size_t(bin)] * volume);
				}
			}
		}
	});
}

std::optional<StarSample> SampleStar(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, int bin,
	double fluxThreshold, double u)
{
	struct Entry
	{
		double     cumulative;
		double     count;
		double     x;
		double     y;
		glm::dvec3 transmittance;
		int        population;
	};
	thread_local std::vector<Entry> entries;
	entries.clear();
	double total = 0.0;
	Walk(scene, params, scene.skyToScene * dirSky, HUGE_VAL, [&](Step const& step)
	{
		for (Contribution const& part : step.parts)
		{
			PopulationTable const& population = scene.populations[size_t(part.population)];
			double const x = step.a;
			double const y = std::min(step.b, ResolvedDistance(fluxThreshold, step.transmittance.y, population.luminosity[size_t(bin)]));
			if (!(y > x))
			{
				continue;
			}
			double const count = part.density * population.fraction[size_t(bin)] * CubeSpan(x, y);
			if (!(count > 0.0))
			{
				continue;
			}
			total += count;
			entries.push_back({total, count, x, y, step.transmittance, part.population});
		}
	});
	if (!(total > 0.0))
	{
		return std::nullopt;
	}
	double const target = u * total;
	auto it = std::partition_point(entries.begin(), entries.end(), [target](Entry const& e) { return e.cumulative < target; });
	if (it == entries.end())
	{
		it = entries.end() - 1;
	}
	double const fraction = std::clamp((target - (it->cumulative - it->count)) / it->count, 0.0, 1.0);
	double const x3 = it->x * it->x * it->x;
	double const y3 = it->y * it->y * it->y;
	return StarSample{std::cbrt(x3 + fraction * (y3 - x3)), it->transmittance, it->population};
}

}
