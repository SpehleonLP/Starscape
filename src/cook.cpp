#include "cook.h"
#include "parallel.h"
#include <kre_starscape.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <format>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

constexpr double kTwoPi = 6.28318530717958647692;
constexpr uint32_t kCalibrationFaceSize = 64;

struct TexelCoord
{
	uint32_t face;
	uint32_t x;
	uint32_t y;
};

TexelCoord Decompose(uint32_t faceSize, uint32_t texel)
{
	return {texel / (faceSize * faceSize), texel % faceSize, (texel / faceSize) % faceSize};
}

glm::dvec3 FaceDirection(uint32_t face, double s, double t)
{
	std::array<double, 3> const d = kre_starscape::FaceToDir(face, s, t);
	return glm::normalize(glm::dvec3(d[0], d[1], d[2]));
}

}

glm::dvec3 TexelCentre(uint32_t faceSize, uint32_t texel)
{
	TexelCoord const c = Decompose(faceSize, texel);
	return FaceDirection(c.face, (c.x + 0.5) / faceSize, (c.y + 0.5) / faceSize);
}

uint32_t SamplePoisson(double lambda, RandomStream& rng)
{
	if (!(lambda > 0.0))
	{
		return 0;
	}
	if (lambda < 30.0)
	{
		// Knuth: multiply uniforms until the product drops below e^-λ.
		double const limit = std::exp(-lambda);
		uint32_t k = 0;
		double product = rng.NextDouble();
		while (product > limit)
		{
			++k;
			product *= rng.NextDouble();
		}
		return k;
	}
	// For large λ a normal approximation is within a fraction of a percent, and only the
	// densest texels near a core ever get here.
	double const u1 = std::max(rng.NextDouble(), 1e-300);
	double const u2 = rng.NextDouble();
	double const z = std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
	// Converting a double outside uint32_t's range is undefined, and an infinite λ would make the sum NaN, so
	// λ is capped (below the type's maximum, leaving room for the normal's spread) before the cast.
	double const capped = std::min(lambda, 4.2e9);
	return uint32_t(std::clamp(std::floor(capped + std::sqrt(capped) * z + 0.5), 0.0, 4294967295.0));
}

glm::dvec3 SampleTexelDirection(uint32_t faceSize, uint32_t texel, RandomStream& rng)
{
	TexelCoord const c = Decompose(faceSize, texel);
	// Solid angle per unit face area is (1 + u² + v²)^-3/2; rejection against its maximum over
	// the texel makes directions uniform in solid angle rather than in face area.
	auto jacobian = [](double u, double v) { return std::pow(1.0 + u * u + v * v, -1.5); };
	double const step = 2.0 / faceSize;
	double const u0 = -1.0 + c.x * step;
	double const v0 = -1.0 + c.y * step;
	double const peak = std::max({jacobian(u0, v0), jacobian(u0 + step, v0), jacobian(u0, v0 + step), jacobian(u0 + step, v0 + step)});
	for (;;)
	{
		double const s = (c.x + rng.NextDouble()) / faceSize;
		double const t = (c.y + rng.NextDouble()) / faceSize;
		if (rng.NextDouble() * peak <= jacobian(2.0 * s - 1.0, 2.0 * t - 1.0))
		{
			return FaceDirection(c.face, s, t);
		}
	}
}

FluxHistogram CalibrationHistogram(PreparedScene const& scene, MarchParams const& params, uint32_t faceSize, unsigned threads)
{
	// One partial histogram per texel row, merged in row order, keeps the floating-point sum
	// independent of the schedule; a per-texel histogram would cost 512 doubles per texel.
	uint32_t const rows = 6 * faceSize;
	std::vector<FluxHistogram> partial(rows);
	ParallelFor(rows, threads, [&](uint32_t row)
	{
		for (uint32_t x = 0; x < faceSize; ++x)
		{
			uint32_t const texel = row * faceSize + x;
			TexelCoord const c = Decompose(faceSize, texel);
			HistogramRay(scene, params, TexelCentre(faceSize, texel), kre_starscape::TexelSolidAngle(faceSize, c.x, c.y), partial[row]);
		}
	}, {});
	FluxHistogram total;
	for (FluxHistogram const& row : partial)
	{
		total.Merge(row);
	}
	return total;
}

Expected<Cooked> Cook(Scene const& scene, CookOptions const& options)
{
	Quality const quality = options.quality.value_or(scene.output.quality);
	MarchParams const params = ParamsFor(quality);
	PreparedScene const prepared = PrepareScene(scene);
	uint32_t const size = scene.output.faceSize;
	uint32_t const texels = 6 * size * size;
	uint64_t const seed = scene.output.seed;
	auto report = [&](std::string_view phase, uint32_t done, uint32_t total)
	{
		if (options.progress)
		{
			options.progress(phase, double(done) / double(total));
		}
	};

	Cooked cooked;
	cooked.faceSize = size;

	// 1. The flux threshold: given, or the one that leaves `target_stars` resolved.
	if (scene.output.limitingMagnitude)
	{
		cooked.fluxThreshold = std::pow(10.0, -0.4 * *scene.output.limitingMagnitude);
	}
	else
	{
		MarchParams calibration = params;
		calibration.samplesPerAxis = 1;
		report("calibrate", 0, 1);
		FluxHistogram const histogram = CalibrationHistogram(prepared, calibration, kCalibrationFaceSize, options.threads);
		std::optional<double> const threshold = histogram.FluxForCount(double(scene.output.targetStars.value_or(0)));
		cooked.stats.budgetUnreachable = !threshold.has_value();
		cooked.fluxThreshold = threshold.value_or(std::pow(10.0, FluxHistogram::kLogMin));
		report("calibrate", 1, 1);
	}

	// 2. Diffuse radiance and the per-bin point budget, one texel per slot.
	cooked.diffuse.assign(texels, glm::dvec3(0.0));
	std::vector<std::array<float, kBins>> budget(texels);
	int const perAxis = params.samplesPerAxis;
	double const weight = 1.0 / double(perAxis * perAxis);
	ParallelFor(texels, options.threads, [&](uint32_t texel)
	{
		TexelCoord const c = Decompose(size, texel);
		double const sampleSolidAngle = kre_starscape::TexelSolidAngle(size, c.x, c.y) * weight;
		RandomStream rng(seed, Purpose::Jitter, texel, 0, 0);
		RayResult result;
		for (int j = 0; j < perAxis; ++j)
		{
			for (int i = 0; i < perAxis; ++i)
			{
				double const s = (c.x + (i + rng.NextDouble()) / perAxis) / size;
				double const t = (c.y + (j + rng.NextDouble()) / perAxis) / size;
				MarchRay(prepared, params, FaceDirection(c.face, s, t), cooked.fluxThreshold, sampleSolidAngle, weight, result);
			}
		}
		cooked.diffuse[texel] = result.radiance;
		for (int bin = 0; bin < kBins; ++bin)
		{
			budget[texel][size_t(bin)] = float(result.budget[size_t(bin)]);
		}
	}, [&](uint32_t done) { report("march", done, texels); });

	// The budget is known before any star is drawn, so an absurd threshold fails here rather than hours later.
	double expectedStars = 0.0;
	double largestBin = 0.0;
	for (std::array<float, kBins> const& texelBudget : budget)
	{
		for (float const count : texelBudget)
		{
			expectedStars += double(count);
			largestBin = std::max(largestBin, double(count));
		}
	}
	double const magnitude = -2.5 * std::log10(cooked.fluxThreshold);
	if (!(expectedStars <= kMaxStars))
	{
		return InputError(std::format("the resolution threshold (magnitude {:.2f}) resolves about {:.3g} stars, over the {:.0f} limit; lower limiting_magnitude or target_stars",
			magnitude, expectedStars, kMaxStars));
	}
	if (!(largestBin <= kMaxStarsPerTexelBin))
	{
		return InputError(std::format("the resolution threshold (magnitude {:.2f}) puts about {:.3g} stars of one luminosity bin in one texel, over the {:.0f} per-texel limit; lower limiting_magnitude or target_stars, or raise face_size",
			magnitude, largestBin, kMaxStarsPerTexelBin));
	}

	// 3. Stars: Poisson counts from the budget, each placed by re-marching its own direction.
	std::vector<std::vector<CookedStar>> perTexel(texels);
	std::atomic<uint64_t> dropped{0};
	ParallelFor(texels, options.threads, [&](uint32_t texel)
	{
		for (int bin = 0; bin < kBins; ++bin)
		{
			RandomStream poisson(seed, Purpose::Poisson, texel, uint32_t(bin), 0);
			uint32_t const count = SamplePoisson(budget[texel][size_t(bin)], poisson);
			for (uint32_t index = 0; index < count; ++index)
			{
				RandomStream rng(seed, Purpose::Star, texel, uint32_t(bin), index);
				glm::dvec3 const dir = SampleTexelDirection(size, texel, rng);
				std::optional<StarSample> const sample = SampleStar(prepared, params, dir, bin, cooked.fluxThreshold, rng.NextDouble());
				if (!sample)
				{
					++dropped;
					continue;
				}
				PopulationTable const& population = prepared.populations[size_t(sample->population)];
				double const atObserver = population.luminosity[size_t(bin)] * 100.0 * kFluxAt10pc / (sample->distance * sample->distance);
				glm::dvec3 flux = atObserver * sample->transmittance * population.chroma[size_t(bin)];
				// Interpolating the distance inside a march step can leave a star a hair under the
				// threshold it was classified by; it was counted as a point, so it stays one, at the threshold.
				double const visual = atObserver * sample->transmittance.y;
				if (visual > 0.0 && visual < cooked.fluxThreshold)
				{
					flux *= cooked.fluxThreshold / visual;
				}
				perTexel[texel].push_back({dir, flux});
			}
		}
	}, [&](uint32_t done) { report("stars", done, texels); });

	for (std::vector<CookedStar>& stars : perTexel)
	{
		cooked.stars.insert(cooked.stars.end(), stars.begin(), stars.end());
	}
	std::stable_sort(cooked.stars.begin(), cooked.stars.end(), [](CookedStar const& a, CookedStar const& b)
	{
		return Luminance(a.flux) > Luminance(b.flux);
	});

	for (uint32_t texel = 0; texel < texels; ++texel)
	{
		TexelCoord const c = Decompose(size, texel);
		cooked.stats.diffuseFlux += Luminance(cooked.diffuse[texel]) * kre_starscape::TexelSolidAngle(size, c.x, c.y);
	}
	for (CookedStar const& star : cooked.stars)
	{
		cooked.stats.pointFlux += Luminance(star.flux);
	}
	cooked.stats.starCount = cooked.stars.size();
	cooked.stats.droppedStars = dropped.load();
	return cooked;
}

}
