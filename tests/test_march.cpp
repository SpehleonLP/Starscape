#include "march.h"
#include "test_scenes.h"
#include <doctest.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <glm/geometric.hpp>

using namespace starscape;
using namespace starscape::test;

TEST_CASE("Ray intervals clip, merge and skip")
{
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 12);
	scene.galaxies.push_back(UniformBall(1.0, 10.0, population, {100.0, 0.0, 0.0}));
	scene.galaxies.push_back(UniformBall(1.0, 10.0, population, {115.0, 0.0, 0.0}));   // overlaps the first
	scene.galaxies.push_back(UniformBall(1.0, 5.0, population, {0.0, 0.0, 0.0}));      // contains the observer
	PreparedScene const prepared = PrepareScene(scene);
	std::vector<Interval> intervals;

	RayIntervals(prepared, {1.0, 0.0, 0.0}, HUGE_VAL, intervals);
	REQUIRE(intervals.size() == 2);
	CHECK(intervals[0].begin == doctest::Approx(scene.observer.sMin));
	CHECK(intervals[0].end == doctest::Approx(5.0));
	CHECK(intervals[1].begin == doctest::Approx(90.0));
	CHECK(intervals[1].end == doctest::Approx(125.0));

	RayIntervals(prepared, {1.0, 0.0, 0.0}, 100.0, intervals);
	REQUIRE(intervals.size() == 2);
	CHECK(intervals[1].end == doctest::Approx(100.0));

	RayIntervals(prepared, {-1.0, 0.0, 0.0}, HUGE_VAL, intervals);
	REQUIRE(intervals.size() == 1);   // only the sphere around the observer

	// Passes 9.95 pc from the first sphere's centre (a short chord) and 11.4 pc from the second's (a miss).
	RayIntervals(prepared, glm::normalize(glm::dvec3(100.0, 10.0, 0.0)), HUGE_VAL, intervals);
	REQUIRE(intervals.size() == 2);
	CHECK(intervals[1].end - intervals[1].begin < 3.0);
}

TEST_CASE("Counts in a uniform medium follow the Euclidean F^-3/2 law")
{
	// Mutation: dropping the s² in the point budget, or depositing once per step in the histogram reddens this test.
	Scene scene = BaseScene();
	int const bin = 12;
	int const population = AddSingleBinPopulation(scene, bin);
	double const density = 0.01;
	scene.galaxies.push_back(UniformBall(density, 2000.0, population, {0.0, 0.0, 0.0}));
	PreparedScene const prepared = PrepareScene(scene);
	FluxHistogram histogram;
	// The medium is isotropic, so one ray standing for the whole sphere gives every count.
	HistogramRay(prepared, ParamsFor(Quality::Quick), {0.0, 0.0, 1.0}, 4.0 * kPi, histogram);

	double const luminosity = BinLuminosity(bin);
	auto fluxAt = [&](double s) { return luminosity * 100.0 * kFluxAt10pc / (s * s); };
	std::vector<double> logF;
	std::vector<double> logN;
	for (double s = 5.0; s <= 500.0; s *= 1.5)
	{
		double const expected = density * 4.0 / 3.0 * kPi * (s * s * s - std::pow(scene.observer.sMin, 3.0));
		double const counted = histogram.CountAbove(fluxAt(s));
		CHECK(counted == doctest::Approx(expected).epsilon(0.03));
		logF.push_back(std::log10(fluxAt(s)));
		logN.push_back(std::log10(counted));
	}
	double meanF = 0.0;
	double meanN = 0.0;
	for (size_t i = 0; i < logF.size(); ++i)
	{
		meanF += logF[i] / double(logF.size());
		meanN += logN[i] / double(logN.size());
	}
	double covariance = 0.0;
	double variance = 0.0;
	for (size_t i = 0; i < logF.size(); ++i)
	{
		covariance += (logF[i] - meanF) * (logN[i] - meanN);
		variance += (logF[i] - meanF) * (logF[i] - meanF);
	}
	CHECK(covariance / variance == doctest::Approx(-1.5).epsilon(0.05 / 1.5));
}

TEST_CASE("The walk through a uniform dust ball dims each channel by exp(-k tau)")
{
	// Mutation: applying one channel's extinction to all channels reddens this test.
	Scene scene = BaseScene();
	double const alpha = 0.002;
	double const radius = 1000.0;
	scene.galaxies.push_back(UniformDust(alpha, radius));
	PreparedScene const dusty = PrepareScene(scene);
	glm::dvec3 const transmittance = TransmittanceTo(dusty, ParamsFor(Quality::Quick), {0.0, 1.0, 0.0}, 1.0e7);
	double const tau = alpha * (radius - scene.observer.sMin) / 1.0857;
	CHECK(transmittance.x == doctest::Approx(std::exp(-0.75 * tau)).epsilon(1e-9));
	CHECK(transmittance.y == doctest::Approx(std::exp(-1.00 * tau)).epsilon(1e-9));
	CHECK(transmittance.z == doctest::Approx(std::exp(-1.32 * tau)).epsilon(1e-9));
}

namespace
{

// One population's light and stars along a ray through a medium of constant density and transmittance, stepped
// as the march steps it: the shared AccumulatePart, compiled here in double.
RayResult StepThrough(PopulationTable const& population, double density, glm::dvec3 transmittance, double fluxThreshold,
	double begin, double end)
{
	MarchParams const params = ParamsFor(Quality::Final);
	RayResult result;
	double s = begin;
	while (s < end)
	{
		double const ds = std::min(std::clamp(params.alpha * s, params.dsMin, params.dsMax), end - s);
		AccumulatePart(population, density, fluxThreshold, s, ds, transmittance, 1.0, 1.0, result);
		s += ds;
	}
	return result;
}

}

TEST_CASE("A step splits its light between points and diffuse without loss")
{
	// Mutation: replacing (ds - w) with 0 in the straddling bins' diffuse radiance (step.glsl.h) loses the part of
	// each bin's light beyond its resolved distance and reddens the sum; the cook-level energy test misses it
	// (second order).
	PreparedScene const prepared = PrepareScene(MilkyWayScene({8200.0, 0.0, 20.0}));
	glm::dvec3 const transmittance(0.8, 0.7, 0.55);
	for (PopulationTable const& population : prepared.populations)
	{
		RayResult const allDiffuse = StepThrough(population, 0.1, transmittance, HUGE_VAL, 0.05, 5000.0);
		REQUIRE(Luminance(allDiffuse.radiance) > 0.0);
		CHECK(allDiffuse.resolvedRadiance == glm::dvec3(0.0));
		for (double const threshold : {1e-2, 1e-4, 1e-6})
		{
			CAPTURE(threshold);
			RayResult const split = StepThrough(population, 0.1, transmittance, threshold, 0.05, 5000.0);
			CHECK(Luminance(split.radiance) < Luminance(allDiffuse.radiance));
			CHECK(Luminance(split.resolvedRadiance) > 0.0);
			double pointCount = 0.0;
			for (double const count : split.budget)
			{
				pointCount += count;
			}
			CHECK(pointCount > 0.0);
			glm::dvec3 const total = split.radiance + split.resolvedRadiance;
			CHECK(total.x == doctest::Approx(allDiffuse.radiance.x).epsilon(1e-9));
			CHECK(total.y == doctest::Approx(allDiffuse.radiance.y).epsilon(1e-9));
			CHECK(total.z == doctest::Approx(allDiffuse.radiance.z).epsilon(1e-9));
		}
	}
}

TEST_CASE("A step's resolved stars fill the cube of their resolved distance")
{
	// Mutation: taking the straddling bins' budget over the whole step (CubeSpanFrom(a, ds)) instead of out to the
	// resolved distance reddens this test.
	Scene scene = BaseScene();
	int const bin = 12;
	int const population = AddSingleBinPopulation(scene, bin);
	PreparedScene const prepared = PrepareScene(scene);
	PopulationTable const& table = prepared.populations[size_t(population)];
	double const density = 0.02;
	double const tv = 0.6;
	double const luminosity = table.luminosity[size_t(bin)];
	double const threshold = luminosity * 100.0 * kFluxAt10pc * tv / (150.0 * 150.0);   // resolved out to 150 pc
	RayResult const result = StepThrough(table, density, glm::dvec3(tv), threshold, 0.05, 400.0);
	CHECK(result.budget[size_t(bin)] == doctest::Approx(density * (std::pow(150.0, 3.0) - std::pow(0.05, 3.0)) / 3.0).epsilon(1e-9));
	glm::dvec3 const light = density * 100.0 * kFluxAt10pc * luminosity * tv * table.chroma[size_t(bin)];
	CHECK(Luminance(result.radiance) == doctest::Approx(Luminance(light) * (400.0 - 150.0)).epsilon(1e-9));
	CHECK(Luminance(result.resolvedRadiance) == doctest::Approx(Luminance(light) * (150.0 - 0.05)).epsilon(1e-9));
}

TEST_CASE("SampleStar picks a resolved distance in proportion to the count")
{
	Scene scene = BaseScene();
	int const bin = 12;
	int const population = AddSingleBinPopulation(scene, bin);
	scene.galaxies.push_back(UniformBall(0.01, 2000.0, population, {0.0, 0.0, 0.0}));
	PreparedScene const prepared = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	double const luminosity = BinLuminosity(bin);
	double const threshold = luminosity * 100.0 * kFluxAt10pc / (100.0 * 100.0);   // resolved out to 100 pc
	// Count ∝ s², so the median distance is 100 × 2^{-1/3}.
	std::optional<StarSample> const median = SampleStar(prepared, params, {0.0, 0.0, 1.0}, bin, threshold, 0.5);
	REQUIRE(median.has_value());
	CHECK(median->distance == doctest::Approx(100.0 * std::pow(0.5, 1.0 / 3.0)).epsilon(0.02));
	CHECK(median->population == population);
	std::optional<StarSample> const farthest = SampleStar(prepared, params, {0.0, 0.0, 1.0}, bin, threshold, 0.999999);
	REQUIRE(farthest.has_value());
	CHECK(farthest->distance <= 100.0 * 1.06);
	CHECK_FALSE(SampleStar(prepared, params, {0.0, 0.0, 1.0}, bin, HUGE_VAL, 0.5).has_value());
}

TEST_CASE("Star saturation moves colour, not light or counts")
{
	// Mutation: PrepareScene not passing starSaturation to Tabulate reddens the spread check.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	PreparedScene const physical = PrepareScene(scene);
	scene.style.starSaturation = 4.0;
	PreparedScene const vivid = PrepareScene(scene);
	// Grey transmittance: only then is luminance exactly kept.
	glm::dvec3 const grey(0.7);
	RayResult a;
	RayResult b;
	for (size_t p = 0; p < physical.populations.size(); ++p)
	{
		RayResult const one = StepThrough(physical.populations[p], 0.1, grey, 1e-4, 0.05, 5000.0);
		RayResult const other = StepThrough(vivid.populations[p], 0.1, grey, 1e-4, 0.05, 5000.0);
		a.radiance += one.radiance;
		b.radiance += other.radiance;
		for (int bin = 0; bin < kBins; ++bin)
		{
			a.budget[size_t(bin)] += one.budget[size_t(bin)];
			b.budget[size_t(bin)] += other.budget[size_t(bin)];
		}
	}
	CHECK(Luminance(b.radiance) == doctest::Approx(Luminance(a.radiance)).epsilon(1e-12));
	CHECK(b.budget == a.budget);
	auto spread = [](glm::dvec3 c) { return (std::max({c.x, c.y, c.z}) - std::min({c.x, c.y, c.z})) / Luminance(c); };
	CHECK(spread(b.radiance) > 2.0 * spread(a.radiance));
}

TEST_CASE("Gas leaves the calibration's star counts alone")
{
	// Mutation: PrepareScene not setting youngLight reddens the first check; any light term that leaked into a part's
	// density would move the histogram. The GPU cook test "Gas adds diffuse light and leaves the stars alone" holds the
	// light itself.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	Galaxy& galaxy = scene.galaxies[0];
	galaxy.dust.bubbleCenter = {8200.0, 0.0, 20.0};
	galaxy.dust.bubbleRadius = 300.0;
	galaxy.dust.bubbleHeight = 800.0;
	PreparedScene const plain = PrepareScene(scene);
	galaxy.disk.gas.glow = 1.5;
	galaxy.disk.gas.palette = {glm::dvec3(1.0, 0.18, 0.32), glm::dvec3(0.25, 1.0, 0.85)};
	galaxy.disk.gas.paletteCount = 2;
	galaxy.disk.gas.largeGlow = 1.5;
	galaxy.dust.bubbleRim = 4.0;
	galaxy.dust.reflection = 2.0;
	PreparedScene const glowing = PrepareScene(scene);
	CHECK(glowing.galaxies[0].galaxy.disk.gas.youngLight == doctest::Approx(Luminance(glowing.populations[0].lightBelow[kBins])));
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir = glm::normalize(glm::dvec3(-1.0, 0.1, 0.02));
	FluxHistogram ha;
	FluxHistogram hb;
	HistogramRay(plain, params, dir, 1.0, ha);
	HistogramRay(glowing, params, dir, 1.0, hb);
	CHECK(ha.count == hb.count);
}
