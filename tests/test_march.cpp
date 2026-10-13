#include "march.h"
#include "test_scenes.h"
#include <doctest.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <glm/geometric.hpp>

using namespace starscape;
using namespace starscape::test;

namespace
{

// Total sky flux of everything within `halfAngle` of +x (sky frame = scene frame here):
// a polar grid fine enough that the integrand is smooth inside each cell.
glm::dvec3 FluxInCone(PreparedScene const& prepared, MarchParams const& params, double halfAngle)
{
	constexpr int thetaSteps = 200;
	constexpr int phiSteps = 64;
	glm::dvec3 total(0.0);
	double const dTheta = halfAngle / thetaSteps;
	double const dPhi = 2.0 * kPi / phiSteps;
	for (int i = 0; i < thetaSteps; ++i)
	{
		double const theta = (i + 0.5) * dTheta;
		for (int j = 0; j < phiSteps; ++j)
		{
			double const phi = (j + 0.5) * dPhi;
			glm::dvec3 const dir(std::cos(theta), std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi));
			RayResult result;
			MarchRay(prepared, params, dir, HUGE_VAL, 0.0, 1.0, result);
			total += result.radiance * (std::sin(theta) * dTheta * dPhi);
		}
	}
	return total;
}

}

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

TEST_CASE("A distant dust-free galaxy delivers L (10/d)² C in total")
{
	// Mutation: restoring an unbounded density term reddens this test.
	Scene scene = BaseScene();
	int const bin = 14;
	int const population = AddSingleBinPopulation(scene, bin);
	double const density = 0.05;
	double const scale = 5000.0;
	double const cutoff = 20000.0;
	double const distance = 1.0e6;
	Galaxy galaxy;
	galaxy.position = {distance, 0.0, 0.0};
	galaxy.halo = {.enabled = true, .density = density, .scale = scale, .cutoff = cutoff, .population = population};
	scene.galaxies.push_back(galaxy);
	PreparedScene const prepared = PrepareScene(scene);

	// ∫₀^R 4π r² (1 + r/a)^-3.5 dr = 4π a³ [−2x^-½ + (4/3)x^-3/2 − 0.4x^-5/2]₁^{1+R/a}.
	auto antiderivative = [](double x) { return -2.0 * std::pow(x, -0.5) + 4.0 / 3.0 * std::pow(x, -1.5) - 0.4 * std::pow(x, -2.5); };
	double const stars = density * 4.0 * kPi * scale * scale * scale * (antiderivative(1.0 + cutoff / scale) - antiderivative(1.0));
	double const expected = stars * BinLuminosity(bin) * 100.0 * kFluxAt10pc / (distance * distance);

	glm::dvec3 const flux = FluxInCone(prepared, ParamsFor(Quality::Quick), std::asin(cutoff / distance) * 1.01);
	CHECK(Luminance(flux) == doctest::Approx(expected).epsilon(0.01));
}

TEST_CASE("A uniform dust ball dims each channel by exp(-k tau)")
{
	// Mutation: applying one channel's extinction to all channels reddens this test.
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 14);
	Galaxy source;
	source.position = {1.0e6, 0.0, 0.0};
	source.halo = {.enabled = true, .density = 0.05, .scale = 5000.0, .cutoff = 20000.0, .population = population};
	scene.galaxies.push_back(source);
	MarchParams const params = ParamsFor(Quality::Quick);
	double const halfAngle = std::asin(0.02) * 1.01;
	glm::dvec3 const clear = FluxInCone(PrepareScene(scene), params, halfAngle);

	double const alpha = 0.002;
	double const radius = 1000.0;
	scene.galaxies.push_back(UniformDust(alpha, radius));
	PreparedScene const dusty = PrepareScene(scene);
	glm::dvec3 const dimmed = FluxInCone(dusty, params, halfAngle);

	double const tau = alpha * (radius - scene.observer.sMin) / 1.0857;
	CHECK(dimmed.x / clear.x == doctest::Approx(std::exp(-0.75 * tau)).epsilon(0.005).scale(0.0));
	CHECK(dimmed.y / clear.y == doctest::Approx(std::exp(-1.00 * tau)).epsilon(0.005).scale(0.0));
	CHECK(dimmed.z / clear.z == doctest::Approx(std::exp(-1.32 * tau)).epsilon(0.005).scale(0.0));
	CHECK(TransmittanceTo(dusty, params, {0.0, 1.0, 0.0}, 1.0e7).y == doctest::Approx(std::exp(-tau)).epsilon(0.005).scale(0.0));
}

TEST_CASE("Points and diffuse split one ray's light without loss")
{
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	PreparedScene const prepared = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir = glm::normalize(glm::dvec3(-1.0, 0.1, 0.02));
	RayResult allDiffuse;
	MarchRay(prepared, params, dir, HUGE_VAL, 1.0, 1.0, allDiffuse);
	for (double const threshold : {1e-2, 1e-4, 1e-6})
	{
		RayResult split;
		MarchRay(prepared, params, dir, threshold, 1.0, 1.0, split);
		CHECK(Luminance(split.radiance) < Luminance(allDiffuse.radiance));
		double pointCount = 0.0;
		for (double const count : split.budget)
		{
			pointCount += count;
		}
		CHECK(pointCount > 0.0);
	}
}

TEST_CASE("Diffuse plus resolved light of one ray is the same at every threshold")
{
	// Mutation: replacing (step.b - c) with 0 in the straddling bins' diffuse radiance (march.cpp) loses the part of
	// each bin's light beyond its resolved distance and reddens this; the totals-only tests miss it (second order).
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	PreparedScene const prepared = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir = glm::normalize(glm::dvec3(-1.0, 0.1, 0.02));
	RayResult allDiffuse;
	MarchRay(prepared, params, dir, HUGE_VAL, 1.0, 1.0, allDiffuse);
	CHECK(allDiffuse.resolvedRadiance == glm::dvec3(0.0));
	for (double const threshold : {1e-2, 1e-4, 1e-6})
	{
		RayResult split;
		MarchRay(prepared, params, dir, threshold, 1.0, 1.0, split);
		CHECK(Luminance(split.resolvedRadiance) > 0.0);
		CHECK(Luminance(split.radiance) > 0.0);
		glm::dvec3 const total = split.radiance + split.resolvedRadiance;
		CHECK(total.x == doctest::Approx(allDiffuse.radiance.x).epsilon(1e-9));
		CHECK(total.y == doctest::Approx(allDiffuse.radiance.y).epsilon(1e-9));
		CHECK(total.z == doctest::Approx(allDiffuse.radiance.z).epsilon(1e-9));
	}
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
	scene.galaxies[0].dust.enabled = false;   // grey transmittance: only then is luminance exactly kept
	PreparedScene const physical = PrepareScene(scene);
	scene.style.starSaturation = 4.0;
	PreparedScene const vivid = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir = glm::normalize(glm::dvec3(-1.0, 0.1, 0.02));
	RayResult a;
	RayResult b;
	MarchRay(physical, params, dir, 1e-4, 1.0, 1.0, a);
	MarchRay(vivid, params, dir, 1e-4, 1.0, 1.0, b);
	CHECK(Luminance(b.radiance) == doctest::Approx(Luminance(a.radiance)).epsilon(1e-12));
	CHECK(b.budget == a.budget);
	auto spread = [](glm::dvec3 c) { return (std::max({c.x, c.y, c.z}) - std::min({c.x, c.y, c.z})) / Luminance(c); };
	CHECK(spread(b.radiance) > 2.0 * spread(a.radiance));
}
