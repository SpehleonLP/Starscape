#include "galaxy.h"
#include "presets.h"
#include <doctest.h>
#include <cmath>

using namespace starscape;

namespace
{
constexpr double kPi = 3.14159265358979323846;

double TotalDensity(GalaxySample const& sample)
{
	double total = 0.0;
	for (int i = 0; i < sample.count; ++i)
	{
		total += sample.parts[size_t(i)].density;
	}
	return total;
}
}

TEST_CASE("Arm proximity peaks on the arm and falls as a Gaussian of arc length")
{
	Arms const arms{.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 2.0, .clumping = 0.0, .phaseDeg = 10.0, .r0 = 1000.0};
	double const radius = 5000.0;
	double const onArm = 10.0 * kPi / 180.0 + std::log(radius / arms.r0) / std::tan(arms.pitchDeg * kPi / 180.0);
	CHECK(ArmProximity(arms, radius, onArm, 0.0) == doctest::Approx(1.0));
	CHECK(ArmProximity(arms, radius, onArm + kPi, 0.0) == doctest::Approx(1.0));   // the second arm
	CHECK(ArmProximity(arms, radius, onArm + arms.width / radius, 0.0) == doctest::Approx(std::exp(-0.5)));
	// A positive offset moves the peak toward the inner edge, which lies at positive azimuth.
	CHECK(ArmProximity(arms, radius, onArm + 150.0 / radius, 150.0) == doctest::Approx(1.0));
	CHECK(ArmProximity(Arms{}, radius, onArm, 0.0) == 0.0);
}

TEST_CASE("Evaluation at the centre and on the axis is finite")
{
	std::vector<Population> const populations = PresetPopulations();
	for (char const* name : {"milky-way", "small-spiral", "elliptical", "dwarf"})
	{
		std::optional<Galaxy> const galaxy = PresetGalaxy(name, populations);
		REQUIRE(galaxy.has_value());
		for (glm::dvec3 const p : {glm::dvec3(0.0), glm::dvec3(0.0, 0.0, 300.0), glm::dvec3(1e-9, 0.0, 0.0)})
		{
			GalaxySample sample;
			EvaluateLocal(*galaxy, p, sample);
			CHECK(std::isfinite(TotalDensity(sample)));
			CHECK(std::isfinite(sample.extinction));
		}
	}
}

TEST_CASE("The bulge integrates to the Hernquist total")
{
	Galaxy galaxy;
	galaxy.bulge = {.enabled = true, .density = 1000.0, .scale = 500.0, .axisY = 1.0, .axisZ = 1.0, .barAngleDeg = 0.0, .cutoff = 5000.0, .population = 0};
	// ∫₀^{10a} 4π r² n_b a / (r (r + a)³) dr = 4π n_b × 100 / 242.
	double const expected = 4.0 * kPi * galaxy.bulge.density * 100.0 / 242.0;
	double total = 0.0;
	constexpr int steps = 200000;
	double const dr = galaxy.bulge.cutoff / steps;
	for (int i = 0; i < steps; ++i)
	{
		double const r = (i + 0.5) * dr;
		GalaxySample sample;
		EvaluateLocal(galaxy, {r, 0.0, 0.0}, sample);
		total += 4.0 * kPi * r * r * TotalDensity(sample) * dr;
	}
	CHECK(total == doctest::Approx(expected).epsilon(0.005));
}

TEST_CASE("Disks fall off as sech² in height and exponentially in radius")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 2.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0, .arms = {}, .population = 0};
	GalaxySample mid;
	GalaxySample high;
	GalaxySample far;
	EvaluateLocal(galaxy, {1000.0, 0.0, 0.0}, mid);
	EvaluateLocal(galaxy, {1000.0, 0.0, 600.0}, high);
	EvaluateLocal(galaxy, {4000.0, 0.0, 0.0}, far);
	CHECK(TotalDensity(high) / TotalDensity(mid) == doctest::Approx(Sech2(1.0)));
	CHECK(TotalDensity(far) / TotalDensity(mid) == doctest::Approx(std::exp(-1.0)));
	GalaxySample outside;
	EvaluateLocal(galaxy, {30001.0, 0.0, 0.0}, outside);
	CHECK(outside.count == 0);
}

TEST_CASE("Arms split the thin disk between its old and young populations")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0,
		.arms = {.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 3.0, .clumping = 0.0, .phaseDeg = 0.0, .r0 = 1000.0},
		.population = 1, .armPopulation = 0, .armYoungFraction = 0.1};
	double const radius = 5000.0;
	double const onArm = std::log(radius / 1000.0) / std::tan(15.0 * kPi / 180.0);
	GalaxySample sample;
	EvaluateLocal(galaxy, {radius * std::cos(onArm), radius * std::sin(onArm), 0.0}, sample);
	REQUIRE(sample.count == 2);
	CHECK(sample.parts[1].population == 0);
	CHECK(sample.parts[1].density / TotalDensity(sample) == doctest::Approx(0.1));
	CHECK(TotalDensity(sample) == doctest::Approx(3.0 * std::exp(-radius / 3000.0)));
}

TEST_CASE("Rotation maps scene positions into the galaxy frame")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0, .arms = {}, .population = 0};
	galaxy.position = {100.0, 0.0, 0.0};
	galaxy.rotationDeg = {90.0, 0.0, 0.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	GalaxySample rotated;
	GalaxySample reference;
	Evaluate(prepared, {100.0, 0.0, 2000.0}, rotated);   // in the rotated disk's plane
	EvaluateLocal(galaxy, {0.0, 2000.0, 0.0}, reference);
	CHECK(TotalDensity(rotated) == doctest::Approx(TotalDensity(reference)));
}

TEST_CASE("Dust is a uniform field when its scales are huge")
{
	Galaxy galaxy;
	galaxy.dust = {.enabled = true, .extinction = 0.002, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1000.0, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0};
	GalaxySample sample;
	EvaluateLocal(galaxy, {300.0, -200.0, 100.0}, sample);
	CHECK(sample.extinction == doctest::Approx(0.002));
	CHECK(BoundRadius(galaxy) == 1000.0);
}

TEST_CASE("Galaxy presets resolve their populations")
{
	std::vector<Population> const populations = PresetPopulations();
	std::optional<Galaxy> const milkyWay = PresetGalaxy("milky-way", populations);
	REQUIRE(milkyWay.has_value());
	CHECK(milkyWay->disk.population == 1);
	CHECK(milkyWay->disk.armPopulation == 0);
	CHECK(milkyWay->bulge.population == 2);
	CHECK(milkyWay->halo.population == 3);
	CHECK(BoundRadius(*milkyWay) == 30000.0);
	CHECK_FALSE(PresetGalaxy("andromeda", populations).has_value());
}
