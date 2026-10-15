#include "galaxy.h"
#include "population.h"
#include "presets.h"
#include "rng.h"
#include <doctest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <glm/geometric.hpp>

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
	// Mutation: removing the axis guard reddens this test.
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
	// Mutation: a wrong Hernquist power reddens this test.
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

// Mutation: give the arm population the disk's scale height and the off-plane ratio check fails.
TEST_CASE("The arm population keeps its own, thinner scale height")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0,
		.arms = {.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 3.0, .clumping = 0.0, .phaseDeg = 0.0, .r0 = 1000.0},
		.population = 1, .armPopulation = 0, .armYoungFraction = 0.1, .armScaleHeight = 90.0};
	double const radius = 5000.0;
	double const onArm = std::log(radius / 1000.0) / std::tan(15.0 * kPi / 180.0);
	GalaxySample plane;
	GalaxySample above;
	EvaluateLocal(galaxy, {radius * std::cos(onArm), radius * std::sin(onArm), 0.0}, plane);
	EvaluateLocal(galaxy, {radius * std::cos(onArm), radius * std::sin(onArm), 400.0}, above);
	REQUIRE(plane.count == 2);
	REQUIRE(above.count == 2);
	CHECK(above.parts[0].density / plane.parts[0].density == doctest::Approx(Sech2(400.0 / 600.0)));
	CHECK(above.parts[1].density / plane.parts[1].density == doctest::Approx(Sech2(400.0 / 180.0)));
}

// Mutation: drop youngFraction from the young share and the between-arms check fails.
TEST_CASE("Young stars sit between the arms too")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0,
		.arms = {.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 3.0, .clumping = 0.0, .phaseDeg = 0.0, .r0 = 1000.0},
		.population = 1, .armPopulation = 0, .armYoungFraction = 0.1, .armScaleHeight = 90.0, .youngFraction = 0.01};
	double const radius = 5000.0;
	// Half an arm period from the crest: ArmProximity is ~e^-120 there.
	double const between = std::log(radius / 1000.0) / std::tan(15.0 * kPi / 180.0) + kPi / 2.0;
	GalaxySample sample;
	EvaluateLocal(galaxy, {radius * std::cos(between), radius * std::sin(between), 0.0}, sample);
	REQUIRE(sample.count == 2);
	CHECK(sample.parts[1].density / TotalDensity(sample) == doctest::Approx(0.01));
}

TEST_CASE("Rotation maps scene positions into the galaxy frame")
{
	// Mutation: dropping the transpose or reordering the axes reddens this test.
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0,
		.arms = {.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 3.0, .clumping = 0.0, .phaseDeg = 0.0, .r0 = 1000.0},
		.population = 0};
	galaxy.bulge = {.enabled = true, .density = 1000.0, .scale = 500.0, .axisY = 0.5, .axisZ = 0.3, .barAngleDeg = 30.0, .cutoff = 5000.0, .population = 0};
	galaxy.position = {100.0, -50.0, 25.0};
	galaxy.rotationDeg = {30.0, 50.0, 20.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);

	// Built independently of PrepareGalaxy: the documented order is Z, then Y, then X.
	auto rx = [](double a) { return glm::dmat3(1, 0, 0, 0, std::cos(a), std::sin(a), 0, -std::sin(a), std::cos(a)); };
	auto ry = [](double a) { return glm::dmat3(std::cos(a), 0, -std::sin(a), 0, 1, 0, std::sin(a), 0, std::cos(a)); };
	auto rz = [](double a) { return glm::dmat3(std::cos(a), std::sin(a), 0, -std::sin(a), std::cos(a), 0, 0, 0, 1); };
	double const d = kPi / 180.0;
	glm::dmat3 const toScene = rx(30.0 * d) * ry(50.0 * d) * rz(20.0 * d);

	for (glm::dvec3 const local : {glm::dvec3(1500.0, 700.0, 100.0), glm::dvec3(-800.0, 1200.0, -150.0), glm::dvec3(300.0, -400.0, 200.0), glm::dvec3(4000.0, 2500.0, 50.0)})
	{
		GalaxySample rotated;
		GalaxySample reference;
		Evaluate(prepared, galaxy.position + toScene * local, rotated);
		EvaluateLocal(galaxy, local, reference);
		CHECK(TotalDensity(reference) > 0.0);
		CHECK(TotalDensity(rotated) == doctest::Approx(TotalDensity(reference)));
	}
}

TEST_CASE("Two disks that split by arm population fill six parts")
{
	Galaxy galaxy;
	Arms const arms{.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 3.0, .clumping = 0.0, .phaseDeg = 0.0, .r0 = 1000.0};
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0, .arms = arms,
		.population = 1, .armPopulation = 0, .armYoungFraction = 0.1};
	galaxy.thickDisk = {.enabled = true, .density = 0.1, .scaleLength = 3600.0, .scaleHeight = 900.0, .cutoff = 30000.0, .arms = arms,
		.population = 1, .armPopulation = 0, .armYoungFraction = 0.1};
	galaxy.bulge = {.enabled = true, .density = 1000.0, .scale = 500.0, .axisY = 1.0, .axisZ = 1.0, .barAngleDeg = 0.0, .cutoff = 30000.0, .population = 2};
	galaxy.halo = {.enabled = true, .density = 0.1, .scale = 1000.0, .cutoff = 30000.0, .population = 3};
	double const radius = 5000.0;
	double const onArm = std::log(radius / 1000.0) / std::tan(15.0 * kPi / 180.0);
	GalaxySample sample;
	EvaluateLocal(galaxy, {radius * std::cos(onArm), radius * std::sin(onArm), 0.0}, sample);
	CHECK(sample.count == 6);
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

TEST_CASE("Dust filaments move extinction around without changing its mean")
{
	// Mutation: dropping the division by filamentMean leaves the mean near 0.65 α₀ at σ 2.
	Galaxy galaxy;
	galaxy.dust = {.enabled = true, .extinction = 0.002, .scaleLength = 1e12, .scaleHeight = 1e4, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 2.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	RandomStream rng(5, Purpose::Test, 0, 0, 0);
	double sum = 0.0;
	double lowest = 1.0;
	int const count = 100000;
	for (int i = 0; i < count; ++i)
	{
		// A slab eight noise cells thick: lattice noise is not stationary within a cell, so any one
		// plane reads a few percent off the volume mean the normalisation holds.
		glm::dvec3 const p(rng.NextDouble() * 1e7 - 5e6, rng.NextDouble() * 1e7 - 5e6, rng.NextDouble() * 4e4 - 2e4);
		GalaxySample sample;
		Evaluate(prepared, p, sample);
		double const field = sample.extinction / (0.002 * Sech2(p.z / 2e4));
		sum += field;
		lowest = std::min(lowest, field);
	}
	CHECK(std::abs(sum / count - 1.0) < 0.03);
	CHECK(lowest < 0.1);   // the structure is real: some sight lines nearly clear
}

TEST_CASE("Associations gather young stars without changing how many there are")
{
	// Mutation: dropping the division by associationMean moves the mean by 16% at σ 2.
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .arms = {}, .population = 1, .armPopulation = 0,
		.youngFraction = 0.01, .associationSize = 150.0, .associationClumping = 2.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	RandomStream rng(7, Purpose::Test, 0, 0, 0);
	double young = 0.0;
	double old = 0.0;
	double densest = 0.0;
	int split = 0;
	int const count = 200000;
	for (int i = 0; i < count; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		Evaluate(prepared, p, sample);
		split += sample.count == 2 ? 1 : 0;
		old += sample.parts[0].density;
		young += sample.parts[1].density;
		densest = std::max(densest, sample.parts[1].density);
	}
	REQUIRE(split == count);
	CHECK(std::abs(young / count / 0.01 - 1.0) < 0.03);
	CHECK(old / count == doctest::Approx(0.99));   // associations add young stars, they do not convert old ones
	CHECK(densest > 0.1);                          // the structure is real: groups ten times the mean density
}

TEST_CASE("A dust bubble empties its core and leaves the outside alone")
{
	// Mutation: dropping the bubble leaves the core at full extinction.
	Galaxy galaxy;
	galaxy.dust = {.enabled = true, .extinction = 0.002, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0,
		.bubbleCenter = {100.0, 0.0, 0.0}, .bubbleRadius = 300.0, .bubbleHeight = 800.0};
	Galaxy open = galaxy;
	open.dust.bubbleRadius = 0.0;
	auto extinction = [](Galaxy const& g, glm::dvec3 p)
	{
		GalaxySample sample;
		EvaluateLocal(g, p, sample);
		return sample.extinction;
	};
	CHECK(extinction(galaxy, {100.0, 0.0, 0.0}) == 0.0);
	CHECK(extinction(galaxy, {100.0, 0.0, 500.0}) == 0.0);   // 0.625 of the vertical semi-axis: inside the chimney
	CHECK(extinction(galaxy, {100.0, 250.0, 0.0}) > 0.0);    // 0.83 of the radius: on the ramp
	CHECK(extinction(galaxy, {100.0, 250.0, 0.0}) < 0.5 * extinction(open, {100.0, 250.0, 0.0}));
	for (glm::dvec3 const outside : {glm::dvec3(100.0, 0.0, 801.0), glm::dvec3(401.0, 0.0, 0.0)})
	{
		CHECK(extinction(open, outside) > 0.0);
		CHECK(extinction(galaxy, outside) == extinction(open, outside));
	}
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

namespace
{

// A disk uniform over any test box, all young stars lighting gas; dust off unless a test turns it on.
Galaxy GasDisk(double glow)
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .arms = {}, .population = 1, .armPopulation = 0,
		.youngFraction = 0.01};
	galaxy.disk.gas.glow = glow;
	galaxy.disk.gas.youngLight = 1.0;   // PrepareScene fills this from the population table
	return galaxy;
}

}

TEST_CASE("Gas glows at glow times the young stars' light, on average")
{
	// Mutation: dropping the division by textureMean leaves the ratio at ⟨G⟩ ≈ 0.88 of glow at filament 2.
	Galaxy galaxy = GasDisk(2.0);
	// Dust present only for its filaments; scale height 400 gives the Milky Way's 200 pc filament wavelength.
	galaxy.dust = {.enabled = true, .extinction = 0.001, .scaleLength = 1e12, .scaleHeight = 400.0, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 2.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	RandomStream rng(11, Purpose::Test, 0, 0, 0);
	double gas = 0.0;
	double young = 0.0;
	int const count = 200000;
	for (int i = 0; i < count; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		Evaluate(prepared, p, sample);
		REQUIRE(sample.count == 2);
		gas += Luminance(sample.gas);
		young += sample.parts[1].density * 100.0 * kFluxAt10pc;
	}
	CHECK(gas / young == doctest::Approx(2.0).epsilon(0.03));
}

TEST_CASE("Gas shares the dust's bubble")
{
	// Mutation: leaving the bubble factor off the gas lights the cavity.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.dust = {.enabled = true, .extinction = 0.002, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0,
		.bubbleCenter = {100.0, 0.0, 0.0}, .bubbleRadius = 300.0, .bubbleHeight = 800.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	GalaxySample inside;
	Evaluate(prepared, {100.0, 0.0, 0.0}, inside);
	REQUIRE(inside.count == 2);   // young stars are there to light it
	CHECK(inside.gas == glm::dvec3(0.0));
	GalaxySample outside;
	Evaluate(prepared, {100.0, 0.0, 900.0}, outside);
	CHECK(Luminance(outside.gas) > 0.0);
}

TEST_CASE("Gas colours are normalised and the densest clumps take the core colour")
{
	// Mutation: swapping the mix's endpoints makes pure-core samples the majority and reddens the ratio check.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.disk.gas.colour = {2.0, 0.0, 0.0};
	galaxy.disk.gas.coreColour = {0.0, 0.0, 3.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	CHECK(Luminance(prepared.galaxy.disk.gas.colour) == doctest::Approx(1.0));
	CHECK(Luminance(prepared.galaxy.disk.gas.coreColour) == doctest::Approx(1.0));
	CHECK(prepared.galaxy.disk.gas.colour.y == 0.0);
	CHECK(prepared.galaxy.disk.gas.coreColour.x == 0.0);
	RandomStream rng(13, Purpose::Test, 0, 0, 0);
	int red = 0;
	int blue = 0;
	for (int i = 0; i < 20000; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		Evaluate(prepared, p, sample);
		REQUIRE(sample.gas.y == 0.0);
		red += sample.gas.x > 0.0 && sample.gas.z == 0.0 ? 1 : 0;
		blue += sample.gas.z > 0.0 && sample.gas.x == 0.0 ? 1 : 0;
	}
	// Clumps are LN(1.2): below 1 for ~73% of points, above 3 for ~6.5%.
	CHECK(blue > 0);
	CHECK(red > 5 * blue);
}

TEST_CASE("Without glow there is no gas")
{
	GalaxySample sample;
	EvaluateLocal(GasDisk(0.0), {300.0, -200.0, 100.0}, sample);
	CHECK(sample.gas == glm::dvec3(0.0));
}

TEST_CASE("Disabled dust leaves the gas texture alone")
{
	// Mutation: dropping the dust.enabled guard in GasTexture lets a disabled dust's filaments shape the gas
	// and its measured mean; the attenuation test cannot see that, since it compares ratios.
	Galaxy const plain = GasDisk(1.0);
	Galaxy ghost = GasDisk(1.0);
	ghost.dust = {.enabled = false, .extinction = 0.002, .scaleLength = 1e12, .scaleHeight = 400.0, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 3.0};
	PreparedGalaxy const a = PrepareGalaxy(plain);
	PreparedGalaxy const b = PrepareGalaxy(ghost);
	CHECK(b.galaxy.disk.gas.textureMean == a.galaxy.disk.gas.textureMean);
	for (glm::dvec3 const p : {glm::dvec3(300.0, -200.0, 100.0), glm::dvec3(-4100.0, 900.0, -30.0), glm::dvec3(12.5, 7000.0, 260.0)})
	{
		GalaxySample sa;
		GalaxySample sb;
		Evaluate(a, p, sa);
		Evaluate(b, p, sb);
		CHECK(Luminance(sa.gas) > 0.0);
		CHECK(sb.gas == sa.gas);
	}
}

TEST_CASE("A palette colours the gas by a large-scale field without changing its light")
{
	// Mutation: dropping the per-entry normalisation in NormaliseColours leaves entries off Y = 1.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.disk.gas.palette = {glm::dvec3(2.0, 0.0, 0.0), glm::dvec3(0.0, 3.0, 0.0), glm::dvec3(0.0, 0.0, 4.0), glm::dvec3(1.0, 1.0, 0.0)};
	galaxy.disk.gas.paletteCount = 4;
	galaxy.disk.gas.paletteWavelength = 1500.0;
	// A core colour equal to nothing in the palette would hide which entry a point took; make it one.
	galaxy.disk.gas.coreColour = galaxy.disk.gas.palette[0];
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	for (int i = 0; i < 4; ++i)
	{
		CHECK(Luminance(prepared.galaxy.disk.gas.palette[size_t(i)]) == doctest::Approx(1.0).epsilon(1e-12));
	}
	Galaxy plain = galaxy;
	plain.disk.gas.paletteCount = 0;
	PreparedGalaxy const reference = PrepareGalaxy(plain);
	RandomStream rng(17, Purpose::Test, 0, 0, 0);
	std::array<int, 4> nearest{};
	int const count = 100000;
	for (int i = 0; i < count; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		GalaxySample without;
		Evaluate(prepared, p, sample);
		Evaluate(reference, p, without);
		REQUIRE(Luminance(without.gas) > 0.0);
		// Colour never changes the amount of light: the palette's gas is exactly as bright as the plain gas.
		REQUIRE(Luminance(sample.gas) == doctest::Approx(Luminance(without.gas)).epsilon(1e-12));
		glm::dvec3 const colour = sample.gas / Luminance(sample.gas);
		int best = 0;
		for (int e = 1; e < 4; ++e)
		{
			if (glm::distance(colour, prepared.galaxy.disk.gas.palette[size_t(e)]) < glm::distance(colour, prepared.galaxy.disk.gas.palette[size_t(best)]))
			{
				best = e;
			}
		}
		++nearest[size_t(best)];
	}
	// Every entry appears. The spec asked for 10% each, but with 0.9 entries per σ the entry two steps
	// from the field's mean is reached only in the tails (|z| > 1.1), so the floor is 1%. Task 7 corrects
	// the spec's row when it folds in.
	for (int e = 0; e < 4; ++e)
	{
		CHECK(nearest[size_t(e)] > count / 100);
	}
}

TEST_CASE("A three-colour palette wraps from its last entry to its first")
{
	// Mutation: wrapping with `% kMaxPalette` instead of `% paletteCount` blends into a zeroed entry, and
	// the gas's luminance drops below one.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.disk.gas.palette = {glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(0.0, 1.0, 0.0), glm::dvec3(0.0, 0.0, 1.0)};
	galaxy.disk.gas.paletteCount = 3;
	// A green core: any core mixing then shows as green, so red with blue and no green can only be the wrap.
	galaxy.disk.gas.coreColour = {0.0, 1.0, 0.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	Galaxy plain = galaxy;
	plain.disk.gas.paletteCount = 0;
	PreparedGalaxy const reference = PrepareGalaxy(plain);
	RandomStream rng(43, Purpose::Test, 0, 0, 0);
	bool blendedLastIntoFirst = false;
	for (int i = 0; i < 100000; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		GalaxySample without;
		Evaluate(prepared, p, sample);
		Evaluate(reference, p, without);
		REQUIRE(Luminance(without.gas) > 0.0);
		// Colour never changes the amount of light: the palette's gas is exactly as bright as the plain gas.
		REQUIRE(Luminance(sample.gas) == doctest::Approx(Luminance(without.gas)).epsilon(1e-12));
		blendedLastIntoFirst = blendedLastIntoFirst || (sample.gas.x > 0.0 && sample.gas.z > 0.0 && sample.gas.y == 0.0);
	}
	CHECK(blendedLastIntoFirst);
}

TEST_CASE("The thick disk never takes the thin disk's palette")
{
	// Mutation: passing galaxy.disk.gas instead of the thick disk's own gas into its gas block lights it;
	// a thick disk alone must stay dark.
	Galaxy galaxy;
	galaxy.thickDisk = {.enabled = true, .density = 1.0, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .arms = {}, .population = 1,
		.armPopulation = 1, .youngFraction = 0.5};
	galaxy.disk.gas.glow = 1.0;
	galaxy.disk.gas.youngLight = 1.0;
	galaxy.disk.gas.paletteCount = 2;
	galaxy.disk.gas.palette = {glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0)};
	GalaxySample sample;
	EvaluateLocal(PrepareGalaxy(galaxy).galaxy, {300.0, -200.0, 100.0}, sample);
	// Without young stars the gas block never runs and the mutation could not show.
	CHECK(sample.count >= 2);
	CHECK(sample.gas == glm::dvec3(0.0));
}

TEST_CASE("Large-scale gas glows at large_glow times the young stars' light, on average")
{
	// Mutation: lighting the term by the disk's whole midplane density instead of its young density puts
	// the ratio at large_glow / young_fraction. The second wavelength spans few noise cells across the box,
	// which only a mean measured in noise space keeps calibrated.
	for (double const wavelength : {100.0, 20000.0})
	{
		Galaxy galaxy = GasDisk(0.0);
		galaxy.disk.gas.largeGlow = 2.0;
		galaxy.disk.gas.largeWavelength = wavelength;
		PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
		RandomStream rng(19, Purpose::Test, 0, 0, 0);
		double gas = 0.0;
		double young = 0.0;
		for (int i = 0; i < 200000; ++i)
		{
			glm::dvec3 const p(rng.NextDouble() * 1e8 - 5e7, rng.NextDouble() * 1e8 - 5e7, rng.NextDouble() * 1e8 - 5e7);
			GalaxySample sample;
			Evaluate(prepared, p, sample);
			REQUIRE(sample.count == 2);
			REQUIRE(std::isfinite(Luminance(sample.gas)));
			gas += Luminance(sample.gas);
			young += sample.parts[1].density * 100.0 * kFluxAt10pc;
		}
		CHECK(gas / young == doctest::Approx(2.0).epsilon(0.03));
	}
}

TEST_CASE("The large-scale field's mean is measured in noise space and divided out")
{
	// Mutations: replacing LargeGasMean's walk with a constant 1 reddens the first three checks; dropping
	// the division by largeMean in EvaluateLocal reddens the last. At σ 1 the mean is 1.007, too close to
	// one for the calibration test's 3% to see, so the division is pinned directly.
	Galaxy galaxy = GasDisk(0.0);
	galaxy.disk.gas.largeGlow = 1.0;
	PreparedGalaxy const a = PrepareGalaxy(galaxy);
	CHECK(a.galaxy.disk.gas.largeMean == doctest::Approx(1.006899).epsilon(1e-5));   // measured at seed 1, 2^18 points
	galaxy.seed = 2;
	CHECK(PrepareGalaxy(galaxy).galaxy.disk.gas.largeMean != a.galaxy.disk.gas.largeMean);
	galaxy.seed = 1;
	galaxy.disk.gas.largeWavelength = 5000.0;
	CHECK(PrepareGalaxy(galaxy).galaxy.disk.gas.largeMean == a.galaxy.disk.gas.largeMean);   // noise space: no wavelength
	PreparedGalaxy doubled = a;
	doubled.galaxy.disk.gas.largeMean *= 2.0;
	glm::dvec3 const p(300.0, -200.0, 100.0);
	GalaxySample once;
	GalaxySample halved;
	Evaluate(a, p, once);
	Evaluate(doubled, p, halved);
	REQUIRE(Luminance(once.gas) > 0.0);
	CHECK(Luminance(halved.gas) == doctest::Approx(0.5 * Luminance(once.gas)).epsilon(1e-12));
}

namespace
{

// GasDisk with a bubble lit on its rim; the dust is real but faint, so the rim is the only gas.
Galaxy RimDisk(double rim)
{
	Galaxy galaxy = GasDisk(0.0);
	galaxy.dust = {.enabled = true, .extinction = 1e-6, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0,
		.bubbleCenter = {0.0, 0.0, 0.0}, .bubbleRadius = 300.0, .bubbleHeight = 800.0};
	galaxy.dust.bubbleRim = rim;
	return galaxy;
}

double RimLight(PreparedGalaxy const& prepared, glm::dvec3 p)
{
	GalaxySample sample;
	Evaluate(prepared, p, sample);
	return Luminance(sample.gas);
}

}

TEST_CASE("The rim glows on a shell just inside the bubble's edge, and nowhere else")
{
	// Mutation: centring the shell at q = 1 moves the peak out of [0.83, 0.87].
	PreparedGalaxy const prepared = PrepareGalaxy(RimDisk(4.0));
	RandomStream rng(23, Purpose::Test, 0, 0, 0);
	for (int i = 0; i < 20000; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		// q = 1 / |(d.x/300, d.y/300, d.z/800)| at unit scale; pick q away from the shell, then scale to it.
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		for (double const q : {0.2, 0.6, 0.64, 1.06, 1.5, 3.0})
		{
			CHECK(RimLight(prepared, d * (q / unit)) == 0.0);
		}
	}
	// Along radial lines through lit patches, the brightest point usually sits on the shell. The wisps
	// swing the light by an order of magnitude over 30 pc, about the shell's own width, so a line's peak
	// can ride a wisp off the shell; the averaged test below pins the profile itself.
	int lines = 0;
	int onShell = 0;
	for (int i = 0; i < 2000 && lines < 50; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		double bestQ = 0.0;
		double best = 0.0;
		for (double q = 0.66; q <= 1.04; q += 0.002)
		{
			double const light = RimLight(prepared, d * (q / unit));
			if (light > best)
			{
				best = light;
				bestQ = q;
			}
		}
		if (best > 0.0)
		{
			++lines;
			onShell += (bestQ >= 0.80 && bestQ <= 0.90) ? 1 : 0;
			CHECK(bestQ >= 0.66);
			CHECK(bestQ <= 1.04);
		}
	}
	CHECK(lines == 50);
	CHECK(onShell >= 20);
}

TEST_CASE("Averaged over the wall, the rim peaks at q = 0.85")
{
	// Mutation: centring the shell at q = 1 moves the averaged peak out of [0.83, 0.87].
	PreparedGalaxy const prepared = PrepareGalaxy(RimDisk(4.0));
	RandomStream rng(29, Purpose::Test, 0, 0, 0);
	std::array<double, 61> profile{};   // q from 0.70 to 1.00 in steps of 0.005
	for (int i = 0; i < 4000; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		for (size_t k = 0; k < profile.size(); ++k)
		{
			profile[k] += RimLight(prepared, d * ((0.70 + 0.005 * double(k)) / unit));
		}
	}
	size_t const peak = size_t(std::max_element(profile.begin(), profile.end()) - profile.begin());
	CHECK(0.70 + 0.005 * double(peak) >= 0.83);
	CHECK(0.70 + 0.005 * double(peak) <= 0.87);
}

TEST_CASE("The rim lights patches of the wall, not all of it")
{
	// Mutation: dropping the - 0.5 offset in the rim's mask lights most of the wall.
	PreparedGalaxy const prepared = PrepareGalaxy(RimDisk(4.0));
	RandomStream rng(31, Purpose::Test, 0, 0, 0);
	int dark = 0;
	int const count = 100000;
	for (int i = 0; i < count; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		dark += RimLight(prepared, d * (0.85 / unit)) == 0.0 ? 1 : 0;
	}
	CHECK(double(dark) / count > 0.60);
	CHECK(double(dark) / count < 0.80);
}

TEST_CASE("Without dust there is no rim")
{
	// Mutation: dropping the dust.enabled guard in RimGlow lights a bubble whose dust is gone.
	Galaxy galaxy = RimDisk(4.0);
	galaxy.dust.enabled = false;
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	RandomStream rng(37, Purpose::Test, 0, 0, 0);
	for (int i = 0; i < 2000; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		CHECK(RimLight(prepared, d * (0.85 / unit)) == 0.0);
	}
}

namespace
{

// GasDisk under uniform dust that reflects; no other gas.
Galaxy ReflectingDisk(double extinction)
{
	Galaxy galaxy = GasDisk(0.0);
	galaxy.dust = {.enabled = true, .extinction = extinction, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0};
	galaxy.dust.reflection = 2.0;
	galaxy.dust.reflectionColour = {0.9, 1.24, 2.0};   // twice the default: normalisation must bring it to Y = 1
	return galaxy;
}

}

TEST_CASE("Reflection tracks the dust and the young stars")
{
	// Mutation: using the dust's filament field without α₀ (out.extinction / dust.extinction) breaks the doubling.
	glm::dvec3 const p(300.0, -200.0, 100.0);
	GalaxySample once;
	GalaxySample twice;
	Evaluate(PrepareGalaxy(ReflectingDisk(0.001)), p, once);
	Evaluate(PrepareGalaxy(ReflectingDisk(0.002)), p, twice);
	REQUIRE(Luminance(once.gas) > 0.0);
	CHECK(Luminance(twice.gas) == doctest::Approx(2.0 * Luminance(once.gas)).epsilon(1e-12));
	// The dust's own exp(-R/1e12) falloff (R ~ 360 pc) is a 4e-10 relative gap to the young density, hence 1e-8.
	// At 1 mag/kpc the reflected light is `reflection` times the young stars' own.
	CHECK(Luminance(once.gas) == doctest::Approx(2.0 * once.parts[1].density * 100.0 * kFluxAt10pc).epsilon(1e-8));
	Galaxy dark = ReflectingDisk(0.001);
	dark.dust.enabled = false;
	GalaxySample none;
	Evaluate(PrepareGalaxy(dark), p, none);
	CHECK(none.gas == glm::dvec3(0.0));
}

TEST_CASE("Reflected light takes its normalised colour")
{
	// Mutation: skipping the reflection colour's normalisation doubles the luminance.
	PreparedGalaxy const prepared = PrepareGalaxy(ReflectingDisk(0.001));
	CHECK(Luminance(prepared.galaxy.dust.reflectionColour) == doctest::Approx(1.0).epsilon(1e-12));
	GalaxySample sample;
	Evaluate(prepared, {300.0, -200.0, 100.0}, sample);
	glm::dvec3 const expected = prepared.galaxy.dust.reflectionColour * Luminance(sample.gas);
	CHECK(sample.gas.x == doctest::Approx(expected.x).epsilon(1e-12));
	CHECK(sample.gas.y == doctest::Approx(expected.y).epsilon(1e-12));
	CHECK(sample.gas.z == doctest::Approx(expected.z).epsilon(1e-12));
}
