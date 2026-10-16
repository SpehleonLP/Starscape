#include "galaxy.h"
#include "gpu/device.h"
#include "gpu/pack.h"
#include "gpu/probe.h"
#include "march.h"
#include "presets.h"
#include "rng.h"
#include "scene.h"
#include <doctest.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <glm/geometric.hpp>

using namespace starscape;
using namespace starscape::gpu;

namespace
{

// stars/pc³ far above float's smallest normal (1.2e-38) and far below anything that lights a pixel.
constexpr double kSignificantDensity = 1e-30;

// Float precision, not a slip (owner ruling 2026-10-16): the GPU gets each point rounded to float, and the dust
// log-normal magnifies that rounding about 26-fold, so a perfect float GPU still scores 2.9e-4 on milky-way. The
// small spiral's arm phase (~10 rad) needs ~3e-7 that float cannot hold; it measured p99.9 1.5e-3, max 2.0e-3.
// The cook oracle (per-texel 1% at p99) is the gate on what this costs in the image.
constexpr double kFieldTolerance = 2e-3;

}

TEST_CASE("Packed structs have the sizes and member offsets GLSL's scalar layout gives them")
{
	// The numbers come from glslangValidator reflection of `layout(scalar) buffer { Galaxy g[]; }` and
	// `{ PopulationTable t[]; }` (spirv-dis: OpDecorate ArrayStride, OpMemberDecorate Offset), 2026-10-16.
	// Mutation: dropping a member from PackGas/PackDust moves the size off 624; swapping density and
	// scaleLength in PackDisk, or extinction and scaleLength in PackDust, moves a sentinel.
	// Sentinels sit at the reflected offsets so a same-size swap is caught, not only a dropped member.
	constexpr size_t kGalaxyStride = 624;
	constexpr size_t kPopulationStride = 780;
	Galaxy galaxy;
	galaxy.seed = 0x0123456789ABCDEFull;
	galaxy.disk.population = 0x11223344;
	galaxy.thickDisk.population = 0x55667788;
	galaxy.bulge.population = 0x0A0B0C0D;
	galaxy.halo.population = 0x0E0F1011;
	galaxy.dust.reflection = 2.5;
	galaxy.disk.density = 0.75;
	galaxy.disk.scaleLength = 1.25;
	galaxy.dust.extinction = 3.5;
	galaxy.dust.scaleLength = 4.5;
	ScalarWriter one;
	PackGalaxy(one, galaxy);
	CHECK(one.Bytes().size() == kGalaxyStride);
	ScalarWriter two;
	PackGalaxy(two, galaxy);
	PackGalaxy(two, galaxy);
	CHECK(two.Bytes().size() == 2 * kGalaxyStride);
	auto word = [&](size_t offset)
	{
		uint32_t value;
		std::memcpy(&value, one.Bytes().data() + offset, sizeof value);
		return value;
	};
	uint64_t seed;
	std::memcpy(&seed, one.Bytes().data() + 24, sizeof seed);
	CHECK(seed == galaxy.seed);
	CHECK(word(32 + 48) == 0x11223344u);    // disk.population
	CHECK(word(264 + 48) == 0x55667788u);   // thickDisk.population
	CHECK(word(496 + 28) == 0x0A0B0C0Du);   // bulge.population
	CHECK(word(528 + 16) == 0x0E0F1011u);   // halo.population
	float reflection;
	std::memcpy(&reflection, one.Bytes().data() + 548 + 60, sizeof reflection);
	CHECK(reflection == 2.5f);
	auto real = [&](size_t offset)
	{
		float value;
		std::memcpy(&value, one.Bytes().data() + offset, sizeof value);
		return value;
	};
	CHECK(real(32 + 4) == 0.75f);    // disk.density
	CHECK(real(32 + 8) == 1.25f);    // disk.scaleLength
	CHECK(real(548 + 4) == 3.5f);    // dust.extinction
	CHECK(real(548 + 8) == 4.5f);    // dust.scaleLength

	ScalarWriter table;
	PopulationTable population;
	population.chroma[1] = glm::dvec3(0.5, 0.25, 0.125);
	population.fraction[1] = 0.375;
	population.luminosity[1] = 6.5;
	population.lightBelow[1] = glm::dvec3(7.5, 8.5, 9.5);
	PackPopulation(table, population);
	CHECK(table.Bytes().size() == kPopulationStride);
	// Element 1 of each array at the reflected offsets (chroma@0, fraction@288, luminosity@384,
	// lightBelow@480): swapping two of PackPopulation's loops moves them.
	auto entry = [&](size_t offset)
	{
		float value;
		std::memcpy(&value, table.Bytes().data() + offset, sizeof value);
		return value;
	};
	CHECK(entry(12 + 4) == 0.25f);
	CHECK(entry(288 + 4) == 0.375f);
	CHECK(entry(384 + 4) == 6.5f);
	CHECK(entry(480 + 12 + 8) == 9.5f);
}

TEST_CASE("The GPU's lattice hash and Philox equal the CPU's, bit for bit")
{
	// Mutation: truncating HashMix to 32 bits in GLSL (casting h to uint before the shifts) reddens the hash
	// check; dropping one Philox round reddens the other.
	auto device = OpenDevice({.allowSoftware = true});
	REQUIRE(device.has_value());
	RandomStream rng(3, Purpose::Test, 0, 0, 0);
	std::vector<glm::i64vec3> cells(100000);
	for (glm::i64vec3& c : cells)
	{
		c = glm::i64vec3(int64_t(rng.NextU32()) - 2000000000, int64_t(rng.NextU32() % 40000) - 20000, int64_t(rng.NextU32() % 7) - 3);
	}
	auto hashes = ProbeLatticeHash(**device, cells, 0xD057D057D057D057ull ^ 1ull);
	REQUIRE(hashes.has_value());
	for (size_t i = 0; i < cells.size(); ++i)
	{
		REQUIRE((*hashes)[i] == LatticeHash(cells[i].x, cells[i].y, cells[i].z, 0xD057D057D057D057ull ^ 1ull));
	}
	std::vector<std::array<uint32_t, 4>> counters(100000);
	for (size_t i = 0; i < counters.size(); ++i)
	{
		counters[i] = {uint32_t(i), 0u, (1u << 24), uint32_t(i % 16)};
	}
	auto blocks = ProbePhilox(**device, counters, 0x0123456789ABCDEFull);
	REQUIRE(blocks.has_value());
	for (size_t i = 0; i < counters.size(); ++i)
	{
		REQUIRE((*blocks)[i] == Philox4x32_10(counters[i], {0x89ABCDEFu, 0x01234567u}));
	}
}

TEST_CASE("The GPU evaluates every galaxy of the stylized scene within float precision of the CPU")
{
	// Mutation: evaluating in the scene frame instead of the galaxy's local frame (adding the galaxy's
	// position back before EvaluateLocal) blows the distant spiral's error far past kFieldTolerance.
	auto device = OpenDevice({.allowSoftware = true});
	REQUIRE(device.has_value());
	Expected<Scene> const scene = LoadScene(STARSCAPE_SOURCE_DIR "/scenes/stylized-milky-way.toml");
	REQUIRE(scene.has_value());
	PreparedScene const prepared = PrepareScene(*scene);
	for (PreparedGalaxy const& galaxy : prepared.galaxies)
	{
		CAPTURE(galaxy.galaxy.preset);
		RandomStream rng(5, Purpose::Test, 0, 0, 0);
		std::vector<glm::dvec3> points;
		while (points.size() < 100000)
		{
			glm::dvec3 const p = (glm::dvec3(rng.NextDouble(), rng.NextDouble(), rng.NextDouble()) * 2.0 - 1.0) * galaxy.boundRadius;
			if (glm::dot(p, p) <= galaxy.boundRadius * galaxy.boundRadius)
			{
				// Flatten half the points into the disk, where the structure is.
				points.push_back(points.size() % 2 == 0 ? p : glm::dvec3(p.x, p.y, p.z * 0.01));
			}
		}
		auto gpu = ProbeEvaluate(**device, galaxy.galaxy, points);
		REQUIRE(gpu.has_value());
		std::vector<double> errors;
		for (size_t i = 0; i < points.size(); ++i)
		{
			GalaxySample cpu;
			EvaluateLocal(galaxy.galaxy, points[i], cpu);
			GalaxySample const& g = (*gpu)[i];
			// A density below float's range (the young stars' Sech2 tail reaches 1e-65) is zero on the GPU, so Push
			// drops the part there. Such a part is no star in any march; the parts that matter must agree in
			// number, order and population.
			auto significant = [](GalaxySample const& sample)
			{
				std::vector<Contribution> kept;
				for (int k = 0; k < sample.count; ++k)
				{
					if (sample.parts[size_t(k)].density > kSignificantDensity)
					{
						kept.push_back(sample.parts[size_t(k)]);
					}
				}
				return kept;
			};
			std::vector<Contribution> const cpuParts = significant(cpu);
			std::vector<Contribution> const gpuParts = significant(g);
			REQUIRE(gpuParts.size() == cpuParts.size());
			auto relative = [](double a, double b) { return std::abs(a - b) / std::max(std::abs(b), 1e-30); };
			double worst = relative(g.extinction, cpu.extinction);
			worst = std::max(worst, relative(Luminance(g.gas), Luminance(cpu.gas)));
			for (size_t k = 0; k < cpuParts.size(); ++k)
			{
				REQUIRE(gpuParts[k].population == cpuParts[k].population);
				worst = std::max(worst, relative(gpuParts[k].density, cpuParts[k].density));
			}
			errors.push_back(worst);
		}
		std::sort(errors.begin(), errors.end());
		double const p999 = errors[size_t(0.999 * double(errors.size()))];
		MESSAGE("p99.9 relative error " << p999 << ", max " << errors.back() << " for " << galaxy.galaxy.preset);
		CHECK(p999 < kFieldTolerance);
	}
}
