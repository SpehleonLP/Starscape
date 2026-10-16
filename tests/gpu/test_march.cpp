#include "cook.h"
#include "gpu/device.h"
#include "gpu/gpu_march.h"
#include "scene.h"
#include "test_scenes.h"
#include <doctest.h>
#include <kre_starscape.h>
#include <algorithm>
#include <cmath>

using namespace starscape;

namespace
{

struct Comparison
{
	double fluxError;      // relative, whole sky
	double radianceP99;    // relative per-texel radiance error, 99th percentile over texels with light
	double budgetP99;      // relative per-texel-bin budget error, 99th percentile over non-zero bins
	double budgetError;    // relative, whole sky: also counts bins the CPU leaves at zero, which budgetP99 skips
	uint64_t cpuStars;
	uint64_t gpuStars;
};

Comparison CompareCooks(Scene const& scene)
{
	CookOptions cpu;
	cpu.threads = 8;
	cpu.device = DeviceChoice::Cpu;
	cpu.keepBudget = true;
	CookOptions gpu = cpu;
	gpu.device = DeviceChoice::Gpu;
	Expected<Cooked> const a = Cook(scene, cpu);
	Expected<Cooked> const b = Cook(scene, gpu);
	REQUIRE(a.has_value());
	REQUIRE(b.has_value());
	REQUIRE(b->device.starts_with("GPU"));
	Comparison result{};
	result.fluxError = std::abs(b->stats.diffuseFlux - a->stats.diffuseFlux) / a->stats.diffuseFlux;
	std::vector<double> radiance;
	for (size_t i = 0; i < a->diffuse.size(); ++i)
	{
		double const reference = Luminance(a->diffuse[i]);
		if (reference > 0.0)
		{
			radiance.push_back(std::abs(Luminance(b->diffuse[i]) - reference) / reference);
		}
	}
	std::sort(radiance.begin(), radiance.end());
	result.radianceP99 = radiance.empty() ? 0.0 : radiance[size_t(0.99 * double(radiance.size()))];
	std::vector<double> budget;
	double cpuBudget = 0.0;
	double gpuBudget = 0.0;
	for (size_t i = 0; i < a->budget.size(); ++i)
	{
		for (size_t bin = 0; bin < size_t(kBins); ++bin)
		{
			double const reference = a->budget[i][bin];
			cpuBudget += reference;
			gpuBudget += double(b->budget[i][bin]);
			if (reference > 0.0)
			{
				budget.push_back(std::abs(double(b->budget[i][bin]) - reference) / reference);
			}
		}
	}
	std::sort(budget.begin(), budget.end());
	result.budgetP99 = budget.empty() ? 0.0 : budget[size_t(0.99 * double(budget.size()))];
	REQUIRE(cpuBudget > 0.0);
	result.budgetError = std::abs(gpuBudget - cpuBudget) / cpuBudget;
	result.cpuStars = a->stats.starCount;
	result.gpuStars = b->stats.starCount;
	MESSAGE("flux ", result.fluxError, " radiance p99 ", result.radianceP99, " budget p99 ", result.budgetP99, " budget total ",
		result.budgetError, " stars ", result.cpuStars, " vs ", result.gpuStars);
	return result;
}

// The spec's cook-oracle tolerances. The whole-sky budget is held to the flux's 0.1%: both are sums over the
// sky, where per-texel noise averages out.
void CheckTolerances(Comparison const& c)
{
	CHECK(c.fluxError < 1e-3);
	CHECK(c.radianceP99 < 1e-2);
	CHECK(c.budgetP99 < 1e-2);
	CHECK(c.budgetError < 1e-3);
	double const sigma = std::sqrt(double(c.cpuStars));
	CHECK(std::abs(double(c.gpuStars) - double(c.cpuStars)) < 3.0 * sigma + 1.0);
}

}

TEST_CASE("The GPU cooks the stylized scene within tolerance of the CPU")
{
	// Mutation: marching from each step's start (mid = t) instead of its midpoint in march.comp shifts every
	// sample by half a step and fails the flux and radiance checks.
	Expected<Scene> scene = LoadScene(STARSCAPE_SOURCE_DIR "/scenes/stylized-milky-way.toml");
	REQUIRE(scene.has_value());
	scene->output.faceSize = 64;
	scene->output.quality = Quality::Final;
	CheckTolerances(CompareCooks(*scene));
}

TEST_CASE("A face smaller than a workgroup cooks and matches")
{
	// Mutation: dropping the edge test in march.comp lets invocations past the face edge add into other texels' sums. A
	// uniform ball gives every texel nearly the same radiance, so the budget check is the one that reddens.
	Scene scene = test::BaseScene();
	int const population = test::AddSingleBinPopulation(scene, 12);
	scene.galaxies.push_back(test::UniformBall(0.02, 400.0, population, {0.0, 0.0, 0.0}));
	scene.output.faceSize = 4;
	scene.output.quality = Quality::Final;
	scene.output.limitingMagnitude = 4.0;
	CheckTolerances(CompareCooks(scene));
}

TEST_CASE("More galaxies than the GPU takes: auto uses the CPU, gpu refuses")
{
	// Mutation: dropping the galaxy cap in OpenGpuMarcher lets a 17-galaxy scene reach a shader whose
	// arrays hold 16.
	Scene scene = test::BaseScene();
	int const population = test::AddSingleBinPopulation(scene, 12);
	for (int g = 0; g < 17; ++g)
	{
		scene.galaxies.push_back(test::UniformBall(0.001, 100.0, population, {1000.0 * g, 5000.0, 0.0}));
	}
	CookOptions options;
	options.device = DeviceChoice::Gpu;
	Expected<Cooked> const refused = Cook(scene, options);
	REQUIRE_FALSE(refused.has_value());
	CHECK(refused.error().message.find("at most 16 galaxies") != std::string::npos);
	options.device = DeviceChoice::Auto;
	std::string logged;
	options.log = [&](std::string_view line) { logged += line; };
	Expected<Cooked> const cpu = Cook(scene, options);
	REQUIRE(cpu.has_value());
	CHECK(cpu->device.starts_with("CPU"));
	CHECK(logged.find("at most 16 galaxies") != std::string::npos);
}

namespace
{

// The progress reports of one pass over a face-64 sky, one per dispatched tile, with tiles sized toward `tileSeconds`.
std::vector<double> TileReports(double tileSeconds)
{
	Scene scene = test::BaseScene();
	int const population = test::AddSingleBinPopulation(scene, 12);
	scene.galaxies.push_back(test::UniformBall(0.02, 400.0, population, {0.0, 0.0, 0.0}));
	constexpr uint32_t size = 64;
	PreparedScene const prepared = PrepareScene(scene);
	std::vector<double> reports;
	Expected<std::unique_ptr<PassMarcher>> marcher =
		OpenGpuMarcher(prepared, ParamsFor(Quality::Quick), size, [&](double done) { reports.push_back(done); }, tileSeconds);
	REQUIRE(marcher.has_value());
	uint32_t const texels = 6 * size * size;
	PassSums sums{std::vector<glm::dvec3>(texels, glm::dvec3(0.0)), std::vector<std::array<float, kBins>>(texels)};
	PassSpec const spec{size, 1, {0, 0}, 0.025, scene.output.seed};
	REQUIRE((*marcher)->MarchPass(spec, sums).has_value());
	return reports;
}

}

TEST_CASE("Tiles double after fast dispatches and drop to a row after slow ones")
{
	// Mutation: never updating the marcher's rows keeps every tile at 8 rows, 48 a pass, and reddens both counts;
	// sizing that ignores the measured seconds doubles every time and reddens the slow count.
	// A target no dispatch can miss: rows go 8, 16, 32, then the face's last 8, and each later face is one tile.
	std::vector<double> const fast = TileReports(1.0e6);
	CHECK(fast.size() == 4 + 5);
	REQUIRE_FALSE(fast.empty());
	CHECK(fast.back() == 1.0);
	// A target every dispatch misses: after the first 8 rows, one row a tile.
	std::vector<double> const slow = TileReports(1.0e-9);
	CHECK(slow.size() == 1 + 56 + 5 * 64);
	REQUIRE_FALSE(slow.empty());
	CHECK(slow.back() == 1.0);
}
