#include "cook.h"
#include "gpu/device.h"
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
	uint64_t cpuStars;
	uint64_t gpuStars;
};

Comparison CompareCooks(Scene const& scene)
{
	CookOptions cpu;
	cpu.threads = 8;
	cpu.device = DeviceChoice::Cpu;
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
	for (size_t i = 0; i < a->budget.size(); ++i)
	{
		for (size_t bin = 0; bin < size_t(kBins); ++bin)
		{
			double const reference = a->budget[i][bin];
			if (reference > 0.0)
			{
				budget.push_back(std::abs(double(b->budget[i][bin]) - reference) / reference);
			}
		}
	}
	std::sort(budget.begin(), budget.end());
	result.budgetP99 = budget.empty() ? 0.0 : budget[size_t(0.99 * double(budget.size()))];
	result.cpuStars = a->stats.starCount;
	result.gpuStars = b->stats.starCount;
	return result;
}

}

TEST_CASE("The GPU cooks the stylized scene within tolerance of the CPU")
{
	// Mutation: marching from each step's start (sNow) instead of its midpoint in march.comp shifts every
	// sample by half a step and fails the flux and radiance checks.
	Expected<Scene> scene = LoadScene(STARSCAPE_SOURCE_DIR "/scenes/stylized-milky-way.toml");
	REQUIRE(scene.has_value());
	scene->output.faceSize = 64;
	scene->output.quality = Quality::Final;
	Comparison const c = CompareCooks(*scene);
	MESSAGE("flux ", c.fluxError, " radiance p99 ", c.radianceP99, " budget p99 ", c.budgetP99, " stars ", c.cpuStars, " vs ", c.gpuStars);
	CHECK(c.fluxError < 1e-3);
	CHECK(c.radianceP99 < 1e-2);
	CHECK(c.budgetP99 < 1e-2);
	double const sigma = std::sqrt(double(c.cpuStars));
	CHECK(std::abs(double(c.gpuStars) - double(c.cpuStars)) < 3.0 * sigma + 1.0);
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
	Comparison const c = CompareCooks(scene);
	MESSAGE("flux ", c.fluxError, " radiance p99 ", c.radianceP99, " budget p99 ", c.budgetP99, " stars ", c.cpuStars, " vs ", c.gpuStars);
	CHECK(c.fluxError < 1e-3);
	CHECK(c.radianceP99 < 1e-2);
	CHECK(c.budgetP99 < 1e-2);
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
