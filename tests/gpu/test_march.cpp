// The GPU march held to answers known by construction, plus its tiling, its refusals and its determinism.
#include "cook.h"
#include "gpu/gpu_march.h"
#include "scene.h"
#include "test_scenes.h"
#include "writer.h"
#include <doctest.h>
#include <kre_starscape.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <string>
#include <thread>

using namespace starscape;
using namespace starscape::test;

namespace
{

CookOptions Options()
{
	CookOptions options;
	options.threads = std::max(1u, std::thread::hardware_concurrency());
	return options;
}

// A uniform single-bin ball around the observer: every ray crosses the same medium, so each texel's mean radiance
// and budget have closed forms whatever the jitter.
Scene BallScene(uint32_t faceSize)
{
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 12);
	scene.galaxies.push_back(UniformBall(0.02, 400.0, population, {0.0, 0.0, 0.0}));
	scene.output.faceSize = faceSize;
	scene.output.quality = Quality::Final;
	scene.output.limitingMagnitude = 4.0;
	return scene;
}

// Holds every texel of `cooked` to the ball's radiance n·100·F₁₀·L·C·(R − r) and budget n·Ω·(r³ − sMin³)/3, where r
// is the distance out to which the bin is resolved. Float sums on the GPU, so 10⁻⁴.
void CheckBall(Scene const& scene, Cooked const& cooked)
{
	int const bin = 12;
	PreparedScene const prepared = PrepareScene(scene);
	PopulationTable const& table = prepared.populations.back();
	double const density = 0.02;
	double const radius = 400.0;
	double const sMin = scene.observer.sMin;
	double const luminosity = table.luminosity[size_t(bin)];
	double const resolved = std::sqrt(luminosity * 100.0 * kFluxAt10pc / cooked.fluxThreshold);
	REQUIRE(resolved > sMin);
	REQUIRE(resolved < radius);
	glm::dvec3 const radiance = density * 100.0 * kFluxAt10pc * luminosity * (radius - resolved) * table.chroma[size_t(bin)];
	double const budgetPerSteradian = density * (resolved * resolved * resolved - sMin * sMin * sMin) / 3.0;
	uint32_t const size = cooked.faceSize;
	REQUIRE(cooked.budget.size() == 6 * size * size);
	int wrong = 0;
	for (uint32_t texel = 0; texel < 6 * size * size; ++texel)
	{
		double const omega = kre_starscape::TexelSolidAngle(size, texel % size, (texel / size) % size);
		double const budget = budgetPerSteradian * omega;
		bool const right = std::abs(cooked.diffuse[texel].x - radiance.x) <= 1e-4 * radiance.x
			&& std::abs(cooked.diffuse[texel].y - radiance.y) <= 1e-4 * radiance.y
			&& std::abs(cooked.diffuse[texel].z - radiance.z) <= 1e-4 * radiance.z
			&& std::abs(double(cooked.budget[texel][size_t(bin)]) - budget) <= 1e-4 * budget;
		if (!right && wrong++ < 3)
		{
			MESSAGE("pass ", cooked.passes, " texel ", texel, ": radiance y ", cooked.diffuse[texel].y, " vs ", radiance.y, ", budget ",
				cooked.budget[texel][size_t(bin)], " vs ", budget);
		}
	}
	CHECK(wrong == 0);
}

// Σ radiance × texel solid angle per channel: the whole sky's flux. A bright limiting magnitude keeps every star
// diffuse, so this is all the light.
glm::dvec3 SkyFlux(Scene const& base)
{
	Scene scene = base;
	scene.output.limitingMagnitude = -30.0;
	Expected<Cooked> const cooked = Cook(scene, Options());
	REQUIRE(cooked.has_value());
	REQUIRE(cooked->stars.empty());
	uint32_t const size = cooked->faceSize;
	glm::dvec3 total(0.0);
	for (uint32_t texel = 0; texel < 6 * size * size; ++texel)
	{
		total += cooked->diffuse[texel] * kre_starscape::TexelSolidAngle(size, texel % size, (texel / size) % size);
	}
	return total;
}

// A uniform dust ball around the observer puts the same column in front of every direction, so whatever lies beyond
// it is dimmed by exactly exp(-k_c τ) in each channel, however coarsely the sky is sampled.
void CheckDimmedByDustBall(glm::dvec3 clear, glm::dvec3 dimmed, double alpha, double radius, double sMin)
{
	double const tau = alpha * (radius - sMin) / 1.0857;
	CHECK(dimmed.x / clear.x == doctest::Approx(std::exp(-0.75 * tau)).epsilon(0.005).scale(0.0));
	CHECK(dimmed.y / clear.y == doctest::Approx(std::exp(-1.00 * tau)).epsilon(0.005).scale(0.0));
	CHECK(dimmed.z / clear.z == doctest::Approx(std::exp(-1.32 * tau)).epsilon(0.005).scale(0.0));
}

// A uniform disk 7 to 13 kpc away, wholly beyond the dust ball, with dust disabled so the gas's filament factor is
// 1; or, with `dust`, faint dust and a bubble at its centre, so the rim and the reflection have something to light.
// The bubble is 3.4° across, wide enough for a face-64 sky's rays to cross it.
Scene NearDisk(bool dust)
{
	Scene scene = BaseScene();
	scene.output.faceSize = 64;
	scene.output.quality = Quality::Final;
	int const population = AddSingleBinPopulation(scene, 14);
	Galaxy source;
	source.position = {1.0e4, 0.0, 0.0};
	source.disk = {.enabled = true, .density = 0.05, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 3000.0, .arms = {},
		.population = population, .armPopulation = population, .youngFraction = 0.5};
	if (dust)
	{
		source.dust = {.enabled = true, .extinction = 1e-6, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 3000.0, .laneOffset = 0.0,
			.contrast = 1.0, .filament = 0.0, .bubbleCenter = {0.0, 0.0, 0.0}, .bubbleRadius = 300.0, .bubbleHeight = 800.0};
	}
	scene.galaxies.push_back(source);
	return scene;
}

}

TEST_CASE("A uniform ball around the observer cooks to its closed form after every pass")
{
	// Mutation: scaling the sums by S² instead of S²/n leaves every partial snapshot n/S² too dark; dropping the
	// scale on the budget alone leaves the radiance right and the budget a sum, which the budget check reddens.
	Scene const scene = BallScene(8);
	CookOptions options = Options();
	options.keepBudget = true;
	int calls = 0;
	options.onPass = [&](Cooked const& cooked, int pass, int passes) -> Expected<void>
	{
		CHECK(pass == calls + 1);
		CHECK(passes == 16);
		CHECK(cooked.passes == pass);
		CheckBall(scene, cooked);
		++calls;
		return {};
	};
	Expected<Cooked> const cooked = Cook(scene, options);
	REQUIRE(cooked.has_value());
	CHECK(calls == 16);
	CHECK(cooked->device.starts_with("GPU "));
}

TEST_CASE("A face smaller than a workgroup cooks to the closed form")
{
	// Mutation: dropping the edge test in march.comp lets invocations past the face edge add into other texels' sums.
	Scene const scene = BallScene(4);
	CookOptions options = Options();
	options.keepBudget = true;
	Expected<Cooked> const cooked = Cook(scene, options);
	REQUIRE(cooked.has_value());
	CheckBall(scene, *cooked);
}

TEST_CASE("A halo peaked on the observer integrates to its closed form along every ray")
{
	// Mutation: sampling each step at its start (mid = t in march.comp) instead of its midpoint biases every step by
	// half its width times the density's slope, about 1% here, where the midpoint rule's error is a few 1e-4.
	Scene scene = BaseScene();
	int const bin = 12;
	int const population = AddSingleBinPopulation(scene, bin);
	double const density = 0.02;
	double const scale = 10.0;
	double const cutoff = 200.0;
	Galaxy halo;
	halo.halo = {.enabled = true, .density = density, .scale = scale, .cutoff = cutoff, .population = population};
	scene.galaxies.push_back(halo);
	scene.output.faceSize = 4;
	scene.output.quality = Quality::Final;
	scene.output.limitingMagnitude = -30.0;   // every star diffuse
	Expected<Cooked> const cooked = Cook(scene, Options());
	REQUIRE(cooked.has_value());

	// ∫ n₀ (1 + s/a)^-3.5 ds from sMin to the cutoff, in radiance units.
	double const sMin = scene.observer.sMin;
	double const column = density * scale / 2.5 * (std::pow(1.0 + sMin / scale, -2.5) - std::pow(1.0 + cutoff / scale, -2.5));
	PreparedScene const prepared = PrepareScene(scene);
	PopulationTable const& table = prepared.populations.back();
	double const expected = Luminance(100.0 * kFluxAt10pc * table.luminosity[size_t(bin)] * column * table.chroma[size_t(bin)]);
	double worst = 0.0;
	for (glm::dvec3 const& radiance : cooked->diffuse)
	{
		worst = std::max(worst, std::abs(Luminance(radiance) / expected - 1.0));
	}
	MESSAGE("worst relative error over texels ", worst);
	CHECK(worst < 1e-3);
}

TEST_CASE("A distant dust-free galaxy delivers L (10/d)² C in total")
{
	// Mutation: weighting each ray by 1/S instead of 1/S² (march.cpp PassConstants) quadruples the sky's light.
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
	// The galaxy is 2.3° across: texels of a third of a degree, sixteen rays each, sample it to 0.1%. Face 128 is
	// off by 1.6%.
	scene.output.faceSize = 256;
	scene.output.quality = Quality::Final;

	// ∫₀^R 4π r² (1 + r/a)^-3.5 dr = 4π a³ [−2x^-½ + (4/3)x^-3/2 − 0.4x^-5/2]₁^{1+R/a}.
	auto antiderivative = [](double x) { return -2.0 * std::pow(x, -0.5) + 4.0 / 3.0 * std::pow(x, -1.5) - 0.4 * std::pow(x, -2.5); };
	double const stars = density * 4.0 * kPi * scale * scale * scale * (antiderivative(1.0 + cutoff / scale) - antiderivative(1.0));
	double const expected = stars * BinLuminosity(bin) * 100.0 * kFluxAt10pc / (distance * distance);
	glm::dvec3 const flux = SkyFlux(scene);
	MESSAGE("flux ", Luminance(flux), " expected ", expected, " relative error ", Luminance(flux) / expected - 1.0);
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
	scene.output.faceSize = 64;
	scene.output.quality = Quality::Final;
	glm::dvec3 const clear = SkyFlux(scene);
	REQUIRE(Luminance(clear) > 0.0);
	scene.galaxies.push_back(UniformDust(0.002, 1000.0));
	CheckDimmedByDustBall(clear, SkyFlux(scene), 0.002, 1000.0, scene.observer.sMin);
}

TEST_CASE("Gas glow is dimmed by the dust in front of it")
{
	// Mutation: adding the gas without the step's transmittance (march.comp) reddens this test.
	// The gas's own light: the same rays with and without glow, so the starlight cancels.
	auto gasFlux = [](Scene const& base)
	{
		Scene lit = base;
		lit.galaxies[0].disk.gas.glow = 1.0;
		return SkyFlux(lit) - SkyFlux(base);
	};
	Scene scene = NearDisk(false);
	glm::dvec3 const clear = gasFlux(scene);
	REQUIRE(Luminance(clear) > 0.0);
	scene.galaxies.push_back(UniformDust(0.002, 1000.0));
	CheckDimmedByDustBall(clear, gasFlux(scene), 0.002, 1000.0, scene.observer.sMin);
}

TEST_CASE("Large-scale gas, the rim and reflection are dimmed by the dust in front of them")
{
	// Mutation: adding the gas without the step's transmittance (march.comp) leaves every term undimmed, and every
	// section fails. The terms reach the march as one gas sum, so no mutation can single one out there.
	struct Light
	{
		char const* name;
		void (*apply)(Galaxy&);
	};
	Light const lights[] = {
		{"large", [](Galaxy& g) { g.disk.gas.largeGlow = 1.0; }},
		{"rim", [](Galaxy& g) { g.dust.bubbleRim = 4.0; }},
		{"reflection", [](Galaxy& g) { g.dust.reflection = 2.0; }},
	};
	// The unlit sky is shared by every term: four cooks of it would be two cooks too many.
	Scene clearScene = NearDisk(true);
	Scene dustyScene = clearScene;
	dustyScene.galaxies.push_back(UniformDust(0.002, 1000.0));
	glm::dvec3 const clearDark = SkyFlux(clearScene);
	glm::dvec3 const dustyDark = SkyFlux(dustyScene);
	for (Light const& light : lights)
	{
		std::string const name = light.name;
		CAPTURE(name);
		auto lit = [&](Scene const& base)
		{
			Scene scene = base;
			light.apply(scene.galaxies[0]);
			return SkyFlux(scene);
		};
		glm::dvec3 const clear = lit(clearScene) - clearDark;
		REQUIRE(Luminance(clear) > 0.0);
		CheckDimmedByDustBall(clear, lit(dustyScene) - dustyDark, 0.002, 1000.0, clearScene.observer.sMin);
	}
}

TEST_CASE("Gas adds diffuse light and leaves the stars alone")
{
	// Mutation: adding the gas into the resolved light instead of the diffuse radiance (march.comp) leaves the sky as
	// dark as without it. Every term is on at once, around a bubble both scenes share, so only the light terms differ.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	scene.output.limitingMagnitude = 4.0;
	Galaxy& galaxy = scene.galaxies[0];
	galaxy.dust.bubbleCenter = {8200.0, 0.0, 20.0};
	galaxy.dust.bubbleRadius = 300.0;
	galaxy.dust.bubbleHeight = 800.0;
	CookOptions options = Options();
	options.keepBudget = true;
	Expected<Cooked> const plain = Cook(scene, options);
	galaxy.disk.gas.glow = 1.5;
	galaxy.disk.gas.palette = {glm::dvec3(1.0, 0.18, 0.32), glm::dvec3(0.25, 1.0, 0.85)};
	galaxy.disk.gas.paletteCount = 2;
	galaxy.disk.gas.largeGlow = 1.5;
	galaxy.dust.bubbleRim = 4.0;
	galaxy.dust.reflection = 2.0;
	Expected<Cooked> const glowing = Cook(scene, options);
	REQUIRE(plain.has_value());
	REQUIRE(glowing.has_value());
	CHECK(glowing->stats.diffuseFlux > 1.01 * plain->stats.diffuseFlux);
	CHECK(glowing->budget == plain->budget);
	REQUIRE(plain->stars.size() > 0);
	CHECK(glowing->stars.size() == plain->stars.size());
}

TEST_CASE("More galaxies than the GPU takes are refused, naming the cap")
{
	// Mutation: dropping the galaxy cap in OpenGpuMarcher lets a 17-galaxy scene reach a shader whose arrays hold 16.
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 12);
	for (int g = 0; g < 17; ++g)
	{
		scene.galaxies.push_back(UniformBall(0.001, 100.0, population, {1000.0 * g, 5000.0, 0.0}));
	}
	Expected<Cooked> const refused = Cook(scene, Options());
	REQUIRE_FALSE(refused.has_value());
	CHECK(refused.error().kind == Error::Kind::Input);
	CHECK(refused.error().message.find("at most 16 galaxies") != std::string::npos);
}

TEST_CASE("A cook is the same bytes twice in one process and in a separate CLI run")
{
	// Mutation: mixing the clock into the march's jitter seed (march.cpp PassConstants) makes every run differ.
	// The CLI's file is removed first and written here, so the comparison never sees bytes from an older binary.
	char const* const scenePath = STARSCAPE_SOURCE_DIR "/tests/data/determinism.toml";
	std::filesystem::path const cliOut = STARSCAPE_CLI_OUT "/determinism-gate.kre_starscape";
	std::filesystem::remove(cliOut);
	std::string const command = std::format("\"{}\" cook \"{}\" -o \"{}\" > /dev/null", STARSCAPE_CLI, scenePath, cliOut.string());
	REQUIRE(std::system(command.c_str()) == 0);
	Expected<std::vector<uint8_t>> const cli = ReadFileBytes(cliOut);
	REQUIRE(cli.has_value());

	Expected<Scene> const scene = LoadScene(scenePath);
	REQUIRE(scene.has_value());
	Expected<Cooked> const first = Cook(*scene, Options());
	Expected<Cooked> const second = Cook(*scene, Options());
	REQUIRE(first.has_value());
	REQUIRE(second.has_value());
	REQUIRE(first->stars.size() > 0);
	std::vector<uint8_t> const a = BuildFile(*first, scene->hash);
	std::vector<uint8_t> const b = BuildFile(*second, scene->hash);
	CHECK(a == b);
	CHECK(*cli == a);
}

namespace
{

// The progress reports of one pass over a face-64 sky, one per dispatched tile, with tiles sized toward `tileSeconds`.
std::vector<double> TileReports(double tileSeconds)
{
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 12);
	scene.galaxies.push_back(UniformBall(0.02, 400.0, population, {0.0, 0.0, 0.0}));
	constexpr uint32_t size = 64;
	PreparedScene const prepared = PrepareScene(scene);
	std::vector<double> reports;
	Expected<std::unique_ptr<PassMarcher>> marcher =
		OpenGpuMarcher(prepared, ParamsFor(Quality::Quick), size, [&](double done) { reports.push_back(done); }, tileSeconds);
	REQUIRE(marcher.has_value());
	uint32_t const texels = 6 * size * size;
	PassSums sums{std::vector<glm::vec3>(texels, glm::vec3(0.0f)), std::vector<std::array<float, kBins>>(texels)};
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
