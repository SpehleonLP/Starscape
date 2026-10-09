#include "cook.h"
#include "preview.h"
#include "test_scenes.h"
#include "writer.h"
#include <kre_starscape.h>
#include <doctest.h>
#include <algorithm>
#include <cmath>
#include <thread>

using namespace starscape;
using namespace starscape::test;

namespace
{

unsigned Threads()
{
	return std::max(1u, std::thread::hardware_concurrency());
}

double TotalFlux(Cooked const& cooked)
{
	return cooked.stats.pointFlux + cooked.stats.diffuseFlux;
}

double AngleDegrees(glm::dvec3 a, glm::dvec3 b)
{
	return std::acos(std::clamp(glm::dot(glm::normalize(a), glm::normalize(b)), -1.0, 1.0)) * 180.0 / kPi;
}

// A compact, centrally peaked lump 300 pc away along `dir` (sky frame = scene frame for the default observer),
// no wider than 2 degrees. Final quality: sixteen rays per texel, so the lump cannot fall between them.
Scene LumpScene(glm::dvec3 dir, double limitingMagnitude)
{
	Scene scene = BaseScene();
	scene.output.faceSize = 32;
	scene.output.quality = Quality::Final;
	scene.output.limitingMagnitude = limitingMagnitude;
	int const population = AddSingleBinPopulation(scene, 14);
	Galaxy galaxy;
	galaxy.position = 300.0 * dir;
	galaxy.halo = {.enabled = true, .density = 5.0, .scale = 1.5, .cutoff = 10.0, .population = population};
	scene.galaxies.push_back(galaxy);
	return scene;
}

// Brightest pixel of a preview, as the direction PixelDirection assigns it.
glm::dvec3 BrightestPixelDirection(Image const& image)
{
	size_t best = 0;
	float bestValue = -1.0f;
	for (size_t i = 0; i + 2 < image.rgb.size(); i += 3)
	{
		float const value = image.rgb[i] + image.rgb[i + 1] + image.rgb[i + 2];
		if (value > bestValue)
		{
			bestValue = value;
			best = i / 3;
		}
	}
	return PixelDirection(image.width, double(best % image.width) + 0.5, double(best / image.width) + 0.5);
}

double SkyLatitude(glm::dvec3 dir)
{
	return std::asin(std::clamp(dir.z, -1.0, 1.0)) * 180.0 / kPi;
}

}

TEST_CASE("Poisson draws have the right mean and variance")
{
	for (double const lambda : {0.3, 5.0, 100.0})
	{
		RandomStream rng(9, Purpose::Test, uint32_t(lambda * 10), 0, 0);
		constexpr int n = 200000;
		double sum = 0.0;
		double sumSquares = 0.0;
		for (int i = 0; i < n; ++i)
		{
			double const k = SamplePoisson(lambda, rng);
			sum += k;
			sumSquares += k * k;
		}
		double const mean = sum / n;
		CHECK(mean == doctest::Approx(lambda).epsilon(0.02).scale(0.0));
		CHECK(sumSquares / n - mean * mean == doctest::Approx(lambda).epsilon(0.03).scale(0.0));
	}
	RandomStream rng(9, Purpose::Test, 0, 0, 1);
	CHECK(SamplePoisson(0.0, rng) == 0);
}

TEST_CASE("Texel directions stay in their texel and are uniform in solid angle")
{
	constexpr uint32_t size = 2;
	uint32_t const texel = 3 * size * size + 0;   // face NegY, a corner texel: the strongest Jacobian
	RandomStream rng(5, Purpose::Test, 0, 0, 0);
	// Split the texel into quarters and compare hit fractions with the quarters' solid angles.
	constexpr uint32_t fine = 4;
	double hits[2][2] = {};
	constexpr int n = 100000;
	for (int i = 0; i < n; ++i)
	{
		glm::dvec3 const dir = SampleTexelDirection(size, texel, rng);
		kre_starscape::FaceCoord const coord = kre_starscape::DirToFace(dir.x, dir.y, dir.z);
		REQUIRE(coord.face == kre_starscape::NegY);
		REQUIRE(coord.s < 0.5);
		REQUIRE(coord.t < 0.5);
		hits[int(coord.t * fine) % 2][int(coord.s * fine) % 2] += 1.0;
	}
	double const whole = kre_starscape::TexelSolidAngle(size, 0, 0);
	for (uint32_t y = 0; y < 2; ++y)
	{
		for (uint32_t x = 0; x < 2; ++x)
		{
			double const share = kre_starscape::TexelSolidAngle(fine, x, y) / whole;
			CHECK(hits[y][x] / n == doctest::Approx(share).epsilon(0.02).scale(0.0));
		}
	}
}

TEST_CASE("Calibration is deterministic across thread counts")
{
	Scene const scene = MilkyWayScene({8200.0, 0.0, 20.0});
	PreparedScene const prepared = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	FluxHistogram const one = CalibrationHistogram(prepared, params, 16, 1);
	FluxHistogram const many = CalibrationHistogram(prepared, params, 16, 8);
	CHECK(one.count == many.count);
}

TEST_CASE("Energy is conserved whatever the star budget")
{
	// Mutation: double-counting resolved light reddens this test.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	scene.output.faceSize = 32;
	scene.output.limitingMagnitude.reset();
	std::vector<Cooked> cooks;
	for (uint32_t const target : {500u, 4000u, 20000u})
	{
		scene.output.targetStars = target;
		cooks.push_back(Cook(scene, {.threads = Threads()}).value());
	}
	CHECK(double(cooks[1].stars.size()) == doctest::Approx(4000.0).epsilon(0.25));
	for (Cooked const& cooked : cooks)
	{
		double pointVariance = 0.0;
		for (CookedStar const& star : cooked.stars)
		{
			pointVariance += Luminance(star.flux) * Luminance(star.flux);
		}
		// Compound-Poisson noise of the points, plus RGB9E5-scale slack for the diffuse sum.
		double const tolerance = 3.0 * std::sqrt(pointVariance) + 0.005 * TotalFlux(cooks[0]);
		CHECK(std::abs(TotalFlux(cooked) - TotalFlux(cooks[0])) <= tolerance);
	}
	CHECK(cooks[2].stats.pointFlux > cooks[0].stats.pointFlux);
	CHECK(cooks[2].stats.diffuseFlux < cooks[0].stats.diffuseFlux);
}

TEST_CASE("From inside the disk the band is diffuse and the points are nearly uniform")
{
	// Mutation: replacing the sech² disk with a flat one reddens this test.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});   // spin axis = galactic pole, so sky latitude = galactic latitude
	scene.output.faceSize = 64;
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 20000;
	Cooked const cooked = Cook(scene, {.threads = Threads()}).value();

	double planeRadiance = 0.0;
	double planeAngle = 0.0;
	double poleRadiance = 0.0;
	double poleAngle = 0.0;
	uint32_t const size = cooked.faceSize;
	for (uint32_t texel = 0; texel < 6 * size * size; ++texel)
	{
		double const latitude = std::abs(SkyLatitude(TexelCentre(size, texel)));
		double const omega = kre_starscape::TexelSolidAngle(size, texel % size, (texel / size) % size);
		if (latitude < 5.0)
		{
			planeRadiance += Luminance(cooked.diffuse[texel]) * omega;
			planeAngle += omega;
		}
		else if (latitude > 85.0)
		{
			poleRadiance += Luminance(cooked.diffuse[texel]) * omega;
			poleAngle += omega;
		}
	}
	// The in-plane column is dust-limited (~1 kpc) while the polar column is 2h ~ 600 pc plus the thick disk,
	// so the diffuse contrast is 2-4, not an order of magnitude (spec: plane/pole >= 2).
	CHECK((planeRadiance / planeAngle) >= 2.0 * (poleRadiance / poleAngle));

	double planeStars = 0.0;
	double highStars = 0.0;
	for (CookedStar const& star : cooked.stars)
	{
		double const latitude = std::abs(SkyLatitude(star.dir));
		planeStars += latitude < 10.0 ? 1.0 : 0.0;
		highStars += latitude > 60.0 ? 1.0 : 0.0;
	}
	double const planeBand = 4.0 * kPi * std::sin(10.0 * kPi / 180.0);
	double const highCaps = 4.0 * kPi * (1.0 - std::sin(60.0 * kPi / 180.0));
	double const ratio = (planeStars / planeBand) / (highStars / highCaps);
	// Resolved stars are the nearby ones, whose density is only mildly flattened (spec: ratio in [1.15, 4]).
	CHECK(ratio >= 1.15);
	CHECK(ratio <= 4.0);
}

TEST_CASE("Observer at the core")
{
	Scene scene = MilkyWayScene({0.0, 0.0, 0.0});
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 2000;
	Cooked const cooked = Cook(scene, {.threads = Threads()}).value();
	CHECK(cooked.stars.size() >= 1000);
	CHECK(cooked.stars.size() <= 3000);
	for (glm::dvec3 const& radiance : cooked.diffuse)
	{
		REQUIRE(std::isfinite(radiance.x + radiance.y + radiance.z));
	}
}

TEST_CASE("Unreachable budget")
{
	Scene scene = BaseScene();
	Galaxy tiny = *PresetGalaxy("dwarf", scene.populations);
	tiny.bulge.density = 160.0;   // ≈ 580 stars in all
	// Near enough (14° across) that the coarse calibration texels cannot miss it, far enough
	// that the threshold floor still resolves every one of its stars.
	tiny.position = {20000.0, 0.0, 0.0};
	scene.galaxies.push_back(tiny);
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 4000;
	Cooked const cooked = Cook(scene, {.threads = Threads()}).value();
	CHECK(cooked.stats.budgetUnreachable);
	CHECK(cooked.stars.size() > 0);
	CHECK(cooked.stars.size() < 2000);
	CHECK(cooked.stats.diffuseFlux < 1e-3 * cooked.stats.pointFlux);
}

TEST_CASE("A source in a known direction lands in that direction in the file and the preview")
{
	// Mutation: swapping x and y in Decompose (cook.cpp) transposes every stored face and reddens the diffuse
	// check; each other test derives its expected direction through the same Decompose, so none could notice.
	glm::dvec3 const dir = glm::normalize(glm::dvec3(0.3, 0.8, 0.52));   // off every axis and diagonal

	// Everything diffuse: the lump is a smudge in one or two texels.
	Cooked const smudge = Cook(LumpScene(dir, -30.0), {.threads = Threads()}).value();
	CHECK(smudge.stars.empty());
	std::vector<uint8_t> const smudgeBytes = BuildFile(smudge, 0);
	kre_starscape::View smudgeView;
	REQUIRE(kre_starscape::Parse(smudgeBytes, smudgeView) == kre_starscape::Result::Ok);
	uint32_t const size = smudgeView.header.diffuseFaceSize;
	double const texelDegrees = 90.0 / size;
	uint32_t bestTexel = 0;
	uint32_t bestValue = 0;
	for (uint32_t texel = 0; texel < 6 * size * size; ++texel)
	{
		uint32_t const face = texel / (size * size);
		uint32_t const packed = kre_starscape::DiffuseTexel(smudgeView, 0, face, texel % size, (texel / size) % size);
		auto const v = kre_starscape::DecodeRGB9E5(packed);
		uint32_t const value = uint32_t(1000.0f * (v[0] + v[1] + v[2]));
		if (value > bestValue)
		{
			bestValue = value;
			bestTexel = texel;
		}
	}
	REQUIRE(bestValue > 0);
	// The file's own face table, not Decompose, says which texel holds the direction.
	kre_starscape::FaceCoord const expected = kre_starscape::DirToFace(dir.x, dir.y, dir.z);
	CHECK(bestTexel / (size * size) == expected.face);
	CHECK(std::abs(int(bestTexel % size) - int(expected.s * size)) <= 1);
	CHECK(std::abs(int((bestTexel / size) % size) - int(expected.t * size)) <= 1);

	Image const diffusePreview = RenderPreview(smudgeView, 512, 0.0, Threads());
	CHECK(AngleDegrees(BrightestPixelDirection(diffusePreview), dir) < 1.5 * texelDegrees);

	// Everything a point: the stars sit on the direction, in the file and in the splats.
	Cooked const stars = Cook(LumpScene(dir, 30.0), {.threads = Threads()}).value();
	REQUIRE(stars.stars.size() > 10);
	std::vector<uint8_t> const starBytes = BuildFile(stars, 0);
	kre_starscape::View starView;
	REQUIRE(kre_starscape::Parse(starBytes, starView) == kre_starscape::Result::Ok);
	REQUIRE(starView.header.starCount == stars.stars.size());
	for (uint32_t i = 0; i < starView.header.starCount; ++i)
	{
		kre_starscape::Star const star = kre_starscape::StarAt(starView, i);
		CHECK(AngleDegrees({star.dir[0], star.dir[1], star.dir[2]}, dir) < 2.5);
	}
	Image const starPreview = RenderPreview(starView, 1024, 0.0, Threads());
	CHECK(AngleDegrees(BrightestPixelDirection(starPreview), dir) < 2.5);
}

TEST_CASE("A limiting magnitude that resolves too many stars is rejected up front")
{
	Scene scene = BaseScene();
	scene.output.limitingMagnitude = 25.0;
	int const population = AddSingleBinPopulation(scene, 12);
	scene.galaxies.push_back(UniformBall(0.01, 2000.0, population, {0.0, 0.0, 0.0}));   // about 3e8 stars
	Expected<Cooked> const cooked = Cook(scene, {.threads = Threads()});
	REQUIRE_FALSE(cooked.has_value());
	CHECK(cooked.error().kind == Error::Kind::Input);
	CHECK(cooked.error().message.find("magnitude 25.00") != std::string::npos);
	CHECK(cooked.error().message.find("stars") != std::string::npos);
}

TEST_CASE("Poisson draws far past uint32 stay near the largest count instead of overflowing")
{
	// Mutation: removing the clamp before the cast in SamplePoisson makes this undefined (UBSan reddens it).
	RandomStream rng(9, Purpose::Test, 1, 0, 0);
	for (double const lambda : {1.5e10, double(HUGE_VAL)})
	{
		for (int i = 0; i < 20; ++i)
		{
			uint32_t const k = SamplePoisson(lambda, rng);
			CHECK(k > 4000000000u);
		}
	}
}
