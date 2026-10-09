#include "cook.h"
#include "test_scenes.h"
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
		CHECK(mean == doctest::Approx(lambda).epsilon(0.02));
		CHECK(sumSquares / n - mean * mean == doctest::Approx(lambda).epsilon(0.03));
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
			CHECK(hits[y][x] / n == doctest::Approx(share).epsilon(0.02));
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
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	scene.output.faceSize = 32;
	scene.output.limitingMagnitude.reset();
	std::vector<Cooked> cooks;
	for (uint32_t const target : {500u, 4000u, 20000u})
	{
		scene.output.targetStars = target;
		cooks.push_back(Cook(scene, {.threads = Threads()}));
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
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});   // spin axis = galactic pole, so sky latitude = galactic latitude
	scene.output.faceSize = 64;
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 20000;
	Cooked const cooked = Cook(scene, {.threads = Threads()});

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
	// so the diffuse contrast is a few, not an order of magnitude (spec: plane/pole >= 2.5).
	CHECK((planeRadiance / planeAngle) >= 2.5 * (poleRadiance / poleAngle));

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
	Cooked const cooked = Cook(scene, {.threads = Threads()});
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
	Cooked const cooked = Cook(scene, {.threads = Threads()});
	CHECK(cooked.stats.budgetUnreachable);
	CHECK(cooked.stars.size() > 0);
	CHECK(cooked.stars.size() < 2000);
	CHECK(cooked.stats.diffuseFlux < 1e-3 * cooked.stats.pointFlux);
}
