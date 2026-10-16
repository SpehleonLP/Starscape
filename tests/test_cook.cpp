#include "cook.h"
#include "rng.h"
#include "test_scenes.h"
#include <kre_starscape.h>
#include <doctest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string_view>
#include <set>

using namespace starscape;
using namespace starscape::test;

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

TEST_CASE("Passes visit every stratum once and every prefix stays spread")
{
	// Mutation: qy = p₁ without the XOR with p₀ makes passes 0 and 1 share a row.
	std::set<std::pair<int, int>> seen;
	for (int pass = 0; pass < 16; ++pass)
	{
		Stratum const s = PassStratum(4, pass);
		CHECK(s.x >= 0);
		CHECK(s.x < 4);
		CHECK(s.y >= 0);
		CHECK(s.y < 4);
		seen.insert({s.x, s.y});
	}
	CHECK(seen.size() == 16);
	Stratum const first = PassStratum(4, 0);
	Stratum const second = PassStratum(4, 1);
	CHECK(first.x / 2 != second.x / 2);
	CHECK(first.y / 2 != second.y / 2);
	std::array<int, 4> quadrants{};
	for (int pass = 0; pass < 8; ++pass)
	{
		Stratum const s = PassStratum(4, pass);
		++quadrants[size_t((s.y / 2) * 2 + s.x / 2)];
		if (pass == 3)
		{
			CHECK(quadrants == std::array<int, 4>{1, 1, 1, 1});
		}
	}
	CHECK(quadrants == std::array<int, 4>{2, 2, 2, 2});
	CHECK(PassStratum(1, 0).x == 0);
	CHECK(PassStratum(1, 0).y == 0);
}

TEST_CASE("Tiles grow toward a quarter second and shrink after a slow one")
{
	// Mutation: sizing by seconds / target instead of target / seconds grows tiles after a slow dispatch.
	CHECK(NextTileRows(8, 0.0625, 512) == 16);    // 4× too fast, but growth is capped at doubling
	CHECK(NextTileRows(32, 4.0, 512) == 2);
	CHECK(NextTileRows(8, 0.0, 512) == 16);     // too fast to time: double, never divide by zero
	CHECK(NextTileRows(400, 0.01, 512) == 512);   // never past the face
	CHECK(NextTileRows(1, 100.0, 512) == 1);      // never below one row
}

TEST_CASE("Without a GPU a cook refuses before the calibration starts")
{
	// Mutation: opening the GPU march after the calibration (cook.cpp) spends the calibration first and reddens the
	// progress count. STARSCAPE_GPU_DISABLE makes every machine one without a GPU, so this needs no device.
	struct DisableGpu
	{
		DisableGpu() { setenv("STARSCAPE_GPU_DISABLE", "1", 1); }
		~DisableGpu() { unsetenv("STARSCAPE_GPU_DISABLE"); }
	} const disabled;
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 300;   // a calibrated threshold, so the calibration is the next thing a cook does
	int calibrateReports = 0;
	CookOptions options;
	options.progress = [&](std::string_view phase, double)
	{
		calibrateReports += phase == "calibrate" ? 1 : 0;
	};
	Expected<Cooked> const cooked = Cook(scene, options);
	REQUIRE_FALSE(cooked.has_value());
	CHECK(cooked.error().kind == Error::Kind::Input);
	CHECK(cooked.error().message.find("no GPU to march on: the GPU is disabled by STARSCAPE_GPU_DISABLE") != std::string::npos);
	CHECK(calibrateReports == 0);
}
