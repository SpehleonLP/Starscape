#include "writer.h"
#include "test_scenes.h"
#include <kre_starscape.h>
#include <doctest.h>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <unistd.h>

using namespace starscape;
using namespace starscape::test;

namespace
{

Cooked Synthetic(uint32_t faceSize)
{
	Cooked cooked;
	cooked.faceSize = faceSize;
	cooked.fluxThreshold = 0.004;
	RandomStream rng(2, Purpose::Test, 0, 0, 0);
	for (uint32_t i = 0; i < 6 * faceSize * faceSize; ++i)
	{
		cooked.diffuse.push_back({rng.NextDouble() * 100.0, rng.NextDouble() * 100.0, rng.NextDouble() * 100.0});
	}
	cooked.stars.push_back({{0.0, 0.0, 1.0}, {3.0, 2.0, 1.0}});
	cooked.stars.push_back({{0.6, 0.8, 0.0}, {0.01, 0.02, 0.03}});
	return cooked;
}

glm::dvec3 Decoded(kre_starscape::View const& view, uint32_t mip, uint32_t face, uint32_t x, uint32_t y)
{
	auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::DiffuseTexel(view, mip, face, x, y));
	return glm::dvec3(v[0], v[1], v[2]) / double(view.header.valueScale);
}

void ExpectClose(glm::dvec3 actual, glm::dvec3 expected)
{
	double const tolerance = std::max({expected.x, expected.y, expected.z}) / 256.0;
	CHECK(std::abs(actual.x - expected.x) <= tolerance);
	CHECK(std::abs(actual.y - expected.y) <= tolerance);
	CHECK(std::abs(actual.z - expected.z) <= tolerance);
}

}

TEST_CASE("Mips are solid-angle-weighted means")
{
	// Mutation: unweighted mips reddens this test.
	std::vector<glm::dvec3> level0(6 * 4 * 4, glm::dvec3(0.0));
	level0[0] = {4.0, 0.0, 0.0};   // face 0, texel (0,0): a corner of the face
	auto const levels = BuildMips(4, level0);
	REQUIRE(levels.size() == 3);
	CHECK(levels[1].size() == 6 * 2 * 2);
	double const weights = kre_starscape::TexelSolidAngle(4, 0, 0) + kre_starscape::TexelSolidAngle(4, 1, 0)
		+ kre_starscape::TexelSolidAngle(4, 0, 1) + kre_starscape::TexelSolidAngle(4, 1, 1);
	CHECK(levels[1][0].x == doctest::Approx(4.0 * kre_starscape::TexelSolidAngle(4, 0, 0) / weights));
	CHECK(levels[2].size() == 6);
}

TEST_CASE("A cooked sky round-trips through the file")
{
	Cooked const cooked = Synthetic(4);
	std::vector<uint8_t> const bytes = BuildFile(cooked, 0x1234);
	kre_starscape::View view;
	REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
	CHECK(view.header.sceneHash == 0x1234);
	CHECK(view.header.diffuseFaceSize == 4);
	CHECK(view.header.starCount == 2);
	CHECK(view.header.fluxThreshold == doctest::Approx(0.004));
	for (uint32_t texel = 0; texel < 6 * 16; ++texel)
	{
		ExpectClose(Decoded(view, 0, texel / 16, texel % 4, (texel / 4) % 4), cooked.diffuse[texel]);
	}
	auto const levels = BuildMips(4, cooked.diffuse);
	ExpectClose(Decoded(view, 1, 5, 1, 1), levels[1][5 * 4 + 1 * 2 + 1]);
	// The 1x1 level, against a weighted mean computed here rather than by BuildMips.
	for (uint32_t face = 0; face < 6; ++face)
	{
		glm::dvec3 sum(0.0);
		double weight = 0.0;
		for (uint32_t y = 0; y < 4; ++y)
		{
			for (uint32_t x = 0; x < 4; ++x)
			{
				double const omega = kre_starscape::TexelSolidAngle(4, x, y);
				sum += omega * cooked.diffuse[(face * 4 + y) * 4 + x];
				weight += omega;
			}
		}
		ExpectClose(Decoded(view, 2, face, 0, 0), sum / weight);
	}
	for (uint32_t i = 0; i < cooked.stars.size(); ++i)
	{
		kre_starscape::Star const star = kre_starscape::StarAt(view, i);
		CHECK(star.dir[0] == float(cooked.stars[i].dir.x));
		CHECK(star.dir[1] == float(cooked.stars[i].dir.y));
		CHECK(star.dir[2] == float(cooked.stars[i].dir.z));
		auto const flux = kre_starscape::DecodeRGB9E5(star.fluxRGB9E5);
		ExpectClose(glm::dvec3(flux[0], flux[1], flux[2]) / double(view.header.valueScale), cooked.stars[i].flux);
	}
}

TEST_CASE("A dark sky keeps its precision")
{
	// Nine decades of radiance (30 stops) next to a bright star: each texel must decode
	// within one part in 256 of itself, however far below the brightest value it sits.
	Cooked cooked;
	cooked.faceSize = 4;
	for (uint32_t i = 0; i < 6 * 16; ++i)
	{
		double const value = std::pow(10.0, -4.0 + 9.0 * double(i) / 95.0);
		cooked.diffuse.push_back({value, 0.5 * value, 0.25 * value});
	}
	cooked.stars.push_back({{1.0, 0.0, 0.0}, {100.0, 100.0, 100.0}});
	std::vector<uint8_t> const bytes = BuildFile(cooked, 0);
	kre_starscape::View view;
	REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
	for (uint32_t texel = 0; texel < 6 * 16; ++texel)
	{
		glm::dvec3 const decoded = Decoded(view, 0, texel / 16, texel % 4, (texel / 4) % 4);
		CHECK(decoded.x == doctest::Approx(cooked.diffuse[texel].x).epsilon(1.0 / 256.0).scale(0.0));
	}
}

TEST_CASE("The value scale is a power of two that keeps the maximum in range")
{
	// Mutation: the unclamped exponent cast reddens this test.
	for (double const maxValue : {1e-9, 0.37, 1.0, 400.0, 1e7})
	{
		float const scale = ChooseValueScale(maxValue);
		int exponent = 0;
		CHECK(std::frexp(scale, &exponent) == 0.5f);
		CHECK(maxValue * scale <= 32768.0);
		CHECK(maxValue * scale > 16384.0 * 0.999);
	}
	CHECK(ChooseValueScale(0.0) == 1.0f);
	CHECK(ChooseValueScale(std::numeric_limits<double>::infinity()) == 1.0f);
	CHECK(ChooseValueScale(std::numeric_limits<double>::quiet_NaN()) == 1.0f);
	CHECK(ChooseValueScale(1e-310) == std::ldexp(1.0f, 100));
	CHECK(ChooseValueScale(1e300) == std::ldexp(1.0f, -100));
}

TEST_CASE("Cooks are byte-identical on any thread count")
{
	// Mutation: seeding draws by thread index reddens this test.
	Scene const scene = MilkyWayScene({8200.0, 0.0, 20.0});   // limiting magnitude 6, face 16, quick
	std::vector<uint8_t> const one = BuildFile(Cook(scene, {.threads = 1}).value(), scene.hash);
	std::vector<uint8_t> const many = BuildFile(Cook(scene, {.threads = 8}).value(), scene.hash);
	CHECK(one == many);
}

TEST_CASE("Writes land atomically")
{
	// Mutation: leaving the .partial file behind reddens this test.
	std::filesystem::path const dir = std::filesystem::temp_directory_path() / std::format("starscape-writer-test-{}", getpid());
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	std::vector<uint8_t> const bytes = {1, 2, 3, 4};
	REQUIRE(WriteFileAtomic(dir / "sky.kre_starscape", bytes).has_value());
	Expected<std::vector<uint8_t>> const read = ReadFileBytes(dir / "sky.kre_starscape");
	REQUIRE(read.has_value());
	CHECK(*read == bytes);

	SUBCASE("failed rename leaves nothing")
	{
		std::filesystem::create_directories(dir / "occupied.kre_starscape" / "child");
		Expected<void> const result = WriteFileAtomic(dir / "occupied.kre_starscape", bytes);
		REQUIRE_FALSE(result.has_value());
		CHECK(result.error().kind == Error::Kind::Io);
		for (auto const& entry : std::filesystem::directory_iterator(dir))
		{
			CHECK(entry.path().filename().string().find(".partial") == std::string::npos);
		}
	}
	SUBCASE("missing directory")
	{
		Expected<void> const result = WriteFileAtomic(dir / "no" / "such" / "dir" / "sky.kre_starscape", bytes);
		REQUIRE_FALSE(result.has_value());
		CHECK(result.error().kind == Error::Kind::Io);
	}
	CHECK(ReadFileBytes(dir / "absent.kre_starscape").error().kind == Error::Kind::Io);
	Expected<std::vector<uint8_t>> const directory = ReadFileBytes(dir);
	REQUIRE_FALSE(directory.has_value());
	CHECK(directory.error().kind == Error::Kind::Io);
	std::filesystem::remove_all(dir);
}
