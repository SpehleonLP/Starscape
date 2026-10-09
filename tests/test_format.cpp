#include <kre_starscape.h>
#include <doctest.h>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

using namespace kre_starscape;

namespace
{

std::vector<uint8_t> ValidFile(uint32_t faceSize = 2, uint32_t starCount = 3)
{
	Header h{};
	h.magic           = kMagic;
	h.version         = kVersion;
	h.diffuseEncoding = uint32_t(Encoding::RGB9E5);
	h.diffuseFaceSize = faceSize;
	h.diffuseMipCount = MipCountFor(faceSize);
	h.starCount       = starCount;
	h.diffuseOffset   = sizeof(Header);
	h.diffuseBytes    = DiffuseBytesFor(faceSize);
	h.starsOffset     = h.diffuseOffset + h.diffuseBytes;
	h.fluxThreshold   = 0.01f;
	h.valueScale      = 1.0f;
	std::vector<uint8_t> bytes(size_t(h.starsOffset + uint64_t(starCount) * sizeof(Star)));
	std::memcpy(bytes.data(), &h, sizeof h);
	for (size_t i = sizeof h; i < bytes.size(); ++i)
	{
		bytes[i] = uint8_t(i * 31u);
	}
	return bytes;
}

template<class Mutate>
Result ParseMutated(Mutate&& mutate)
{
	std::vector<uint8_t> bytes = ValidFile();
	Header h;
	std::memcpy(&h, bytes.data(), sizeof h);
	mutate(h);
	std::memcpy(bytes.data(), &h, sizeof h);
	View view;
	return Parse(bytes, view);
}

}

TEST_CASE("A well-formed file parses into spans over the caller's bytes")
{
	std::vector<uint8_t> const bytes = ValidFile();
	View view;
	REQUIRE(Parse(bytes, view) == Result::Ok);
	CHECK(view.diffuse.size() == 120);   // 6 faces × (2×2 + 1×1) texels × 4 bytes
	CHECK(view.stars.size() == 48);
	CHECK(view.diffuse.data() == bytes.data() + 64);
	CHECK(DiffuseBytesFor(256) == 2097144);
	CHECK(MipCountFor(256) == 9);
	CHECK(MipCountFor(1) == 1);
}

TEST_CASE("Every truncation is rejected")
{
	std::vector<uint8_t> const bytes = ValidFile();
	for (size_t length = 0; length < bytes.size(); ++length)
	{
		View view;
		CHECK(Parse(std::span(bytes.data(), length), view) != Result::Ok);
	}
}

TEST_CASE("Each header defect is named")
{
	// Mutation: dropping the region-overlap check reddens this test.
	CHECK(ParseMutated([](Header& h) { h.magic ^= 1u; }) == Result::BadMagic);
	CHECK(ParseMutated([](Header& h) { h.version = 2; }) == Result::BadVersion);
	CHECK(ParseMutated([](Header& h) { h.diffuseEncoding = uint32_t(Encoding::BC6H); }) == Result::BadEncoding);
	CHECK(ParseMutated([](Header& h) { h.diffuseFaceSize = 3; }) == Result::BadFaceSize);
	CHECK(ParseMutated([](Header& h) { h.diffuseFaceSize = 0; }) == Result::BadFaceSize);
	CHECK(ParseMutated([](Header& h) { h.diffuseMipCount = 1; }) == Result::BadMipCount);
	CHECK(ParseMutated([](Header& h) { h.valueScale = 0.0f; }) == Result::BadValueScale);
	CHECK(ParseMutated([](Header& h) { h.valueScale = std::numeric_limits<float>::quiet_NaN(); }) == Result::BadValueScale);
	CHECK(ParseMutated([](Header& h) { h.diffuseBytes += 4; }) == Result::DiffuseSizeMismatch);
	CHECK(ParseMutated([](Header& h) { h.diffuseOffset = 8; }) == Result::DiffuseOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.diffuseOffset = 1000; }) == Result::DiffuseOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.starsOffset = ~uint64_t(0) - 8; }) == Result::StarsOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.starCount = ~0u; }) == Result::StarsOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.starsOffset = h.diffuseOffset + 4; }) == Result::RegionsOverlap);
}

TEST_CASE("Accessors read the documented byte positions")
{
	std::vector<uint8_t> bytes = ValidFile();
	uint32_t const texelValue = 0x12345678u;
	// mip 1, face NegY (3), texel (0,0): mip 0 holds 6 × 2 × 2 texels.
	size_t const texelOffset = 64 + size_t(MipOffset(2, 1)) + 3u * 4u;
	CHECK(MipOffset(2, 1) == 96);
	std::memcpy(bytes.data() + texelOffset, &texelValue, 4);
	Star const star{{0.0f, 0.6f, 0.8f}, 0xABCDEF01u};
	std::memcpy(bytes.data() + 64 + 120 + 2 * sizeof(Star), &star, sizeof star);
	View view;
	REQUIRE(Parse(bytes, view) == Result::Ok);
	CHECK(DiffuseTexel(view, 1, NegY, 0, 0) == texelValue);
	Star const read = StarAt(view, 2);
	CHECK(read.dir[2] == 0.8f);
	CHECK(read.fluxRGB9E5 == 0xABCDEF01u);
}

TEST_CASE("RGB9E5 encodes the reference values")
{
	// Mutation: a wrong exponent bias reddens this test.
	CHECK(EncodeRGB9E5(1.0f, 1.0f, 1.0f) == 0x84020100u);
	CHECK(EncodeRGB9E5(kRGB9E5Max, kRGB9E5Max, kRGB9E5Max) == 0xFFFFFFFFu);
	CHECK(EncodeRGB9E5(1.0e9f, 0.0f, 0.0f) == EncodeRGB9E5(kRGB9E5Max, 0.0f, 0.0f));
	CHECK(EncodeRGB9E5(0.0f, 0.0f, 0.0f) == 0u);
	CHECK(EncodeRGB9E5(-1.0f, -2.0f, -3.0f) == 0u);
	CHECK(EncodeRGB9E5(std::numeric_limits<float>::quiet_NaN(), 0.0f, 0.0f) == 0u);
	auto const decoded = DecodeRGB9E5(EncodeRGB9E5(1.0f, 0.5f, 0.25f));
	CHECK(decoded[0] == 1.0f);
	CHECK(decoded[1] == 0.5f);
	CHECK(decoded[2] == 0.25f);
}

TEST_CASE("RGB9E5 round trip error stays within half a quantum (<= largest channel / 511.5)")
{
	constexpr int n = 1000000;
	for (int i = 0; i < n; ++i)
	{
		float const big = std::ldexp(1.0f, -14) * std::pow(65408.0f / std::ldexp(1.0f, -14), float(i) / n);
		float const r = big;
		float const g = big * float((i % 97) + 1) / 97.0f;
		float const b = big * float((i % 13)) / 13.0f;
		auto const d = DecodeRGB9E5(EncodeRGB9E5(r, g, b));
		float const tolerance = big / 511.5f + 1e-30f;   // a rounding bump makes the quantum 2^k with big >= 511.5 * 2^(k-1)
		REQUIRE(std::abs(d[0] - r) <= tolerance);
		REQUIRE(std::abs(d[1] - g) <= tolerance);
		REQUIRE(std::abs(d[2] - b) <= tolerance);
	}
}

TEST_CASE("Cube faces follow the Vulkan selection table")
{
	FaceCoord c = DirToFace(1.0, 0.0, 0.0);
	CHECK(c.face == PosX);
	CHECK(c.s == 0.5);
	CHECK(c.t == 0.5);
	c = DirToFace(1.0, 0.5, 0.0);   // +y is up: the upper half of the face
	CHECK(c.face == PosX);
	CHECK(c.t == 0.25);
	c = DirToFace(0.5, 0.0, 1.0);
	CHECK(c.face == PosZ);
	CHECK(c.s == 0.75);
	c = DirToFace(0.0, -1.0, 0.0);
	CHECK(c.face == NegY);
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (double s = 0.05; s < 1.0; s += 0.1)
		{
			for (double t = 0.05; t < 1.0; t += 0.1)
			{
				auto const d = FaceToDir(face, s, t);
				FaceCoord const back = DirToFace(d[0], d[1], d[2]);
				REQUIRE(back.face == face);
				REQUIRE(back.s == doctest::Approx(s));
				REQUIRE(back.t == doctest::Approx(t));
			}
		}
	}
}

TEST_CASE("Texel solid angles cover the sphere")
{
	constexpr uint32_t size = 8;
	double total = 0.0;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t y = 0; y < size; ++y)
		{
			for (uint32_t x = 0; x < size; ++x)
			{
				total += TexelSolidAngle(size, x, y);
			}
		}
	}
	CHECK(total == doctest::Approx(4.0 * 3.14159265358979323846).epsilon(1e-12));
	CHECK(TexelSolidAngle(size, 3, 3) > TexelSolidAngle(size, 0, 0));
}
