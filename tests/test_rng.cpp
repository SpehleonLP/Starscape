#include "rng.h"
#include <doctest.h>

using namespace starscape;

TEST_CASE("Philox4x32-10 matches the Random123 known-answer vectors")
{
	// Mutation: a wrong Philox multiplier reddens this test.
	CHECK(Philox4x32_10({0u, 0u, 0u, 0u}, {0u, 0u}) == PhiloxCounter{0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu, 0x9b00dbd8u});
	CHECK(Philox4x32_10({~0u, ~0u, ~0u, ~0u}, {~0u, ~0u}) == PhiloxCounter{0x408f276du, 0x41c83b0eu, 0xa20bc7c6u, 0x6d5451fdu});
	CHECK(Philox4x32_10({0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u}, {0xa4093822u, 0x299f31d0u})
		== PhiloxCounter{0xd16cfe09u, 0x94fdccebu, 0x5001e420u, 0x24126ea1u});
}

TEST_CASE("RandomStream is a pure function of its key")
{
	RandomStream a(7, Purpose::Star, 1, 2, 3);
	RandomStream b(7, Purpose::Star, 1, 2, 3);
	for (int i = 0; i < 100; ++i)
	{
		CHECK(a.NextU32() == b.NextU32());
	}
	uint32_t const reference = RandomStream(7, Purpose::Star, 1, 2, 3).NextU32();
	CHECK(RandomStream(7, Purpose::Star, 1, 2, 4).NextU32() != reference);
	CHECK(RandomStream(7, Purpose::Poisson, 1, 2, 3).NextU32() != reference);
	CHECK(RandomStream(8, Purpose::Star, 1, 2, 3).NextU32() != reference);
}

TEST_CASE("RandomStream doubles are uniform in [0, 1)")
{
	RandomStream rng(1, Purpose::Test, 0, 0, 0);
	constexpr int n = 200000;
	double sum = 0.0;
	double sumSquares = 0.0;
	for (int i = 0; i < n; ++i)
	{
		double const u = rng.NextDouble();
		REQUIRE(u >= 0.0);
		REQUIRE(u < 1.0);
		sum += u;
		sumSquares += u * u;
	}
	CHECK(sum / n == doctest::Approx(0.5).epsilon(0.01));
	CHECK(sumSquares / n == doctest::Approx(1.0 / 3.0).epsilon(0.01));
}
