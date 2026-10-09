#pragma once
#include <array>
#include <cstddef>
#include <cstdint>

namespace starscape
{

using PhiloxCounter = std::array<uint32_t, 4>;
using PhiloxKey     = std::array<uint32_t, 2>;

// Philox4x32-10 (Salmon et al., SC11): a counter-based generator, so a draw is a pure
// function of its key and counter and never of which thread reached it first.
constexpr PhiloxCounter Philox4x32_10(PhiloxCounter counter, PhiloxKey key)
{
	constexpr uint32_t kMultiplier0 = 0xD2511F53u;
	constexpr uint32_t kMultiplier1 = 0xCD9E8D57u;
	constexpr uint32_t kWeyl0       = 0x9E3779B9u;
	constexpr uint32_t kWeyl1       = 0xBB67AE85u;
	for (int round = 0; round < 10; ++round)
	{
		uint64_t const product0 = uint64_t(kMultiplier0) * counter[0];
		uint64_t const product1 = uint64_t(kMultiplier1) * counter[2];
		counter = {
			uint32_t(product1 >> 32) ^ counter[1] ^ key[0],
			uint32_t(product1),
			uint32_t(product0 >> 32) ^ counter[3] ^ key[1],
			uint32_t(product0)};
		key[0] += kWeyl0;
		key[1] += kWeyl1;
	}
	return counter;
}

enum class Purpose : uint32_t
{
	Jitter  = 1,
	Poisson = 2,
	Star    = 3,
	Test    = 4,
};

class RandomStream
{
public:
	// `c` must fit in 24 bits: the purpose occupies the top byte of its counter word.
	RandomStream(uint64_t seed, Purpose purpose, uint32_t a, uint32_t b, uint32_t c)
	: m_key{uint32_t(seed), uint32_t(seed >> 32)}
	, m_counter{a, b, (uint32_t(purpose) << 24) | (c & 0x00FFFFFFu), 0u}
	{
	}

	uint32_t NextU32()
	{
		if (m_used == 4)
		{
			m_block = Philox4x32_10(m_counter, m_key);
			++m_counter[3];
			m_used = 0;
		}
		return m_block[size_t(m_used++)];
	}

	// Uniform in [0, 1) with 53 random bits.
	double NextDouble()
	{
		uint64_t const high = NextU32() >> 5;
		uint64_t const low  = NextU32() >> 6;
		return double(high * 67108864u + low) * (1.0 / 9007199254740992.0);
	}

private:
	PhiloxKey     m_key;
	PhiloxCounter m_counter;
	PhiloxCounter m_block{};
	int           m_used = 4;
};

}
