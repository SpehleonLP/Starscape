// Perlin gradient noise and its octave sums, shared by the CPU and the GPU.
#ifndef STARSCAPE_FIELD_NOISE_GLSL_H
#define STARSCAPE_FIELD_NOISE_GLSL_H

#include "common.glsl.h"

#ifdef __cplusplus
namespace starscape
{
#endif

FIELD_INLINE uint64_t HashMix(uint64_t h, uint64_t v)
{
	h ^= v + 0x9E3779B97F4A7C15ul + (h << 6) + (h >> 2);
	h *= 0xBF58476D1CE4E5B9ul;
	h ^= h >> 31;
	return h;
}

// 64-bit integer mixing: unlike a sin() hash it is exact at any lattice coordinate. The GPU keeps the
// 64-bit form because every approved sky's noise was tuned on it.
FIELD_INLINE uint LatticeHash(int64_t x, int64_t y, int64_t z, uint64_t seed)
{
	uint64_t h = seed * 0x9E3779B97F4A7C15ul + 0xD1B54A32D192ED03ul;
	h = HashMix(h, uint64_t(x));
	h = HashMix(h, uint64_t(y));
	h = HashMix(h, uint64_t(z));
	return uint(h >> 32);
}

FIELD_INLINE real Fade(real t)
{
	return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

FIELD_INLINE real Lerp(real a, real b, real t)
{
	return a + (b - a) * t;
}

// The twelve edge gradients, as g[0] fx + g[1] fy + g[2] fz in that order.
FIELD_INLINE real Corner(int64_t ix, int64_t iy, int64_t iz, real fx, real fy, real fz, uint64_t seed)
{
	const uint g = LatticeHash(ix, iy, iz, seed) % 12u;
	real gx = 0.0;
	real gy = 0.0;
	real gz = 0.0;
	if (g < 4u)
	{
		gx = (g == 0u || g == 2u) ? 1.0 : -1.0;
		gy = (g < 2u) ? 1.0 : -1.0;
	}
	else if (g < 8u)
	{
		gx = (g == 4u || g == 6u) ? 1.0 : -1.0;
		gz = (g < 6u) ? 1.0 : -1.0;
	}
	else
	{
		gy = (g == 8u || g == 10u) ? 1.0 : -1.0;
		gz = (g < 10u) ? 1.0 : -1.0;
	}
	return gx * fx + gy * fy + gz * fz;
}

// Gradient noise is zero on its lattice, so its variance sags on lattice planes, and doubling keeps
// z = 0 on one at every octave: a galaxy's midplane, where most of its dust is, would be the least
// structured place in it. An irrational shift per octave moves every octave off the lattice.
const real3 kOctaveShift = real3(0.5497004779, 0.6710436067, 0.8191725134);

FIELD_INLINE real GradientNoise(real3 p, uint64_t seed)
{
	const real x0 = floor(p.x);
	const real y0 = floor(p.y);
	const real z0 = floor(p.z);
	const int64_t ix = int64_t(x0);
	const int64_t iy = int64_t(y0);
	const int64_t iz = int64_t(z0);
	const real fx = p.x - x0;
	const real fy = p.y - y0;
	const real fz = p.z - z0;
	const real u = Fade(fx);
	const real v = Fade(fy);
	const real w = Fade(fz);
	const real c000 = Corner(ix, iy, iz, fx, fy, fz, seed);
	const real c100 = Corner(ix + 1, iy, iz, fx - 1.0, fy, fz, seed);
	const real c010 = Corner(ix, iy + 1, iz, fx, fy - 1.0, fz, seed);
	const real c110 = Corner(ix + 1, iy + 1, iz, fx - 1.0, fy - 1.0, fz, seed);
	const real c001 = Corner(ix, iy, iz + 1, fx, fy, fz - 1.0, seed);
	const real c101 = Corner(ix + 1, iy, iz + 1, fx - 1.0, fy, fz - 1.0, seed);
	const real c011 = Corner(ix, iy + 1, iz + 1, fx, fy - 1.0, fz - 1.0, seed);
	const real c111 = Corner(ix + 1, iy + 1, iz + 1, fx - 1.0, fy - 1.0, fz - 1.0, seed);
	return Lerp(
		Lerp(Lerp(c000, c100, u), Lerp(c010, c110, u), v),
		Lerp(Lerp(c001, c101, u), Lerp(c011, c111, u), v),
		w);
}

// `octaves` octaves at doubling frequency and halving amplitude, divided by the amplitude sum.
FIELD_INLINE real Fbm(real3 p, int octaves, uint64_t seed)
{
	real sum = 0.0;
	real amplitude = 1.0;
	real norm = 0.0;
	p += kOctaveShift;
	for (int octave = 0; octave < octaves; ++octave)
	{
		sum += amplitude * GradientNoise(p, seed + uint64_t(octave) * 0x632BE59BD9B4E019ul);
		norm += amplitude;
		amplitude *= 0.5;
		p = p * 2.0 + kOctaveShift;
	}
	return sum / norm;
}

// Squared, so the ridges where the noise crosses zero are thin and the dust they carry is filaments.
FIELD_INLINE real RidgedFbm(real3 p, int octaves, uint64_t seed)
{
	real sum = 0.0;
	real amplitude = 1.0;
	real norm = 0.0;
	p += kOctaveShift;
	for (int octave = 0; octave < octaves; ++octave)
	{
		const real n = GradientNoise(p, seed + uint64_t(octave) * 0x632BE59BD9B4E019ul);
		sum += amplitude * ((1.0 - abs(n)) * (1.0 - abs(n)));
		norm += amplitude;
		amplitude *= 0.5;
		p = p * 2.0 + kOctaveShift;
	}
	return sum / norm;
}

#ifdef __cplusplus
}
#endif

#endif
