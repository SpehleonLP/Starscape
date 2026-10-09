#include "noise.h"
#include <cmath>
#include <initializer_list>

namespace starscape
{
namespace
{

// 64-bit integer mixing: unlike a sin() hash it is exact at any lattice coordinate.
uint32_t LatticeHash(int64_t x, int64_t y, int64_t z, uint64_t seed)
{
	uint64_t h = seed * 0x9E3779B97F4A7C15ull + 0xD1B54A32D192ED03ull;
	for (uint64_t const v : {uint64_t(x), uint64_t(y), uint64_t(z)})
	{
		h ^= v + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
		h *= 0xBF58476D1CE4E5B9ull;
		h ^= h >> 31;
	}
	return uint32_t(h >> 32);
}

constexpr double kGradients[12][3] = {
	{1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
	{1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
	{0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}};

double Fade(double t)
{
	return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

double Lerp(double a, double b, double t)
{
	return a + (b - a) * t;
}

double Corner(int64_t ix, int64_t iy, int64_t iz, double fx, double fy, double fz, uint64_t seed)
{
	double const* g = kGradients[LatticeHash(ix, iy, iz, seed) % 12u];
	return g[0] * fx + g[1] * fy + g[2] * fz;
}

// Gradient noise is zero on its lattice, so its variance sags on lattice planes, and doubling keeps
// z = 0 on one at every octave: a galaxy's midplane, where most of its dust is, would be the least
// structured place in it. An irrational shift per octave moves every octave off the lattice.
constexpr glm::dvec3 kOctaveShift(0.5497004779, 0.6710436067, 0.8191725134);

}

double GradientNoise(glm::dvec3 p, uint64_t seed)
{
	double const x0 = std::floor(p.x);
	double const y0 = std::floor(p.y);
	double const z0 = std::floor(p.z);
	int64_t const ix = int64_t(x0);
	int64_t const iy = int64_t(y0);
	int64_t const iz = int64_t(z0);
	double const fx = p.x - x0;
	double const fy = p.y - y0;
	double const fz = p.z - z0;
	double const u = Fade(fx);
	double const v = Fade(fy);
	double const w = Fade(fz);

	double const c000 = Corner(ix, iy, iz, fx, fy, fz, seed);
	double const c100 = Corner(ix + 1, iy, iz, fx - 1.0, fy, fz, seed);
	double const c010 = Corner(ix, iy + 1, iz, fx, fy - 1.0, fz, seed);
	double const c110 = Corner(ix + 1, iy + 1, iz, fx - 1.0, fy - 1.0, fz, seed);
	double const c001 = Corner(ix, iy, iz + 1, fx, fy, fz - 1.0, seed);
	double const c101 = Corner(ix + 1, iy, iz + 1, fx - 1.0, fy, fz - 1.0, seed);
	double const c011 = Corner(ix, iy + 1, iz + 1, fx, fy - 1.0, fz - 1.0, seed);
	double const c111 = Corner(ix + 1, iy + 1, iz + 1, fx - 1.0, fy - 1.0, fz - 1.0, seed);

	return Lerp(
		Lerp(Lerp(c000, c100, u), Lerp(c010, c110, u), v),
		Lerp(Lerp(c001, c101, u), Lerp(c011, c111, u), v),
		w);
}

namespace
{

template<class Shape>
double Octaves(glm::dvec3 p, int octaves, uint64_t seed, Shape shape)
{
	double sum = 0.0;
	double amplitude = 1.0;
	double norm = 0.0;
	p += kOctaveShift;
	for (int octave = 0; octave < octaves; ++octave)
	{
		sum += amplitude * shape(GradientNoise(p, seed + uint64_t(octave) * 0x632BE59BD9B4E019ull));
		norm += amplitude;
		amplitude *= 0.5;
		p = p * 2.0 + kOctaveShift;
	}
	return sum / norm;
}

}

double Fbm(glm::dvec3 p, int octaves, uint64_t seed)
{
	return Octaves(p, octaves, seed, [](double noise) { return noise; });
}

// Squared, so the ridges where the noise crosses zero are thin and the dust they carry is filaments.
double RidgedFbm(glm::dvec3 p, int octaves, uint64_t seed)
{
	return Octaves(p, octaves, seed, [](double noise) { return (1.0 - std::abs(noise)) * (1.0 - std::abs(noise)); });
}

}
