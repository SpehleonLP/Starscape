#include "noise.h"
#include "rng.h"
#include <doctest.h>
#include <algorithm>
#include <cmath>

using namespace starscape;

TEST_CASE("Gradient noise is zero on the lattice and deterministic")
{
	CHECK(GradientNoise({3.0, -2.0, 7.0}, 5) == 0.0);
	glm::dvec3 const p(0.3, 1.7, -2.2);
	CHECK(GradientNoise(p, 5) == GradientNoise(p, 5));
	CHECK(GradientNoise(p, 5) != GradientNoise(p, 6));
}

TEST_CASE("Fbm stays bounded with a near-zero mean")
{
	RandomStream rng(3, Purpose::Test, 0, 0, 0);
	constexpr int n = 100000;
	double sum = 0.0;
	double extreme = 0.0;
	for (int i = 0; i < n; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1000.0 - 500.0, rng.NextDouble() * 1000.0 - 500.0, rng.NextDouble() * 1000.0 - 500.0);
		double const v = Fbm(p, 6, 11);
		sum += v;
		extreme = std::max(extreme, std::abs(v));
	}
	CHECK(extreme <= 1.1);
	CHECK(extreme > 0.3);
	CHECK(std::abs(sum / n) < 0.02);
}

// The earlier experiment hashed with fract(sin(dot(p, k)) * 43758) in single precision,
// which collapses at coordinates like 250 000. Noise here must still vary there.
TEST_CASE("Noise keeps its variance at galactic coordinates")
{
	RandomStream rng(4, Purpose::Test, 0, 0, 0);
	constexpr int n = 2000;
	double sum = 0.0;
	double sumSquares = 0.0;
	for (int i = 0; i < n; ++i)
	{
		glm::dvec3 const p(250000.0 + rng.NextDouble() * 10.0, 1.0e6 + rng.NextDouble() * 10.0, rng.NextDouble() * 10.0);
		double const v = Fbm(p, 4, 9);
		sum += v;
		sumSquares += v * v;
	}
	double const mean = sum / n;
	CHECK(sumSquares / n - mean * mean > 0.01);
}

// A galaxy's midplane is z = 0; on an unshifted lattice it is a lattice plane at every octave, where
// gradient noise sags toward zero. Unshifted, fbm's variance there drops to 0.021 and ridged fbm's mean
// rises to 0.659, 0.22σ: at the presets' σ 2 that is 1.55 times the calibrated dust along the plane.
TEST_CASE("Noise keeps its statistics on the coordinate planes")
{
	RandomStream rng(6, Purpose::Test, 0, 0, 0);
	constexpr int n = 100000;
	double sumSquares = 0.0;
	double ridged = 0.0;
	for (int i = 0; i < n; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1000.0 - 500.0, rng.NextDouble() * 1000.0 - 500.0, 0.0);
		double const v = Fbm(p, 6, 11);
		sumSquares += v * v;
		ridged += RidgedFbm(p, 6, 11);
	}
	// Over the volume: variance 0.0252, ridged mean 0.6305.
	CHECK(sumSquares / n > 0.023);
	CHECK(std::abs(ridged / n - 0.6305) < 0.01);
}
