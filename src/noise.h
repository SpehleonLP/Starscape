#pragma once
#include <cstdint>
#include <glm/vec3.hpp>

namespace starscape
{

// Perlin gradient noise: roughly in [-1, 1], exactly zero on integer lattice points.
double GradientNoise(glm::dvec3 p, uint64_t seed);

// `octaves` octaves of GradientNoise at doubling frequency and halving amplitude,
// divided by the amplitude sum so the result stays roughly in [-1, 1].
double Fbm(glm::dvec3 p, int octaves, uint64_t seed);

}
