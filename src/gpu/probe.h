#pragma once
#include "common.h"
#include "field/galaxy.glsl.h"
#include "gpu/device.h"
#include <array>
#include <cstdint>
#include <span>
#include <vector>
#include <glm/vec3.hpp>
#include <glm/gtc/type_precision.hpp>

namespace starscape::gpu
{

// Test-only: one-shot kernels that run a piece of the shared field model on the GPU so a test can hold it
// against the host's compilation of the same code. They allocate and free their own buffers and pipeline on every call.
Expected<std::vector<uint32_t>> ProbeLatticeHash(Device& device, std::span<glm::i64vec3 const> cells, uint64_t seed);
Expected<std::vector<std::array<uint32_t, 4>>> ProbePhilox(Device& device, std::span<std::array<uint32_t, 4> const> counters, uint64_t key);
// One march ray's jitter: stratum (stratumX, stratumY) of perAxis² in texel (x, y) of `face`.
struct JitterQuery
{
	uint32_t faceSize;
	uint32_t face;
	uint32_t x;
	uint32_t y;
	uint32_t stratumX;
	uint32_t stratumY;
	uint32_t perAxis;
};
// The sky-frame direction march.comp fires for each query under `seed`, from the march's own jitter function.
Expected<std::vector<std::array<float, 3>>> ProbeJitterDirection(Device& device, std::span<JitterQuery const> queries, uint64_t seed);
// GalaxySample at each point, evaluated by the GPU in the galaxy's local frame.
Expected<std::vector<GalaxySample>> ProbeEvaluate(Device& device, Galaxy const& galaxy, std::span<glm::dvec3 const> points);

}
