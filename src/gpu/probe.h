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
// against the CPU's. They allocate and free their own buffers and pipeline on every call.
Expected<std::vector<uint32_t>> ProbeLatticeHash(Device& device, std::span<glm::i64vec3 const> cells, uint64_t seed);
Expected<std::vector<std::array<uint32_t, 4>>> ProbePhilox(Device& device, std::span<std::array<uint32_t, 4> const> counters, uint64_t key);
// GalaxySample at each point, evaluated by the GPU in the galaxy's local frame.
Expected<std::vector<GalaxySample>> ProbeEvaluate(Device& device, Galaxy const& galaxy, std::span<glm::dvec3 const> points);

}
