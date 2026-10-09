#pragma once
#include "common.h"
#include "cook.h"
#include <cstdint>
#include <filesystem>
#include <span>
#include <vector>

namespace starscape
{

// levels[0] is `level0`; each further level halves the face size down to 1×1.
std::vector<std::vector<glm::dvec3>> BuildMips(uint32_t faceSize, std::vector<glm::dvec3> const& level0);

// Power of two placing `maxValue` in (16384, 32768], half the RGB9E5 maximum.
float ChooseValueScale(double maxValue);

std::vector<uint8_t> BuildFile(Cooked const& cooked, uint64_t sceneHash);

// Writes beside the target and renames into place, so a failure never leaves a half-written file.
Expected<void> WriteFileAtomic(std::filesystem::path const& path, std::span<uint8_t const> bytes);
Expected<std::vector<uint8_t>> ReadFileBytes(std::filesystem::path const& path);

}
