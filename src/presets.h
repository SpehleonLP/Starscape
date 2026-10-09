#pragma once
#include "galaxy.h"
#include "population.h"
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace starscape
{

// In order: young-disk, old-disk, bulge, halo. Scenes index populations, so the order is stable.
std::vector<Population> PresetPopulations();

// Every field filled; a scene overrides any of them. Population indices resolve against `populations`.
std::optional<Galaxy> PresetGalaxy(std::string_view name, std::span<Population const> populations);

}
