#pragma once
#include "population.h"
#include <vector>

namespace starscape
{

// In order: young-disk, old-disk, bulge, halo. Scenes index populations, so the order is stable.
std::vector<Population> PresetPopulations();

}
