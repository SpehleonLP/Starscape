#pragma once
#include "common.h"
#include "cook.h"
#include <cstdint>
#include <functional>
#include <memory>

namespace starscape
{

// The GPU march, when this build has one and a device qualifies; otherwise an InputError saying why.
// `progress` receives the fraction of the pass done, as the CpuMarcher's does, once per dispatched tile.
// `tileSeconds` is the dispatch duration tiles are sized toward; tests force it to pin the sizing.
Expected<std::unique_ptr<PassMarcher>> OpenGpuMarcher(PreparedScene const& prepared, MarchParams const& params, uint32_t faceSize,
	std::function<void(double)> progress, double tileSeconds = kTileSeconds);

}
