#pragma once
#include "common.h"
#include "cook.h"
#include <cstdint>
#include <functional>
#include <memory>

namespace starscape
{

// The GPU march, or an InputError saying why there is none: no qualifying device ("no GPU to march on: " and what
// the machine lacks), too many galaxies, or a face too large for the device's buffers.
// `progress` receives the fraction of the pass done, once per dispatched tile.
// `tileSeconds` is the dispatch duration tiles are sized toward; tests force it to pin the sizing.
Expected<std::unique_ptr<PassMarcher>> OpenGpuMarcher(PreparedScene const& prepared, MarchParams const& params, uint32_t faceSize,
	std::function<void(double)>&& progress, double tileSeconds = kTileSeconds);

}
