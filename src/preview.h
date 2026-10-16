#pragma once
#include <kre_starscape.h>
#include <cstdint>
#include <vector>
#include <glm/vec3.hpp>

namespace starscape
{

struct Image
{
	uint32_t           width = 0;
	uint32_t           height = 0;
	std::vector<float> rgb;   // row-major, top row first, physical radiance × 2^EV
};

// Equirectangular: column → longitude −π…π left to right, row → latitude +π/2…−π/2 top to bottom.
double PixelSolidAngle(uint32_t width, uint32_t row);
glm::dvec3 PixelDirection(uint32_t width, double column, double row);   // continuous pixel coordinates

// Bilinear and seamless: taps beyond a face edge are fetched from the neighbouring face.
glm::dvec3 SampleDiffuse(kre_starscape::View const& view, glm::dvec3 dir);

Image RenderPreview(kre_starscape::View const& view, uint32_t width, double exposureEV, unsigned threads);

std::vector<uint8_t> EncodeRadiance(Image const& image);

}
