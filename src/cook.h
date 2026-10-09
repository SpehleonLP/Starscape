#pragma once
#include "march.h"
#include "rng.h"
#include "scene.h"
#include <cstdint>
#include <functional>
#include <optional>
#include <string_view>
#include <vector>

namespace starscape
{

struct CookedStar
{
	glm::dvec3 dir;    // sky frame, unit
	glm::dvec3 flux;   // mag-0 units, linear Rec.709, physical
};

struct CookStats
{
	double   pointFlux = 0.0;      // Σ Luminance(star flux)
	double   diffuseFlux = 0.0;    // Σ Luminance(radiance) × texel solid angle
	uint64_t starCount = 0;
	uint64_t droppedStars = 0;     // drawn, but their direction held no resolved stars on re-march
	bool     budgetUnreachable = false;
};

struct Cooked
{
	uint32_t                faceSize = 0;
	std::vector<glm::dvec3> diffuse;   // mip 0, texel order, physical radiance
	std::vector<CookedStar> stars;     // brightest first
	double                  fluxThreshold = 0.0;
	CookStats               stats;
};

struct CookOptions
{
	unsigned                                                 threads = 1;
	std::optional<Quality>                                   quality{};    // overrides the scene's
	std::function<void(std::string_view phase, double done)> progress{};
};

Cooked Cook(Scene const& scene, CookOptions const& options);

FluxHistogram CalibrationHistogram(PreparedScene const& scene, MarchParams const& params, uint32_t faceSize, unsigned threads);
uint32_t SamplePoisson(double lambda, RandomStream& rng);
glm::dvec3 SampleTexelDirection(uint32_t faceSize, uint32_t texel, RandomStream& rng);
glm::dvec3 TexelCentre(uint32_t faceSize, uint32_t texel);

}
