#pragma once
#include "common.h"
#include "galaxy.h"
#include "population.h"
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string_view>
#include <vector>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace starscape
{

enum class Quality
{
	Final,
	Quick,
};

struct Output
{
	uint32_t                faceSize = 256;
	std::optional<uint32_t> targetStars;
	std::optional<double>   limitingMagnitude;
	Quality                 quality = Quality::Final;
	uint64_t                seed = 1;   // star draws and jitter; galaxies carry their own noise seeds
	uint32_t                previewWidth = 8192;
	bool                    previewMarch = false;
};

struct Observer
{
	glm::dvec3 position{0.0};
	glm::dvec3 spinAxis{0.0, 0.0, 1.0};       // sky-frame +z
	glm::dvec3 zeroMeridian{1.0, 0.0, 0.0};   // sky-frame +x after removing its spin-axis component
	double     sMin = 0.05;                   // parsecs; nearer stars belong to the observer's own system
};

// Scene-wide departures from the physical sky, for a stylised look.
struct Style
{
	double starSaturation = 1.0;   // 1: physical star colours; above 1 pushes them away from white
};

struct Scene
{
	Output                  output;
	Observer                observer;
	Style                   style;
	std::vector<Galaxy>     galaxies;
	std::vector<Population> populations;   // presets first; tests may append their own
	uint64_t                hash = 0;      // FNV-1a 64 of the scene text
};

Expected<Scene> ParseScene(std::string_view text, std::string_view sourceName);
Expected<Scene> LoadScene(std::filesystem::path const& path);
glm::dmat3 SkyToScene(Observer const& observer);

}
