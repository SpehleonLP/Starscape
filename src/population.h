#pragma once
#include <array>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <glm/vec3.hpp>
#include "field/step.glsl.h"

namespace starscape
{

inline constexpr double kLogLuminosityMin  = -4.0;   // log10 L_V of the lowest bin edge
inline constexpr double kLogLuminositySpan = 10.0;   // decades covered by the bins

double BinLuminosity(int bin);   // geometric centre of the bin, L_V☉
double BinEdge(int edge);        // edge in [0, kBins]

struct Population
{
	std::string                   name;
	std::array<double, kBins>     fraction{};      // number fraction of the population's stars per bin, sums to 1
	std::array<double, kBins>     temperature{};   // representative effective temperature per bin, kelvin
};

glm::dvec3 BlackbodyChroma(double kelvin);
// `saturation` scales each bin's colour away from white (1: physical); every bin keeps its luminance.
PopulationTable Tabulate(Population const& population, double saturation = 1.0);
std::optional<int> FindPopulation(std::span<Population const> populations, std::string_view name);

}
