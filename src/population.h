#pragma once
#include <array>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <glm/vec3.hpp>

namespace starscape
{

inline constexpr int    kBins              = 24;
inline constexpr double kLogLuminosityMin  = -4.0;   // log10 L_V of the lowest bin edge
inline constexpr double kLogLuminositySpan = 10.0;   // decades covered by the bins

// Flux of a 1 L_V☉ source at 10 pc in magnitude-0 units: 10^(−0.4 × 4.83), the Sun's M_V.
inline constexpr double kFluxAt10pc = 0.01169499;

double BinLuminosity(int bin);   // geometric centre of the bin, L_V☉
double BinEdge(int edge);        // edge in [0, kBins]
double Luminance(glm::dvec3 rgb);   // Rec.709 Y; for a star's flux this is the V-band proxy

struct Population
{
	std::string                   name;
	std::array<double, kBins>     fraction{};      // number fraction of the population's stars per bin, sums to 1
	std::array<double, kBins>     temperature{};   // representative effective temperature per bin, kelvin
};

// What the march reads per step, derived once per population.
struct PopulationTable
{
	std::array<glm::dvec3, kBins>     chroma{};       // linear Rec.709, luminance 1
	std::array<double, kBins>         fraction{};
	std::array<double, kBins>         luminosity{};
	// lightBelow[i] = Σ_{j<i} fraction_j × luminosity_j × chroma_j: the diffuse light of every
	// bin below a resolution split, so a step costs one lookup instead of a loop over bins.
	std::array<glm::dvec3, kBins + 1> lightBelow{};
};

glm::dvec3 BlackbodyChroma(double kelvin);
PopulationTable Tabulate(Population const& population);
std::optional<int> FindPopulation(std::span<Population const> populations, std::string_view name);

}
