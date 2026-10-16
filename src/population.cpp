#include "population.h"
#include <cmath>
#include <cstddef>
#include <glm/common.hpp>

namespace starscape
{

double BinLuminosity(int bin)
{
	return std::pow(10.0, kLogLuminosityMin + kLogLuminositySpan * (bin + 0.5) / kBins);
}

double BinEdge(int edge)
{
	return std::pow(10.0, kLogLuminosityMin + kLogLuminositySpan * edge / kBins);
}

namespace
{

// Wyman, Sloan & Shirley (2013), "Simple Analytic Approximations to the CIE XYZ Color Matching
// Functions": a piecewise-Gaussian fit to the CIE 1931 2° observer, no 471-row table needed.
double Lobe(double nm, double mu, double sigmaBelow, double sigmaAbove)
{
	double const t = (nm - mu) / (nm < mu ? sigmaBelow : sigmaAbove);
	return std::exp(-0.5 * t * t);
}

glm::dvec3 ColourMatch(double nm)
{
	double const x = 1.056 * Lobe(nm, 599.8, 37.9, 31.0) + 0.362 * Lobe(nm, 442.0, 16.0, 26.7) - 0.065 * Lobe(nm, 501.1, 20.4, 26.2);
	double const y = 0.821 * Lobe(nm, 568.8, 46.9, 40.5) + 0.286 * Lobe(nm, 530.9, 16.3, 31.1);
	double const z = 1.217 * Lobe(nm, 437.0, 11.8, 36.0) + 0.681 * Lobe(nm, 459.0, 26.0, 13.8);
	return {x, y, z};
}

// Planck's law up to a constant factor; the chroma normalisation removes the constant.
double Planck(double nm, double kelvin)
{
	constexpr double kSecondRadiation = 1.438776877e-2;   // h c / k, metre kelvin
	double const metres = nm * 1e-9;
	return 1.0 / (std::pow(metres, 5) * std::expm1(kSecondRadiation / (metres * kelvin)));
}

// Chroma has luminance 1, so white is (1, 1, 1) and Y(1 + s(c - 1)) = 1: inside the gamut the push keeps
// each bin's light by construction. Clipping a negative channel to zero breaks that, so renormalise. The
// clipped vector keeps a channel >= 1 (the one that was >= 1 before), so its luminance stays positive.
glm::dvec3 Saturate(glm::dvec3 chroma, double saturation)
{
	if (saturation == 1.0)
	{
		return chroma;
	}
	glm::dvec3 const pushed = glm::max(glm::dvec3(1.0) + saturation * (chroma - glm::dvec3(1.0)), glm::dvec3(0.0));
	return pushed / Luminance(pushed);
}

}

glm::dvec3 BlackbodyChroma(double kelvin)
{
	glm::dvec3 xyz(0.0);
	for (int nm = 360; nm <= 830; ++nm)
	{
		xyz += Planck(nm, kelvin) * ColourMatch(nm);
	}
	glm::dvec3 rgb(
		3.2404542 * xyz.x - 1.5371385 * xyz.y - 0.4985314 * xyz.z,
		-0.9692660 * xyz.x + 1.8760108 * xyz.y + 0.0415560 * xyz.z,
		0.0556434 * xyz.x - 0.2040259 * xyz.y + 1.0572252 * xyz.z);
	// The hottest and coolest blackbodies sit just outside the Rec.709 gamut.
	rgb = glm::max(rgb, glm::dvec3(0.0));
	return rgb / Luminance(rgb);
}

PopulationTable Tabulate(Population const& population, double saturation)
{
	PopulationTable table;
	for (int bin = 0; bin < kBins; ++bin)
	{
		size_t const i = size_t(bin);
		table.fraction[i]   = population.fraction[i];
		table.luminosity[i] = BinLuminosity(bin);
		table.chroma[i]     = Saturate(BlackbodyChroma(population.temperature[i]), saturation);
		table.lightBelow[i + 1] = table.lightBelow[i] + table.fraction[i] * table.luminosity[i] * table.chroma[i];
	}
	return table;
}

std::optional<int> FindPopulation(std::span<Population const> populations, std::string_view name)
{
	for (size_t i = 0; i < populations.size(); ++i)
	{
		if (populations[i].name == name)
		{
			return int(i);
		}
	}
	return std::nullopt;
}

}
