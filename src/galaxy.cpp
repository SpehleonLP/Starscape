#include "galaxy.h"
#include "noise.h"
#include "population.h"
#include <algorithm>
#include <cmath>
#include <glm/common.hpp>
#include <glm/geometric.hpp>
#include <glm/matrix.hpp>

namespace starscape
{
namespace
{

// A log-normal of bounded noise does not average one at large σ, so the mean is measured and divided
// out: structure moves stars and dust around without changing how much there is.
template<class Field>
double SpatialMean(Field field)
{
	// An R3 low-discrepancy walk over a box many noise cells wide, so the samples are independent.
	constexpr int kSamples = 1 << 18;
	double sum = 0.0;
	for (int i = 1; i <= kSamples; ++i)
	{
		sum += field(glm::fract(double(i) * glm::dvec3(0.8191725134, 0.6710436067, 0.5497004779)) * 4096.0);
	}
	return sum / kSamples;
}

double DustFieldMean(Dust const& dust, uint64_t seed)
{
	if (dust.filament == 0.0)
	{
		return 1.0;
	}
	return SpatialMean([&](glm::dvec3 p) { return DustField(dust, p, 1.0, seed); });
}

double AssociationMean(Disk const& disk, uint64_t seed)
{
	if (disk.associationSize <= 0.0 || disk.associationClumping == 0.0)
	{
		return 1.0;
	}
	return SpatialMean([&](glm::dvec3 p) { return AssociationField(disk, p, 1.0, seed); });
}

// Needs the dust's filamentMean already measured.
double GasTextureMean(Galaxy const& galaxy, Disk const& disk)
{
	if (disk.gas.glow <= 0.0)
	{
		return 1.0;
	}
	return SpatialMean([&](glm::dvec3 p) { return GasTexture(galaxy, p, GasClumps(galaxy, p)); });
}

void NormaliseColours(Gas& gas)
{
	gas.colour /= Luminance(gas.colour);
	gas.coreColour /= Luminance(gas.coreColour);
	for (int i = 0; i < gas.paletteCount; ++i)
	{
		gas.palette[i] /= Luminance(gas.palette[i]);
	}
}

glm::dmat3 RotationX(double angle)
{
	double const c = std::cos(angle);
	double const s = std::sin(angle);
	return glm::dmat3(1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c);
}

glm::dmat3 RotationY(double angle)
{
	double const c = std::cos(angle);
	double const s = std::sin(angle);
	return glm::dmat3(c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c);
}

glm::dmat3 RotationZ(double angle)
{
	double const c = std::cos(angle);
	double const s = std::sin(angle);
	return glm::dmat3(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0);
}

}

double BoundRadius(Galaxy const& galaxy)
{
	double radius = 0.0;
	if (galaxy.disk.enabled)
	{
		radius = std::max(radius, galaxy.disk.cutoff);
	}
	if (galaxy.thickDisk.enabled)
	{
		radius = std::max(radius, galaxy.thickDisk.cutoff);
	}
	if (galaxy.bulge.enabled)
	{
		radius = std::max(radius, galaxy.bulge.cutoff);
	}
	if (galaxy.halo.enabled)
	{
		radius = std::max(radius, galaxy.halo.cutoff);
	}
	if (galaxy.dust.enabled)
	{
		radius = std::max(radius, galaxy.dust.cutoff);
	}
	return radius;
}

// Measured in noise space, so the mean holds at any wavelength: a box in parsecs would span a single
// noise cell at the largest wavelengths and measure nothing.
double LargeGasMean(Gas const& gas, uint64_t seed)
{
	if (gas.largeGlow <= 0.0)
	{
		return 1.0;
	}
	return SpatialMean([&](glm::dvec3 q) { return LargeGasField(q, 1.0, seed); });
}

PreparedGalaxy PrepareGalaxy(Galaxy const& galaxy)
{
	PreparedGalaxy prepared;
	prepared.galaxy = galaxy;
	prepared.galaxy.dust.filamentMean = DustFieldMean(galaxy.dust, galaxy.seed);
	prepared.galaxy.disk.associationMean = AssociationMean(galaxy.disk, galaxy.seed);
	prepared.galaxy.thickDisk.associationMean = AssociationMean(galaxy.thickDisk, galaxy.seed);
	prepared.galaxy.disk.gas.textureMean = GasTextureMean(prepared.galaxy, prepared.galaxy.disk);
	prepared.galaxy.disk.gas.largeMean = LargeGasMean(prepared.galaxy.disk.gas, prepared.galaxy.seed);
	NormaliseColours(prepared.galaxy.disk.gas);
	prepared.galaxy.dust.reflectionColour /= Luminance(prepared.galaxy.dust.reflectionColour);
	glm::dmat3 const toScene = RotationX(Radians(galaxy.rotationDeg.x)) * RotationY(Radians(galaxy.rotationDeg.y)) * RotationZ(Radians(galaxy.rotationDeg.z));
	prepared.toLocal = glm::transpose(toScene);
	prepared.position = galaxy.position;
	prepared.boundRadius = BoundRadius(galaxy);
	return prepared;
}

}
