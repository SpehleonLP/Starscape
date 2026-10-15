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

constexpr double kPi = 3.14159265358979323846;

// Salts so the dust filaments do not trace the same noise as the arm clumps.
constexpr uint64_t kDustNoiseSalt = 0xD057D057D057D057ull;

// Moments measured over 4·10⁵ random points. Fbm: σ 0.167 at 4 octaves, 0.159 at 6. Ridged fbm at 6 octaves.
constexpr double kFbmSigma = 0.163;
constexpr double kRidgedMean = 0.6307;
constexpr double kRidgedSigma = 0.1339;
constexpr int kDustOctaves = 6;
constexpr int kAssociationOctaves = 3;
constexpr uint64_t kAssociationNoiseSalt = 0xA550C1A7E5ull;

// The gas texture's noise. Fixed rather than scene keys: the owner-approved look was tuned on exactly these.
constexpr uint64_t kGasClumpSalt = 0x6A5ull;
constexpr uint64_t kGasWispSalt = 0x6A6ull;
constexpr double kGasClumpWavelength = 40.0;   // parsecs: a bright nebula's core
constexpr double kGasWispWavelength = 25.0;
constexpr int kGasOctaves = 5;
// The palette's hue field. Fixed rather than a scene key: the owner-approved nebula look was tuned on it.
constexpr uint64_t kPaletteSalt = 0x9A1E77Eull;
constexpr uint64_t kLargeGasSalt = 0xB16B16ull;
constexpr int kLargeGasOctaves = 4;
// The rim's patches and wisps. Fixed rather than scene keys: the owner-approved nebula look was tuned on them.
constexpr uint64_t kRimMaskSalt = 0x81B1ull;
constexpr uint64_t kRimWispSalt = 0x81B0ull;

// A log-normal field whose logarithm has standard deviation `sigma`, given z of unit variance. Holding
// the mean near one lets structure move stars and dust around without changing how much there is,
// so star counts and the dust's magnitudes-per-kiloparsec calibration survive any setting. For fbm the
// mean stays within 2% of one up to σ 1.5.
double LogNormal(double sigma, double z)
{
	return std::exp(sigma * z - 0.5 * sigma * sigma);
}

// Dust is filamentary, so its field rides ridged noise; that noise is bounded above, which pulls the
// log-normal's mean below one (0.65 at σ 2). Measured here rather than tabulated, so any σ stays calibrated.
double DustField(Dust const& dust, glm::dvec3 noiseSpace, uint64_t seed)
{
	return LogNormal(dust.filament, (RidgedFbm(noiseSpace, kDustOctaves, seed ^ kDustNoiseSalt) - kRidgedMean) / kRidgedSigma);
}

// Young stars are born in OB associations tens to a couple of hundred parsecs across, so the bright
// stars, which are mostly young, gather into groups like Orion's while the old field stays smooth.
double AssociationField(Disk const& disk, glm::dvec3 noiseSpace, uint64_t seed)
{
	return LogNormal(disk.associationClumping, Fbm(noiseSpace, kAssociationOctaves, seed ^ kAssociationNoiseSalt) / kFbmSigma);
}

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
	return SpatialMean([&](glm::dvec3 p) { return DustField(dust, p, seed); });
}

double AssociationMean(Disk const& disk, uint64_t seed)
{
	if (disk.associationSize <= 0.0 || disk.associationClumping == 0.0)
	{
		return 1.0;
	}
	return SpatialMean([&](glm::dvec3 p) { return AssociationField(disk, p, seed); });
}

double DustWavelength(Galaxy const& galaxy)
{
	return 0.5 * (galaxy.disk.arms.count > 0 ? galaxy.disk.arms.width : galaxy.dust.scaleHeight);
}

// Log-normal like the clumps: dust is filaments with clear holes between them, not a modulated sheet.
// Mean one over space: the dust's filaments with the extinction's normalisation.
double DustFilaments(Galaxy const& galaxy, glm::dvec3 p)
{
	return DustField(galaxy.dust, p / DustWavelength(galaxy), galaxy.seed) / galaxy.dust.filamentMean;
}

// The bubble's normalised radius: 1 on its ellipsoid. Only meaningful when it has a radius.
double BubbleRadius(Dust const& dust, glm::dvec3 p)
{
	glm::dvec3 const d = p - dust.bubbleCenter;
	return std::sqrt((d.x * d.x + d.y * d.y) / (dust.bubbleRadius * dust.bubbleRadius) + d.z * d.z / (dust.bubbleHeight * dust.bubbleHeight));
}

// Empty within 70% of the bubble's semi-axes, ramping smoothly to one at its edge; one without a bubble.
// Without the cavity, filaments within a few hundred parsecs of an observer in the plane blot the high
// latitudes, where the real sky's are faint cirrus.
double BubbleFactor(Dust const& dust, glm::dvec3 p)
{
	if (dust.bubbleRadius <= 0.0)
	{
		return 1.0;
	}
	double const t = std::clamp((BubbleRadius(dust, p) - 0.7) / 0.3, 0.0, 1.0);
	return t * t * (3.0 - 2.0 * t);
}

// Mean one: the soft clumps that also pick the gas colour.
double GasClumps(Galaxy const& galaxy, glm::dvec3 p)
{
	return LogNormal(1.2, Fbm(p / kGasClumpWavelength, kGasOctaves, galaxy.seed ^ kGasClumpSalt) / kFbmSigma);
}

// G before its mean is divided out: clumps plus thin bright wisps, gathered onto the dust's filaments
// because gas and dust are one medium.
double GasTexture(Galaxy const& galaxy, glm::dvec3 p, double clumps)
{
	double const wisps = 6.0 * std::pow(std::max(0.0, RidgedFbm(p / kGasWispWavelength, kGasOctaves, galaxy.seed ^ kGasWispSalt)), 6.0);
	double const filaments = galaxy.dust.enabled ? DustFilaments(galaxy, p) : 1.0;
	return (0.6 * clumps + wisps) * std::sqrt(filaments);
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

// A one-σ change in the hue field moves about one palette entry, so a sky shows every entry without
// flickering between them. Entries have luminance 1, so any blend of two does too.
glm::dvec3 PaletteColour(Galaxy const& galaxy, Gas const& gas, glm::dvec3 p)
{
	if (gas.paletteCount == 0)
	{
		return gas.colour;
	}
	double const count = double(gas.paletteCount);
	double const z = Fbm(p / gas.paletteWavelength, 3, galaxy.seed ^ kPaletteSalt) / kFbmSigma;
	double u = std::fmod(0.9 * z, count);
	if (u < 0.0)
	{
		u += count;
	}
	double const index = std::floor(u);
	double f = u - index;
	f = f * f * (3.0 - 2.0 * f);
	int const i = int(index) % gas.paletteCount;
	return glm::mix(gas.palette[size_t(i)], gas.palette[size_t((i + 1) % gas.paletteCount)], f);
}

void NormaliseColours(Gas& gas)
{
	gas.colour /= Luminance(gas.colour);
	gas.coreColour /= Luminance(gas.coreColour);
	for (int i = 0; i < gas.paletteCount; ++i)
	{
		gas.palette[size_t(i)] /= Luminance(gas.palette[size_t(i)]);
	}
}

double Radians(double degrees)
{
	return degrees * kPi / 180.0;
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

// 4e^{-2|x|} / (1 + e^{-2|x|})² equals 1/cosh²(x) without overflowing for large |x|.
double Sech2(double x)
{
	double const e = std::exp(-2.0 * std::abs(x));
	return 4.0 * e / ((1.0 + e) * (1.0 + e));
}

double ArmProximity(Arms const& arms, double radius, double azimuth, double offset)
{
	// ln(R / r0) diverges on the axis, where arms have no meaning anyway.
	if (arms.count <= 0 || radius < 1e-6)
	{
		return 0.0;
	}
	double const period = 2.0 * kPi / arms.count;
	double const phase = Radians(arms.phaseDeg) + std::log(radius / arms.r0) / std::tan(Radians(arms.pitchDeg));
	double delta = azimuth - phase;
	delta -= period * std::floor(delta / period + 0.5);
	// Arc length to the nearest arm. Positive delta is the inner edge: the arm through that
	// azimuth crosses it at a larger radius.
	double const distance = delta * radius - offset;
	return std::exp(-distance * distance / (2.0 * arms.width * arms.width));
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

double LargeGasField(glm::dvec3 noiseSpace, uint64_t seed)
{
	return LogNormal(1.0, Fbm(noiseSpace, kLargeGasOctaves, seed ^ kLargeGasSalt) / kFbmSigma);
}

// Measured in noise space, so the mean holds at any wavelength: a box in parsecs would span a single
// noise cell at the largest wavelengths and measure nothing.
double LargeGasMean(Gas const& gas, uint64_t seed)
{
	if (gas.largeGlow <= 0.0)
	{
		return 1.0;
	}
	return SpatialMean([&](glm::dvec3 q) { return LargeGasField(q, seed); });
}

// A shell just inside the dust's ramp, so most of it sits in front of the wall's dust. Not confined to the
// young disk's thickness: the bubble opens into chimneys, so the wall wraps the whole sky.
glm::dvec3 RimGlow(Galaxy const& galaxy, double radius, glm::dvec3 p)
{
	Dust const& dust = galaxy.dust;
	Disk const& disk = galaxy.disk;
	if (dust.bubbleRim <= 0.0 || !dust.enabled || dust.bubbleRadius <= 0.0 || !disk.enabled)
	{
		return glm::dvec3(0.0);
	}
	double const offset = BubbleRadius(dust, p) - 0.85;
	// Cut where the Gaussian is below e^-25, so the term costs nothing away from the wall.
	if (std::abs(offset) > 0.2)
	{
		return glm::dvec3(0.0);
	}
	// Only patches of the wall glow: the eye finds arcs where a lit patch's edge runs along the line of sight.
	double const mask = std::clamp(Fbm(p / 150.0, 3, galaxy.seed ^ kRimMaskSalt) / kFbmSigma - 0.5, 0.0, 1.0);
	if (mask <= 0.0)
	{
		return glm::dvec3(0.0);
	}
	double const shell = std::exp(-(offset / 0.04) * (offset / 0.04));
	double const wisps = 0.3 + 6.0 * std::pow(std::max(0.0, RidgedFbm(p / 30.0, 5, galaxy.seed ^ kRimWispSalt)), 6.0);
	// The young density the disk would have at this radius, without arms, associations or its thickness.
	double const reference = disk.density * std::exp(-radius / disk.scaleLength) * (disk.youngFraction + disk.armYoungFraction);
	double const light = dust.bubbleRim * disk.gas.youngLight * reference * shell * mask * wisps;
	// Four times the palette's frequency gives the rim its own, finer colour changes.
	return (light * 100.0 * kFluxAt10pc) * PaletteColour(galaxy, disk.gas, 4.0 * p);
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
	glm::dmat3 const toScene = RotationX(Radians(galaxy.rotationDeg.x)) * RotationY(Radians(galaxy.rotationDeg.y)) * RotationZ(Radians(galaxy.rotationDeg.z));
	prepared.toLocal = glm::transpose(toScene);
	prepared.position = galaxy.position;
	prepared.boundRadius = BoundRadius(galaxy);
	return prepared;
}

void EvaluateLocal(Galaxy const& galaxy, glm::dvec3 p, GalaxySample& out)
{
	out.count = 0;
	out.extinction = 0.0;
	out.gas = glm::dvec3(0.0);
	double const radius = std::hypot(p.x, p.y);
	double const r = glm::length(p);
	double const azimuth = std::atan2(p.y, p.x);

	auto push = [&out](double density, int population)
	{
		if (density > 0.0 && population >= 0)
		{
			out.parts[size_t(out.count++)] = {density, population};
		}
	};

	auto addDisk = [&](Disk const& disk)
	{
		if (!disk.enabled || r > disk.cutoff)
		{
			return;
		}
		double const crest = ArmProximity(disk.arms, radius, azimuth, 0.0);
		double modulation = 1.0;
		double clumps = 1.0;
		if (disk.arms.count > 0)
		{
			// Log-normal, so clumping opens real gaps between star clouds instead of a ripple on the crest.
			clumps = LogNormal(disk.arms.clumping, Fbm(p / disk.arms.width, 4, galaxy.seed) / kFbmSigma);
			modulation = std::max(0.0, 1.0 + (disk.arms.contrast - 1.0) * crest * clumps);
		}
		double const midplane = disk.density * std::exp(-radius / disk.scaleLength) * modulation;
		// Young stars concentrate on arms but are not confined to them: the Sun sits between arms,
		// yet half the brightest stars in Earth's sky are young B stars.
		double const young = disk.armPopulation >= 0 ? std::min(1.0, (disk.youngFraction + disk.armYoungFraction * crest) * clumps) : 0.0;
		push(midplane * (1.0 - young) * Sech2(p.z / (2.0 * disk.scaleHeight)), disk.population);
		// Young stars have not yet been heated out of the plane: sharing the old disk's thickness turns
		// each nearby arm into a wall whose tangent draws a hard edge from horizon to zenith.
		double const youngHeight = disk.armScaleHeight > 0.0 ? disk.armScaleHeight : disk.scaleHeight;
		// Associations add young stars where they are rather than taking a share of the old ones, so the
		// field multiplies the young density, not the share capped at one.
		double associations = 1.0;
		if (disk.associationSize > 0.0)
		{
			associations = AssociationField(disk, p / disk.associationSize, galaxy.seed) / disk.associationMean;
		}
		double const youngDensity = midplane * young * associations * Sech2(p.z / (2.0 * youngHeight));
		push(youngDensity, disk.armPopulation);
		if (youngDensity > 0.0 && (disk.gas.glow > 0.0 || disk.gas.largeGlow > 0.0))
		{
			glm::dvec3 const base = PaletteColour(galaxy, disk.gas, p);
			double const bubble = BubbleFactor(galaxy.dust, p);
			if (disk.gas.glow > 0.0)
			{
				double const clumps = GasClumps(galaxy, p);
				double const texture = GasTexture(galaxy, p, clumps) / disk.gas.textureMean;
				// Dense clumps turn from the base colour to the core colour; both have luminance 1, so the mix keeps the light.
				double const t = std::clamp((clumps - 1.0) / 2.0, 0.0, 1.0);
				glm::dvec3 const colour = glm::mix(base, disk.gas.coreColour, t * t * (3.0 - 2.0 * t));
				double const light = disk.gas.glow * disk.gas.youngLight * youngDensity * texture * bubble;
				out.gas += (light * 100.0 * kFluxAt10pc) * colour;
			}
			if (disk.gas.largeGlow > 0.0)
			{
				double const field = LargeGasField(p / disk.gas.largeWavelength, galaxy.seed) / disk.gas.largeMean;
				double const light = disk.gas.largeGlow * disk.gas.youngLight * youngDensity * field * bubble;
				out.gas += (light * 100.0 * kFluxAt10pc) * base;
			}
		}
	};
	addDisk(galaxy.disk);
	addDisk(galaxy.thickDisk);
	if (r <= galaxy.disk.cutoff)
	{
		out.gas += RimGlow(galaxy, radius, p);
	}

	if (galaxy.bulge.enabled && r <= galaxy.bulge.cutoff)
	{
		Bulge const& bulge = galaxy.bulge;
		double const c = std::cos(Radians(bulge.barAngleDeg));
		double const s = std::sin(Radians(bulge.barAngleDeg));
		double const bx = c * p.x + s * p.y;
		double const by = (-s * p.x + c * p.y) / bulge.axisY;
		double const bz = p.z / bulge.axisZ;
		// The Hernquist cusp diverges at r' = 0; a 1 pc floor keeps an observer at the core finite.
		double const rb = std::max(1.0, std::sqrt(bx * bx + by * by + bz * bz));
		double const a = bulge.scale;
		push(bulge.density * a / (rb * (rb + a) * (rb + a) * (rb + a)), bulge.population);
	}

	if (galaxy.halo.enabled && r <= galaxy.halo.cutoff)
	{
		push(galaxy.halo.density * std::pow(1.0 + r / galaxy.halo.scale, -3.5), galaxy.halo.population);
	}

	if (galaxy.dust.enabled && r <= galaxy.dust.cutoff)
	{
		Dust const& dust = galaxy.dust;
		Arms const& arms = galaxy.disk.arms;
		double const lane = ArmProximity(arms, radius, azimuth, dust.laneOffset * arms.width);
		double const laneModulation = 1.0 + (dust.contrast - 1.0) * lane;
		double const filaments = DustFilaments(galaxy, p);
		out.extinction = dust.extinction * std::exp(-radius / dust.scaleLength) * Sech2(p.z / (2.0 * dust.scaleHeight)) * laneModulation * filaments
			* BubbleFactor(dust, p);
	}
}

}
