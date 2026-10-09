#include "galaxy.h"
#include "noise.h"
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

double DustFieldMean(Dust const& dust, uint64_t seed)
{
	if (dust.filament == 0.0)
	{
		return 1.0;
	}
	// An R3 low-discrepancy walk over a box many noise cells wide, so the samples are independent.
	constexpr int kSamples = 1 << 16;
	double sum = 0.0;
	for (int i = 1; i <= kSamples; ++i)
	{
		glm::dvec3 const u = glm::fract(double(i) * glm::dvec3(0.8191725134, 0.6710436067, 0.5497004779));
		sum += DustField(dust, u * 4096.0, seed);
	}
	return sum / kSamples;
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

PreparedGalaxy PrepareGalaxy(Galaxy const& galaxy)
{
	PreparedGalaxy prepared;
	prepared.galaxy = galaxy;
	prepared.galaxy.dust.filamentMean = DustFieldMean(galaxy.dust, galaxy.seed);
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
		push(midplane * young * Sech2(p.z / (2.0 * youngHeight)), disk.armPopulation);
	};
	addDisk(galaxy.disk);
	addDisk(galaxy.thickDisk);

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
		double const wavelength = 0.5 * (arms.count > 0 ? arms.width : dust.scaleHeight);
		// Log-normal like the clumps: dust is filaments with clear holes between them, not a modulated sheet.
		double const filaments = DustField(dust, p / wavelength, galaxy.seed) / dust.filamentMean;
		out.extinction = dust.extinction * std::exp(-radius / dust.scaleLength) * Sech2(p.z / (2.0 * dust.scaleHeight)) * laneModulation * filaments;
	}
}

}
