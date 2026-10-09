#include "galaxy.h"
#include "noise.h"
#include <algorithm>
#include <cmath>
#include <glm/geometric.hpp>
#include <glm/matrix.hpp>

namespace starscape
{
namespace
{

constexpr double kPi = 3.14159265358979323846;

// Salts so the dust filaments do not trace the same noise as the arm clumps.
constexpr uint64_t kDustNoiseSalt = 0xD057D057D057D057ull;

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
		if (disk.arms.count > 0)
		{
			double const clumps = 1.0 + disk.arms.clumping * Fbm(p / disk.arms.width, 4, galaxy.seed);
			modulation = std::max(0.0, 1.0 + (disk.arms.contrast - 1.0) * crest * clumps);
		}
		double const n = disk.density * std::exp(-radius / disk.scaleLength) * Sech2(p.z / (2.0 * disk.scaleHeight)) * modulation;
		double const young = disk.armPopulation >= 0 ? disk.armYoungFraction * crest : 0.0;
		push(n * (1.0 - young), disk.population);
		push(n * young, disk.armPopulation);
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
		double const filaments = std::max(0.0, 1.0 + dust.filament * Fbm(p / wavelength, 6, galaxy.seed ^ kDustNoiseSalt));
		out.extinction = dust.extinction * std::exp(-radius / dust.scaleLength) * Sech2(p.z / (2.0 * dust.scaleHeight)) * laneModulation * filaments;
	}
}

}
