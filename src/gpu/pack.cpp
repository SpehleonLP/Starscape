#include "gpu/pack.h"

namespace starscape::gpu
{
namespace
{

void PackArms(ScalarWriter& w, Arms const& arms)
{
	w.Int(arms.count);
	w.Float(arms.pitchDeg);
	w.Float(arms.width);
	w.Float(arms.contrast);
	w.Float(arms.clumping);
	w.Float(arms.phaseDeg);
	w.Float(arms.r0);
}

void PackGas(ScalarWriter& w, Gas const& gas)
{
	w.Float(gas.glow);
	w.Vec3(gas.colour);
	w.Vec3(gas.coreColour);
	// The whole array whatever paletteCount says: the GLSL array has a fixed length.
	for (glm::dvec3 const& colour : gas.palette)
	{
		w.Vec3(colour);
	}
	w.Int(gas.paletteCount);
	w.Float(gas.paletteWavelength);
	w.Float(gas.largeGlow);
	w.Float(gas.largeWavelength);
	w.Float(gas.largeMean);
	w.Float(gas.textureMean);
	w.Float(gas.youngLight);
}

void PackDisk(ScalarWriter& w, Disk const& disk)
{
	w.Bool(disk.enabled);
	w.Float(disk.density);
	w.Float(disk.scaleLength);
	w.Float(disk.scaleHeight);
	w.Float(disk.cutoff);
	PackArms(w, disk.arms);
	w.Int(disk.population);
	w.Int(disk.armPopulation);
	w.Float(disk.armYoungFraction);
	w.Float(disk.armScaleHeight);
	w.Float(disk.youngFraction);
	w.Float(disk.associationSize);
	w.Float(disk.associationClumping);
	w.Float(disk.associationMean);
	PackGas(w, disk.gas);
}

void PackBulge(ScalarWriter& w, Bulge const& bulge)
{
	w.Bool(bulge.enabled);
	w.Float(bulge.density);
	w.Float(bulge.scale);
	w.Float(bulge.axisY);
	w.Float(bulge.axisZ);
	w.Float(bulge.barAngleDeg);
	w.Float(bulge.cutoff);
	w.Int(bulge.population);
}

void PackHalo(ScalarWriter& w, Halo const& halo)
{
	w.Bool(halo.enabled);
	w.Float(halo.density);
	w.Float(halo.scale);
	w.Float(halo.cutoff);
	w.Int(halo.population);
}

void PackDust(ScalarWriter& w, Dust const& dust)
{
	w.Bool(dust.enabled);
	w.Float(dust.extinction);
	w.Float(dust.scaleLength);
	w.Float(dust.scaleHeight);
	w.Float(dust.cutoff);
	w.Float(dust.laneOffset);
	w.Float(dust.contrast);
	w.Float(dust.filament);
	w.Vec3(dust.bubbleCenter);
	w.Float(dust.bubbleRadius);
	w.Float(dust.bubbleHeight);
	w.Float(dust.filamentMean);
	w.Float(dust.bubbleRim);
	w.Float(dust.reflection);
	w.Vec3(dust.reflectionColour);
}

}

void ScalarWriter::Put(void const* data, size_t size)
{
	Align(size);
	uint8_t const* bytes = static_cast<uint8_t const*>(data);
	m_bytes.insert(m_bytes.end(), bytes, bytes + size);
}

void ScalarWriter::Align(size_t alignment)
{
	m_bytes.resize((m_bytes.size() + alignment - 1) / alignment * alignment, 0);
}

void ScalarWriter::Float(double value)
{
	float const narrowed = float(value);
	Put(&narrowed, sizeof narrowed);
}

void ScalarWriter::Int(int32_t value)
{
	Put(&value, sizeof value);
}

void ScalarWriter::Uint(uint32_t value)
{
	Put(&value, sizeof value);
}

void ScalarWriter::Bool(bool value)
{
	Uint(value ? 1u : 0u);
}

void ScalarWriter::U64(uint64_t value)
{
	Put(&value, sizeof value);
}

void ScalarWriter::Vec3(glm::dvec3 value)
{
	Float(value.x);
	Float(value.y);
	Float(value.z);
}

void ScalarWriter::Mat3(glm::dmat3 value)
{
	for (int column = 0; column < 3; ++column)
	{
		Vec3(value[column]);
	}
}

void PackGalaxy(ScalarWriter& w, Galaxy const& galaxy)
{
	w.Align(kGalaxyAlignment);
	w.Vec3(galaxy.position);
	w.Vec3(galaxy.rotationDeg);
	w.U64(galaxy.seed);
	PackDisk(w, galaxy.disk);
	PackDisk(w, galaxy.thickDisk);
	PackBulge(w, galaxy.bulge);
	PackHalo(w, galaxy.halo);
	PackDust(w, galaxy.dust);
	w.Align(kGalaxyAlignment);
}

void PackPopulation(ScalarWriter& w, PopulationTable const& table)
{
	for (glm::dvec3 const& chroma : table.chroma)
	{
		w.Vec3(chroma);
	}
	for (double const fraction : table.fraction)
	{
		w.Float(fraction);
	}
	for (double const luminosity : table.luminosity)
	{
		w.Float(luminosity);
	}
	for (glm::dvec3 const& light : table.lightBelow)
	{
		w.Vec3(light);
	}
}

}
