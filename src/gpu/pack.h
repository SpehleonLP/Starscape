#pragma once
#include "field/step.glsl.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace starscape::gpu
{

// Writes values in GLSL's scalar block layout: each scalar aligned to its own size, a vector to its
// component's, a struct to its largest member's, arrays of structs at a stride rounded to that alignment.
class ScalarWriter
{
public:
	void Float(double value);
	void Int(int32_t value);
	void Uint(uint32_t value);
	void Bool(bool value);   // 4 bytes, as GLSL stores bool in a buffer
	void U64(uint64_t value);
	void Vec3(glm::dvec3 value);
	void Mat3(glm::dmat3 value);   // column-major, three tight vec3 columns
	void Align(size_t alignment);
	std::vector<uint8_t> const& Bytes() const { return m_bytes; }

private:
	void Put(void const* data, size_t size);

	std::vector<uint8_t> m_bytes;
};

// The GLSL struct Galaxy, member by member. Ends padded to kGalaxyAlignment, so consecutive galaxies pack at
// the array stride.
void PackGalaxy(ScalarWriter& w, Galaxy const& galaxy);
void PackPopulation(ScalarWriter& w, PopulationTable const& table);
inline constexpr size_t kGalaxyAlignment = 8;   // Galaxy holds a uint64_t seed

}
