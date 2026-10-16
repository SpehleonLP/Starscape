// Shared between C++ (double precision, glm) and GLSL (float). Everything in src/field/ is written once and
// compiled by both, so the GPU march and the host's calibration and star re-march evaluate the same model. Code here obeys both
// languages: no lambdas, templates, std:: calls, auto or recursion; `const real x`, never `real const x`;
// and no name GLSL reserves or uses as a built-in (out, sample, step, distance, length, texture, ...).
#ifndef STARSCAPE_FIELD_COMMON_GLSL_H
#define STARSCAPE_FIELD_COMMON_GLSL_H

#ifdef __cplusplus

#include <array>
#include <cmath>
#include <cstdint>
#include <glm/common.hpp>
#include <glm/geometric.hpp>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#define IN(T) T const&
#define INOUT(T) T&
#define ARRAY(T, name, N) std::array<T, N> name
#define ZEROED {}
#define DEFAULT(value) = value
#define FIELD_INLINE inline
// Namespace-scope constants are read by inline functions in many TUs: inline constexpr gives them one definition.
#define FIELD_CONST inline constexpr

// 64-bit constants are written with a `ul` suffix, the one spelling both languages accept: glslang lexes a
// macro argument before pasting, so an unsuffixed 64-bit hex literal fails there, and `ull` is not GLSL.
static_assert(sizeof(unsigned long) == 8, "shared field code spells 64-bit constants with a ul suffix");

namespace starscape
{

using real = double;
using real3 = glm::dvec3;
using real3x3 = glm::dmat3;
using uint = uint32_t;

// GLSL's built-ins under their GLSL names, so shared code reads the same in both languages.
using std::abs;
using std::cos;
using std::exp;
using std::floor;
using std::log;
using std::pow;
using std::sin;
using std::sqrt;
using std::tan;
using glm::clamp;
using glm::dot;
using glm::length;
using glm::max;
using glm::min;
using glm::mix;

inline double atan(double y, double x)
{
	return std::atan2(y, x);
}

// C's fmod, exact in double on the host; GLSL spells it out below.
inline double FMod(double x, double y)
{
	return std::fmod(x, y);
}

// The disc radius of (x, y), exact in double on the host; GLSL has no hypot.
inline double Hypot(double x, double y)
{
	return std::hypot(x, y);
}

inline constexpr double kHugeReal = HUGE_VAL;

}

#else

#define IN(T) in T
#define INOUT(T) inout T
#define ARRAY(T, name, N) T name[N]
#define ZEROED
#define DEFAULT(value)
#define FIELD_INLINE
#define FIELD_CONST const
#define real float
#define real3 vec3
#define real3x3 mat3

real FMod(real x, real y)
{
	return x - y * trunc(x / y);
}

real Hypot(real x, real y)
{
	return sqrt(x * x + y * y);
}

// Large enough that no luminosity reaches it; GLSL has no portable constant infinity.
FIELD_CONST real kHugeReal = 3.0e38;

#endif

#ifdef __cplusplus
namespace starscape
{
#endif

// Moments measured over 4·10⁵ random points. Fbm: σ 0.167 at 4 octaves, 0.159 at 6.
FIELD_CONST real kFbmSigma = 0.163;

FIELD_INLINE real Luminance(real3 rgb)
{
	return 0.2126 * rgb.x + 0.7152 * rgb.y + 0.0722 * rgb.z;
}

// 4e^{-2|x|} / (1 + e^{-2|x|})² equals 1/cosh²(x) without overflowing for large |x|.
FIELD_INLINE real Sech2(real x)
{
	const real e = exp(-2.0 * abs(x));
	return 4.0 * e / ((1.0 + e) * (1.0 + e));
}

// A log-normal field whose logarithm has standard deviation `sigma`, given z of unit variance. Holding
// the mean near one lets structure move stars and dust around without changing how much there is,
// so star counts and the dust's magnitudes-per-kiloparsec calibration survive any setting. For fbm the
// mean stays within 2% of one up to σ 1.5.
FIELD_INLINE real LogNormal(real sigma, real z)
{
	return exp(sigma * z - 0.5 * sigma * sigma);
}

#ifdef __cplusplus
}
#endif

#endif
