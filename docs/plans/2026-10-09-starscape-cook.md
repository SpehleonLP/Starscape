# Starscape Cook Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `starscape`, an offline C++ CLI that cooks a `.kre_starscape` night sky (RGB9E5 diffuse cubemap + point stars) and an equirectangular `.hdr` preview from a TOML scene of procedural galaxies.

**Architecture:** A static library `starscape-core` (scene parsing, galaxy model, ray march, star drawing, file writer, preview) plus a thin CLI and a doctest test binary. One ray march per cubemap texel accumulates both the diffuse radiance and a per-luminosity-bin point budget, split by a single apparent-flux threshold; stars are Poisson-drawn from the budget and placed by re-marching. The file format lives in a std-only public header `include/kre_starscape.h` that the engine will include unchanged.

**Tech Stack:** C++23 (GCC 13), CMake ≥ 3.25 + Ninja, vendored glm 1.0.1, toml++ 3.4.0, doctest 2.4.11.

**Spec:** `docs/spec.md` (read it first; every formula here is from it).

## Global Constraints

- Repository: `/mnt/Passport/Worktrees/starscape` (a worktree of the Starscape repo on branch `starscape-v1`; the main checkout at `/mnt/Passport/Kreatures-Engine/Tooling/Starscape` stays on `main` and is not touched until Task 11). All paths below are relative to it unless absolute.
- No dependency on `Spehleon/lib` or anything private. Third-party code is vendored under `third_party/` with versions in `third_party/VERSIONS`.
- C++23, `-Wall -Wextra -Wpedantic` clean. Tabs, Allman braces, **every** control-flow body braced, even one statement.
- Containers, strings, `std::function` and structs larger than a couple of registers are passed `const&` (read) or `&&` (consume), never by value. `glm` vectors, spans and enums go by value.
- Errors are returned as values (`Expected<T>`), printed once in `main`. No exceptions, no `abort`, no asserts on input.
- Comments say *why*, not *what*.
- Determinism: every random draw comes from `RandomStream` keyed on (seed, purpose, a, b, c); results are written to per-index slots and reduced in index order. Output is byte-identical on any thread count.
- Units: parsecs; luminosity in L_V☉; flux in magnitude-0 units; `kFluxAt10pc = 0.01169499` (= 10^(−0.4·4.83)).
- Build (outside the engine's jobserver wrapper, which covers engine trees only — keep it small and low-priority):
  - configure: `cmake -S /mnt/Passport/Worktrees/starscape -B /mnt/Passport/Worktrees/starscape/Build -G Ninja`
  - build: `nice -n 19 ionice -c 3 cmake --build /mnt/Passport/Worktrees/starscape/Build -j 4`
  - one test: `/mnt/Passport/Worktrees/starscape/Build/starscape-tests --test-case='<name>'`
  - all: `ctest --test-dir /mnt/Passport/Worktrees/starscape/Build --output-on-failure`
- Every test that names a mutation must be shown red under that mutation once, then the mutation reverted (`touch` the file afterwards so Ninja rebuilds it). Record "mutation verified: <what>" in the commit message.
- Never weaken a failing physics test to make it pass. If a model-level assertion (inside-disk ratios, energy) fails, stop and report the numbers.

## Review Focus

1. Observer inside a galaxy core (bulge cusp, density spikes near `s_min`): the cook must finish with finite values and roughly the requested star count. → Task 7 test "observer at the core".
2. Star budget larger than the scene holds: the cook must finish, flag the budget unreachable, and make every star a point. → Task 7 test "unreachable budget".
3. Observer on a galaxy's axis or at its centre (R = 0 in the arm phase `ln(R/R0)`): densities must be finite. → Task 4 test "evaluation at the centre is finite".
4. Rays that start inside, graze, or cross several bounding spheres: intervals merged, no stretch marched twice. → Task 6 test "ray intervals".
5. An output path whose final rename fails: I/O error, no `.partial` file left behind. → Task 8 test "failed rename leaves nothing".

---

### Task 1: Scaffold, Philox RNG and gradient noise

**Files:**
- Create: `.gitignore`, `CMakeLists.txt`, `third_party/VERSIONS`, `third_party/glm/…`, `third_party/toml.hpp`, `third_party/doctest.h`
- Create: `src/common.h`, `src/rng.h`, `src/noise.h`, `src/noise.cpp`
- Test: `tests/main.cpp`, `tests/test_rng.cpp`, `tests/test_noise.cpp`

**Interfaces:**
- Produces: `starscape::Error{Kind kind; std::string message}`, `Expected<T>`, `InputError(std::string&&)`, `IoError(std::string&&)`, macros `STARSCAPE_TRY(expr)`, `STARSCAPE_TRY_ASSIGN(name, expr)`; `Philox4x32_10(PhiloxCounter, PhiloxKey) -> PhiloxCounter`; `enum class Purpose { Jitter = 1, Poisson = 2, Star = 3, Test = 4 }`; `RandomStream(uint64_t seed, Purpose, uint32_t a, uint32_t b, uint32_t c)` with `NextU32()`, `NextDouble()`; `double GradientNoise(glm::dvec3, uint64_t seed)`, `double Fbm(glm::dvec3, int octaves, uint64_t seed)`.

- [ ] **Step 1: Vendor third-party headers**

Download into a fresh scratch directory (untrusted data), then copy the headers in:

```bash
mkdir -p /tmp/claude-1000/starscape-vendor && cd /tmp/claude-1000/starscape-vendor && curl -fsSLo toml.hpp https://raw.githubusercontent.com/marzer/tomlplusplus/v3.4.0/toml.hpp && curl -fsSLo doctest.h https://raw.githubusercontent.com/doctest/doctest/v2.4.11/doctest/doctest.h && curl -fsSLo glm.zip https://github.com/g-truc/glm/releases/download/1.0.1/glm-1.0.1-light.zip && curl -fsSLo glm-copying.txt https://raw.githubusercontent.com/g-truc/glm/1.0.1/copying.txt && unzip -q glm.zip -d glm-extract && sha256sum toml.hpp doctest.h glm.zip glm-copying.txt
```

Copy `toml.hpp` and `doctest.h` to `third_party/`; copy the extracted `glm` directory so that `third_party/glm/glm.hpp` exists; copy `glm-copying.txt` to `third_party/glm/copying.txt`. Write `third_party/VERSIONS` listing each file, its URL and the sha256 printed above:

```
toml.hpp      tomlplusplus v3.4.0   https://raw.githubusercontent.com/marzer/tomlplusplus/v3.4.0/toml.hpp   sha256 <hash>
doctest.h     doctest v2.4.11       https://raw.githubusercontent.com/doctest/doctest/v2.4.11/doctest/doctest.h   sha256 <hash>
glm/          glm 1.0.1 (light)     https://github.com/g-truc/glm/releases/download/1.0.1/glm-1.0.1-light.zip   sha256 <hash>
```

- [ ] **Step 2: Write `.gitignore` and `CMakeLists.txt`**

`.gitignore`:
```
Build/
```

`CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.25)
project(Starscape VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

# Every source except the CLI entry point; the tests link the same library.
set(STARSCAPE_SOURCES
	src/noise.cpp
)

add_library(starscape-core STATIC ${STARSCAPE_SOURCES})
target_include_directories(starscape-core PUBLIC include src)
target_include_directories(starscape-core SYSTEM PUBLIC third_party)
target_compile_definitions(starscape-core PUBLIC TOML_EXCEPTIONS=0)
if(NOT MSVC)
	target_compile_options(starscape-core PUBLIC -Wall -Wextra -Wpedantic)
endif()
target_link_libraries(starscape-core PUBLIC Threads::Threads)

enable_testing()
set(STARSCAPE_TESTS
	tests/main.cpp
	tests/test_rng.cpp
	tests/test_noise.cpp
)
add_executable(starscape-tests ${STARSCAPE_TESTS})
target_link_libraries(starscape-tests PRIVATE starscape-core)
add_test(NAME unit COMMAND starscape-tests)
```

Later tasks add their `.cpp` files to `STARSCAPE_SOURCES` and their test files to `STARSCAPE_TESTS`.

- [ ] **Step 3: Write `src/common.h`**

```cpp
#pragma once
#include <expected>
#include <string>
#include <utility>

namespace starscape
{

struct Error
{
	// Input = bad scene or file contents (exit 1); Io = the filesystem refused (exit 2).
	enum class Kind { Input, Io };
	Kind        kind = Kind::Input;
	std::string message;
};

template<class T>
using Expected = std::expected<T, Error>;

inline std::unexpected<Error> InputError(std::string&& message)
{
	return std::unexpected(Error{Error::Kind::Input, std::move(message)});
}

inline std::unexpected<Error> IoError(std::string&& message)
{
	return std::unexpected(Error{Error::Kind::Io, std::move(message)});
}

}

#define STARSCAPE_TRY(expr) \
	do \
	{ \
		if (auto starscapeTry_ = (expr); !starscapeTry_) \
		{ \
			return std::unexpected(std::move(starscapeTry_.error())); \
		} \
	} while (0)

#define STARSCAPE_TRY_ASSIGN(name, expr) \
	auto name##Result_ = (expr); \
	if (!name##Result_) \
	{ \
		return std::unexpected(std::move(name##Result_.error())); \
	} \
	auto name = std::move(*name##Result_)
```

- [ ] **Step 4: Write the failing tests**

`tests/main.cpp`:
```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
```

`tests/test_rng.cpp`:
```cpp
#include "rng.h"
#include <doctest.h>

using namespace starscape;

TEST_CASE("Philox4x32-10 matches the Random123 known-answer vectors")
{
	CHECK(Philox4x32_10({0u, 0u, 0u, 0u}, {0u, 0u}) == PhiloxCounter{0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu, 0x9b00dbd8u});
	CHECK(Philox4x32_10({~0u, ~0u, ~0u, ~0u}, {~0u, ~0u}) == PhiloxCounter{0x408f276du, 0x41c83b0eu, 0xa20bc7c6u, 0x6d5451fdu});
	CHECK(Philox4x32_10({0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u}, {0xa4093822u, 0x299f31d0u})
		== PhiloxCounter{0xd16cfe09u, 0x94fdccebu, 0x5001e420u, 0x24126ea1u});
}

TEST_CASE("RandomStream is a pure function of its key")
{
	RandomStream a(7, Purpose::Star, 1, 2, 3);
	RandomStream b(7, Purpose::Star, 1, 2, 3);
	for (int i = 0; i < 100; ++i)
	{
		CHECK(a.NextU32() == b.NextU32());
	}
	uint32_t const reference = RandomStream(7, Purpose::Star, 1, 2, 3).NextU32();
	CHECK(RandomStream(7, Purpose::Star, 1, 2, 4).NextU32() != reference);
	CHECK(RandomStream(7, Purpose::Poisson, 1, 2, 3).NextU32() != reference);
	CHECK(RandomStream(8, Purpose::Star, 1, 2, 3).NextU32() != reference);
}

TEST_CASE("RandomStream doubles are uniform in [0, 1)")
{
	RandomStream rng(1, Purpose::Test, 0, 0, 0);
	constexpr int n = 200000;
	double sum = 0.0;
	double sumSquares = 0.0;
	for (int i = 0; i < n; ++i)
	{
		double const u = rng.NextDouble();
		REQUIRE(u >= 0.0);
		REQUIRE(u < 1.0);
		sum += u;
		sumSquares += u * u;
	}
	CHECK(sum / n == doctest::Approx(0.5).epsilon(0.01));
	CHECK(sumSquares / n == doctest::Approx(1.0 / 3.0).epsilon(0.01));
}
```

`tests/test_noise.cpp`:
```cpp
#include "noise.h"
#include "rng.h"
#include <doctest.h>
#include <algorithm>
#include <cmath>

using namespace starscape;

TEST_CASE("Gradient noise is zero on the lattice and deterministic")
{
	CHECK(GradientNoise({3.0, -2.0, 7.0}, 5) == 0.0);
	glm::dvec3 const p(0.3, 1.7, -2.2);
	CHECK(GradientNoise(p, 5) == GradientNoise(p, 5));
	CHECK(GradientNoise(p, 5) != GradientNoise(p, 6));
}

TEST_CASE("Fbm stays bounded with a near-zero mean")
{
	RandomStream rng(3, Purpose::Test, 0, 0, 0);
	constexpr int n = 100000;
	double sum = 0.0;
	double extreme = 0.0;
	for (int i = 0; i < n; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1000.0 - 500.0, rng.NextDouble() * 1000.0 - 500.0, rng.NextDouble() * 1000.0 - 500.0);
		double const v = Fbm(p, 6, 11);
		sum += v;
		extreme = std::max(extreme, std::abs(v));
	}
	CHECK(extreme <= 1.1);
	CHECK(extreme > 0.3);
	CHECK(std::abs(sum / n) < 0.02);
}

// The earlier experiment hashed with fract(sin(dot(p, k)) * 43758) in single precision,
// which collapses at coordinates like 250 000. Noise here must still vary there.
TEST_CASE("Noise keeps its variance at galactic coordinates")
{
	RandomStream rng(4, Purpose::Test, 0, 0, 0);
	constexpr int n = 2000;
	double sum = 0.0;
	double sumSquares = 0.0;
	for (int i = 0; i < n; ++i)
	{
		glm::dvec3 const p(250000.0 + rng.NextDouble() * 10.0, 1.0e6 + rng.NextDouble() * 10.0, rng.NextDouble() * 10.0);
		double const v = Fbm(p, 4, 9);
		sum += v;
		sumSquares += v * v;
	}
	double const mean = sum / n;
	CHECK(sumSquares / n - mean * mean > 0.01);
}
```

- [ ] **Step 5: Configure and build to verify the tests fail**

Run the configure and build commands from Global Constraints.
Expected: compile errors (`rng.h`, `noise.h` not found).

- [ ] **Step 6: Write `src/rng.h`**

```cpp
#pragma once
#include <array>
#include <cstdint>

namespace starscape
{

using PhiloxCounter = std::array<uint32_t, 4>;
using PhiloxKey     = std::array<uint32_t, 2>;

// Philox4x32-10 (Salmon et al., SC11): a counter-based generator, so a draw is a pure
// function of its key and counter and never of which thread reached it first.
constexpr PhiloxCounter Philox4x32_10(PhiloxCounter counter, PhiloxKey key)
{
	constexpr uint32_t kMultiplier0 = 0xD2511F53u;
	constexpr uint32_t kMultiplier1 = 0xCD9E8D57u;
	constexpr uint32_t kWeyl0       = 0x9E3779B9u;
	constexpr uint32_t kWeyl1       = 0xBB67AE85u;
	for (int round = 0; round < 10; ++round)
	{
		uint64_t const product0 = uint64_t(kMultiplier0) * counter[0];
		uint64_t const product1 = uint64_t(kMultiplier1) * counter[2];
		counter = {
			uint32_t(product1 >> 32) ^ counter[1] ^ key[0],
			uint32_t(product1),
			uint32_t(product0 >> 32) ^ counter[3] ^ key[1],
			uint32_t(product0)};
		key[0] += kWeyl0;
		key[1] += kWeyl1;
	}
	return counter;
}

enum class Purpose : uint32_t
{
	Jitter  = 1,
	Poisson = 2,
	Star    = 3,
	Test    = 4,
};

class RandomStream
{
public:
	// `c` must fit in 24 bits: the purpose occupies the top byte of its counter word.
	RandomStream(uint64_t seed, Purpose purpose, uint32_t a, uint32_t b, uint32_t c)
	: m_key{uint32_t(seed), uint32_t(seed >> 32)}
	, m_counter{a, b, (uint32_t(purpose) << 24) | (c & 0x00FFFFFFu), 0u}
	{
	}

	uint32_t NextU32()
	{
		if (m_used == 4)
		{
			m_block = Philox4x32_10(m_counter, m_key);
			++m_counter[3];
			m_used = 0;
		}
		return m_block[size_t(m_used++)];
	}

	// Uniform in [0, 1) with 53 random bits.
	double NextDouble()
	{
		uint64_t const high = NextU32() >> 5;
		uint64_t const low  = NextU32() >> 6;
		return double(high * 67108864u + low) * (1.0 / 9007199254740992.0);
	}

private:
	PhiloxKey     m_key;
	PhiloxCounter m_counter;
	PhiloxCounter m_block{};
	int           m_used = 4;
};

}
```

- [ ] **Step 7: Write `src/noise.h` and `src/noise.cpp`**

`src/noise.h`:
```cpp
#pragma once
#include <cstdint>
#include <glm/vec3.hpp>

namespace starscape
{

// Perlin gradient noise: roughly in [-1, 1], exactly zero on integer lattice points.
double GradientNoise(glm::dvec3 p, uint64_t seed);

// `octaves` octaves of GradientNoise at doubling frequency and halving amplitude,
// divided by the amplitude sum so the result stays roughly in [-1, 1].
double Fbm(glm::dvec3 p, int octaves, uint64_t seed);

}
```

`src/noise.cpp`:
```cpp
#include "noise.h"
#include <cmath>

namespace starscape
{
namespace
{

// 64-bit integer mixing: unlike a sin() hash it is exact at any lattice coordinate.
uint32_t LatticeHash(int64_t x, int64_t y, int64_t z, uint64_t seed)
{
	uint64_t h = seed * 0x9E3779B97F4A7C15ull + 0xD1B54A32D192ED03ull;
	for (uint64_t const v : {uint64_t(x), uint64_t(y), uint64_t(z)})
	{
		h ^= v + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
		h *= 0xBF58476D1CE4E5B9ull;
		h ^= h >> 31;
	}
	return uint32_t(h >> 32);
}

constexpr double kGradients[12][3] = {
	{1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
	{1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
	{0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}};

double Fade(double t)
{
	return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

double Lerp(double a, double b, double t)
{
	return a + (b - a) * t;
}

double Corner(int64_t ix, int64_t iy, int64_t iz, double fx, double fy, double fz, uint64_t seed)
{
	double const* g = kGradients[LatticeHash(ix, iy, iz, seed) % 12u];
	return g[0] * fx + g[1] * fy + g[2] * fz;
}

}

double GradientNoise(glm::dvec3 p, uint64_t seed)
{
	double const x0 = std::floor(p.x);
	double const y0 = std::floor(p.y);
	double const z0 = std::floor(p.z);
	int64_t const ix = int64_t(x0);
	int64_t const iy = int64_t(y0);
	int64_t const iz = int64_t(z0);
	double const fx = p.x - x0;
	double const fy = p.y - y0;
	double const fz = p.z - z0;
	double const u = Fade(fx);
	double const v = Fade(fy);
	double const w = Fade(fz);

	double const c000 = Corner(ix, iy, iz, fx, fy, fz, seed);
	double const c100 = Corner(ix + 1, iy, iz, fx - 1.0, fy, fz, seed);
	double const c010 = Corner(ix, iy + 1, iz, fx, fy - 1.0, fz, seed);
	double const c110 = Corner(ix + 1, iy + 1, iz, fx - 1.0, fy - 1.0, fz, seed);
	double const c001 = Corner(ix, iy, iz + 1, fx, fy, fz - 1.0, seed);
	double const c101 = Corner(ix + 1, iy, iz + 1, fx - 1.0, fy, fz - 1.0, seed);
	double const c011 = Corner(ix, iy + 1, iz + 1, fx, fy - 1.0, fz - 1.0, seed);
	double const c111 = Corner(ix + 1, iy + 1, iz + 1, fx - 1.0, fy - 1.0, fz - 1.0, seed);

	return Lerp(
		Lerp(Lerp(c000, c100, u), Lerp(c010, c110, u), v),
		Lerp(Lerp(c001, c101, u), Lerp(c011, c111, u), v),
		w);
}

double Fbm(glm::dvec3 p, int octaves, uint64_t seed)
{
	double sum = 0.0;
	double amplitude = 1.0;
	double norm = 0.0;
	for (int octave = 0; octave < octaves; ++octave)
	{
		sum += amplitude * GradientNoise(p, seed + uint64_t(octave) * 0x632BE59BD9B4E019ull);
		norm += amplitude;
		amplitude *= 0.5;
		p *= 2.0;
	}
	return sum / norm;
}

}
```

- [ ] **Step 8: Build and run the tests**

Run build, then `ctest --test-dir /mnt/Passport/Worktrees/starscape/Build --output-on-failure`.
Expected: all 6 test cases pass.

- [ ] **Step 9: Run the reddening mutation**

Change `kMultiplier1` to `0xCD9E8D58u`; rebuild; the Philox test must fail. Revert, `touch src/rng.h`, rebuild, green.

- [ ] **Step 10: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Scaffold: CMake, vendored deps, Philox RNG, gradient noise

Mutation verified: wrong Philox multiplier reddens the KAT test.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Public format header `kre_starscape.h`

**Files:**
- Create: `include/kre_starscape.h`
- Test: `tests/test_format.cpp` (add to `STARSCAPE_TESTS`)

**Interfaces:**
- Produces (namespace `kre_starscape`): `kMagic`, `kVersion`, `enum class Encoding : uint32_t { RGB9E5 = 1, BC6H = 2 }`, `struct Header` (64 bytes), `struct Star` (16 bytes), `kRGB9E5Max`, `uint32_t EncodeRGB9E5(float r, float g, float b)`, `std::array<float,3> DecodeRGB9E5(uint32_t)`, `enum Face : uint32_t { PosX, NegX, PosY, NegY, PosZ, NegZ }`, `struct FaceCoord { uint32_t face; double s, t; }`, `FaceCoord DirToFace(double x, double y, double z)`, `std::array<double,3> FaceToDir(uint32_t face, double s, double t)`, `double TexelSolidAngle(uint32_t faceSize, uint32_t x, uint32_t y)`, `uint32_t MipCountFor(uint32_t faceSize)`, `uint64_t DiffuseBytesFor(uint32_t faceSize)`, `uint64_t MipOffset(uint32_t faceSize, uint32_t mip)`, `enum class Result`, `std::string_view Describe(Result)`, `struct View { Header header; std::span<uint8_t const> diffuse, stars; }`, `Result Parse(std::span<uint8_t const>, View&)`, `uint32_t DiffuseTexel(View const&, uint32_t mip, uint32_t face, uint32_t x, uint32_t y)`, `Star StarAt(View const&, uint32_t index)`.

Header layout (spec amendment: the spec's `pad[4]` became `valueScale`; stored values = physical × valueScale, so readers divide):

| offset | field |
|---|---|
| 0 | `uint32_t magic` |
| 4 | `uint32_t version` |
| 8 | `uint64_t sceneHash` |
| 16 | `uint32_t diffuseEncoding` |
| 20 | `uint32_t diffuseFaceSize` |
| 24 | `uint32_t diffuseMipCount` |
| 28 | `uint32_t starCount` |
| 32 | `uint64_t diffuseOffset` |
| 40 | `uint64_t diffuseBytes` |
| 48 | `uint64_t starsOffset` |
| 56 | `float fluxThreshold` |
| 60 | `float valueScale` |

- [ ] **Step 1: Write the failing tests**

`tests/test_format.cpp`:
```cpp
#include <kre_starscape.h>
#include <doctest.h>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

using namespace kre_starscape;

namespace
{

std::vector<uint8_t> ValidFile(uint32_t faceSize = 2, uint32_t starCount = 3)
{
	Header h{};
	h.magic           = kMagic;
	h.version         = kVersion;
	h.diffuseEncoding = uint32_t(Encoding::RGB9E5);
	h.diffuseFaceSize = faceSize;
	h.diffuseMipCount = MipCountFor(faceSize);
	h.starCount       = starCount;
	h.diffuseOffset   = sizeof(Header);
	h.diffuseBytes    = DiffuseBytesFor(faceSize);
	h.starsOffset     = h.diffuseOffset + h.diffuseBytes;
	h.fluxThreshold   = 0.01f;
	h.valueScale      = 1.0f;
	std::vector<uint8_t> bytes(size_t(h.starsOffset + uint64_t(starCount) * sizeof(Star)));
	std::memcpy(bytes.data(), &h, sizeof h);
	for (size_t i = sizeof h; i < bytes.size(); ++i)
	{
		bytes[i] = uint8_t(i * 31u);
	}
	return bytes;
}

template<class Mutate>
Result ParseMutated(Mutate&& mutate)
{
	std::vector<uint8_t> bytes = ValidFile();
	Header h;
	std::memcpy(&h, bytes.data(), sizeof h);
	mutate(h);
	std::memcpy(bytes.data(), &h, sizeof h);
	View view;
	return Parse(bytes, view);
}

}

TEST_CASE("A well-formed file parses into spans over the caller's bytes")
{
	std::vector<uint8_t> const bytes = ValidFile();
	View view;
	REQUIRE(Parse(bytes, view) == Result::Ok);
	CHECK(view.diffuse.size() == 120);   // 6 faces × (2×2 + 1×1) texels × 4 bytes
	CHECK(view.stars.size() == 48);
	CHECK(view.diffuse.data() == bytes.data() + 64);
	CHECK(DiffuseBytesFor(256) == 2097144);
	CHECK(MipCountFor(256) == 9);
	CHECK(MipCountFor(1) == 1);
}

TEST_CASE("Every truncation is rejected")
{
	std::vector<uint8_t> const bytes = ValidFile();
	for (size_t length = 0; length < bytes.size(); ++length)
	{
		View view;
		CHECK(Parse(std::span(bytes.data(), length), view) != Result::Ok);
	}
}

TEST_CASE("Each header defect is named")
{
	CHECK(ParseMutated([](Header& h) { h.magic ^= 1u; }) == Result::BadMagic);
	CHECK(ParseMutated([](Header& h) { h.version = 2; }) == Result::BadVersion);
	CHECK(ParseMutated([](Header& h) { h.diffuseEncoding = uint32_t(Encoding::BC6H); }) == Result::BadEncoding);
	CHECK(ParseMutated([](Header& h) { h.diffuseFaceSize = 3; }) == Result::BadFaceSize);
	CHECK(ParseMutated([](Header& h) { h.diffuseFaceSize = 0; }) == Result::BadFaceSize);
	CHECK(ParseMutated([](Header& h) { h.diffuseMipCount = 1; }) == Result::BadMipCount);
	CHECK(ParseMutated([](Header& h) { h.valueScale = 0.0f; }) == Result::BadValueScale);
	CHECK(ParseMutated([](Header& h) { h.valueScale = std::numeric_limits<float>::quiet_NaN(); }) == Result::BadValueScale);
	CHECK(ParseMutated([](Header& h) { h.diffuseBytes += 4; }) == Result::DiffuseSizeMismatch);
	CHECK(ParseMutated([](Header& h) { h.diffuseOffset = 8; }) == Result::DiffuseOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.diffuseOffset = 1000; }) == Result::DiffuseOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.starsOffset = ~uint64_t(0) - 8; }) == Result::StarsOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.starCount = ~0u; }) == Result::StarsOutOfBounds);
	CHECK(ParseMutated([](Header& h) { h.starsOffset = h.diffuseOffset + 4; }) == Result::RegionsOverlap);
}

TEST_CASE("Accessors read the documented byte positions")
{
	std::vector<uint8_t> bytes = ValidFile();
	uint32_t const texelValue = 0x12345678u;
	// mip 1, face NegY (3), texel (0,0): mip 0 holds 6 × 2 × 2 texels.
	size_t const texelOffset = 64 + size_t(MipOffset(2, 1)) + 3u * 4u;
	CHECK(MipOffset(2, 1) == 96);
	std::memcpy(bytes.data() + texelOffset, &texelValue, 4);
	Star const star{{0.0f, 0.6f, 0.8f}, 0xABCDEF01u};
	std::memcpy(bytes.data() + 64 + 120 + 2 * sizeof(Star), &star, sizeof star);
	View view;
	REQUIRE(Parse(bytes, view) == Result::Ok);
	CHECK(DiffuseTexel(view, 1, NegY, 0, 0) == texelValue);
	Star const read = StarAt(view, 2);
	CHECK(read.dir[2] == 0.8f);
	CHECK(read.fluxRGB9E5 == 0xABCDEF01u);
}

TEST_CASE("RGB9E5 encodes the reference values")
{
	CHECK(EncodeRGB9E5(1.0f, 1.0f, 1.0f) == 0x84020100u);
	CHECK(EncodeRGB9E5(kRGB9E5Max, kRGB9E5Max, kRGB9E5Max) == 0xFFFFFFFFu);
	CHECK(EncodeRGB9E5(1.0e9f, 0.0f, 0.0f) == EncodeRGB9E5(kRGB9E5Max, 0.0f, 0.0f));
	CHECK(EncodeRGB9E5(0.0f, 0.0f, 0.0f) == 0u);
	CHECK(EncodeRGB9E5(-1.0f, -2.0f, -3.0f) == 0u);
	CHECK(EncodeRGB9E5(std::numeric_limits<float>::quiet_NaN(), 0.0f, 0.0f) == 0u);
	auto const decoded = DecodeRGB9E5(EncodeRGB9E5(1.0f, 0.5f, 0.25f));
	CHECK(decoded[0] == 1.0f);
	CHECK(decoded[1] == 0.5f);
	CHECK(decoded[2] == 0.25f);
}

TEST_CASE("RGB9E5 round trip error stays within half a mantissa step of the largest channel")
{
	constexpr int n = 1000000;
	for (int i = 0; i < n; ++i)
	{
		float const big = std::ldexp(1.0f, -14) * std::pow(65408.0f / std::ldexp(1.0f, -14), float(i) / n);
		float const r = big;
		float const g = big * float((i % 97) + 1) / 97.0f;
		float const b = big * float((i % 13)) / 13.0f;
		auto const d = DecodeRGB9E5(EncodeRGB9E5(r, g, b));
		float const tolerance = big / 512.0f + 1e-30f;
		REQUIRE(std::abs(d[0] - r) <= tolerance);
		REQUIRE(std::abs(d[1] - g) <= tolerance);
		REQUIRE(std::abs(d[2] - b) <= tolerance);
	}
}

TEST_CASE("Cube faces follow the Vulkan selection table")
{
	FaceCoord c = DirToFace(1.0, 0.0, 0.0);
	CHECK(c.face == PosX);
	CHECK(c.s == 0.5);
	CHECK(c.t == 0.5);
	c = DirToFace(1.0, 0.5, 0.0);   // +y is up: the upper half of the face
	CHECK(c.face == PosX);
	CHECK(c.t == 0.25);
	c = DirToFace(0.5, 0.0, 1.0);
	CHECK(c.face == PosZ);
	CHECK(c.s == 0.75);
	c = DirToFace(0.0, -1.0, 0.0);
	CHECK(c.face == NegY);
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (double s = 0.05; s < 1.0; s += 0.1)
		{
			for (double t = 0.05; t < 1.0; t += 0.1)
			{
				auto const d = FaceToDir(face, s, t);
				FaceCoord const back = DirToFace(d[0], d[1], d[2]);
				REQUIRE(back.face == face);
				REQUIRE(back.s == doctest::Approx(s));
				REQUIRE(back.t == doctest::Approx(t));
			}
		}
	}
}

TEST_CASE("Texel solid angles cover the sphere")
{
	constexpr uint32_t size = 8;
	double total = 0.0;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t y = 0; y < size; ++y)
		{
			for (uint32_t x = 0; x < size; ++x)
			{
				total += TexelSolidAngle(size, x, y);
			}
		}
	}
	CHECK(total == doctest::Approx(4.0 * 3.14159265358979323846).epsilon(1e-12));
	CHECK(TexelSolidAngle(size, 3, 3) > TexelSolidAngle(size, 0, 0));
}
```

- [ ] **Step 2: Build to verify the tests fail**

Expected: compile error, `kre_starscape.h` not found.

- [ ] **Step 3: Write `include/kre_starscape.h`**

```cpp
#pragma once
// kre_starscape: the cooked night sky written by Starscape (https://github.com/SpehleonLP/Starscape).
// Header-only and standard-library-only so engines can include it as-is. Layout: docs/spec.md, "File format".
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <string_view>

namespace kre_starscape
{

inline constexpr uint32_t kMagic   = 0x5254534Bu;   // "KSTR" read as little-endian bytes
inline constexpr uint32_t kVersion = 1;
inline constexpr uint32_t kMaxFaceSize = 4096;

enum class Encoding : uint32_t
{
	RGB9E5 = 1,
	BC6H   = 2,   // reserved; no writer or reader yet
};

// Read and written with one memcpy; the static_asserts pin every offset.
struct Header
{
	uint32_t magic;
	uint32_t version;
	uint64_t sceneHash;        // FNV-1a 64 of the scene file bytes, for provenance
	uint32_t diffuseEncoding;
	uint32_t diffuseFaceSize;  // power of two
	uint32_t diffuseMipCount;  // full chain down to 1×1
	uint32_t starCount;
	uint64_t diffuseOffset;
	uint64_t diffuseBytes;
	uint64_t starsOffset;
	float    fluxThreshold;    // F_res, magnitude-0 flux units (physical, not scaled)
	float    valueScale;       // every stored radiance and flux = physical × valueScale
};
static_assert(sizeof(Header) == 64);
static_assert(offsetof(Header, sceneHash) == 8);
static_assert(offsetof(Header, diffuseOffset) == 32);
static_assert(offsetof(Header, starsOffset) == 48);
static_assert(offsetof(Header, valueScale) == 60);

// Uploads as a vertex buffer without repacking.
struct Star
{
	float    dir[3];       // unit vector, sky frame
	uint32_t fluxRGB9E5;   // linear Rec.709 flux × valueScale
};
static_assert(sizeof(Star) == 16);
static_assert(std::endian::native == std::endian::little, "kre_starscape is little-endian on disk");

// ---- RGB9E5 (VK_FORMAT_E5B9G9R9_UFLOAT_PACK32): 9-bit mantissas, shared 5-bit exponent, bias 15.

inline constexpr float kRGB9E5Max = 65408.0f;   // (511 / 512) × 2^16

inline uint32_t EncodeRGB9E5(float r, float g, float b)
{
	// `c > 0` is false for NaN, so NaN and negatives both clamp to zero, as the reference does.
	auto clampChannel = [](float c) { return c > 0.0f ? std::min(c, kRGB9E5Max) : 0.0f; };
	float const rc = clampChannel(r);
	float const gc = clampChannel(g);
	float const bc = clampChannel(b);
	float const maxChannel = std::max({rc, gc, bc});
	if (maxChannel == 0.0f)
	{
		return 0u;
	}
	// frexp gives maxChannel = m × 2^e with m in [0.5, 1), so floor(log2(maxChannel)) = e − 1 exactly,
	// where the reference's floor(log2()) can round the wrong way next to a power of two.
	int exponent = 0;
	std::frexp(maxChannel, &exponent);
	int sharedExponent = std::max(-16, exponent - 1) + 16;
	double scale = std::ldexp(1.0, sharedExponent - 24);
	if (int(std::floor(maxChannel / scale + 0.5)) == 512)
	{
		sharedExponent += 1;
		scale *= 2.0;
	}
	auto quantize = [scale](float c) { return uint32_t(std::floor(c / scale + 0.5)); };
	return quantize(rc) | (quantize(gc) << 9) | (quantize(bc) << 18) | (uint32_t(sharedExponent) << 27);
}

inline std::array<float, 3> DecodeRGB9E5(uint32_t packed)
{
	float const scale = std::ldexp(1.0f, int(packed >> 27) - 24);
	return {float(packed & 511u) * scale, float((packed >> 9) & 511u) * scale, float((packed >> 18) & 511u) * scale};
}

// ---- Cube map faces, Vulkan order and orientation ("Cube Map Face Selection" in the Vulkan spec).

enum Face : uint32_t
{
	PosX,
	NegX,
	PosY,
	NegY,
	PosZ,
	NegZ,
};

struct FaceCoord
{
	uint32_t face;
	double   s;   // [0, 1], left to right
	double   t;   // [0, 1], top to bottom (row 0 is the first row in memory)
};

inline FaceCoord DirToFace(double x, double y, double z)
{
	double const ax = std::abs(x);
	double const ay = std::abs(y);
	double const az = std::abs(z);
	uint32_t face;
	double sc;
	double tc;
	double ma;
	if (ax >= ay && ax >= az)
	{
		ma = ax;
		face = x >= 0.0 ? PosX : NegX;
		sc = x >= 0.0 ? -z : z;
		tc = -y;
	}
	else if (ay >= az)
	{
		ma = ay;
		face = y >= 0.0 ? PosY : NegY;
		sc = x;
		tc = y >= 0.0 ? z : -z;
	}
	else
	{
		ma = az;
		face = z >= 0.0 ? PosZ : NegZ;
		sc = z >= 0.0 ? x : -x;
		tc = -y;
	}
	return {face, 0.5 * (sc / ma + 1.0), 0.5 * (tc / ma + 1.0)};
}

// Unnormalised direction through (s, t) on `face`; s and t may lie outside [0, 1].
inline std::array<double, 3> FaceToDir(uint32_t face, double s, double t)
{
	double const u = 2.0 * s - 1.0;
	double const v = 2.0 * t - 1.0;
	switch (face)
	{
	case PosX:
		return {1.0, -v, -u};
	case NegX:
		return {-1.0, -v, u};
	case PosY:
		return {u, 1.0, v};
	case NegY:
		return {u, -1.0, -v};
	case PosZ:
		return {u, -v, 1.0};
	default:
		return {-u, -v, -1.0};
	}
}

inline double TexelSolidAngle(uint32_t faceSize, uint32_t x, uint32_t y)
{
	auto area = [](double u, double v) { return std::atan2(u * v, std::sqrt(u * u + v * v + 1.0)); };
	double const step = 2.0 / double(faceSize);
	double const u0 = -1.0 + double(x) * step;
	double const v0 = -1.0 + double(y) * step;
	double const u1 = u0 + step;
	double const v1 = v0 + step;
	return area(u0, v0) - area(u0, v1) - area(u1, v0) + area(u1, v1);
}

inline uint32_t MipCountFor(uint32_t faceSize)
{
	return uint32_t(std::bit_width(faceSize));
}

// Bytes from the start of the diffuse payload to the first texel of `mip`.
inline uint64_t MipOffset(uint32_t faceSize, uint32_t mip)
{
	uint64_t offset = 0;
	for (uint32_t level = 0; level < mip; ++level)
	{
		uint64_t const size = faceSize >> level;
		offset += 6u * size * size * 4u;
	}
	return offset;
}

inline uint64_t DiffuseBytesFor(uint32_t faceSize)
{
	return MipOffset(faceSize, MipCountFor(faceSize));
}

// ---- Parsing

enum class Result
{
	Ok,
	TooSmall,
	BadMagic,
	BadVersion,
	BadEncoding,
	BadFaceSize,
	BadMipCount,
	BadValueScale,
	DiffuseSizeMismatch,
	DiffuseOutOfBounds,
	StarsOutOfBounds,
	RegionsOverlap,
};

inline std::string_view Describe(Result result)
{
	switch (result)
	{
	case Result::Ok:                  return "ok";
	case Result::TooSmall:            return "file is smaller than its header";
	case Result::BadMagic:            return "not a kre_starscape file (bad magic)";
	case Result::BadVersion:          return "unsupported kre_starscape version";
	case Result::BadEncoding:         return "unsupported diffuse encoding";
	case Result::BadFaceSize:         return "diffuse face size is not a power of two in [1, 4096]";
	case Result::BadMipCount:         return "mip count does not match the face size";
	case Result::BadValueScale:       return "value scale is not a positive finite number";
	case Result::DiffuseSizeMismatch: return "diffuse byte count does not match the face size";
	case Result::DiffuseOutOfBounds:  return "diffuse payload lies outside the file";
	case Result::StarsOutOfBounds:    return "star array lies outside the file";
	case Result::RegionsOverlap:      return "diffuse payload and star array overlap";
	}
	return "unknown error";
}

struct View
{
	Header                   header{};
	std::span<uint8_t const> diffuse;   // mip-major, then face, then rows top to bottom
	std::span<uint8_t const> stars;     // header.starCount × Star
};

// Never asserts: every field is checked before anything is read through it.
inline Result Parse(std::span<uint8_t const> file, View& out)
{
	if (file.size() < sizeof(Header))
	{
		return Result::TooSmall;
	}
	Header h;
	std::memcpy(&h, file.data(), sizeof h);
	if (h.magic != kMagic)
	{
		return Result::BadMagic;
	}
	if (h.version != kVersion)
	{
		return Result::BadVersion;
	}
	if (h.diffuseEncoding != uint32_t(Encoding::RGB9E5))
	{
		return Result::BadEncoding;
	}
	if (h.diffuseFaceSize == 0 || h.diffuseFaceSize > kMaxFaceSize || !std::has_single_bit(h.diffuseFaceSize))
	{
		return Result::BadFaceSize;
	}
	if (h.diffuseMipCount != MipCountFor(h.diffuseFaceSize))
	{
		return Result::BadMipCount;
	}
	if (!(h.valueScale > 0.0f) || !std::isfinite(h.valueScale))
	{
		return Result::BadValueScale;
	}
	if (h.diffuseBytes != DiffuseBytesFor(h.diffuseFaceSize))
	{
		return Result::DiffuseSizeMismatch;
	}
	uint64_t const size = file.size();
	// Written as `length <= size − offset` so a hostile offset cannot overflow the sum.
	auto inBounds = [size](uint64_t offset, uint64_t length)
	{
		return offset >= sizeof(Header) && offset <= size && length <= size - offset;
	};
	if (!inBounds(h.diffuseOffset, h.diffuseBytes))
	{
		return Result::DiffuseOutOfBounds;
	}
	uint64_t const starBytes = uint64_t(h.starCount) * sizeof(Star);
	if (!inBounds(h.starsOffset, starBytes))
	{
		return Result::StarsOutOfBounds;
	}
	bool const disjoint = h.diffuseOffset + h.diffuseBytes <= h.starsOffset || h.starsOffset + starBytes <= h.diffuseOffset;
	if (starBytes != 0 && !disjoint)
	{
		return Result::RegionsOverlap;
	}
	out.header  = h;
	out.diffuse = file.subspan(size_t(h.diffuseOffset), size_t(h.diffuseBytes));
	out.stars   = file.subspan(size_t(h.starsOffset), size_t(starBytes));
	return Result::Ok;
}

// The payload need not be 4-byte aligned in the caller's buffer, so reads go through memcpy.
inline uint32_t DiffuseTexel(View const& view, uint32_t mip, uint32_t face, uint32_t x, uint32_t y)
{
	uint32_t const size = view.header.diffuseFaceSize >> mip;
	uint64_t const offset = MipOffset(view.header.diffuseFaceSize, mip) + ((uint64_t(face) * size + y) * size + x) * 4u;
	uint32_t value;
	std::memcpy(&value, view.diffuse.data() + offset, 4);
	return value;
}

inline Star StarAt(View const& view, uint32_t index)
{
	Star star;
	std::memcpy(&star, view.stars.data() + size_t(index) * sizeof(Star), sizeof star);
	return star;
}

}
```

- [ ] **Step 4: Build and run the tests**

Expected: all format test cases pass.

- [ ] **Step 5: Run the reddening mutations**

(a) Delete the `RegionsOverlap` check → "Each header defect is named" fails. (b) Change `sharedExponent - 24` to `- 23` in `EncodeRGB9E5` → "RGB9E5 encodes the reference values" fails. Revert both, `touch include/kre_starscape.h`, rebuild, green.

- [ ] **Step 6: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "kre_starscape.h: format structs, RGB9E5, cube faces, validating parser

Mutation verified: dropped overlap check and wrong exponent bias both redden.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Stellar populations and blackbody colour

**Files:**
- Create: `src/population.h`, `src/population.cpp`, `src/presets.h`, `src/presets.cpp`
- Test: `tests/test_population.cpp`

**Interfaces:**
- Consumes: nothing from earlier tasks except `glm`.
- Produces: `kBins = 24`, `kLogLuminosityMin = -4.0`, `kLogLuminositySpan = 10.0`, `kFluxAt10pc = 0.01169499`, `double BinLuminosity(int bin)`, `double BinEdge(int edge)`, `double Luminance(glm::dvec3 rgb)`, `struct Population { std::string name; std::array<double,kBins> fraction, temperature; }`, `struct PopulationTable { std::array<glm::dvec3,kBins> chroma; std::array<double,kBins> fraction, luminosity; std::array<glm::dvec3,kBins+1> lightBelow; }`, `glm::dvec3 BlackbodyChroma(double kelvin)`, `PopulationTable Tabulate(Population const&)`, `std::optional<int> FindPopulation(std::span<Population const>, std::string_view)`, `std::vector<Population> PresetPopulations()` (order: young-disk, old-disk, bulge, halo).

- [ ] **Step 1: Write the failing tests**

`tests/test_population.cpp`:
```cpp
#include "population.h"
#include "presets.h"
#include <doctest.h>
#include <cmath>
#include <numeric>

using namespace starscape;

TEST_CASE("Luminosity bins are log-spaced with geometric centres")
{
	CHECK(BinEdge(0) == doctest::Approx(1e-4));
	CHECK(BinEdge(kBins) == doctest::Approx(1e6));
	for (int bin = 0; bin < kBins; ++bin)
	{
		CHECK(BinLuminosity(bin) == doctest::Approx(std::sqrt(BinEdge(bin) * BinEdge(bin + 1))));
	}
}

TEST_CASE("Preset populations are normalised and keep their derived mean luminosity")
{
	std::vector<Population> const presets = PresetPopulations();
	REQUIRE(presets.size() == 4);
	CHECK(presets[0].name == "young-disk");
	CHECK(presets[1].name == "old-disk");
	CHECK(presets[2].name == "bulge");
	CHECK(presets[3].name == "halo");
	// Means from the Kroupa-IMF derivation recorded in presets.cpp.
	double const expectedMean[] = {29.42, 1.124, 1.295, 1.295};
	for (size_t i = 0; i < presets.size(); ++i)
	{
		double const sum = std::accumulate(presets[i].fraction.begin(), presets[i].fraction.end(), 0.0);
		CHECK(sum == doctest::Approx(1.0).epsilon(1e-5));
		PopulationTable const table = Tabulate(presets[i]);
		CHECK(Luminance(table.lightBelow[kBins]) == doctest::Approx(expectedMean[i]).epsilon(0.01));
		CHECK(Luminance(table.lightBelow[0]) == 0.0);
	}
}

TEST_CASE("Blackbody chroma has unit luminance and the right hue order")
{
	for (double kelvin : {2700.0, 5800.0, 6500.0, 20000.0, 60000.0})
	{
		CHECK(Luminance(BlackbodyChroma(kelvin)) == doctest::Approx(1.0));
	}
	glm::dvec3 const white = BlackbodyChroma(6500.0);
	CHECK(white.x == doctest::Approx(1.0).epsilon(0.2));
	CHECK(white.y == doctest::Approx(1.0).epsilon(0.2));
	CHECK(white.z == doctest::Approx(1.0).epsilon(0.2));
	glm::dvec3 const red = BlackbodyChroma(3000.0);
	CHECK(red.x > red.y);
	CHECK(red.y > red.z);
	glm::dvec3 const blue = BlackbodyChroma(20000.0);
	CHECK(blue.z > blue.y);
	CHECK(blue.y > blue.x);
}

TEST_CASE("Populations are found by name")
{
	std::vector<Population> const presets = PresetPopulations();
	CHECK(FindPopulation(presets, "bulge") == 2);
	CHECK_FALSE(FindPopulation(presets, "nebula").has_value());
}
```

- [ ] **Step 2: Build to verify the tests fail**

Add `src/population.cpp src/presets.cpp` to `STARSCAPE_SOURCES` and `tests/test_population.cpp` to `STARSCAPE_TESTS`. Expected: compile errors.

- [ ] **Step 3: Write `src/population.h` and `src/population.cpp`**

`src/population.h`:
```cpp
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
```

`src/population.cpp`:
```cpp
#include "population.h"
#include <cmath>
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

double Luminance(glm::dvec3 rgb)
{
	return 0.2126 * rgb.x + 0.7152 * rgb.y + 0.0722 * rgb.z;
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

PopulationTable Tabulate(Population const& population)
{
	PopulationTable table;
	for (int bin = 0; bin < kBins; ++bin)
	{
		size_t const i = size_t(bin);
		table.fraction[i]   = population.fraction[i];
		table.luminosity[i] = BinLuminosity(bin);
		table.chroma[i]     = BlackbodyChroma(population.temperature[i]);
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
```

- [ ] **Step 4: Write `src/presets.h` and the population half of `src/presets.cpp`**

`src/presets.h`:
```cpp
#pragma once
#include "population.h"
#include <vector>

namespace starscape
{

// In order: young-disk, old-disk, bulge, halo. Scenes index populations, so the order is stable.
std::vector<Population> PresetPopulations();

}
```

`src/presets.cpp`:
```cpp
#include "presets.h"

namespace starscape
{

// Derivation (offline script, reproduced here so the numbers can be regenerated):
// - Kroupa (2001) IMF, ξ(m) ∝ m^-1.3 for 0.08–0.5 M☉ and m^-2.3 above, sampled 0.08–150 M☉.
// - Main sequence: L_bol = 0.23 m^2.3 (m < 0.43), m^4 (< 2), 1.4 m^3.5 (< 55), 32000 m;
//   R = m^0.8 (m < 1) or m^0.57; T = 5772 (L / R²)^0.25.
// - L_V = L_bol × (blackbody flux under the CIE ȳ curve / σT⁴), normalised so the Sun has L_V = 1.
// - A turnoff mass removes stars above it; stars within 4 % above the turnoff become giants
//   (the post-main-sequence lifetime is ~10 % of the main-sequence one, lifetime ∝ m^-2.5),
//   spread log-uniformly over L_V 10–1000 at T 4900 → 4000 K. Heavier stars are remnants and,
//   like every star below the lowest edge, are counted in bin 0.
// - young-disk: no turnoff. old-disk: 1.0 M☉. bulge and halo: 0.85 M☉.
// - T per bin is the light-weighted mean; empty bins repeat the previous bin's T.
// Mean L_V per star: young-disk 29.42, old-disk 1.124, bulge 1.295, halo 1.295.
std::vector<Population> PresetPopulations()
{
	Population young{"young-disk",
		{2.041409e-01, 1.388415e-01, 1.299330e-01, 1.212253e-01, 1.121052e-01, 5.238831e-02, 4.934281e-02, 4.202055e-02,
		 3.494578e-02, 2.863004e-02, 2.334727e-02, 1.869960e-02, 1.551932e-02, 1.111291e-02, 7.276815e-03, 4.573593e-03,
		 2.717156e-03, 1.546605e-03, 8.465797e-04, 4.846047e-04, 3.021698e-04, 0.0, 0.0, 0.0},
		{2680, 2842, 2999, 3170, 3354, 3639, 4026, 4495, 5060, 5764, 6780, 8146,
		 9957, 12263, 15408, 19718, 25664, 33876, 45214, 61254, 65806, 65806, 65806, 65806}};

	Population old{"old-disk",
		{2.963755e-01, 1.388415e-01, 1.299330e-01, 1.212253e-01, 1.121052e-01, 5.238831e-02, 4.934281e-02, 4.202055e-02,
		 3.494578e-02, 1.789629e-02, 0.0, 0.0, 1.022090e-03, 1.034405e-03, 1.022090e-03, 1.022090e-03,
		 8.250609e-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2680, 2842, 2999, 3170, 3354, 3639, 4026, 4495, 5060, 5567, 5567, 5567,
		 4792, 4604, 4416, 4229, 4066, 4066, 4066, 4066, 4066, 4066, 4066, 4066}};

	Population bulge{"bulge",
		{3.183034e-01, 1.388415e-01, 1.299330e-01, 1.212253e-01, 1.121052e-01, 5.238831e-02, 4.934281e-02, 4.202055e-02,
		 3.005734e-02, 0.0, 0.0, 0.0, 1.199895e-03, 1.214352e-03, 1.199895e-03, 1.199895e-03,
		 9.685900e-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2680, 2842, 2999, 3170, 3354, 3639, 4026, 4495, 5001, 5001, 5001, 5001,
		 4792, 4604, 4416, 4229, 4066, 4066, 4066, 4066, 4066, 4066, 4066, 4066}};

	Population halo = bulge;
	halo.name = "halo";

	return {std::move(young), std::move(old), std::move(bulge), std::move(halo)};
}

}
```

- [ ] **Step 5: Build and run the tests**

Expected: all population test cases pass.

- [ ] **Step 6: Run the reddening mutation**

Change young-disk bin 20's fraction `3.021698e-04` to `3.021698e-03` → the mean-luminosity check fails. Revert, `touch src/presets.cpp`, rebuild, green.

- [ ] **Step 7: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Stellar populations: Kroupa-derived bin tables, blackbody chroma

Mutation verified: a mistyped bin fraction reddens the mean-luminosity check.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Galaxy model — components, arms, dust, presets

**Files:**
- Create: `src/galaxy.h`, `src/galaxy.cpp`
- Modify: `src/presets.h`, `src/presets.cpp` (add `PresetGalaxy`)
- Test: `tests/test_galaxy.cpp`

**Interfaces:**
- Consumes: `Fbm` (Task 1), `Population`, `FindPopulation` (Task 3).
- Produces: `struct Arms`, `struct Disk`, `struct Bulge`, `struct Halo`, `struct Dust`, `struct Galaxy`, `struct Contribution { double density; int population; }`, `struct GalaxySample { std::array<Contribution,5> parts; int count; double extinction; }`, `struct PreparedGalaxy { Galaxy galaxy; glm::dmat3 toLocal; glm::dvec3 position; double boundRadius; }`, `PreparedGalaxy PrepareGalaxy(Galaxy const&)`, `void EvaluateLocal(Galaxy const&, glm::dvec3 local, GalaxySample&)`, `void Evaluate(PreparedGalaxy const&, glm::dvec3 scenePosition, GalaxySample&)`, `double ArmProximity(Arms const&, double radius, double azimuth, double offset)`, `double BoundRadius(Galaxy const&)`, `double Sech2(double)`, `std::optional<Galaxy> PresetGalaxy(std::string_view, std::span<Population const>)`.

Spec refinements implemented here (already folded into `docs/spec.md`): every component has an explicit `cutoff` in parsecs; the thin disk's young population weight is `armYoungFraction × a` (default 0.05) rather than `a`; the bulge radius is floored at 1 pc; dust lanes use the thin disk's arm geometry; dust noise wavelength is half the arm width (half the dust scale height when there are no arms).

- [ ] **Step 1: Write the failing tests**

`tests/test_galaxy.cpp`:
```cpp
#include "galaxy.h"
#include "presets.h"
#include <doctest.h>
#include <cmath>

using namespace starscape;

namespace
{
constexpr double kPi = 3.14159265358979323846;

double TotalDensity(GalaxySample const& sample)
{
	double total = 0.0;
	for (int i = 0; i < sample.count; ++i)
	{
		total += sample.parts[size_t(i)].density;
	}
	return total;
}
}

TEST_CASE("Arm proximity peaks on the arm and falls as a Gaussian of arc length")
{
	Arms const arms{.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 2.0, .clumping = 0.0, .phaseDeg = 10.0, .r0 = 1000.0};
	double const radius = 5000.0;
	double const onArm = 10.0 * kPi / 180.0 + std::log(radius / arms.r0) / std::tan(arms.pitchDeg * kPi / 180.0);
	CHECK(ArmProximity(arms, radius, onArm, 0.0) == doctest::Approx(1.0));
	CHECK(ArmProximity(arms, radius, onArm + kPi, 0.0) == doctest::Approx(1.0));   // the second arm
	CHECK(ArmProximity(arms, radius, onArm + arms.width / radius, 0.0) == doctest::Approx(std::exp(-0.5)));
	// A positive offset moves the peak toward the inner edge, which lies at positive azimuth.
	CHECK(ArmProximity(arms, radius, onArm + 150.0 / radius, 150.0) == doctest::Approx(1.0));
	CHECK(ArmProximity(Arms{}, radius, onArm, 0.0) == 0.0);
}

TEST_CASE("Evaluation at the centre and on the axis is finite")
{
	std::vector<Population> const populations = PresetPopulations();
	for (char const* name : {"milky-way", "small-spiral", "elliptical", "dwarf"})
	{
		std::optional<Galaxy> const galaxy = PresetGalaxy(name, populations);
		REQUIRE(galaxy.has_value());
		for (glm::dvec3 const p : {glm::dvec3(0.0), glm::dvec3(0.0, 0.0, 300.0), glm::dvec3(1e-9, 0.0, 0.0)})
		{
			GalaxySample sample;
			EvaluateLocal(*galaxy, p, sample);
			CHECK(std::isfinite(TotalDensity(sample)));
			CHECK(std::isfinite(sample.extinction));
		}
	}
}

TEST_CASE("The bulge integrates to the Hernquist total")
{
	Galaxy galaxy;
	galaxy.bulge = {.enabled = true, .density = 1000.0, .scale = 500.0, .axisY = 1.0, .axisZ = 1.0, .barAngleDeg = 0.0, .cutoff = 5000.0, .population = 0};
	// ∫₀^{10a} 4π r² n_b a / (r (r + a)³) dr = 4π n_b × 100 / 242.
	double const expected = 4.0 * kPi * galaxy.bulge.density * 100.0 / 242.0;
	double total = 0.0;
	constexpr int steps = 200000;
	double const dr = galaxy.bulge.cutoff / steps;
	for (int i = 0; i < steps; ++i)
	{
		double const r = (i + 0.5) * dr;
		GalaxySample sample;
		EvaluateLocal(galaxy, {r, 0.0, 0.0}, sample);
		total += 4.0 * kPi * r * r * TotalDensity(sample) * dr;
	}
	CHECK(total == doctest::Approx(expected).epsilon(0.005));
}

TEST_CASE("Disks fall off as sech² in height and exponentially in radius")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 2.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0, .arms = {}, .population = 0};
	GalaxySample mid;
	GalaxySample high;
	GalaxySample far;
	EvaluateLocal(galaxy, {1000.0, 0.0, 0.0}, mid);
	EvaluateLocal(galaxy, {1000.0, 0.0, 600.0}, high);
	EvaluateLocal(galaxy, {4000.0, 0.0, 0.0}, far);
	CHECK(TotalDensity(high) / TotalDensity(mid) == doctest::Approx(Sech2(1.0)));
	CHECK(TotalDensity(far) / TotalDensity(mid) == doctest::Approx(std::exp(-1.0)));
	GalaxySample outside;
	EvaluateLocal(galaxy, {30001.0, 0.0, 0.0}, outside);
	CHECK(outside.count == 0);
}

TEST_CASE("Arms split the thin disk between its old and young populations")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0,
		.arms = {.count = 2, .pitchDeg = 15.0, .width = 300.0, .contrast = 3.0, .clumping = 0.0, .phaseDeg = 0.0, .r0 = 1000.0},
		.population = 1, .armPopulation = 0, .armYoungFraction = 0.1};
	double const radius = 5000.0;
	double const onArm = std::log(radius / 1000.0) / std::tan(15.0 * kPi / 180.0);
	GalaxySample sample;
	EvaluateLocal(galaxy, {radius * std::cos(onArm), radius * std::sin(onArm), 0.0}, sample);
	REQUIRE(sample.count == 2);
	CHECK(sample.parts[1].population == 0);
	CHECK(sample.parts[1].density / TotalDensity(sample) == doctest::Approx(0.1));
	CHECK(TotalDensity(sample) == doctest::Approx(3.0 * std::exp(-radius / 3000.0)));
}

TEST_CASE("Rotation maps scene positions into the galaxy frame")
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 3000.0, .scaleHeight = 300.0, .cutoff = 30000.0, .arms = {}, .population = 0};
	galaxy.position = {100.0, 0.0, 0.0};
	galaxy.rotationDeg = {90.0, 0.0, 0.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	GalaxySample rotated;
	GalaxySample reference;
	Evaluate(prepared, {100.0, 0.0, 2000.0}, rotated);   // in the rotated disk's plane
	EvaluateLocal(galaxy, {0.0, 2000.0, 0.0}, reference);
	CHECK(TotalDensity(rotated) == doctest::Approx(TotalDensity(reference)));
}

TEST_CASE("Dust is a uniform field when its scales are huge")
{
	Galaxy galaxy;
	galaxy.dust = {.enabled = true, .extinction = 0.002, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1000.0, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0};
	GalaxySample sample;
	EvaluateLocal(galaxy, {300.0, -200.0, 100.0}, sample);
	CHECK(sample.extinction == doctest::Approx(0.002));
	CHECK(BoundRadius(galaxy) == 1000.0);
}

TEST_CASE("Galaxy presets resolve their populations")
{
	std::vector<Population> const populations = PresetPopulations();
	std::optional<Galaxy> const milkyWay = PresetGalaxy("milky-way", populations);
	REQUIRE(milkyWay.has_value());
	CHECK(milkyWay->disk.population == 1);
	CHECK(milkyWay->disk.armPopulation == 0);
	CHECK(milkyWay->bulge.population == 2);
	CHECK(milkyWay->halo.population == 3);
	CHECK(BoundRadius(*milkyWay) == 30000.0);
	CHECK_FALSE(PresetGalaxy("andromeda", populations).has_value());
}
```

- [ ] **Step 2: Build to verify the tests fail**

Add `src/galaxy.cpp` and `tests/test_galaxy.cpp` to the lists. Expected: compile errors.

- [ ] **Step 3: Write `src/galaxy.h`**

```cpp
#pragma once
#include "population.h"
#include <array>
#include <cstdint>
#include <string>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace starscape
{

// Logarithmic spiral arms: arm k sits at azimuth phase + 2πk/count + ln(R / r0) / tan(pitch).
struct Arms
{
	int    count     = 0;
	double pitchDeg  = 12.0;
	double width     = 400.0;   // Gaussian σ across the arm, parsecs of arc length
	double contrast  = 1.0;     // density multiplier on the arm's crest
	double clumping  = 0.0;     // fbm amplitude breaking arms into clumps and spurs
	double phaseDeg  = 0.0;
	double r0        = 1000.0;
};

// Field order matters: presets and tests use designated initialisers.
struct Disk
{
	bool   enabled          = false;
	double density          = 0.0;      // stars/pc³ at the centre of the midplane
	double scaleLength      = 1.0;
	double scaleHeight      = 1.0;
	double cutoff           = 0.0;      // spherical radius beyond which the component is zero
	Arms   arms;
	int    population       = -1;
	int    armPopulation    = -1;       // -1: arms modulate density only
	double armYoungFraction = 0.0;      // share of armPopulation on an arm's crest
};

// Hernquist profile n_b a / (r' (r' + a)³) over a radius r' stretched along a bar.
struct Bulge
{
	bool   enabled     = false;
	double density     = 0.0;   // n_b; total stars = 2π n_b before the cutoff
	double scale       = 1.0;   // a
	double axisY       = 1.0;   // bar axis ratios relative to its long axis
	double axisZ       = 1.0;
	double barAngleDeg = 0.0;
	double cutoff      = 0.0;
	int    population  = -1;
};

struct Halo
{
	bool   enabled    = false;
	double density    = 0.0;   // n_h in n_h (1 + r / scale)^-3.5
	double scale      = 1.0;
	double cutoff     = 0.0;
	int    population = -1;
};

struct Dust
{
	bool   enabled     = false;
	double extinction  = 0.0;   // α₀, visual magnitudes per parsec at the centre of the midplane
	double scaleLength = 1.0;
	double scaleHeight = 1.0;
	double cutoff      = 0.0;
	double laneOffset  = 0.0;   // lane shift toward each arm's inner edge, in arm widths
	double contrast    = 1.0;
	double filament    = 0.0;
};

struct Galaxy
{
	std::string preset;
	glm::dvec3  position{0.0};      // scene frame, parsecs
	glm::dvec3  rotationDeg{0.0};   // Euler angles applied Z, then Y, then X
	uint64_t    seed = 1;
	Disk        disk;
	Disk        thickDisk;
	Bulge       bulge;
	Halo        halo;
	Dust        dust;
};

struct Contribution
{
	double density    = 0.0;   // stars/pc³
	int    population = -1;
};

// Thin disk (old + young), thick disk, bulge, halo: at most five parts per galaxy.
struct GalaxySample
{
	std::array<Contribution, 5> parts{};
	int                         count      = 0;
	double                      extinction = 0.0;   // α_V, magnitudes per parsec
};

struct PreparedGalaxy
{
	Galaxy     galaxy;
	glm::dmat3 toLocal{1.0};
	glm::dvec3 position{0.0};
	double     boundRadius = 0.0;
};

double Sech2(double x);
double ArmProximity(Arms const& arms, double radius, double azimuth, double offset);
double BoundRadius(Galaxy const& galaxy);
PreparedGalaxy PrepareGalaxy(Galaxy const& galaxy);
void EvaluateLocal(Galaxy const& galaxy, glm::dvec3 local, GalaxySample& out);

inline void Evaluate(PreparedGalaxy const& prepared, glm::dvec3 scenePosition, GalaxySample& out)
{
	EvaluateLocal(prepared.galaxy, prepared.toLocal * (scenePosition - prepared.position), out);
}

}
```

- [ ] **Step 4: Write `src/galaxy.cpp`**

```cpp
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
```

- [ ] **Step 5: Add `PresetGalaxy`**

In `src/presets.h` add `#include "galaxy.h"`, `#include <optional>`, `#include <span>`, `#include <string_view>` and declare:
```cpp
// Every field filled; a scene overrides any of them. Population indices resolve against `populations`.
std::optional<Galaxy> PresetGalaxy(std::string_view name, std::span<Population const> populations);
```

Append to `src/presets.cpp` (inside the namespace):
```cpp
// Milky Way numbers: thin disk 0.1 stars/pc³ at the Sun's 8.2 kpc (n₀ = 0.1 e^{8200/2600});
// thick disk ~6 % of that locally; bulge ≈ 2 × 10¹⁰ stars (n_b = N / 2π) with a 25° bar;
// halo ≈ 10⁻⁴ stars/pc³ locally; dust ≈ 1 mag/kpc locally (α₀ = 0.001 e^{8200/3200}).
std::optional<Galaxy> PresetGalaxy(std::string_view name, std::span<Population const> populations)
{
	auto population = [&](std::string_view populationName) { return FindPopulation(populations, populationName).value_or(-1); };
	Galaxy galaxy;
	galaxy.preset = std::string(name);
	if (name == "milky-way")
	{
		galaxy.disk = {.enabled = true, .density = 2.34, .scaleLength = 2600.0, .scaleHeight = 300.0, .cutoff = 18000.0,
			.arms = {.count = 4, .pitchDeg = 12.0, .width = 400.0, .contrast = 2.5, .clumping = 0.6, .phaseDeg = 0.0, .r0 = 1000.0},
			.population = population("old-disk"), .armPopulation = population("young-disk"), .armYoungFraction = 0.05};
		galaxy.thickDisk = {.enabled = true, .density = 0.058, .scaleLength = 3600.0, .scaleHeight = 900.0, .cutoff = 22000.0,
			.arms = {}, .population = population("old-disk")};
		galaxy.bulge = {.enabled = true, .density = 3.2e9, .scale = 600.0, .axisY = 0.5, .axisZ = 0.4, .barAngleDeg = 25.0,
			.cutoff = 6000.0, .population = population("bulge")};
		galaxy.halo = {.enabled = true, .density = 0.24, .scale = 1000.0, .cutoff = 30000.0, .population = population("halo")};
		galaxy.dust = {.enabled = true, .extinction = 0.013, .scaleLength = 3200.0, .scaleHeight = 100.0, .cutoff = 18000.0,
			.laneOffset = 0.15, .contrast = 3.0, .filament = 0.7};
		return galaxy;
	}
	if (name == "small-spiral")
	{
		galaxy.disk = {.enabled = true, .density = 1.5, .scaleLength = 1500.0, .scaleHeight = 250.0, .cutoff = 10000.0,
			.arms = {.count = 2, .pitchDeg = 18.0, .width = 300.0, .contrast = 3.0, .clumping = 0.7, .phaseDeg = 0.0, .r0 = 600.0},
			.population = population("old-disk"), .armPopulation = population("young-disk"), .armYoungFraction = 0.08};
		galaxy.bulge = {.enabled = true, .density = 3.2e8, .scale = 300.0, .axisY = 1.0, .axisZ = 1.0, .barAngleDeg = 0.0,
			.cutoff = 3000.0, .population = population("bulge")};
		galaxy.halo = {.enabled = true, .density = 0.05, .scale = 800.0, .cutoff = 15000.0, .population = population("halo")};
		galaxy.dust = {.enabled = true, .extinction = 0.008, .scaleLength = 1800.0, .scaleHeight = 90.0, .cutoff = 10000.0,
			.laneOffset = 0.15, .contrast = 3.0, .filament = 0.7};
		return galaxy;
	}
	if (name == "elliptical")
	{
		galaxy.bulge = {.enabled = true, .density = 1.6e10, .scale = 3000.0, .axisY = 0.8, .axisZ = 0.6, .barAngleDeg = 0.0,
			.cutoff = 30000.0, .population = population("bulge")};
		galaxy.halo = {.enabled = true, .density = 0.5, .scale = 3000.0, .cutoff = 60000.0, .population = population("halo")};
		return galaxy;
	}
	if (name == "dwarf")
	{
		galaxy.bulge = {.enabled = true, .density = 1.6e7, .scale = 500.0, .axisY = 1.0, .axisZ = 0.7, .barAngleDeg = 0.0,
			.cutoff = 5000.0, .population = population("old-disk")};
		return galaxy;
	}
	return std::nullopt;
}
```

- [ ] **Step 6: Build and run the tests**

Expected: all galaxy test cases pass.

- [ ] **Step 7: Run the reddening mutations**

(a) In the bulge, square `(rb + a)` instead of cubing it → "The bulge integrates to the Hernquist total" fails. (b) Remove the `radius < 1e-6` guard in `ArmProximity` → "Evaluation at the centre and on the axis is finite" fails. Revert, `touch src/galaxy.cpp`, rebuild, green.

- [ ] **Step 8: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Galaxy model: disks with log-spiral arms, Hernquist bar, halo, dust lanes, presets

Mutation verified: wrong Hernquist power and a missing axis guard both redden.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Scene file (TOML) and sky frame

**Files:**
- Create: `src/scene.h`, `src/scene.cpp`
- Test: `tests/test_scene.cpp`

**Interfaces:**
- Consumes: `Galaxy` and friends, `PresetGalaxy`, `PresetPopulations`, `FindPopulation`, `Expected`, `InputError`, `IoError`, `STARSCAPE_TRY`, `STARSCAPE_TRY_ASSIGN`.
- Produces: `enum class Quality { Final, Quick }`, `struct Output { uint32_t faceSize; std::optional<uint32_t> targetStars; std::optional<double> limitingMagnitude; Quality quality; uint64_t seed; uint32_t previewWidth; bool previewMarch; }`, `struct Observer { glm::dvec3 position, spinAxis, zeroMeridian; double sMin; }`, `struct Scene { Output output; Observer observer; std::vector<Galaxy> galaxies; std::vector<Population> populations; uint64_t hash; }`, `Expected<Scene> ParseScene(std::string_view text, std::string_view sourceName)`, `Expected<Scene> LoadScene(std::filesystem::path const&)`, `glm::dmat3 SkyToScene(Observer const&)` (columns are the sky x, y, z axes in the scene frame).

Scene keys (spec "Scene (TOML)", plus refinements folded into the spec): `[output]` `face_size`, `target_stars` | `limiting_magnitude`, `quality`, `seed`, `preview = {width, march}`; `[observer]` `position` (required), `spin_axis`, `zero_meridian`, `s_min`; `[[galaxy]]` `preset` (required), `position`, `rotation`, `seed`, and component tables `disk` (+ `arms` inline table, `population`, `arm_population`, `arm_young_fraction`), `thick_disk`, `bulge` (`density`, `scale`, `axis_y`, `axis_z`, `bar_angle_deg`), `halo`, `dust` (`extinction`, `lane_offset`, `contrast`, `filament`); every component also takes `enabled`, `cutoff` and its scale keys. **A component table that is present enables that component unless it says `enabled = false`.**

- [ ] **Step 1: Write the failing tests**

`tests/test_scene.cpp`:
```cpp
#include "scene.h"
#include <doctest.h>
#include <glm/geometric.hpp>
#include <string>

using namespace starscape;

namespace
{

constexpr char const* kMinimal = R"(
[output]
face_size = 64
target_stars = 1000

[observer]
position = [8200.0, 0.0, 20.0]

[[galaxy]]
preset = "milky-way"

  [galaxy.disk]
  scale_height = 250
)";

std::string ErrorOf(std::string_view text)
{
	Expected<Scene> const scene = ParseScene(text, "test.toml");
	REQUIRE_FALSE(scene.has_value());
	return scene.error().message;
}

}

TEST_CASE("A minimal scene takes its preset and applies overrides")
{
	Expected<Scene> const scene = ParseScene(kMinimal, "test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->output.faceSize == 64);
	CHECK(scene->output.targetStars == 1000u);
	CHECK_FALSE(scene->output.limitingMagnitude.has_value());
	CHECK(scene->output.quality == Quality::Final);
	REQUIRE(scene->galaxies.size() == 1);
	CHECK(scene->galaxies[0].disk.scaleLength == 2600.0);   // from the preset
	CHECK(scene->galaxies[0].disk.scaleHeight == 250.0);    // overridden
	CHECK(scene->observer.position.x == 8200.0);
	CHECK(scene->populations.size() == 4);
	CHECK(scene->hash != 0);
	CHECK(scene->hash == ParseScene(kMinimal, "other.toml")->hash);
}

TEST_CASE("Unknown keys are errors with a position")
{
	std::string const message = ErrorOf("[output]\nface_sise = 256\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n");
	CHECK(message.find("test.toml:2:") != std::string::npos);
	CHECK(message.find("output.face_sise") != std::string::npos);
	CHECK(message.find("unknown key") != std::string::npos);
}

TEST_CASE("Exactly one of target_stars and limiting_magnitude")
{
	CHECK(ErrorOf("[output]\ntarget_stars = 10\nlimiting_magnitude = 6.0\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("exactly one") != std::string::npos);
	CHECK(ErrorOf("[output]\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("exactly one") != std::string::npos);
}

TEST_CASE("Values are range- and type-checked")
{
	char const* tail = "[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n";
	CHECK(ErrorOf(std::string("[output]\ntarget_stars = 10\nface_size = 300\n") + tail).find("power of two") != std::string::npos);
	CHECK(ErrorOf(std::string("[output]\ntarget_stars = 10\nface_size = \"big\"\n") + tail).find("expected an integer") != std::string::npos);
	CHECK(ErrorOf(std::string("[output]\ntarget_stars = 10\nquality = \"best\"\n") + tail).find("quality") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"andromeda\"\n").find("unknown preset") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n[galaxy.bulge]\npopulation = \"nebula\"\n").find("unknown population") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\nspin_axis = [1.0, 0.0, 0.0]\nzero_meridian = [2.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("parallel") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n").find("galaxy") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[[galaxy]]\npreset = \"dwarf\"\n").find("observer") != std::string::npos);
}

TEST_CASE("A component table enables its component")
{
	Expected<Scene> const scene = ParseScene(
		"[output]\nlimiting_magnitude = 6.0\n[observer]\nposition = [0.0, 0.0, 0.0]\n"
		"[[galaxy]]\npreset = \"dwarf\"\n[galaxy.thick_disk]\ndensity = 0.1\nscale_length = 1000\nscale_height = 300\ncutoff = 5000\npopulation = \"old-disk\"\n"
		"[galaxy.halo]\nenabled = false\n",
		"test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->galaxies[0].thickDisk.enabled);
	CHECK(scene->galaxies[0].thickDisk.population == 1);
	CHECK_FALSE(scene->galaxies[0].halo.enabled);
}

TEST_CASE("The sky frame is orthonormal with +z on the spin axis")
{
	Observer observer;
	observer.spinAxis = {0.0, -0.46, 0.89};
	observer.zeroMeridian = {1.0, 0.3, 0.0};
	glm::dmat3 const sky = SkyToScene(observer);
	CHECK(glm::length(sky[0]) == doctest::Approx(1.0));
	CHECK(glm::length(sky[1]) == doctest::Approx(1.0));
	CHECK(glm::dot(sky[0], sky[1]) == doctest::Approx(0.0));
	CHECK(glm::dot(sky[0], sky[2]) == doctest::Approx(0.0));
	CHECK(glm::dot(sky[2], glm::normalize(observer.spinAxis)) == doctest::Approx(1.0));
	CHECK(glm::dot(glm::cross(sky[0], sky[1]), sky[2]) == doctest::Approx(1.0));   // right-handed
}
```

- [ ] **Step 2: Build to verify the tests fail**

Add `src/scene.cpp` and `tests/test_scene.cpp` to the lists. Expected: compile errors.

- [ ] **Step 3: Write `src/scene.h`**

```cpp
#pragma once
#include "common.h"
#include "galaxy.h"
#include "population.h"
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string_view>
#include <vector>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace starscape
{

enum class Quality
{
	Final,
	Quick,
};

struct Output
{
	uint32_t                faceSize = 256;
	std::optional<uint32_t> targetStars;
	std::optional<double>   limitingMagnitude;
	Quality                 quality = Quality::Final;
	uint64_t                seed = 1;   // star draws and jitter; galaxies carry their own noise seeds
	uint32_t                previewWidth = 8192;
	bool                    previewMarch = false;
};

struct Observer
{
	glm::dvec3 position{0.0};
	glm::dvec3 spinAxis{0.0, 0.0, 1.0};       // sky-frame +z
	glm::dvec3 zeroMeridian{1.0, 0.0, 0.0};   // sky-frame +x after removing its spin-axis component
	double     sMin = 0.05;                   // parsecs; nearer stars belong to the observer's own system
};

struct Scene
{
	Output                  output;
	Observer                observer;
	std::vector<Galaxy>     galaxies;
	std::vector<Population> populations;   // presets first; tests may append their own
	uint64_t                hash = 0;      // FNV-1a 64 of the scene text
};

Expected<Scene> ParseScene(std::string_view text, std::string_view sourceName);
Expected<Scene> LoadScene(std::filesystem::path const& path);
glm::dmat3 SkyToScene(Observer const& observer);

}
```

- [ ] **Step 4: Write `src/scene.cpp`**

```cpp
#include "scene.h"
#include "presets.h"
#include <toml.hpp>
#include <algorithm>
#include <bit>
#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <sstream>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

uint64_t Fnv1a64(std::string_view bytes)
{
	uint64_t hash = 0xCBF29CE484222325ull;
	for (char const c : bytes)
	{
		hash ^= uint8_t(c);
		hash *= 0x100000001B3ull;
	}
	return hash;
}

constexpr double kHuge = 1e12;

// Reads one TOML table, remembering every key it was asked for so Finish() can reject the rest:
// a misspelt key must not silently fall back to a preset value.
class TableReader
{
public:
	TableReader(toml::table const& table, std::string_view source, std::string&& path)
	: m_table(table)
	, m_source(source)
	, m_path(std::move(path))
	{
	}

	std::string_view Source() const { return m_source; }
	std::string const& Path() const { return m_path; }

	std::unexpected<Error> Fail(toml::source_region const& where, std::string_view key, std::string_view message) const
	{
		std::string name = m_path;
		if (!key.empty())
		{
			name = name.empty() ? std::string(key) : std::format("{}.{}", name, key);
		}
		if (name.empty())
		{
			name = "scene";
		}
		return InputError(std::format("{}:{}:{}: {}: {}", m_source, where.begin.line, where.begin.column, name, message));
	}

	std::unexpected<Error> Fail(std::string_view key, std::string_view message) const
	{
		return Fail(m_table.source(), key, message);
	}

	bool Has(std::string_view key) const { return m_table.contains(key); }

	Expected<void> Number(std::string_view key, double& out, double min, double max)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<double> const value = node->is_boolean() ? std::nullopt : node->value<double>();
		if (!value)
		{
			return Fail(node->source(), key, "expected a number");
		}
		if (!std::isfinite(*value) || *value < min || *value > max)
		{
			return Fail(node->source(), key, std::format("must be in [{}, {}]", min, max));
		}
		out = *value;
		return {};
	}

	Expected<void> Integer(std::string_view key, int64_t& out, int64_t min, int64_t max)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<int64_t> const value = node->value_exact<int64_t>();
		if (!value)
		{
			return Fail(node->source(), key, "expected an integer");
		}
		if (*value < min || *value > max)
		{
			return Fail(node->source(), key, std::format("must be in [{}, {}]", min, max));
		}
		out = *value;
		return {};
	}

	Expected<void> Boolean(std::string_view key, bool& out)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<bool> const value = node->value_exact<bool>();
		if (!value)
		{
			return Fail(node->source(), key, "expected true or false");
		}
		out = *value;
		return {};
	}

	Expected<void> Text(std::string_view key, std::string& out)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<std::string> value = node->value_exact<std::string>();
		if (!value)
		{
			return Fail(node->source(), key, "expected a string");
		}
		out = std::move(*value);
		return {};
	}

	Expected<void> Vector(std::string_view key, glm::dvec3& out)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		toml::array const* array = node->as_array();
		if (array == nullptr || array->size() != 3)
		{
			return Fail(node->source(), key, "expected an array of three numbers");
		}
		for (size_t i = 0; i < 3; ++i)
		{
			toml::node const& element = (*array)[i];
			std::optional<double> const value = element.is_boolean() ? std::nullopt : element.value<double>();
			if (!value || !std::isfinite(*value))
			{
				return Fail(element.source(), key, "expected an array of three numbers");
			}
			out[int(i)] = *value;
		}
		return {};
	}

	Expected<toml::table const*> Table(std::string_view key)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return nullptr;
		}
		if (!node->is_table())
		{
			return Fail(node->source(), key, "expected a table");
		}
		return node->as_table();
	}

	Expected<toml::array const*> Array(std::string_view key)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return nullptr;
		}
		if (!node->is_array())
		{
			return Fail(node->source(), key, "expected an array");
		}
		return node->as_array();
	}

	Expected<void> Finish() const
	{
		for (auto const& [key, node] : m_table)
		{
			if (std::find(m_taken.begin(), m_taken.end(), key.str()) == m_taken.end())
			{
				return Fail(node.source(), key.str(), "unknown key");
			}
		}
		return {};
	}

private:
	toml::node const* Take(std::string_view key)
	{
		m_taken.emplace_back(key);
		return m_table.get(key);
	}

	toml::table const&       m_table;
	std::string_view         m_source;
	std::string              m_path;
	std::vector<std::string> m_taken;
};

Expected<TableReader> ChildReader(TableReader& parent, std::string_view key, toml::table const& table)
{
	std::string path = parent.Path().empty() ? std::string(key) : std::format("{}.{}", parent.Path(), key);
	return TableReader(table, parent.Source(), std::move(path));
}

Expected<void> ReadPopulation(TableReader& reader, std::string_view key, std::span<Population const> populations, int& out)
{
	std::string name;
	STARSCAPE_TRY(reader.Text(key, name));
	if (name.empty())
	{
		return {};
	}
	std::optional<int> const found = FindPopulation(populations, name);
	if (!found)
	{
		return reader.Fail(key, std::format("unknown population '{}'", name));
	}
	out = *found;
	return {};
}

Expected<void> ReadArms(TableReader& parent, Arms& arms)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("arms"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "arms", *table));
	int64_t count = arms.count;
	STARSCAPE_TRY(reader.Integer("count", count, 0, 16));
	arms.count = int(count);
	STARSCAPE_TRY(reader.Number("pitch_deg", arms.pitchDeg, 1.0, 89.0));
	STARSCAPE_TRY(reader.Number("width", arms.width, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("contrast", arms.contrast, 0.0, 1000.0));
	STARSCAPE_TRY(reader.Number("clumping", arms.clumping, 0.0, 10.0));
	STARSCAPE_TRY(reader.Number("phase_deg", arms.phaseDeg, -360.0, 360.0));
	STARSCAPE_TRY(reader.Number("r0", arms.r0, 1e-3, kHuge));
	return reader.Finish();
}

Expected<void> ReadDisk(TableReader& parent, std::string_view key, Disk& disk, std::span<Population const> populations, bool thin)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table(key));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, key, *table));
	disk.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", disk.enabled));
	STARSCAPE_TRY(reader.Number("density", disk.density, 0.0, kHuge));
	STARSCAPE_TRY(reader.Number("scale_length", disk.scaleLength, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("scale_height", disk.scaleHeight, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("cutoff", disk.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(ReadPopulation(reader, "population", populations, disk.population));
	if (thin)
	{
		STARSCAPE_TRY(ReadPopulation(reader, "arm_population", populations, disk.armPopulation));
		STARSCAPE_TRY(reader.Number("arm_young_fraction", disk.armYoungFraction, 0.0, 1.0));
		STARSCAPE_TRY(ReadArms(reader, disk.arms));
	}
	return reader.Finish();
}

Expected<void> ReadBulge(TableReader& parent, Bulge& bulge, std::span<Population const> populations)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("bulge"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "bulge", *table));
	bulge.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", bulge.enabled));
	STARSCAPE_TRY(reader.Number("density", bulge.density, 0.0, kHuge * 1e6));
	STARSCAPE_TRY(reader.Number("scale", bulge.scale, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("axis_y", bulge.axisY, 1e-3, 1.0));
	STARSCAPE_TRY(reader.Number("axis_z", bulge.axisZ, 1e-3, 1.0));
	STARSCAPE_TRY(reader.Number("bar_angle_deg", bulge.barAngleDeg, -180.0, 180.0));
	STARSCAPE_TRY(reader.Number("cutoff", bulge.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(ReadPopulation(reader, "population", populations, bulge.population));
	return reader.Finish();
}

Expected<void> ReadHalo(TableReader& parent, Halo& halo, std::span<Population const> populations)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("halo"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "halo", *table));
	halo.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", halo.enabled));
	STARSCAPE_TRY(reader.Number("density", halo.density, 0.0, kHuge));
	STARSCAPE_TRY(reader.Number("scale", halo.scale, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("cutoff", halo.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(ReadPopulation(reader, "population", populations, halo.population));
	return reader.Finish();
}

Expected<void> ReadDust(TableReader& parent, Dust& dust)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("dust"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "dust", *table));
	dust.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", dust.enabled));
	STARSCAPE_TRY(reader.Number("extinction", dust.extinction, 0.0, 10.0));
	STARSCAPE_TRY(reader.Number("scale_length", dust.scaleLength, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("scale_height", dust.scaleHeight, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("cutoff", dust.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(reader.Number("lane_offset", dust.laneOffset, -5.0, 5.0));
	STARSCAPE_TRY(reader.Number("contrast", dust.contrast, 0.0, 1000.0));
	STARSCAPE_TRY(reader.Number("filament", dust.filament, 0.0, 10.0));
	return reader.Finish();
}

Expected<Galaxy> ReadGalaxy(toml::table const& table, std::string_view source, size_t index, std::span<Population const> populations)
{
	TableReader reader(table, source, std::format("galaxy[{}]", index));
	std::string preset;
	STARSCAPE_TRY(reader.Text("preset", preset));
	if (preset.empty())
	{
		return reader.Fail("preset", "every galaxy needs a preset");
	}
	std::optional<Galaxy> base = PresetGalaxy(preset, populations);
	if (!base)
	{
		return reader.Fail("preset", std::format("unknown preset '{}'", preset));
	}
	Galaxy galaxy = std::move(*base);
	STARSCAPE_TRY(reader.Vector("position", galaxy.position));
	STARSCAPE_TRY(reader.Vector("rotation", galaxy.rotationDeg));
	int64_t seed = int64_t(galaxy.seed);
	STARSCAPE_TRY(reader.Integer("seed", seed, 0, std::numeric_limits<int64_t>::max()));
	galaxy.seed = uint64_t(seed);
	STARSCAPE_TRY(ReadDisk(reader, "disk", galaxy.disk, populations, true));
	STARSCAPE_TRY(ReadDisk(reader, "thick_disk", galaxy.thickDisk, populations, false));
	STARSCAPE_TRY(ReadBulge(reader, galaxy.bulge, populations));
	STARSCAPE_TRY(ReadHalo(reader, galaxy.halo, populations));
	STARSCAPE_TRY(ReadDust(reader, galaxy.dust));
	STARSCAPE_TRY(reader.Finish());
	return galaxy;
}

Expected<void> ReadOutput(TableReader& root, Output& output)
{
	STARSCAPE_TRY_ASSIGN(table, root.Table("output"));
	if (table == nullptr)
	{
		return root.Fail("output", "missing [output] table");
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(root, "output", *table));
	int64_t faceSize = output.faceSize;
	STARSCAPE_TRY(reader.Integer("face_size", faceSize, 1, 4096));
	if (!std::has_single_bit(uint64_t(faceSize)))
	{
		return reader.Fail("face_size", "must be a power of two");
	}
	output.faceSize = uint32_t(faceSize);
	if (reader.Has("target_stars") == reader.Has("limiting_magnitude"))
	{
		return reader.Fail("", "set exactly one of target_stars and limiting_magnitude");
	}
	if (reader.Has("target_stars"))
	{
		int64_t target = 0;
		STARSCAPE_TRY(reader.Integer("target_stars", target, 1, 10000000));
		output.targetStars = uint32_t(target);
	}
	else
	{
		double magnitude = 0.0;
		STARSCAPE_TRY(reader.Number("limiting_magnitude", magnitude, -30.0, 30.0));
		output.limitingMagnitude = magnitude;
	}
	std::string quality;
	STARSCAPE_TRY(reader.Text("quality", quality));
	if (quality == "quick")
	{
		output.quality = Quality::Quick;
	}
	else if (quality.empty() || quality == "final")
	{
		output.quality = Quality::Final;
	}
	else
	{
		return reader.Fail("quality", "must be \"final\" or \"quick\"");
	}
	int64_t seed = int64_t(output.seed);
	STARSCAPE_TRY(reader.Integer("seed", seed, 0, std::numeric_limits<int64_t>::max()));
	output.seed = uint64_t(seed);
	STARSCAPE_TRY_ASSIGN(preview, reader.Table("preview"));
	if (preview != nullptr)
	{
		STARSCAPE_TRY_ASSIGN(previewReader, ChildReader(reader, "preview", *preview));
		int64_t width = output.previewWidth;
		STARSCAPE_TRY(previewReader.Integer("width", width, 16, 65536));
		if (width % 2 != 0)
		{
			return previewReader.Fail("width", "must be even (the preview is width × width/2)");
		}
		output.previewWidth = uint32_t(width);
		STARSCAPE_TRY(previewReader.Boolean("march", output.previewMarch));
		STARSCAPE_TRY(previewReader.Finish());
	}
	return reader.Finish();
}

Expected<void> ReadObserver(TableReader& root, Observer& observer)
{
	STARSCAPE_TRY_ASSIGN(table, root.Table("observer"));
	if (table == nullptr)
	{
		return root.Fail("observer", "missing [observer] table");
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(root, "observer", *table));
	if (!reader.Has("position"))
	{
		return reader.Fail("position", "the observer needs a position");
	}
	STARSCAPE_TRY(reader.Vector("position", observer.position));
	STARSCAPE_TRY(reader.Vector("spin_axis", observer.spinAxis));
	STARSCAPE_TRY(reader.Vector("zero_meridian", observer.zeroMeridian));
	STARSCAPE_TRY(reader.Number("s_min", observer.sMin, 1e-6, 1e6));
	if (glm::length(observer.spinAxis) < 1e-9)
	{
		return reader.Fail("spin_axis", "must not be zero");
	}
	glm::dvec3 const z = glm::normalize(observer.spinAxis);
	glm::dvec3 const across = observer.zeroMeridian - glm::dot(observer.zeroMeridian, z) * z;
	if (glm::length(across) < 1e-6 * std::max(1.0, glm::length(observer.zeroMeridian)))
	{
		return reader.Fail("zero_meridian", "must not be parallel to spin_axis");
	}
	return reader.Finish();
}

}

glm::dmat3 SkyToScene(Observer const& observer)
{
	glm::dvec3 const z = glm::normalize(observer.spinAxis);
	glm::dvec3 const x = glm::normalize(observer.zeroMeridian - glm::dot(observer.zeroMeridian, z) * z);
	glm::dvec3 const y = glm::cross(z, x);
	return glm::dmat3(x, y, z);
}

Expected<Scene> ParseScene(std::string_view text, std::string_view sourceName)
{
	toml::parse_result parsed = toml::parse(text, sourceName);
	if (!parsed)
	{
		toml::parse_error const& error = parsed.error();
		return InputError(std::format("{}:{}:{}: {}", sourceName, error.source().begin.line, error.source().begin.column, error.description()));
	}
	toml::table const& root = parsed.table();

	Scene scene;
	scene.populations = PresetPopulations();
	scene.hash = Fnv1a64(text);

	TableReader reader(root, sourceName, "");
	STARSCAPE_TRY(ReadOutput(reader, scene.output));
	STARSCAPE_TRY(ReadObserver(reader, scene.observer));
	STARSCAPE_TRY_ASSIGN(galaxies, reader.Array("galaxy"));
	if (galaxies == nullptr || galaxies->empty())
	{
		return reader.Fail("galaxy", "the scene needs at least one [[galaxy]]");
	}
	for (size_t i = 0; i < galaxies->size(); ++i)
	{
		toml::table const* table = (*galaxies)[i].as_table();
		if (table == nullptr)
		{
			return reader.Fail((*galaxies)[i].source(), "galaxy", "expected [[galaxy]] tables");
		}
		STARSCAPE_TRY_ASSIGN(galaxy, ReadGalaxy(*table, sourceName, i, scene.populations));
		scene.galaxies.push_back(std::move(galaxy));
	}
	STARSCAPE_TRY(reader.Finish());
	return scene;
}

Expected<Scene> LoadScene(std::filesystem::path const& path)
{
	std::ifstream file(path, std::ios::binary);
	if (!file)
	{
		return IoError(std::format("{}: cannot open scene file", path.string()));
	}
	std::ostringstream text;
	text << file.rdbuf();
	if (file.bad())
	{
		return IoError(std::format("{}: read failed", path.string()));
	}
	return ParseScene(text.str(), path.string());
}

}
```

`toml::parse` keeps a copy of `sourceName` in its source regions, so passing a view is fine. If toml++'s `value<double>()` rejects integer nodes in this version, replace both `value<double>()` calls with a helper that tries `value_exact<double>()` then `value_exact<int64_t>()`; the "A minimal scene…" test (`scale_height = 250`) catches that.

- [ ] **Step 5: Build and run the tests**

Expected: all scene test cases pass.

- [ ] **Step 6: Run the reddening mutation**

Make `Finish()` return `{}` unconditionally → "Unknown keys are errors with a position" fails. Revert, `touch src/scene.cpp`, rebuild, green.

- [ ] **Step 7: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Scene: strict TOML parsing over presets, sky frame

Mutation verified: disabling unknown-key rejection reddens its test.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: The ray march

**Files:**
- Create: `src/march.h`, `src/march.cpp`, `tests/test_scenes.h`
- Test: `tests/test_march.cpp`

**Interfaces:**
- Consumes: `Scene`, `Quality`, `SkyToScene` (Task 5); `PreparedGalaxy`, `PrepareGalaxy`, `Evaluate`, `Contribution`, `GalaxySample` (Task 4); `PopulationTable`, `Tabulate`, `kBins`, `kFluxAt10pc` (Task 3).
- Produces:
  - `struct MarchParams { int samplesPerAxis; double alpha, dsMin, dsMax; }`, `MarchParams ParamsFor(Quality)` — Final {4, 0.02, 0.01, 5}, Quick {1, 0.1, 0.01, 50}.
  - `struct PreparedScene { std::vector<PreparedGalaxy> galaxies; std::vector<PopulationTable> populations; glm::dvec3 origin; glm::dmat3 skyToScene; double sMin; }`, `PreparedScene PrepareScene(Scene const&)`.
  - `struct Interval { double begin, end; }`, `void RayIntervals(PreparedScene const&, glm::dvec3 dirScene, double maxDistance, std::vector<Interval>&)`.
  - `glm::dvec3 TransmittanceForDepth(double magnitudes)`, `glm::dvec3 TransmittanceTo(PreparedScene const&, MarchParams const&, glm::dvec3 dirSky, double distance)`.
  - `struct RayResult { glm::dvec3 radiance; std::array<double,kBins> budget; }`, `void MarchRay(PreparedScene const&, MarchParams const&, glm::dvec3 dirSky, double fluxThreshold, double sampleSolidAngle, double weight, RayResult&)`.
  - `struct FluxHistogram { kBuckets = 512; kLogMin = -20; kLogMax = 10; std::array<double,kBuckets> count; void Add(double flux, double stars); void Merge(FluxHistogram const&); double CountAbove(double flux) const; std::optional<double> FluxForCount(double target) const; }`, `void HistogramRay(PreparedScene const&, MarchParams const&, glm::dvec3 dirSky, double sampleSolidAngle, FluxHistogram&)`.
  - `struct StarSample { double distance; glm::dvec3 transmittance; int population; }`, `std::optional<StarSample> SampleStar(PreparedScene const&, MarchParams const&, glm::dvec3 dirSky, int bin, double fluxThreshold, double u)`.

`dirSky` is a unit vector in the sky frame; every function rotates it to the scene frame with `skyToScene`. `fluxThreshold = HUGE_VAL` makes every star diffuse.

- [ ] **Step 1: Write the shared test scenes**

`tests/test_scenes.h`:
```cpp
#pragma once
#include "presets.h"
#include "scene.h"
#include <format>

namespace starscape::test
{

inline constexpr double kPi = 3.14159265358979323846;

inline Scene BaseScene()
{
	Scene scene;
	scene.populations = PresetPopulations();
	scene.output.faceSize = 16;
	scene.output.quality = Quality::Quick;
	scene.output.limitingMagnitude = 6.0;
	return scene;
}

// Every star of this population sits in one luminosity bin, so counts and fluxes have closed forms.
inline int AddSingleBinPopulation(Scene& scene, int bin, double kelvin = 5800.0)
{
	Population population;
	population.name = std::format("single-{}", bin);
	population.fraction[size_t(bin)] = 1.0;
	population.temperature.fill(kelvin);
	scene.populations.push_back(std::move(population));
	return int(scene.populations.size()) - 1;
}

// (1 + r / 10¹²)^-3.5 is 1 to nine digits, so this is a uniform ball.
inline Galaxy UniformBall(double density, double cutoff, int population, glm::dvec3 position)
{
	Galaxy galaxy;
	galaxy.position = position;
	galaxy.halo = {.enabled = true, .density = density, .scale = 1e12, .cutoff = cutoff, .population = population};
	return galaxy;
}

inline Galaxy UniformDust(double magnitudesPerParsec, double cutoff)
{
	Galaxy galaxy;
	galaxy.dust = {.enabled = true, .extinction = magnitudesPerParsec, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = cutoff,
		.laneOffset = 0.0, .contrast = 1.0, .filament = 0.0};
	return galaxy;
}

inline Scene MilkyWayScene(glm::dvec3 observer)
{
	Scene scene = BaseScene();
	scene.galaxies.push_back(*PresetGalaxy("milky-way", scene.populations));
	scene.observer.position = observer;
	return scene;
}

}
```

- [ ] **Step 2: Write the failing tests**

`tests/test_march.cpp`:
```cpp
#include "march.h"
#include "test_scenes.h"
#include <doctest.h>
#include <cmath>
#include <vector>

using namespace starscape;
using namespace starscape::test;

namespace
{

// Total sky flux of everything within `halfAngle` of +x (sky frame = scene frame here):
// a polar grid fine enough that the integrand is smooth inside each cell.
glm::dvec3 FluxInCone(PreparedScene const& prepared, MarchParams const& params, double halfAngle)
{
	constexpr int thetaSteps = 200;
	constexpr int phiSteps = 64;
	glm::dvec3 total(0.0);
	double const dTheta = halfAngle / thetaSteps;
	double const dPhi = 2.0 * kPi / phiSteps;
	for (int i = 0; i < thetaSteps; ++i)
	{
		double const theta = (i + 0.5) * dTheta;
		for (int j = 0; j < phiSteps; ++j)
		{
			double const phi = (j + 0.5) * dPhi;
			glm::dvec3 const dir(std::cos(theta), std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi));
			RayResult result;
			MarchRay(prepared, params, dir, HUGE_VAL, 0.0, 1.0, result);
			total += result.radiance * (std::sin(theta) * dTheta * dPhi);
		}
	}
	return total;
}

}

TEST_CASE("Ray intervals clip, merge and skip")
{
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 12);
	scene.galaxies.push_back(UniformBall(1.0, 10.0, population, {100.0, 0.0, 0.0}));
	scene.galaxies.push_back(UniformBall(1.0, 10.0, population, {115.0, 0.0, 0.0}));   // overlaps the first
	scene.galaxies.push_back(UniformBall(1.0, 5.0, population, {0.0, 0.0, 0.0}));      // contains the observer
	PreparedScene const prepared = PrepareScene(scene);
	std::vector<Interval> intervals;

	RayIntervals(prepared, {1.0, 0.0, 0.0}, HUGE_VAL, intervals);
	REQUIRE(intervals.size() == 2);
	CHECK(intervals[0].begin == doctest::Approx(scene.observer.sMin));
	CHECK(intervals[0].end == doctest::Approx(5.0));
	CHECK(intervals[1].begin == doctest::Approx(90.0));
	CHECK(intervals[1].end == doctest::Approx(125.0));

	RayIntervals(prepared, {1.0, 0.0, 0.0}, 100.0, intervals);
	REQUIRE(intervals.size() == 2);
	CHECK(intervals[1].end == doctest::Approx(100.0));

	RayIntervals(prepared, {-1.0, 0.0, 0.0}, HUGE_VAL, intervals);
	REQUIRE(intervals.size() == 1);   // only the sphere around the observer

	// Passes 9.95 pc from the first sphere's centre (a short chord) and 11.4 pc from the second's (a miss).
	RayIntervals(prepared, glm::normalize(glm::dvec3(100.0, 10.0, 0.0)), HUGE_VAL, intervals);
	REQUIRE(intervals.size() == 2);
	CHECK(intervals[1].end - intervals[1].begin < 3.0);
}

TEST_CASE("Counts in a uniform medium follow the Euclidean F^-3/2 law")
{
	Scene scene = BaseScene();
	int const bin = 12;
	int const population = AddSingleBinPopulation(scene, bin);
	double const density = 0.01;
	scene.galaxies.push_back(UniformBall(density, 2000.0, population, {0.0, 0.0, 0.0}));
	PreparedScene const prepared = PrepareScene(scene);
	FluxHistogram histogram;
	// The medium is isotropic, so one ray standing for the whole sphere gives every count.
	HistogramRay(prepared, ParamsFor(Quality::Quick), {0.0, 0.0, 1.0}, 4.0 * kPi, histogram);

	double const luminosity = BinLuminosity(bin);
	auto fluxAt = [&](double s) { return luminosity * 100.0 * kFluxAt10pc / (s * s); };
	std::vector<double> logF;
	std::vector<double> logN;
	for (double s = 5.0; s <= 500.0; s *= 1.5)
	{
		double const expected = density * 4.0 / 3.0 * kPi * (s * s * s - std::pow(scene.observer.sMin, 3.0));
		double const counted = histogram.CountAbove(fluxAt(s));
		CHECK(counted == doctest::Approx(expected).epsilon(0.03));
		logF.push_back(std::log10(fluxAt(s)));
		logN.push_back(std::log10(counted));
	}
	double meanF = 0.0;
	double meanN = 0.0;
	for (size_t i = 0; i < logF.size(); ++i)
	{
		meanF += logF[i] / double(logF.size());
		meanN += logN[i] / double(logN.size());
	}
	double covariance = 0.0;
	double variance = 0.0;
	for (size_t i = 0; i < logF.size(); ++i)
	{
		covariance += (logF[i] - meanF) * (logN[i] - meanN);
		variance += (logF[i] - meanF) * (logF[i] - meanF);
	}
	CHECK(covariance / variance == doctest::Approx(-1.5).epsilon(0.05 / 1.5));
}

TEST_CASE("A distant dust-free galaxy delivers L (10/d)² C in total")
{
	Scene scene = BaseScene();
	int const bin = 14;
	int const population = AddSingleBinPopulation(scene, bin);
	double const density = 0.05;
	double const scale = 5000.0;
	double const cutoff = 20000.0;
	double const distance = 1.0e6;
	Galaxy galaxy;
	galaxy.position = {distance, 0.0, 0.0};
	galaxy.halo = {.enabled = true, .density = density, .scale = scale, .cutoff = cutoff, .population = population};
	scene.galaxies.push_back(galaxy);
	PreparedScene const prepared = PrepareScene(scene);

	// ∫₀^R 4π r² (1 + r/a)^-3.5 dr = 4π a³ [−2x^-½ + (4/3)x^-3/2 − 0.4x^-5/2]₁^{1+R/a}.
	auto antiderivative = [](double x) { return -2.0 * std::pow(x, -0.5) + 4.0 / 3.0 * std::pow(x, -1.5) - 0.4 * std::pow(x, -2.5); };
	double const stars = density * 4.0 * kPi * scale * scale * scale * (antiderivative(1.0 + cutoff / scale) - antiderivative(1.0));
	double const expected = stars * BinLuminosity(bin) * 100.0 * kFluxAt10pc / (distance * distance);

	glm::dvec3 const flux = FluxInCone(prepared, ParamsFor(Quality::Quick), std::asin(cutoff / distance) * 1.01);
	CHECK(Luminance(flux) == doctest::Approx(expected).epsilon(0.01));
}

TEST_CASE("A uniform dust ball dims each channel by exp(-k tau)")
{
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 14);
	Galaxy source;
	source.position = {1.0e6, 0.0, 0.0};
	source.halo = {.enabled = true, .density = 0.05, .scale = 5000.0, .cutoff = 20000.0, .population = population};
	scene.galaxies.push_back(source);
	MarchParams const params = ParamsFor(Quality::Quick);
	double const halfAngle = std::asin(0.02) * 1.01;
	glm::dvec3 const clear = FluxInCone(PrepareScene(scene), params, halfAngle);

	double const alpha = 0.002;
	double const radius = 1000.0;
	scene.galaxies.push_back(UniformDust(alpha, radius));
	PreparedScene const dusty = PrepareScene(scene);
	glm::dvec3 const dimmed = FluxInCone(dusty, params, halfAngle);

	double const tau = alpha * (radius - scene.observer.sMin) / 1.0857;
	CHECK(dimmed.x / clear.x == doctest::Approx(std::exp(-0.75 * tau)).epsilon(0.005));
	CHECK(dimmed.y / clear.y == doctest::Approx(std::exp(-1.00 * tau)).epsilon(0.005));
	CHECK(dimmed.z / clear.z == doctest::Approx(std::exp(-1.32 * tau)).epsilon(0.005));
	CHECK(TransmittanceTo(dusty, params, {0.0, 1.0, 0.0}, 1.0e7).y == doctest::Approx(std::exp(-tau)).epsilon(0.005));
}

TEST_CASE("Points and diffuse split one ray's light without loss")
{
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	PreparedScene const prepared = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir = glm::normalize(glm::dvec3(-1.0, 0.1, 0.02));
	RayResult allDiffuse;
	MarchRay(prepared, params, dir, HUGE_VAL, 1.0, 1.0, allDiffuse);
	for (double const threshold : {1e-2, 1e-4, 1e-6})
	{
		RayResult split;
		MarchRay(prepared, params, dir, threshold, 1.0, 1.0, split);
		CHECK(Luminance(split.radiance) < Luminance(allDiffuse.radiance));
		double pointCount = 0.0;
		for (double const count : split.budget)
		{
			pointCount += count;
		}
		CHECK(pointCount > 0.0);
	}
}

TEST_CASE("SampleStar picks a resolved distance in proportion to the count")
{
	Scene scene = BaseScene();
	int const bin = 12;
	int const population = AddSingleBinPopulation(scene, bin);
	scene.galaxies.push_back(UniformBall(0.01, 2000.0, population, {0.0, 0.0, 0.0}));
	PreparedScene const prepared = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	double const luminosity = BinLuminosity(bin);
	double const threshold = luminosity * 100.0 * kFluxAt10pc / (100.0 * 100.0);   // resolved out to 100 pc
	// Count ∝ s², so the median distance is 100 × 2^{-1/3}.
	std::optional<StarSample> const median = SampleStar(prepared, params, {0.0, 0.0, 1.0}, bin, threshold, 0.5);
	REQUIRE(median.has_value());
	CHECK(median->distance == doctest::Approx(100.0 * std::pow(0.5, 1.0 / 3.0)).epsilon(0.02));
	CHECK(median->population == population);
	std::optional<StarSample> const farthest = SampleStar(prepared, params, {0.0, 0.0, 1.0}, bin, threshold, 0.999999);
	REQUIRE(farthest.has_value());
	CHECK(farthest->distance <= 100.0 * 1.06);
	CHECK_FALSE(SampleStar(prepared, params, {0.0, 0.0, 1.0}, bin, HUGE_VAL, 0.5).has_value());
}
```

- [ ] **Step 3: Build to verify the tests fail**

Add `src/march.cpp` and `tests/test_march.cpp` to the lists. Expected: compile errors.

- [ ] **Step 4: Write `src/march.h`**

```cpp
#pragma once
#include "galaxy.h"
#include "population.h"
#include "scene.h"
#include <array>
#include <optional>
#include <vector>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace starscape
{

// Step length ds = clamp(alpha × s, dsMin, dsMax); samplesPerAxis² jittered rays per texel.
struct MarchParams
{
	int    samplesPerAxis = 4;
	double alpha          = 0.02;
	double dsMin          = 0.01;
	double dsMax          = 5.0;
};

MarchParams ParamsFor(Quality quality);

struct PreparedScene
{
	std::vector<PreparedGalaxy>  galaxies;
	std::vector<PopulationTable> populations;
	glm::dvec3                   origin{0.0};
	glm::dmat3                   skyToScene{1.0};
	double                       sMin = 0.05;
};

PreparedScene PrepareScene(Scene const& scene);

struct Interval
{
	double begin;
	double end;
};

// Sorted, merged stretches of the ray inside any galaxy's bounding sphere, clipped to [sMin, maxDistance].
void RayIntervals(PreparedScene const& scene, glm::dvec3 dirScene, double maxDistance, std::vector<Interval>& out);

// Per-channel transmittance behind `magnitudes` of visual extinction.
glm::dvec3 TransmittanceForDepth(double magnitudes);
glm::dvec3 TransmittanceTo(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double distance);

struct RayResult
{
	glm::dvec3                radiance{0.0};   // mag-0 flux per steradian, linear Rec.709
	std::array<double, kBins> budget{};        // expected resolved stars per luminosity bin
};

// Adds `weight` × this ray's diffuse radiance and its point budget over `sampleSolidAngle` to `out`.
void MarchRay(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double fluxThreshold,
	double sampleSolidAngle, double weight, RayResult& out);

// Expected star counts by log10 of apparent flux, for choosing the threshold that meets a star budget.
struct FluxHistogram
{
	static constexpr int    kBuckets = 512;
	static constexpr double kLogMin  = -20.0;
	static constexpr double kLogMax  = 10.0;
	static constexpr double kBucketWidth = (kLogMax - kLogMin) / kBuckets;

	std::array<double, kBuckets> count{};

	void Add(double flux, double stars);
	void Merge(FluxHistogram const& other);
	double CountAbove(double flux) const;
	// Smallest flux whose CountAbove is `target`; nullopt when the whole histogram holds fewer.
	std::optional<double> FluxForCount(double target) const;
};

void HistogramRay(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double sampleSolidAngle,
	FluxHistogram& out);

struct StarSample
{
	double     distance;
	glm::dvec3 transmittance;
	int        population;
};

// Re-marches `dirSky` and inverts the cumulative count of resolved stars of `bin` at quantile `u`.
// nullopt when nothing of that bin is resolved along the ray.
std::optional<StarSample> SampleStar(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, int bin,
	double fluxThreshold, double u);

}
```

- [ ] **Step 5: Write `src/march.cpp`**

```cpp
#include "march.h"
#include <algorithm>
#include <cmath>
#include <span>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

struct Step
{
	double                        s;
	double                        ds;
	glm::dvec3                    transmittance;
	std::span<Contribution const> parts;
};

// Steps every interval at its midpoints, handing each step's contributions and the
// transmittance up to it to `visit`. Returns the total visual extinction in magnitudes.
template<class Visit>
double Walk(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dir, double maxDistance, Visit&& visit)
{
	thread_local std::vector<Interval> intervals;
	thread_local std::vector<Contribution> parts;
	RayIntervals(scene, dir, maxDistance, intervals);
	double depth = 0.0;
	for (Interval const& interval : intervals)
	{
		double s = interval.begin;
		while (s < interval.end)
		{
			double const remaining = interval.end - s;
			double const ds = std::min(std::clamp(params.alpha * s, params.dsMin, params.dsMax), remaining);
			double const mid = s + 0.5 * ds;
			glm::dvec3 const position = scene.origin + mid * dir;
			parts.clear();
			double extinction = 0.0;
			for (PreparedGalaxy const& galaxy : scene.galaxies)
			{
				glm::dvec3 const offset = position - galaxy.position;
				if (glm::dot(offset, offset) > galaxy.boundRadius * galaxy.boundRadius)
				{
					continue;
				}
				GalaxySample sample;
				Evaluate(galaxy, position, sample);
				parts.insert(parts.end(), sample.parts.begin(), sample.parts.begin() + sample.count);
				extinction += sample.extinction;
			}
			visit(Step{mid, ds, TransmittanceForDepth(depth + 0.5 * extinction * ds), parts});
			depth += extinction * ds;
			s = ds == remaining ? interval.end : s + ds;
		}
	}
	return depth;
}

// A bin is resolved at this step when L × T_V × 100 × C / s² ≥ F_res.
double MinResolvedLuminosity(double fluxThreshold, Step const& step)
{
	double const tv = step.transmittance.y;
	if (!(tv > 0.0))
	{
		return HUGE_VAL;
	}
	return fluxThreshold * step.s * step.s / (100.0 * kFluxAt10pc * tv);
}

int FirstResolvedBin(PopulationTable const& population, double minLuminosity)
{
	return int(std::lower_bound(population.luminosity.begin(), population.luminosity.end(), minLuminosity) - population.luminosity.begin());
}

}

MarchParams ParamsFor(Quality quality)
{
	if (quality == Quality::Quick)
	{
		return {.samplesPerAxis = 1, .alpha = 0.1, .dsMin = 0.01, .dsMax = 50.0};
	}
	return {.samplesPerAxis = 4, .alpha = 0.02, .dsMin = 0.01, .dsMax = 5.0};
}

PreparedScene PrepareScene(Scene const& scene)
{
	PreparedScene prepared;
	for (Galaxy const& galaxy : scene.galaxies)
	{
		prepared.galaxies.push_back(PrepareGalaxy(galaxy));
	}
	for (Population const& population : scene.populations)
	{
		prepared.populations.push_back(Tabulate(population));
	}
	prepared.origin = scene.observer.position;
	prepared.skyToScene = SkyToScene(scene.observer);
	prepared.sMin = scene.observer.sMin;
	return prepared;
}

glm::dvec3 TransmittanceForDepth(double magnitudes)
{
	// Channel extinction relative to V approximates Cardelli R_V = 3.1 at 610, 550 and 465 nm;
	// 1.0857 = 2.5 / ln 10 turns magnitudes into optical depth.
	double const tau = magnitudes / 1.0857;
	return {std::exp(-0.75 * tau), std::exp(-1.00 * tau), std::exp(-1.32 * tau)};
}

void RayIntervals(PreparedScene const& scene, glm::dvec3 dir, double maxDistance, std::vector<Interval>& out)
{
	out.clear();
	for (PreparedGalaxy const& galaxy : scene.galaxies)
	{
		glm::dvec3 const offset = scene.origin - galaxy.position;
		double const b = glm::dot(offset, dir);
		double const c = glm::dot(offset, offset) - galaxy.boundRadius * galaxy.boundRadius;
		double const discriminant = b * b - c;
		if (discriminant <= 0.0)
		{
			continue;
		}
		double const root = std::sqrt(discriminant);
		double const begin = std::max(-b - root, scene.sMin);
		double const end = std::min(-b + root, maxDistance);
		if (end > begin)
		{
			out.push_back({begin, end});
		}
	}
	std::sort(out.begin(), out.end(), [](Interval const& a, Interval const& b) { return a.begin < b.begin; });
	// Overlapping spheres share their stretch of ray; merging marches each point once.
	size_t merged = 0;
	for (size_t i = 0; i < out.size(); ++i)
	{
		if (merged > 0 && out[i].begin <= out[merged - 1].end)
		{
			out[merged - 1].end = std::max(out[merged - 1].end, out[i].end);
		}
		else
		{
			out[merged++] = out[i];
		}
	}
	out.resize(merged);
}

glm::dvec3 TransmittanceTo(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double distance)
{
	return TransmittanceForDepth(Walk(scene, params, scene.skyToScene * dirSky, distance, [](Step const&) {}));
}

void MarchRay(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double fluxThreshold,
	double sampleSolidAngle, double weight, RayResult& out)
{
	Walk(scene, params, scene.skyToScene * dirSky, HUGE_VAL, [&](Step const& step)
	{
		double const resolvedFrom = MinResolvedLuminosity(fluxThreshold, step);
		for (Contribution const& part : step.parts)
		{
			PopulationTable const& population = scene.populations[size_t(part.population)];
			int const split = FirstResolvedBin(population, resolvedFrom);
			// n s² dΩ ds stars at flux 100 L C / s² each: the s² cancels, leaving radiance 100 n L C ds.
			out.radiance += (weight * part.density * 100.0 * kFluxAt10pc * step.ds) * (step.transmittance * population.lightBelow[size_t(split)]);
			for (int bin = split; bin < kBins; ++bin)
			{
				out.budget[size_t(bin)] += part.density * population.fraction[size_t(bin)] * step.s * step.s * sampleSolidAngle * step.ds;
			}
		}
	});
}

void FluxHistogram::Add(double flux, double stars)
{
	if (!(flux > 0.0) || !(stars > 0.0))
	{
		return;
	}
	int const bucket = std::clamp(int(std::floor((std::log10(flux) - kLogMin) / kBucketWidth)), 0, kBuckets - 1);
	count[size_t(bucket)] += stars;
}

void FluxHistogram::Merge(FluxHistogram const& other)
{
	for (size_t i = 0; i < count.size(); ++i)
	{
		count[i] += other.count[i];
	}
}

// Within a bucket, counts are taken as uniform in log flux.
double FluxHistogram::CountAbove(double flux) const
{
	double const position = (std::log10(flux) - kLogMin) / kBucketWidth;
	if (position >= kBuckets)
	{
		return 0.0;
	}
	int const bucket = std::max(0, int(std::floor(position)));
	double sum = count[size_t(bucket)] * std::clamp(1.0 - (position - bucket), 0.0, 1.0);
	for (int b = bucket + 1; b < kBuckets; ++b)
	{
		sum += count[size_t(b)];
	}
	return sum;
}

std::optional<double> FluxHistogram::FluxForCount(double target) const
{
	double cumulative = 0.0;
	for (int b = kBuckets - 1; b >= 0; --b)
	{
		double const here = count[size_t(b)];
		if (here > 0.0 && cumulative + here >= target)
		{
			double const fractionAbove = (target - cumulative) / here;
			return std::pow(10.0, kLogMin + (b + 1 - fractionAbove) * kBucketWidth);
		}
		cumulative += here;
	}
	return std::nullopt;
}

void HistogramRay(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double sampleSolidAngle,
	FluxHistogram& out)
{
	Walk(scene, params, scene.skyToScene * dirSky, HUGE_VAL, [&](Step const& step)
	{
		double const perLuminosity = step.transmittance.y * 100.0 * kFluxAt10pc / (step.s * step.s);
		double const volume = step.s * step.s * sampleSolidAngle * step.ds;
		for (Contribution const& part : step.parts)
		{
			PopulationTable const& population = scene.populations[size_t(part.population)];
			for (int bin = 0; bin < kBins; ++bin)
			{
				out.Add(population.luminosity[size_t(bin)] * perLuminosity, part.density * population.fraction[size_t(bin)] * volume);
			}
		}
	});
}

std::optional<StarSample> SampleStar(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, int bin,
	double fluxThreshold, double u)
{
	struct Entry
	{
		double     cumulative;
		double     count;
		double     begin;
		double     ds;
		glm::dvec3 transmittance;
		int        population;
	};
	thread_local std::vector<Entry> entries;
	entries.clear();
	double total = 0.0;
	Walk(scene, params, scene.skyToScene * dirSky, HUGE_VAL, [&](Step const& step)
	{
		double const resolvedFrom = MinResolvedLuminosity(fluxThreshold, step);
		for (Contribution const& part : step.parts)
		{
			PopulationTable const& population = scene.populations[size_t(part.population)];
			if (population.luminosity[size_t(bin)] < resolvedFrom)
			{
				continue;
			}
			double const count = part.density * population.fraction[size_t(bin)] * step.s * step.s * step.ds;
			if (!(count > 0.0))
			{
				continue;
			}
			total += count;
			entries.push_back({total, count, step.s - 0.5 * step.ds, step.ds, step.transmittance, part.population});
		}
	});
	if (!(total > 0.0))
	{
		return std::nullopt;
	}
	double const target = u * total;
	auto it = std::partition_point(entries.begin(), entries.end(), [target](Entry const& e) { return e.cumulative < target; });
	if (it == entries.end())
	{
		it = entries.end() - 1;
	}
	double const fraction = std::clamp((target - (it->cumulative - it->count)) / it->count, 0.0, 1.0);
	return StarSample{it->begin + fraction * it->ds, it->transmittance, it->population};
}

}
```

- [ ] **Step 6: Build and run the tests**

Expected: all march test cases pass. If the Euclidean absolute count misses by a few percent, check the step midpoint rule before touching tolerances: at `alpha = 0.1` the midpoint error in ∫s²ds is ~0.1 %.

- [ ] **Step 7: Run the reddening mutations**

(a) Drop `step.s * step.s` from the budget/histogram volume → the Euclidean test fails (slope −0.5). (b) In `TransmittanceForDepth` use `-1.00 * tau` for all three channels → the dust test fails. (c) Add `+ r / galaxy.halo.scale` to the halo density in `galaxy.cpp` → the distant-galaxy test fails. Revert all, `touch` the files, rebuild, green.

- [ ] **Step 8: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Ray march: intervals, diffuse/point split, flux histogram, star re-march

Mutation verified: dropped s², single-channel extinction and an unbounded density term each redden.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: The cook — calibration, texel march, star drawing

**Files:**
- Create: `src/parallel.h`, `src/cook.h`, `src/cook.cpp`
- Test: `tests/test_cook.cpp`

**Interfaces:**
- Consumes: everything from Task 6; `RandomStream`, `Purpose` (Task 1); `kre_starscape::FaceToDir`, `DirToFace`, `TexelSolidAngle` (Task 2); `Luminance` (Task 3).
- Produces:
  - `void ParallelFor(uint32_t count, unsigned threads, std::function<void(uint32_t)> const& body, std::function<void(uint32_t done)> const& progress)`.
  - `struct CookedStar { glm::dvec3 dir; glm::dvec3 flux; }`, `struct CookStats { double pointFlux, diffuseFlux; uint64_t starCount, droppedStars; bool budgetUnreachable; }`, `struct Cooked { uint32_t faceSize; std::vector<glm::dvec3> diffuse; std::vector<CookedStar> stars; double fluxThreshold; CookStats stats; }`, `struct CookOptions { unsigned threads; std::optional<Quality> quality; std::function<void(std::string_view phase, double fraction)> progress; }`.
  - `Cooked Cook(Scene const&, CookOptions const&)`, `FluxHistogram CalibrationHistogram(PreparedScene const&, MarchParams const&, uint32_t faceSize, unsigned threads)`, `uint32_t SamplePoisson(double lambda, RandomStream&)`, `glm::dvec3 SampleTexelDirection(uint32_t faceSize, uint32_t texel, RandomStream&)`, `glm::dvec3 TexelCentre(uint32_t faceSize, uint32_t texel)`.

Texel index = `face × S² + y × S + x`, matching the file's payload order. `diffuse` holds mip 0 in that order, in physical units. Stars are sorted by `Luminance(flux)`, brightest first, ties kept in texel order.

- [ ] **Step 1: Write the failing tests**

`tests/test_cook.cpp`:
```cpp
#include "cook.h"
#include "test_scenes.h"
#include <kre_starscape.h>
#include <doctest.h>
#include <cmath>
#include <thread>

using namespace starscape;
using namespace starscape::test;

namespace
{

unsigned Threads()
{
	return std::max(1u, std::thread::hardware_concurrency());
}

double TotalFlux(Cooked const& cooked)
{
	return cooked.stats.pointFlux + cooked.stats.diffuseFlux;
}

double SkyLatitude(glm::dvec3 dir)
{
	return std::asin(std::clamp(dir.z, -1.0, 1.0)) * 180.0 / kPi;
}

}

TEST_CASE("Poisson draws have the right mean and variance")
{
	for (double const lambda : {0.3, 5.0, 100.0})
	{
		RandomStream rng(9, Purpose::Test, uint32_t(lambda * 10), 0, 0);
		constexpr int n = 200000;
		double sum = 0.0;
		double sumSquares = 0.0;
		for (int i = 0; i < n; ++i)
		{
			double const k = SamplePoisson(lambda, rng);
			sum += k;
			sumSquares += k * k;
		}
		double const mean = sum / n;
		CHECK(mean == doctest::Approx(lambda).epsilon(0.02));
		CHECK(sumSquares / n - mean * mean == doctest::Approx(lambda).epsilon(0.03));
	}
	RandomStream rng(9, Purpose::Test, 0, 0, 1);
	CHECK(SamplePoisson(0.0, rng) == 0);
}

TEST_CASE("Texel directions stay in their texel and are uniform in solid angle")
{
	constexpr uint32_t size = 2;
	uint32_t const texel = 3 * size * size + 0;   // face NegY, a corner texel: the strongest Jacobian
	RandomStream rng(5, Purpose::Test, 0, 0, 0);
	// Split the texel into quarters and compare hit fractions with the quarters' solid angles.
	constexpr uint32_t fine = 4;
	double hits[2][2] = {};
	constexpr int n = 100000;
	for (int i = 0; i < n; ++i)
	{
		glm::dvec3 const dir = SampleTexelDirection(size, texel, rng);
		kre_starscape::FaceCoord const coord = kre_starscape::DirToFace(dir.x, dir.y, dir.z);
		REQUIRE(coord.face == kre_starscape::NegY);
		REQUIRE(coord.s < 0.5);
		REQUIRE(coord.t < 0.5);
		hits[int(coord.t * fine) % 2][int(coord.s * fine) % 2] += 1.0;
	}
	double const whole = kre_starscape::TexelSolidAngle(size, 0, 0);
	for (uint32_t y = 0; y < 2; ++y)
	{
		for (uint32_t x = 0; x < 2; ++x)
		{
			double const share = kre_starscape::TexelSolidAngle(fine, x, y) / whole;
			CHECK(hits[y][x] / n == doctest::Approx(share).epsilon(0.02));
		}
	}
}

TEST_CASE("Calibration is deterministic across thread counts")
{
	Scene const scene = MilkyWayScene({8200.0, 0.0, 20.0});
	PreparedScene const prepared = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	FluxHistogram const one = CalibrationHistogram(prepared, params, 16, 1);
	FluxHistogram const many = CalibrationHistogram(prepared, params, 16, 8);
	CHECK(one.count == many.count);
}

TEST_CASE("Energy is conserved whatever the star budget")
{
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	scene.output.faceSize = 32;
	scene.output.limitingMagnitude.reset();
	std::vector<Cooked> cooks;
	for (uint32_t const target : {500u, 4000u, 20000u})
	{
		scene.output.targetStars = target;
		cooks.push_back(Cook(scene, {.threads = Threads()}));
	}
	CHECK(double(cooks[1].stars.size()) == doctest::Approx(4000.0).epsilon(0.25));
	for (Cooked const& cooked : cooks)
	{
		double pointVariance = 0.0;
		for (CookedStar const& star : cooked.stars)
		{
			pointVariance += Luminance(star.flux) * Luminance(star.flux);
		}
		// Compound-Poisson noise of the points, plus RGB9E5-scale slack for the diffuse sum.
		double const tolerance = 3.0 * std::sqrt(pointVariance) + 0.005 * TotalFlux(cooks[0]);
		CHECK(std::abs(TotalFlux(cooked) - TotalFlux(cooks[0])) <= tolerance);
	}
	CHECK(cooks[2].stats.pointFlux > cooks[0].stats.pointFlux);
	CHECK(cooks[2].stats.diffuseFlux < cooks[0].stats.diffuseFlux);
}

TEST_CASE("From inside the disk the band is diffuse and the points are nearly uniform")
{
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});   // spin axis = galactic pole, so sky latitude = galactic latitude
	scene.output.faceSize = 64;
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 20000;
	Cooked const cooked = Cook(scene, {.threads = Threads()});

	double planeRadiance = 0.0;
	double planeAngle = 0.0;
	double poleRadiance = 0.0;
	double poleAngle = 0.0;
	uint32_t const size = cooked.faceSize;
	for (uint32_t texel = 0; texel < 6 * size * size; ++texel)
	{
		double const latitude = std::abs(SkyLatitude(TexelCentre(size, texel)));
		double const omega = kre_starscape::TexelSolidAngle(size, texel % size, (texel / size) % size);
		if (latitude < 5.0)
		{
			planeRadiance += Luminance(cooked.diffuse[texel]) * omega;
			planeAngle += omega;
		}
		else if (latitude > 85.0)
		{
			poleRadiance += Luminance(cooked.diffuse[texel]) * omega;
			poleAngle += omega;
		}
	}
	CHECK((planeRadiance / planeAngle) >= 10.0 * (poleRadiance / poleAngle));

	double planeStars = 0.0;
	double highStars = 0.0;
	for (CookedStar const& star : cooked.stars)
	{
		double const latitude = std::abs(SkyLatitude(star.dir));
		planeStars += latitude < 10.0 ? 1.0 : 0.0;
		highStars += latitude > 60.0 ? 1.0 : 0.0;
	}
	double const planeBand = 4.0 * kPi * std::sin(10.0 * kPi / 180.0);
	double const highCaps = 4.0 * kPi * (1.0 - std::sin(60.0 * kPi / 180.0));
	double const ratio = (planeStars / planeBand) / (highStars / highCaps);
	CHECK(ratio >= 1.5);
	CHECK(ratio <= 4.0);
}

TEST_CASE("Observer at the core")
{
	Scene scene = MilkyWayScene({0.0, 0.0, 0.0});
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 2000;
	Cooked const cooked = Cook(scene, {.threads = Threads()});
	CHECK(cooked.stars.size() >= 1000);
	CHECK(cooked.stars.size() <= 3000);
	for (glm::dvec3 const& radiance : cooked.diffuse)
	{
		REQUIRE(std::isfinite(radiance.x + radiance.y + radiance.z));
	}
}

TEST_CASE("Unreachable budget")
{
	Scene scene = BaseScene();
	Galaxy tiny = *PresetGalaxy("dwarf", scene.populations);
	tiny.bulge.density = 160.0;   // ≈ 580 stars in all
	// Near enough (14° across) that the coarse calibration texels cannot miss it, far enough
	// that the threshold floor still resolves every one of its stars.
	tiny.position = {20000.0, 0.0, 0.0};
	scene.galaxies.push_back(tiny);
	scene.output.limitingMagnitude.reset();
	scene.output.targetStars = 4000;
	Cooked const cooked = Cook(scene, {.threads = Threads()});
	CHECK(cooked.stats.budgetUnreachable);
	CHECK(cooked.stars.size() > 0);
	CHECK(cooked.stars.size() < 2000);
	CHECK(cooked.stats.diffuseFlux < 1e-3 * cooked.stats.pointFlux);
}
```

- [ ] **Step 2: Build to verify the tests fail**

Add `src/cook.cpp` and `tests/test_cook.cpp` to the lists. Expected: compile errors.

- [ ] **Step 3: Write `src/parallel.h`**

```cpp
#pragma once
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace starscape
{

// Runs body(i) for every i in [0, count) on `threads` threads, the caller included. The schedule
// is not deterministic, so bodies must write only to per-index slots; callers reduce in index order.
inline void ParallelFor(uint32_t count, unsigned threads, std::function<void(uint32_t)> const& body,
	std::function<void(uint32_t done)> const& progress)
{
	std::atomic<uint32_t> next{0};
	std::atomic<uint32_t> done{0};
	std::mutex progressMutex;
	auto worker = [&]
	{
		for (uint32_t index = next++; index < count; index = next++)
		{
			body(index);
			uint32_t const finished = ++done;
			if (progress)
			{
				std::lock_guard lock(progressMutex);
				progress(finished);
			}
		}
	};
	std::vector<std::jthread> pool;
	for (unsigned t = 1; t < std::max(1u, threads); ++t)
	{
		pool.emplace_back(worker);
	}
	worker();
}

}
```

- [ ] **Step 4: Write `src/cook.h`**

```cpp
#pragma once
#include "march.h"
#include "rng.h"
#include "scene.h"
#include <cstdint>
#include <functional>
#include <optional>
#include <string_view>
#include <vector>

namespace starscape
{

struct CookedStar
{
	glm::dvec3 dir;    // sky frame, unit
	glm::dvec3 flux;   // mag-0 units, linear Rec.709, physical
};

struct CookStats
{
	double   pointFlux = 0.0;      // Σ Luminance(star flux)
	double   diffuseFlux = 0.0;    // Σ Luminance(radiance) × texel solid angle
	uint64_t starCount = 0;
	uint64_t droppedStars = 0;     // drawn, but their direction held no resolved stars on re-march
	bool     budgetUnreachable = false;
};

struct Cooked
{
	uint32_t                faceSize = 0;
	std::vector<glm::dvec3> diffuse;   // mip 0, texel order, physical radiance
	std::vector<CookedStar> stars;     // brightest first
	double                  fluxThreshold = 0.0;
	CookStats               stats;
};

struct CookOptions
{
	unsigned                                                 threads = 1;
	std::optional<Quality>                                   quality;    // overrides the scene's
	std::function<void(std::string_view phase, double done)> progress;
};

Cooked Cook(Scene const& scene, CookOptions const& options);

FluxHistogram CalibrationHistogram(PreparedScene const& scene, MarchParams const& params, uint32_t faceSize, unsigned threads);
uint32_t SamplePoisson(double lambda, RandomStream& rng);
glm::dvec3 SampleTexelDirection(uint32_t faceSize, uint32_t texel, RandomStream& rng);
glm::dvec3 TexelCentre(uint32_t faceSize, uint32_t texel);

}
```

- [ ] **Step 5: Write `src/cook.cpp`**

```cpp
#include "cook.h"
#include "parallel.h"
#include <kre_starscape.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

constexpr double kTwoPi = 6.28318530717958647692;
constexpr uint32_t kCalibrationFaceSize = 64;

struct TexelCoord
{
	uint32_t face;
	uint32_t x;
	uint32_t y;
};

TexelCoord Decompose(uint32_t faceSize, uint32_t texel)
{
	return {texel / (faceSize * faceSize), texel % faceSize, (texel / faceSize) % faceSize};
}

glm::dvec3 FaceDirection(uint32_t face, double s, double t)
{
	std::array<double, 3> const d = kre_starscape::FaceToDir(face, s, t);
	return glm::normalize(glm::dvec3(d[0], d[1], d[2]));
}

}

glm::dvec3 TexelCentre(uint32_t faceSize, uint32_t texel)
{
	TexelCoord const c = Decompose(faceSize, texel);
	return FaceDirection(c.face, (c.x + 0.5) / faceSize, (c.y + 0.5) / faceSize);
}

uint32_t SamplePoisson(double lambda, RandomStream& rng)
{
	if (!(lambda > 0.0))
	{
		return 0;
	}
	if (lambda < 30.0)
	{
		// Knuth: multiply uniforms until the product drops below e^-λ.
		double const limit = std::exp(-lambda);
		uint32_t k = 0;
		double product = rng.NextDouble();
		while (product > limit)
		{
			++k;
			product *= rng.NextDouble();
		}
		return k;
	}
	// For large λ a normal approximation is within a fraction of a percent, and only the
	// densest texels near a core ever get here.
	double const u1 = std::max(rng.NextDouble(), 1e-300);
	double const u2 = rng.NextDouble();
	double const z = std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
	return uint32_t(std::max(0.0, std::floor(lambda + std::sqrt(lambda) * z + 0.5)));
}

glm::dvec3 SampleTexelDirection(uint32_t faceSize, uint32_t texel, RandomStream& rng)
{
	TexelCoord const c = Decompose(faceSize, texel);
	// Solid angle per unit face area is (1 + u² + v²)^-3/2; rejection against its maximum over
	// the texel makes directions uniform in solid angle rather than in face area.
	auto jacobian = [](double u, double v) { return std::pow(1.0 + u * u + v * v, -1.5); };
	double const step = 2.0 / faceSize;
	double const u0 = -1.0 + c.x * step;
	double const v0 = -1.0 + c.y * step;
	double const peak = std::max({jacobian(u0, v0), jacobian(u0 + step, v0), jacobian(u0, v0 + step), jacobian(u0 + step, v0 + step)});
	for (;;)
	{
		double const s = (c.x + rng.NextDouble()) / faceSize;
		double const t = (c.y + rng.NextDouble()) / faceSize;
		if (rng.NextDouble() * peak <= jacobian(2.0 * s - 1.0, 2.0 * t - 1.0))
		{
			return FaceDirection(c.face, s, t);
		}
	}
}

FluxHistogram CalibrationHistogram(PreparedScene const& scene, MarchParams const& params, uint32_t faceSize, unsigned threads)
{
	// One partial histogram per texel row, merged in row order, keeps the floating-point sum
	// independent of the schedule; a per-texel histogram would cost 512 doubles per texel.
	uint32_t const rows = 6 * faceSize;
	std::vector<FluxHistogram> partial(rows);
	ParallelFor(rows, threads, [&](uint32_t row)
	{
		for (uint32_t x = 0; x < faceSize; ++x)
		{
			uint32_t const texel = row * faceSize + x;
			TexelCoord const c = Decompose(faceSize, texel);
			HistogramRay(scene, params, TexelCentre(faceSize, texel), kre_starscape::TexelSolidAngle(faceSize, c.x, c.y), partial[row]);
		}
	}, {});
	FluxHistogram total;
	for (FluxHistogram const& row : partial)
	{
		total.Merge(row);
	}
	return total;
}

Cooked Cook(Scene const& scene, CookOptions const& options)
{
	Quality const quality = options.quality.value_or(scene.output.quality);
	MarchParams const params = ParamsFor(quality);
	PreparedScene const prepared = PrepareScene(scene);
	uint32_t const size = scene.output.faceSize;
	uint32_t const texels = 6 * size * size;
	uint64_t const seed = scene.output.seed;
	auto report = [&](std::string_view phase, uint32_t done, uint32_t total)
	{
		if (options.progress)
		{
			options.progress(phase, double(done) / double(total));
		}
	};

	Cooked cooked;
	cooked.faceSize = size;

	// 1. The flux threshold: given, or the one that leaves `target_stars` resolved.
	if (scene.output.limitingMagnitude)
	{
		cooked.fluxThreshold = std::pow(10.0, -0.4 * *scene.output.limitingMagnitude);
	}
	else
	{
		MarchParams calibration = params;
		calibration.samplesPerAxis = 1;
		report("calibrate", 0, 1);
		FluxHistogram const histogram = CalibrationHistogram(prepared, calibration, kCalibrationFaceSize, options.threads);
		std::optional<double> const threshold = histogram.FluxForCount(double(scene.output.targetStars.value_or(0)));
		cooked.stats.budgetUnreachable = !threshold.has_value();
		cooked.fluxThreshold = threshold.value_or(std::pow(10.0, FluxHistogram::kLogMin));
		report("calibrate", 1, 1);
	}

	// 2. Diffuse radiance and the per-bin point budget, one texel per slot.
	cooked.diffuse.assign(texels, glm::dvec3(0.0));
	std::vector<std::array<float, kBins>> budget(texels);
	int const perAxis = params.samplesPerAxis;
	double const weight = 1.0 / double(perAxis * perAxis);
	ParallelFor(texels, options.threads, [&](uint32_t texel)
	{
		TexelCoord const c = Decompose(size, texel);
		double const sampleSolidAngle = kre_starscape::TexelSolidAngle(size, c.x, c.y) * weight;
		RandomStream rng(seed, Purpose::Jitter, texel, 0, 0);
		RayResult result;
		for (int j = 0; j < perAxis; ++j)
		{
			for (int i = 0; i < perAxis; ++i)
			{
				double const s = (c.x + (i + rng.NextDouble()) / perAxis) / size;
				double const t = (c.y + (j + rng.NextDouble()) / perAxis) / size;
				MarchRay(prepared, params, FaceDirection(c.face, s, t), cooked.fluxThreshold, sampleSolidAngle, weight, result);
			}
		}
		cooked.diffuse[texel] = result.radiance;
		for (int bin = 0; bin < kBins; ++bin)
		{
			budget[texel][size_t(bin)] = float(result.budget[size_t(bin)]);
		}
	}, [&](uint32_t done) { report("march", done, texels); });

	// 3. Stars: Poisson counts from the budget, each placed by re-marching its own direction.
	std::vector<std::vector<CookedStar>> perTexel(texels);
	std::atomic<uint64_t> dropped{0};
	ParallelFor(texels, options.threads, [&](uint32_t texel)
	{
		for (int bin = 0; bin < kBins; ++bin)
		{
			RandomStream poisson(seed, Purpose::Poisson, texel, uint32_t(bin), 0);
			uint32_t const count = SamplePoisson(budget[texel][size_t(bin)], poisson);
			for (uint32_t index = 0; index < count; ++index)
			{
				RandomStream rng(seed, Purpose::Star, texel, uint32_t(bin), index);
				glm::dvec3 const dir = SampleTexelDirection(size, texel, rng);
				std::optional<StarSample> const sample = SampleStar(prepared, params, dir, bin, cooked.fluxThreshold, rng.NextDouble());
				if (!sample)
				{
					++dropped;
					continue;
				}
				PopulationTable const& population = prepared.populations[size_t(sample->population)];
				double const atObserver = population.luminosity[size_t(bin)] * 100.0 * kFluxAt10pc / (sample->distance * sample->distance);
				glm::dvec3 flux = atObserver * sample->transmittance * population.chroma[size_t(bin)];
				// Interpolating the distance inside a march step can leave a star a hair under the
				// threshold it was classified by; it was counted as a point, so it stays one, at the threshold.
				double const visual = atObserver * sample->transmittance.y;
				if (visual > 0.0 && visual < cooked.fluxThreshold)
				{
					flux *= cooked.fluxThreshold / visual;
				}
				perTexel[texel].push_back({dir, flux});
			}
		}
	}, [&](uint32_t done) { report("stars", done, texels); });

	for (std::vector<CookedStar>& stars : perTexel)
	{
		cooked.stars.insert(cooked.stars.end(), stars.begin(), stars.end());
	}
	std::stable_sort(cooked.stars.begin(), cooked.stars.end(), [](CookedStar const& a, CookedStar const& b)
	{
		return Luminance(a.flux) > Luminance(b.flux);
	});

	for (uint32_t texel = 0; texel < texels; ++texel)
	{
		TexelCoord const c = Decompose(size, texel);
		cooked.stats.diffuseFlux += Luminance(cooked.diffuse[texel]) * kre_starscape::TexelSolidAngle(size, c.x, c.y);
	}
	for (CookedStar const& star : cooked.stars)
	{
		cooked.stats.pointFlux += Luminance(star.flux);
	}
	cooked.stats.starCount = cooked.stars.size();
	cooked.stats.droppedStars = dropped.load();
	return cooked;
}

}
```

`budget` is stored as `float` to keep 393k texels × 24 bins at 37 MB; Poisson draws read it back in texel order, so determinism holds.

- [ ] **Step 6: Build and run the tests**

Expected: all cook test cases pass. The inside-disk and energy tests run full quick cooks (a few seconds each on 12 threads). **If "From inside the disk…" fails on its ratio checks, do not change the bounds or presets: stop and report both measured ratios to the coordinator** — they are model-level results the owner decides on.

- [ ] **Step 7: Run the reddening mutations**

(a) In `MarchRay`, replace `lightBelow[size_t(split)]` with `lightBelow[size_t(kBins)]` (resolved stars counted in both layers) → the energy test fails. (b) Replace `Sech2(p.z / (2.0 * disk.scaleHeight))` with `1.0` in `galaxy.cpp` → the inside-disk diffuse check fails. Revert, `touch`, rebuild, green.

- [ ] **Step 8: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Cook: threshold calibration, parallel texel march, Poisson star drawing

Mutation verified: double-counting resolved light breaks energy conservation;
a flat disk loses the band.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Writer — mips, value scale, file assembly, atomic write

**Files:**
- Create: `src/writer.h`, `src/writer.cpp`
- Test: `tests/test_writer.cpp`

**Interfaces:**
- Consumes: `Cooked`, `CookedStar`, `Cook` (Task 7); everything in `kre_starscape.h` (Task 2); `Expected`, `IoError` (Task 1).
- Produces: `std::vector<std::vector<glm::dvec3>> BuildMips(uint32_t faceSize, std::vector<glm::dvec3> const& level0)`, `float ChooseValueScale(double maxValue)`, `std::vector<uint8_t> BuildFile(Cooked const&, uint64_t sceneHash)`, `Expected<void> WriteFileAtomic(std::filesystem::path const&, std::span<uint8_t const>)`, `Expected<std::vector<uint8_t>> ReadFileBytes(std::filesystem::path const&)`.

Why `valueScale` matters (owner question, 2026-10-09): a night sky's difficulty is the dark end, not the bright one. RGB9E5 keeps 9 mantissa bits *relative to each texel's own brightest channel* over about 32 stops (1.5 × 10⁻⁵ to 65408; below that the mantissa loses bits), so a dust lane 7 stops below the band loses no relative precision. `valueScale` is a power of two that puts the scene's brightest value near the top of that range, so the dark end has every remaining stop. The "dark sky keeps its precision" test below pins this.

- [ ] **Step 1: Write the failing tests**

`tests/test_writer.cpp`:
```cpp
#include "writer.h"
#include "test_scenes.h"
#include <kre_starscape.h>
#include <doctest.h>
#include <cmath>
#include <filesystem>
#include <fstream>

using namespace starscape;
using namespace starscape::test;

namespace
{

Cooked Synthetic(uint32_t faceSize)
{
	Cooked cooked;
	cooked.faceSize = faceSize;
	cooked.fluxThreshold = 0.004;
	RandomStream rng(2, Purpose::Test, 0, 0, 0);
	for (uint32_t i = 0; i < 6 * faceSize * faceSize; ++i)
	{
		cooked.diffuse.push_back({rng.NextDouble() * 100.0, rng.NextDouble() * 100.0, rng.NextDouble() * 100.0});
	}
	cooked.stars.push_back({{0.0, 0.0, 1.0}, {3.0, 2.0, 1.0}});
	cooked.stars.push_back({{0.6, 0.8, 0.0}, {0.01, 0.02, 0.03}});
	return cooked;
}

glm::dvec3 Decoded(kre_starscape::View const& view, uint32_t mip, uint32_t face, uint32_t x, uint32_t y)
{
	auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::DiffuseTexel(view, mip, face, x, y));
	return glm::dvec3(v[0], v[1], v[2]) / double(view.header.valueScale);
}

void ExpectClose(glm::dvec3 actual, glm::dvec3 expected)
{
	double const tolerance = std::max({expected.x, expected.y, expected.z}) / 256.0;
	CHECK(std::abs(actual.x - expected.x) <= tolerance);
	CHECK(std::abs(actual.y - expected.y) <= tolerance);
	CHECK(std::abs(actual.z - expected.z) <= tolerance);
}

}

TEST_CASE("Mips are solid-angle-weighted means")
{
	std::vector<glm::dvec3> level0(6 * 4 * 4, glm::dvec3(0.0));
	level0[0] = {4.0, 0.0, 0.0};   // face 0, texel (0,0): a corner of the face
	auto const levels = BuildMips(4, level0);
	REQUIRE(levels.size() == 3);
	CHECK(levels[1].size() == 6 * 2 * 2);
	double const weights = kre_starscape::TexelSolidAngle(4, 0, 0) + kre_starscape::TexelSolidAngle(4, 1, 0)
		+ kre_starscape::TexelSolidAngle(4, 0, 1) + kre_starscape::TexelSolidAngle(4, 1, 1);
	CHECK(levels[1][0].x == doctest::Approx(4.0 * kre_starscape::TexelSolidAngle(4, 0, 0) / weights));
	CHECK(levels[2].size() == 6);
}

TEST_CASE("A cooked sky round-trips through the file")
{
	Cooked const cooked = Synthetic(4);
	std::vector<uint8_t> const bytes = BuildFile(cooked, 0x1234);
	kre_starscape::View view;
	REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
	CHECK(view.header.sceneHash == 0x1234);
	CHECK(view.header.diffuseFaceSize == 4);
	CHECK(view.header.starCount == 2);
	CHECK(view.header.fluxThreshold == doctest::Approx(0.004));
	for (uint32_t texel = 0; texel < 6 * 16; ++texel)
	{
		ExpectClose(Decoded(view, 0, texel / 16, texel % 4, (texel / 4) % 4), cooked.diffuse[texel]);
	}
	auto const levels = BuildMips(4, cooked.diffuse);
	ExpectClose(Decoded(view, 1, 5, 1, 1), levels[1][5 * 4 + 1 * 2 + 1]);
	kre_starscape::Star const star = kre_starscape::StarAt(view, 0);
	CHECK(star.dir[2] == 1.0f);
	auto const flux = kre_starscape::DecodeRGB9E5(star.fluxRGB9E5);
	ExpectClose(glm::dvec3(flux[0], flux[1], flux[2]) / double(view.header.valueScale), {3.0, 2.0, 1.0});
}

TEST_CASE("A dark sky keeps its precision")
{
	// Nine decades of radiance (30 stops) next to a bright star: each texel must decode
	// within one part in 256 of itself, however far below the brightest value it sits.
	Cooked cooked;
	cooked.faceSize = 4;
	for (uint32_t i = 0; i < 6 * 16; ++i)
	{
		double const value = std::pow(10.0, -4.0 + 9.0 * double(i) / 95.0);
		cooked.diffuse.push_back({value, 0.5 * value, 0.25 * value});
	}
	cooked.stars.push_back({{1.0, 0.0, 0.0}, {100.0, 100.0, 100.0}});
	std::vector<uint8_t> const bytes = BuildFile(cooked, 0);
	kre_starscape::View view;
	REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
	for (uint32_t texel = 0; texel < 6 * 16; ++texel)
	{
		glm::dvec3 const decoded = Decoded(view, 0, texel / 16, texel % 4, (texel / 4) % 4);
		CHECK(decoded.x == doctest::Approx(cooked.diffuse[texel].x).epsilon(1.0 / 256.0));
	}
}

TEST_CASE("The value scale is a power of two that keeps the maximum in range")
{
	for (double const maxValue : {1e-9, 0.37, 1.0, 400.0, 1e7})
	{
		float const scale = ChooseValueScale(maxValue);
		int exponent = 0;
		CHECK(std::frexp(scale, &exponent) == 0.5f);
		CHECK(maxValue * scale <= 32768.0);
		CHECK(maxValue * scale > 16384.0 * 0.999);
	}
	CHECK(ChooseValueScale(0.0) == 1.0f);
}

TEST_CASE("Cooks are byte-identical on any thread count")
{
	Scene const scene = MilkyWayScene({8200.0, 0.0, 20.0});   // limiting magnitude 6, face 16, quick
	std::vector<uint8_t> const one = BuildFile(Cook(scene, {.threads = 1}), scene.hash);
	std::vector<uint8_t> const many = BuildFile(Cook(scene, {.threads = 8}), scene.hash);
	CHECK(one == many);
}

TEST_CASE("Writes land atomically")
{
	std::filesystem::path const dir = std::filesystem::temp_directory_path() / "starscape-writer-test";
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	std::vector<uint8_t> const bytes = {1, 2, 3, 4};
	REQUIRE(WriteFileAtomic(dir / "sky.kre_starscape", bytes).has_value());
	Expected<std::vector<uint8_t>> const read = ReadFileBytes(dir / "sky.kre_starscape");
	REQUIRE(read.has_value());
	CHECK(*read == bytes);

	SUBCASE("failed rename leaves nothing")
	{
		std::filesystem::create_directories(dir / "occupied.kre_starscape" / "child");
		Expected<void> const result = WriteFileAtomic(dir / "occupied.kre_starscape", bytes);
		REQUIRE_FALSE(result.has_value());
		CHECK(result.error().kind == Error::Kind::Io);
		for (auto const& entry : std::filesystem::directory_iterator(dir))
		{
			CHECK(entry.path().filename().string().find(".partial") == std::string::npos);
		}
	}
	SUBCASE("missing directory")
	{
		Expected<void> const result = WriteFileAtomic(dir / "no" / "such" / "dir" / "sky.kre_starscape", bytes);
		REQUIRE_FALSE(result.has_value());
		CHECK(result.error().kind == Error::Kind::Io);
	}
	CHECK(ReadFileBytes(dir / "absent.kre_starscape").error().kind == Error::Kind::Io);
	std::filesystem::remove_all(dir);
}
```

- [ ] **Step 2: Build to verify the tests fail**

Add `src/writer.cpp` and `tests/test_writer.cpp` to the lists. Expected: compile errors.

- [ ] **Step 3: Write `src/writer.h`**

```cpp
#pragma once
#include "common.h"
#include "cook.h"
#include <cstdint>
#include <filesystem>
#include <span>
#include <vector>

namespace starscape
{

// levels[0] is `level0`; each further level halves the face size down to 1×1.
std::vector<std::vector<glm::dvec3>> BuildMips(uint32_t faceSize, std::vector<glm::dvec3> const& level0);

// Power of two placing `maxValue` in (16384, 32768], half the RGB9E5 maximum.
float ChooseValueScale(double maxValue);

std::vector<uint8_t> BuildFile(Cooked const& cooked, uint64_t sceneHash);

// Writes beside the target and renames into place, so a failure never leaves a half-written file.
Expected<void> WriteFileAtomic(std::filesystem::path const& path, std::span<uint8_t const> bytes);
Expected<std::vector<uint8_t>> ReadFileBytes(std::filesystem::path const& path);

}
```

- [ ] **Step 4: Write `src/writer.cpp`**

```cpp
#include "writer.h"
#include <kre_starscape.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <format>
#include <fstream>
#include <system_error>

namespace starscape
{

std::vector<std::vector<glm::dvec3>> BuildMips(uint32_t faceSize, std::vector<glm::dvec3> const& level0)
{
	std::vector<std::vector<glm::dvec3>> levels{level0};
	for (uint32_t size = faceSize; size > 1; size /= 2)
	{
		uint32_t const half = size / 2;
		std::vector<glm::dvec3> coarse(size_t(6) * half * half);
		std::vector<glm::dvec3> const& fine = levels.back();
		for (uint32_t face = 0; face < 6; ++face)
		{
			for (uint32_t y = 0; y < half; ++y)
			{
				for (uint32_t x = 0; x < half; ++x)
				{
					glm::dvec3 sum(0.0);
					double weight = 0.0;
					for (uint32_t dy = 0; dy < 2; ++dy)
					{
						for (uint32_t dx = 0; dx < 2; ++dx)
						{
							uint32_t const fx = 2 * x + dx;
							uint32_t const fy = 2 * y + dy;
							double const w = kre_starscape::TexelSolidAngle(size, fx, fy);
							sum += w * fine[(size_t(face) * size + fy) * size + fx];
							weight += w;
						}
					}
					coarse[(size_t(face) * half + y) * half + x] = sum / weight;
				}
			}
		}
		levels.push_back(std::move(coarse));
	}
	return levels;
}

float ChooseValueScale(double maxValue)
{
	if (!(maxValue > 0.0))
	{
		return 1.0f;
	}
	int const exponent = std::clamp(int(std::floor(std::log2(32768.0 / maxValue))), -100, 100);
	return std::ldexp(1.0f, exponent);
}

std::vector<uint8_t> BuildFile(Cooked const& cooked, uint64_t sceneHash)
{
	std::vector<std::vector<glm::dvec3>> const levels = BuildMips(cooked.faceSize, cooked.diffuse);
	double maxValue = 0.0;
	for (glm::dvec3 const& v : cooked.diffuse)
	{
		maxValue = std::max({maxValue, v.x, v.y, v.z});
	}
	for (CookedStar const& star : cooked.stars)
	{
		maxValue = std::max({maxValue, star.flux.x, star.flux.y, star.flux.z});
	}
	float const scale = ChooseValueScale(maxValue);

	kre_starscape::Header header{};
	header.magic           = kre_starscape::kMagic;
	header.version         = kre_starscape::kVersion;
	header.sceneHash       = sceneHash;
	header.diffuseEncoding = uint32_t(kre_starscape::Encoding::RGB9E5);
	header.diffuseFaceSize = cooked.faceSize;
	header.diffuseMipCount = kre_starscape::MipCountFor(cooked.faceSize);
	header.starCount       = uint32_t(cooked.stars.size());
	header.diffuseOffset   = sizeof(kre_starscape::Header);
	header.diffuseBytes    = kre_starscape::DiffuseBytesFor(cooked.faceSize);
	header.starsOffset     = header.diffuseOffset + header.diffuseBytes;
	header.fluxThreshold   = float(cooked.fluxThreshold);
	header.valueScale      = scale;

	std::vector<uint8_t> bytes(size_t(header.starsOffset + uint64_t(header.starCount) * sizeof(kre_starscape::Star)));
	std::memcpy(bytes.data(), &header, sizeof header);
	uint8_t* cursor = bytes.data() + header.diffuseOffset;
	auto pack = [scale](glm::dvec3 const& v)
	{
		return kre_starscape::EncodeRGB9E5(float(v.x * scale), float(v.y * scale), float(v.z * scale));
	};
	for (std::vector<glm::dvec3> const& level : levels)
	{
		for (glm::dvec3 const& texel : level)
		{
			uint32_t const packed = pack(texel);
			std::memcpy(cursor, &packed, 4);
			cursor += 4;
		}
	}
	for (CookedStar const& star : cooked.stars)
	{
		kre_starscape::Star const record{{float(star.dir.x), float(star.dir.y), float(star.dir.z)}, pack(star.flux)};
		std::memcpy(cursor, &record, sizeof record);
		cursor += sizeof record;
	}
	return bytes;
}

Expected<void> WriteFileAtomic(std::filesystem::path const& path, std::span<uint8_t const> bytes)
{
	std::filesystem::path partial = path;
	partial += ".partial";
	{
		std::ofstream file(partial, std::ios::binary | std::ios::trunc);
		if (!file)
		{
			return IoError(std::format("{}: cannot create", partial.string()));
		}
		file.write(reinterpret_cast<char const*>(bytes.data()), std::streamsize(bytes.size()));
		file.close();
		if (!file)
		{
			std::error_code ignored;
			std::filesystem::remove(partial, ignored);
			return IoError(std::format("{}: write failed", partial.string()));
		}
	}
	std::error_code error;
	std::filesystem::rename(partial, path, error);
	if (error)
	{
		std::error_code ignored;
		std::filesystem::remove(partial, ignored);
		return IoError(std::format("{}: {}", path.string(), error.message()));
	}
	return {};
}

Expected<std::vector<uint8_t>> ReadFileBytes(std::filesystem::path const& path)
{
	std::ifstream file(path, std::ios::binary | std::ios::ate);
	if (!file)
	{
		return IoError(std::format("{}: cannot open", path.string()));
	}
	std::streamsize const size = file.tellg();
	file.seekg(0);
	std::vector<uint8_t> bytes(size_t(std::max<std::streamsize>(size, 0)));
	file.read(reinterpret_cast<char*>(bytes.data()), size);
	if (!file)
	{
		return IoError(std::format("{}: read failed", path.string()));
	}
	return bytes;
}

}
```

- [ ] **Step 5: Build and run the tests**

Expected: all writer test cases pass.

- [ ] **Step 6: Run the reddening mutations**

(a) In `Cook`, seed the jitter stream with `uint64_t(std::hash<std::thread::id>{}(std::this_thread::get_id()))` instead of `seed` → "Cooks are byte-identical…" fails. (b) In `WriteFileAtomic`, delete the `remove(partial…)` after a failed rename → "failed rename leaves nothing" fails. (c) Make `BuildMips` weight every child by 1 → "Mips are solid-angle-weighted means" fails. Revert, `touch`, rebuild, green.

- [ ] **Step 7: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Writer: solid-angle mips, power-of-two value scale, atomic file write

Mutation verified: thread-keyed jitter, a leaked .partial and unweighted mips each redden.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Equirectangular preview and Radiance `.hdr`

**Files:**
- Create: `src/preview.h`, `src/preview.cpp`
- Test: `tests/test_preview.cpp`

**Interfaces:**
- Consumes: `kre_starscape::View`, `DirToFace`, `FaceToDir`, `DiffuseTexel`, `StarAt`, `DecodeRGB9E5` (Task 2); `PreparedScene`, `MarchParams`, `MarchRay` (Task 6); `ParallelFor` (Task 7); `BuildFile`, `Cooked` (Tasks 7–8, tests only).
- Produces: `struct Image { uint32_t width, height; std::vector<float> rgb; }`, `double PixelSolidAngle(uint32_t width, uint32_t row)`, `glm::dvec3 PixelDirection(uint32_t width, double column, double row)`, `glm::dvec3 SampleDiffuse(kre_starscape::View const&, glm::dvec3 dir)`, `Image RenderPreview(kre_starscape::View const&, uint32_t width, double exposureEV, unsigned threads)`, `Image RenderMarchedPreview(PreparedScene const&, MarchParams const&, kre_starscape::View const&, uint32_t width, double exposureEV, unsigned threads)`, `std::vector<uint8_t> EncodeRadiance(Image const&)`.

The preview uses its own Radiance RGBE writer (flat scanlines; every reader accepts them) instead of vendoring `stb_image_write` — a spec refinement already folded into `docs/spec.md`.

- [ ] **Step 1: Write the failing tests**

`tests/test_preview.cpp`:
```cpp
#include "preview.h"
#include "writer.h"
#include <kre_starscape.h>
#include <doctest.h>
#include <cmath>
#include <string>

using namespace starscape;

namespace
{

constexpr double kPi = 3.14159265358979323846;

std::vector<uint8_t> SkyFile(std::vector<glm::dvec3>&& diffuse, uint32_t faceSize, std::vector<CookedStar>&& stars)
{
	Cooked cooked;
	cooked.faceSize = faceSize;
	cooked.diffuse = std::move(diffuse);
	cooked.stars = std::move(stars);
	cooked.fluxThreshold = 0.01;
	return BuildFile(cooked, 0);
}

}

TEST_CASE("A splatted star keeps its exact flux at every latitude")
{
	for (double const latitude : {0.0, 60.0, 85.0})
	{
		double const b = latitude * kPi / 180.0;
		glm::dvec3 const dir(std::cos(b) * std::cos(0.3), std::cos(b) * std::sin(0.3), std::sin(b));
		std::vector<uint8_t> const bytes = SkyFile(std::vector<glm::dvec3>(6, glm::dvec3(0.0)), 1, {{dir, {1.0, 2.0, 3.0}}});
		kre_starscape::View view;
		REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
		Image const image = RenderPreview(view, 512, 0.0, 1);
		glm::dvec3 total(0.0);
		for (uint32_t row = 0; row < image.height; ++row)
		{
			for (uint32_t column = 0; column < image.width; ++column)
			{
				size_t const i = (size_t(row) * image.width + column) * 3;
				total += glm::dvec3(image.rgb[i], image.rgb[i + 1], image.rgb[i + 2]) * PixelSolidAngle(image.width, row);
			}
		}
		CHECK(total.x == doctest::Approx(1.0).epsilon(0.005));
		CHECK(total.y == doctest::Approx(2.0).epsilon(0.005));
		CHECK(total.z == doctest::Approx(3.0).epsilon(0.005));
	}
}

TEST_CASE("Pixel solid angles cover the sphere")
{
	double total = 0.0;
	for (uint32_t row = 0; row < 128; ++row)
	{
		total += 256.0 * PixelSolidAngle(256, row);
	}
	CHECK(total == doctest::Approx(4.0 * kPi));
	glm::dvec3 const top = PixelDirection(256, 128.0, 0.5);
	CHECK(top.z > 0.999);
}

TEST_CASE("The diffuse layer is sampled without seams")
{
	// Radiance 1 + z everywhere: a smooth field, so any face-edge mistake shows as an outlier.
	constexpr uint32_t size = 8;
	std::vector<glm::dvec3> diffuse;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t y = 0; y < size; ++y)
		{
			for (uint32_t x = 0; x < size; ++x)
			{
				auto const d = kre_starscape::FaceToDir(face, (x + 0.5) / size, (y + 0.5) / size);
				double const z = d[2] / std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
				diffuse.push_back(glm::dvec3(1.0 + z));
			}
		}
	}
	std::vector<uint8_t> const bytes = SkyFile(std::move(diffuse), size, {});
	kre_starscape::View view;
	REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
	Image const image = RenderPreview(view, 256, 0.0, 4);
	for (uint32_t row = 0; row < image.height; ++row)
	{
		for (uint32_t column = 0; column < image.width; ++column)
		{
			double const z = PixelDirection(image.width, column + 0.5, row + 0.5).z;
			REQUIRE(image.rgb[(size_t(row) * image.width + column) * 3] == doctest::Approx(1.0 + z).epsilon(0.08));
		}
	}
}

TEST_CASE("Radiance RGBE encoding")
{
	Image image{2, 1, {1.0f, 0.5f, 0.25f, 0.0f, 0.0f, 0.0f}};
	std::vector<uint8_t> const bytes = EncodeRadiance(image);
	std::string const header = "#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y 1 +X 2\n";
	REQUIRE(bytes.size() == header.size() + 8);
	CHECK(std::string(bytes.begin(), bytes.begin() + std::ptrdiff_t(header.size())) == header);
	CHECK(bytes[header.size() + 0] == 128);
	CHECK(bytes[header.size() + 1] == 64);
	CHECK(bytes[header.size() + 2] == 32);
	CHECK(bytes[header.size() + 3] == 129);
	CHECK(bytes[header.size() + 7] == 0);
}
```

- [ ] **Step 2: Build to verify the tests fail**

Add `src/preview.cpp` and `tests/test_preview.cpp` to the lists. Expected: compile errors.

- [ ] **Step 3: Write `src/preview.h`**

```cpp
#pragma once
#include "march.h"
#include <kre_starscape.h>
#include <cstdint>
#include <vector>
#include <glm/vec3.hpp>

namespace starscape
{

struct Image
{
	uint32_t           width = 0;
	uint32_t           height = 0;
	std::vector<float> rgb;   // row-major, top row first, physical radiance × 2^EV
};

// Equirectangular: column → longitude −π…π left to right, row → latitude +π/2…−π/2 top to bottom.
double PixelSolidAngle(uint32_t width, uint32_t row);
glm::dvec3 PixelDirection(uint32_t width, double column, double row);   // continuous pixel coordinates

// Bilinear and seamless: taps beyond a face edge are fetched from the neighbouring face.
glm::dvec3 SampleDiffuse(kre_starscape::View const& view, glm::dvec3 dir);

Image RenderPreview(kre_starscape::View const& view, uint32_t width, double exposureEV, unsigned threads);
// Marches the diffuse layer at the preview's own resolution (slow: one ray per pixel);
// stars still come from `view`, and so does the flux threshold that splits the two.
Image RenderMarchedPreview(PreparedScene const& scene, MarchParams const& params, kre_starscape::View const& view,
	uint32_t width, double exposureEV, unsigned threads);

std::vector<uint8_t> EncodeRadiance(Image const& image);

}
```

- [ ] **Step 4: Write `src/preview.cpp`**

```cpp
#include "preview.h"
#include "parallel.h"
#include <algorithm>
#include <cmath>
#include <format>
#include <string>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

constexpr double kPi = 3.14159265358979323846;
// FWHM 1.5 px: wide enough that a star never aliases, narrow enough to stay a point.
constexpr double kSplatSigma = 1.5 / 2.3548200450309493;

glm::dvec3 Fetch(kre_starscape::View const& view, uint32_t face, int x, int y)
{
	int const size = int(view.header.diffuseFaceSize);
	if (x < 0 || y < 0 || x >= size || y >= size)
	{
		// Off this face: aim through the would-be texel centre and let the face table pick the neighbour.
		std::array<double, 3> const d = kre_starscape::FaceToDir(face, (x + 0.5) / size, (y + 0.5) / size);
		kre_starscape::FaceCoord const c = kre_starscape::DirToFace(d[0], d[1], d[2]);
		face = c.face;
		x = std::clamp(int(c.s * size), 0, size - 1);
		y = std::clamp(int(c.t * size), 0, size - 1);
	}
	auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::DiffuseTexel(view, 0, face, uint32_t(x), uint32_t(y)));
	return glm::dvec3(v[0], v[1], v[2]) / double(view.header.valueScale);
}

void SplatStars(kre_starscape::View const& view, double exposure, Image& image)
{
	int const width = int(image.width);
	int const height = int(image.height);
	for (uint32_t i = 0; i < view.header.starCount; ++i)
	{
		kre_starscape::Star const star = kre_starscape::StarAt(view, i);
		auto const packed = kre_starscape::DecodeRGB9E5(star.fluxRGB9E5);
		glm::dvec3 const flux = glm::dvec3(packed[0], packed[1], packed[2]) * (exposure / double(view.header.valueScale));
		glm::dvec3 const dir = glm::normalize(glm::dvec3(star.dir[0], star.dir[1], star.dir[2]));
		double const longitude = std::atan2(dir.y, dir.x);
		double const latitude = std::asin(std::clamp(dir.z, -1.0, 1.0));
		// Continuous coordinates where integer + 0.5 is a pixel centre.
		double const u = (longitude + kPi) / (2.0 * kPi) * width;
		double const v = (0.5 * kPi - latitude) / kPi * height;

		struct Tap
		{
			int    column;
			int    row;
			double weight;
		};
		std::array<Tap, 49> taps{};
		int count = 0;
		double weightSum = 0.0;
		int const column0 = int(std::floor(u));
		int const row0 = int(std::floor(v));
		for (int dy = -3; dy <= 3; ++dy)
		{
			int const row = row0 + dy;
			if (row < 0 || row >= height)
			{
				continue;
			}
			for (int dx = -3; dx <= 3; ++dx)
			{
				int const column = column0 + dx;
				double const ex = column + 0.5 - u;
				double const ey = row + 0.5 - v;
				double const weight = std::exp(-(ex * ex + ey * ey) / (2.0 * kSplatSigma * kSplatSigma));
				taps[size_t(count++)] = {((column % width) + width) % width, row, weight};
				weightSum += weight;
			}
		}
		// Normalising by the taps actually kept and dividing by each pixel's solid angle makes
		// the star's integrated flux exact at any latitude, poles included.
		for (int t = 0; t < count; ++t)
		{
			Tap const& tap = taps[size_t(t)];
			glm::dvec3 const radiance = flux * (tap.weight / weightSum / PixelSolidAngle(image.width, uint32_t(tap.row)));
			size_t const index = (size_t(tap.row) * image.width + size_t(tap.column)) * 3;
			image.rgb[index]     += float(radiance.x);
			image.rgb[index + 1] += float(radiance.y);
			image.rgb[index + 2] += float(radiance.z);
		}
	}
}

template<class Radiance>
Image Render(uint32_t width, unsigned threads, Radiance&& radiance)
{
	Image image{width, width / 2, {}};
	image.rgb.assign(size_t(image.width) * image.height * 3, 0.0f);
	ParallelFor(image.height, threads, [&](uint32_t row)
	{
		for (uint32_t column = 0; column < image.width; ++column)
		{
			glm::dvec3 const value = radiance(PixelDirection(image.width, column + 0.5, row + 0.5));
			size_t const index = (size_t(row) * image.width + column) * 3;
			image.rgb[index]     = float(value.x);
			image.rgb[index + 1] = float(value.y);
			image.rgb[index + 2] = float(value.z);
		}
	}, {});
	return image;
}

}

double PixelSolidAngle(uint32_t width, uint32_t row)
{
	uint32_t const height = width / 2;
	double const top = 0.5 * kPi - kPi * row / height;
	double const bottom = 0.5 * kPi - kPi * (row + 1) / height;
	return 2.0 * kPi / width * (std::sin(top) - std::sin(bottom));
}

glm::dvec3 PixelDirection(uint32_t width, double column, double row)
{
	uint32_t const height = width / 2;
	double const longitude = 2.0 * kPi * column / width - kPi;
	double const latitude = 0.5 * kPi - kPi * row / height;
	return {std::cos(latitude) * std::cos(longitude), std::cos(latitude) * std::sin(longitude), std::sin(latitude)};
}

glm::dvec3 SampleDiffuse(kre_starscape::View const& view, glm::dvec3 dir)
{
	int const size = int(view.header.diffuseFaceSize);
	kre_starscape::FaceCoord const c = kre_starscape::DirToFace(dir.x, dir.y, dir.z);
	double const px = c.s * size - 0.5;
	double const py = c.t * size - 0.5;
	int const x0 = int(std::floor(px));
	int const y0 = int(std::floor(py));
	double const fx = px - x0;
	double const fy = py - y0;
	return (1.0 - fx) * (1.0 - fy) * Fetch(view, c.face, x0, y0)
		+ fx * (1.0 - fy) * Fetch(view, c.face, x0 + 1, y0)
		+ (1.0 - fx) * fy * Fetch(view, c.face, x0, y0 + 1)
		+ fx * fy * Fetch(view, c.face, x0 + 1, y0 + 1);
}

Image RenderPreview(kre_starscape::View const& view, uint32_t width, double exposureEV, unsigned threads)
{
	double const exposure = std::exp2(exposureEV);
	Image image = Render(width, threads, [&](glm::dvec3 dir) { return exposure * SampleDiffuse(view, dir); });
	SplatStars(view, exposure, image);
	return image;
}

Image RenderMarchedPreview(PreparedScene const& scene, MarchParams const& params, kre_starscape::View const& view,
	uint32_t width, double exposureEV, unsigned threads)
{
	double const exposure = std::exp2(exposureEV);
	double const threshold = double(view.header.fluxThreshold);
	Image image = Render(width, threads, [&](glm::dvec3 dir)
	{
		RayResult result;
		MarchRay(scene, params, dir, threshold, 0.0, 1.0, result);
		return exposure * result.radiance;
	});
	SplatStars(view, exposure, image);
	return image;
}

std::vector<uint8_t> EncodeRadiance(Image const& image)
{
	std::string const header = std::format("#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y {} +X {}\n", image.height, image.width);
	std::vector<uint8_t> bytes(header.begin(), header.end());
	bytes.reserve(bytes.size() + size_t(image.width) * image.height * 4);
	for (size_t i = 0; i + 2 < image.rgb.size(); i += 3)
	{
		float const r = std::max(image.rgb[i], 0.0f);
		float const g = std::max(image.rgb[i + 1], 0.0f);
		float const b = std::max(image.rgb[i + 2], 0.0f);
		float const v = std::max({r, g, b});
		if (!(v >= 1e-32f))
		{
			bytes.insert(bytes.end(), {0, 0, 0, 0});
			continue;
		}
		int exponent = 0;
		float const mantissa = std::frexp(v, &exponent);
		float const scale = mantissa * 256.0f / v;
		bytes.insert(bytes.end(), {uint8_t(r * scale), uint8_t(g * scale), uint8_t(b * scale), uint8_t(exponent + 128)});
	}
	return bytes;
}

}
```

- [ ] **Step 5: Build and run the tests**

Expected: all preview test cases pass.

- [ ] **Step 6: Run the reddening mutations**

(a) Drop `/ PixelSolidAngle(...)` from the splat (multiply instead) → the flux test fails at 60° and 85°. (b) In `Fetch`, replace the neighbour lookup with clamping to the current face → the seam test fails. If (b) stays green, tighten the seam test by sampling 1 + x + z instead, and re-run. Revert, `touch src/preview.cpp`, rebuild, green.

- [ ] **Step 7: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "Preview: seamless cubemap sampling, flux-exact star splats, Radiance writer

Mutation verified: a missing solid-angle divide and clamp-only seams both redden.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 10: CLI, example scene, README

**Files:**
- Create: `src/cli.cpp`, `scenes/milky-way-from-earth.toml`, `tests/data/tiny.toml`, `tests/data/bad-key.toml`, `tests/expect_exit.cmake`, `README.md`
- Modify: `CMakeLists.txt` (executable target, CLI tests)

**Interfaces:**
- Consumes: `LoadScene`, `Scene`, `Quality` (Task 5); `Cook`, `CookOptions`, `Cooked` (Task 7); `BuildFile`, `WriteFileAtomic`, `ReadFileBytes` (Task 8); `RenderPreview`, `RenderMarchedPreview`, `EncodeRadiance` (Task 9); `PrepareScene`, `ParamsFor` (Task 6); `kre_starscape::Parse`, `Describe`, `StarAt`, `DecodeRGB9E5`, `DiffuseTexel`, `TexelSolidAngle` (Task 2); `Luminance` (Task 3).
- Produces: the `starscape` executable with `cook`, `preview`, `inspect`. Exit codes: 0 success, 1 invalid input or usage, 2 I/O failure.

- [ ] **Step 1: Write the test data and CLI tests**

`tests/data/tiny.toml`:
```toml
[output]
face_size = 16
target_stars = 300
quality = "quick"
preview = { width = 256 }

[observer]
position = [8200.0, 0.0, 20.0]

[[galaxy]]
preset = "milky-way"
```

`tests/data/bad-key.toml`:
```toml
[output]
face_size = 16
target_stars = 300
qualty = "quick"

[observer]
position = [0.0, 0.0, 0.0]

[[galaxy]]
preset = "dwarf"
```

`tests/expect_exit.cmake`:
```cmake
# Runs EXE with ARG1..ARG6 and fails unless it exits with EXPECT.
set(arguments "")
foreach(i RANGE 1 6)
	if(DEFINED ARG${i})
		list(APPEND arguments "${ARG${i}}")
	endif()
endforeach()
execute_process(COMMAND "${EXE}" ${arguments} RESULT_VARIABLE code)
if(NOT code EQUAL EXPECT)
	message(FATAL_ERROR "expected exit ${EXPECT}, got ${code}")
endif()
```

Append to `CMakeLists.txt`:
```cmake
add_executable(starscape src/cli.cpp)
target_link_libraries(starscape PRIVATE starscape-core)

set(CLI_OUT ${CMAKE_CURRENT_BINARY_DIR}/cli-test)
file(MAKE_DIRECTORY ${CLI_OUT})
file(WRITE ${CLI_OUT}/garbage.kre_starscape "this is not a sky")
function(starscape_cli_test name expect)
	set(args "")
	set(i 1)
	foreach(arg IN LISTS ARGN)
		list(APPEND args "-DARG${i}=${arg}")
		math(EXPR i "${i} + 1")
	endforeach()
	add_test(NAME ${name} COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:starscape> -DEXPECT=${expect} ${args} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
endfunction()
starscape_cli_test(cli-cook 0 cook ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/tiny.toml -o ${CLI_OUT}/tiny.kre_starscape --preview ${CLI_OUT}/tiny.hdr)
set_tests_properties(cli-cook PROPERTIES FIXTURES_SETUP tiny-sky)
starscape_cli_test(cli-inspect 0 inspect ${CLI_OUT}/tiny.kre_starscape)
starscape_cli_test(cli-preview 0 preview ${CLI_OUT}/tiny.kre_starscape -o ${CLI_OUT}/again.hdr --exposure 2)
set_tests_properties(cli-inspect cli-preview PROPERTIES FIXTURES_REQUIRED tiny-sky)
starscape_cli_test(cli-bad-scene 1 cook ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/bad-key.toml -o ${CLI_OUT}/bad.kre_starscape)
starscape_cli_test(cli-missing-scene 2 cook ${CLI_OUT}/absent.toml -o ${CLI_OUT}/absent.kre_starscape)
starscape_cli_test(cli-garbage-file 1 inspect ${CLI_OUT}/garbage.kre_starscape)
starscape_cli_test(cli-missing-file 2 inspect ${CLI_OUT}/absent.kre_starscape)
starscape_cli_test(cli-usage 1 frobnicate)
```

- [ ] **Step 2: Build to verify the tests fail**

Expected: CMake error, `src/cli.cpp` missing.

- [ ] **Step 3: Write `src/cli.cpp`**

```cpp
#include "cook.h"
#include "march.h"
#include "preview.h"
#include "scene.h"
#include "writer.h"
#include <kre_starscape.h>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

using namespace starscape;

namespace
{

constexpr char const* kUsage =
	"usage:\n"
	"  starscape cook    <scene.toml> -o <out.kre_starscape> [--preview <out.hdr>] [--threads N] [--quick]\n"
	"  starscape preview <in.kre_starscape> -o <out.hdr> [--width N] [--exposure EV] [--threads N]\n"
	"  starscape inspect <in.kre_starscape>\n";

auto const gStart = std::chrono::steady_clock::now();

void Log(std::string_view message)
{
	auto const elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - gStart).count();
	std::printf("[+%02lld:%02lld:%02lld] %.*s\n", (long long)(elapsed / 3600), (long long)(elapsed / 60 % 60), (long long)(elapsed % 60),
		int(message.size()), message.data());
	std::fflush(stdout);
}

int Report(Error const& error)
{
	std::fprintf(stderr, "starscape: %s\n", error.message.c_str());
	return error.kind == Error::Kind::Io ? 2 : 1;
}

int Usage(std::string_view problem)
{
	std::fprintf(stderr, "starscape: %.*s\n%s", int(problem.size()), problem.data(), kUsage);
	return 1;
}

struct Arguments
{
	std::vector<std::string_view> positional;
	std::optional<std::string_view> output;
	std::optional<std::string_view> preview;
	unsigned threads = std::max(1u, std::thread::hardware_concurrency());
	bool quick = false;
	std::optional<uint32_t> width;
	double exposure = 0.0;
};

template<class T>
bool ParseNumber(std::string_view text, T& out)
{
	char const* end = text.data() + text.size();
	auto const [pointer, error] = std::from_chars(text.data(), end, out);
	return error == std::errc() && pointer == end;
}

std::optional<std::string> ParseArguments(int argc, char** argv, Arguments& out)
{
	for (int i = 2; i < argc; ++i)
	{
		std::string_view const arg = argv[i];
		auto value = [&]() -> std::optional<std::string_view>
		{
			if (i + 1 >= argc)
			{
				return std::nullopt;
			}
			return std::string_view(argv[++i]);
		};
		if (arg == "--quick")
		{
			out.quick = true;
		}
		else if (arg == "-o" || arg == "--preview" || arg == "--threads" || arg == "--width" || arg == "--exposure")
		{
			std::optional<std::string_view> const v = value();
			if (!v)
			{
				return std::format("{} needs a value", arg);
			}
			if (arg == "-o")
			{
				out.output = *v;
			}
			else if (arg == "--preview")
			{
				out.preview = *v;
			}
			else if (arg == "--threads")
			{
				if (!ParseNumber(*v, out.threads) || out.threads == 0)
				{
					return "--threads needs a positive integer";
				}
			}
			else if (arg == "--width")
			{
				uint32_t width = 0;
				if (!ParseNumber(*v, width) || width < 16 || width % 2 != 0)
				{
					return "--width needs an even integer ≥ 16";
				}
				out.width = width;
			}
			else
			{
				if (!ParseNumber(*v, out.exposure) || !std::isfinite(out.exposure))
				{
					return "--exposure needs a number";
				}
			}
		}
		else if (arg.starts_with("-"))
		{
			return std::format("unknown option {}", arg);
		}
		else
		{
			out.positional.push_back(arg);
		}
	}
	return std::nullopt;
}

Expected<std::vector<uint8_t>> LoadSky(std::string_view path, kre_starscape::View& view)
{
	STARSCAPE_TRY_ASSIGN(bytes, ReadFileBytes(std::string(path)));
	// `view` points into `bytes`, so the caller keeps the returned buffer alive while using it.
	kre_starscape::Result const result = kre_starscape::Parse(bytes, view);
	if (result != kre_starscape::Result::Ok)
	{
		return InputError(std::format("{}: {}", path, kre_starscape::Describe(result)));
	}
	return bytes;
}

int Cook(Arguments const& args)
{
	if (args.positional.size() != 1 || !args.output)
	{
		return Usage("cook needs one scene file and -o");
	}
	Expected<Scene> const scene = LoadScene(std::string(args.positional[0]));
	if (!scene)
	{
		return Report(scene.error());
	}
	Log(std::format("cooking {} on {} threads", args.positional[0], args.threads));
	int lastPercent = -1;
	std::string lastPhase;
	CookOptions options;
	options.threads = args.threads;
	if (args.quick)
	{
		options.quality = Quality::Quick;
	}
	auto const marchStart = std::chrono::steady_clock::now();
	options.progress = [&](std::string_view phase, double done)
	{
		int const percent = int(done * 100.0);
		if (phase == lastPhase && percent == lastPercent)
		{
			return;
		}
		lastPhase = std::string(phase);
		lastPercent = percent;
		double const seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - marchStart).count();
		std::string eta;
		if (phase == "march" && done > 0.01)
		{
			eta = std::format(", ETA {:.0f} s", seconds / done * (1.0 - done));
		}
		Log(std::format("{} {}%{}", phase, percent, eta));
	};
	Cooked const cooked = starscape::Cook(*scene, options);
	std::vector<uint8_t> const bytes = BuildFile(cooked, scene->hash);
	if (Expected<void> const written = WriteFileAtomic(std::string(*args.output), bytes); !written)
	{
		return Report(written.error());
	}

	double const total = cooked.stats.pointFlux + cooked.stats.diffuseFlux;
	Log(std::format("wrote {} ({} bytes)", *args.output, bytes.size()));
	Log(std::format("threshold F_res = {:.4g} (magnitude {:.2f}); {} stars; points carry {:.2f}% of the light",
		cooked.fluxThreshold, -2.5 * std::log10(cooked.fluxThreshold), cooked.stats.starCount,
		total > 0.0 ? 100.0 * cooked.stats.pointFlux / total : 0.0));
	if (cooked.stats.budgetUnreachable)
	{
		Log("warning: the scene holds fewer stars than target_stars; every star is a point");
	}
	if (cooked.stats.droppedStars > 0)
	{
		Log(std::format("warning: {} drawn stars found nothing resolved on re-march and were dropped", cooked.stats.droppedStars));
	}

	if (args.preview)
	{
		kre_starscape::View view;
		kre_starscape::Parse(bytes, view);
		Image const image = scene->output.previewMarch
			? RenderMarchedPreview(PrepareScene(*scene), ParamsFor(options.quality.value_or(scene->output.quality)), view,
				scene->output.previewWidth, 0.0, args.threads)
			: RenderPreview(view, scene->output.previewWidth, 0.0, args.threads);
		if (Expected<void> const written = WriteFileAtomic(std::string(*args.preview), EncodeRadiance(image)); !written)
		{
			return Report(written.error());
		}
		Log(std::format("wrote preview {} ({}×{})", *args.preview, image.width, image.height));
	}
	return 0;
}

int Preview(Arguments const& args)
{
	if (args.positional.size() != 1 || !args.output)
	{
		return Usage("preview needs one .kre_starscape file and -o");
	}
	kre_starscape::View view;
	Expected<std::vector<uint8_t>> const bytes = LoadSky(args.positional[0], view);
	if (!bytes)
	{
		return Report(bytes.error());
	}
	Image const image = RenderPreview(view, args.width.value_or(8192), args.exposure, args.threads);
	if (Expected<void> const written = WriteFileAtomic(std::string(*args.output), EncodeRadiance(image)); !written)
	{
		return Report(written.error());
	}
	Log(std::format("wrote preview {} ({}×{})", *args.output, image.width, image.height));
	return 0;
}

int Inspect(Arguments const& args)
{
	if (args.positional.size() != 1)
	{
		return Usage("inspect needs one .kre_starscape file");
	}
	kre_starscape::View view;
	Expected<std::vector<uint8_t>> const bytes = LoadSky(args.positional[0], view);
	if (!bytes)
	{
		return Report(bytes.error());
	}
	kre_starscape::Header const& h = view.header;
	double const scale = h.valueScale;
	std::printf("scene hash      %016llx\n", (unsigned long long)h.sceneHash);
	std::printf("diffuse         RGB9E5 cubemap, %u² per face, %u mips, %llu bytes\n", h.diffuseFaceSize, h.diffuseMipCount, (unsigned long long)h.diffuseBytes);
	std::printf("stars           %u\n", h.starCount);
	std::printf("flux threshold  %.4g (magnitude %.2f)\n", h.fluxThreshold, -2.5 * std::log10(h.fluxThreshold));
	std::printf("value scale     %g\n", scale);

	double diffuse = 0.0;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t y = 0; y < h.diffuseFaceSize; ++y)
		{
			for (uint32_t x = 0; x < h.diffuseFaceSize; ++x)
			{
				auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::DiffuseTexel(view, 0, face, x, y));
				diffuse += Luminance(glm::dvec3(v[0], v[1], v[2])) / scale * kre_starscape::TexelSolidAngle(h.diffuseFaceSize, x, y);
			}
		}
	}
	double points = 0.0;
	std::array<uint32_t, 40> histogram{};   // magnitude bins from −10 to +30
	for (uint32_t i = 0; i < h.starCount; ++i)
	{
		auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::StarAt(view, i).fluxRGB9E5);
		double const flux = Luminance(glm::dvec3(v[0], v[1], v[2])) / scale;
		points += flux;
		if (flux > 0.0)
		{
			int const bin = std::clamp(int(std::floor(-2.5 * std::log10(flux))) + 10, 0, 39);
			++histogram[size_t(bin)];
		}
	}
	std::printf("total light     %.6g (points %.2f%%, diffuse %.2f%%)\n", points + diffuse,
		100.0 * points / std::max(points + diffuse, 1e-300), 100.0 * diffuse / std::max(points + diffuse, 1e-300));
	std::printf("stars by magnitude:\n");
	for (size_t bin = 0; bin < histogram.size(); ++bin)
	{
		if (histogram[bin] > 0)
		{
			std::printf("  %+3d .. %+3d  %u\n", int(bin) - 10, int(bin) - 9, histogram[bin]);
		}
	}
	return 0;
}

}

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		return Usage("no command");
	}
	Arguments args;
	if (std::optional<std::string> const problem = ParseArguments(argc, argv, args))
	{
		return Usage(*problem);
	}
	std::string_view const command = argv[1];
	if (command == "cook")
	{
		return Cook(args);
	}
	if (command == "preview")
	{
		return Preview(args);
	}
	if (command == "inspect")
	{
		return Inspect(args);
	}
	return Usage(std::format("unknown command {}", command));
}
```

Note the `tiny.toml` preview: the CLI's `cook` uses `scene.output.previewWidth` (256 there).

- [ ] **Step 4: Write `scenes/milky-way-from-earth.toml`**

```toml
# A Milky-Way-like galaxy seen from a Sun-like position: 8.2 kpc out, 20 pc above the plane,
# with Earth's spin axis (celestial north) tilted ~62.6° from the galactic pole.
[output]
face_size    = 256
target_stars = 4000
quality      = "final"
preview      = { width = 8192, march = false }

[observer]
position      = [8200.0, 0.0, 20.0]
spin_axis     = [0.4838, -0.7471, 0.4560]
zero_meridian = [0.0549, -0.4940, -0.8677]

[[galaxy]]
preset = "milky-way"
seed   = 1

[[galaxy]]
preset   = "small-spiral"
position = [780000.0, 120000.0, -40000.0]
rotation = [77.0, 0.0, 35.0]
seed     = 2
```

(The spin axis is the ICRS north pole (l = 122.93°, b = 27.13°) and the meridian the vernal equinox (l = 96.34°, b = −60.19°). Galactic coordinates have +x toward the centre; here the centre lies at −x from the observer, so with z shared, x and y flip sign. The cook re-orthogonalises the meridian, so rounding is harmless.)

- [ ] **Step 5: Write `README.md`**

````markdown
# Starscape

Cooks a night sky for any vantage point in a scene of procedural galaxies: a diffuse HDR
cubemap of the light you cannot resolve into stars, plus the few thousand stars you can.
Output is a `.kre_starscape` file (format: [`include/kre_starscape.h`](include/kre_starscape.h),
header-only, standard library only) and an equirectangular Radiance `.hdr` preview.

## Build

```
cmake -S . -B Build -G Ninja
cmake --build Build
ctest --test-dir Build --output-on-failure
```

C++23 (GCC 13+). Dependencies are vendored in `third_party/` (see `third_party/VERSIONS`).

## Use

```
starscape cook scenes/milky-way-from-earth.toml -o sky.kre_starscape --preview sky.hdr
starscape cook scene.toml -o quick.kre_starscape --quick     # 1 ray per texel, coarse steps
starscape preview sky.kre_starscape -o bright.hdr --exposure 4
starscape inspect sky.kre_starscape
```

A scene is TOML: `[output]` (face size, `target_stars` or `limiting_magnitude`, quality),
`[observer]` (position in parsecs, spin axis, zero meridian) and one or more `[[galaxy]]`
tables starting from a preset (`milky-way`, `small-spiral`, `elliptical`, `dwarf`) with any
field overridden. Unknown keys are errors. See [`docs/spec.md`](docs/spec.md) for every key,
formula and unit.

## Units

Parsecs; visual luminosity in solar units; star flux in units of a magnitude-0 star
(flux = 10^(−0.4 m)); diffuse radiance in magnitude-0 flux per steradian; linear Rec.709.
Stored values are multiplied by the header's `valueScale`; divide by it to read them.

## How the light is split

One apparent-flux threshold decides, at every step of every ray, whether a luminosity bin's
stars are bright enough to be points. Points and diffuse light add up to the whole sky, so
asking for fewer stars moves light into the texture instead of losing it.

## References

- O. Elfwering, *Procedurally Generating an Artificial Galaxy*, Lund University, 2016 — the starting point.
- L. Hernquist, ApJ 356 (1990) — bulge profile.
- P. Kroupa, MNRAS 322 (2001) — initial mass function.
- J. A. Cardelli, G. C. Clayton, J. S. Mathis, ApJ 345 (1989) — extinction curve.
- P. C. van der Kruit, L. Searle, A&A 95 (1981) — disk profile.
- C. Wyman, P.-P. Sloan, P. Shirley, JCGT 2(2) (2013) — analytic CIE colour matching functions.
- J. K. Salmon et al., SC11 (2011) — Philox.

## License

MIT.
````

- [ ] **Step 6: Build and run every test**

Run build, then `ctest --test-dir /mnt/Passport/Worktrees/starscape/Build --output-on-failure`.
Expected: `unit` and all eight `cli-*` tests pass.

- [ ] **Step 7: Run the reddening mutation**

Make `Report` return 1 for every error kind → `cli-missing-scene` and `cli-missing-file` fail. Revert, `touch src/cli.cpp`, rebuild, green.

- [ ] **Step 8: Commit**

```bash
git -C /mnt/Passport/Worktrees/starscape add -A && git -C /mnt/Passport/Worktrees/starscape commit -m "CLI (cook, preview, inspect), example scene, README

Mutation verified: collapsing exit codes reddens the I/O-failure CLI tests.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 11: First real cook, then Kreatures integration (owner-gated)

**Files:**
- Modify (Kreatures superproject, `/mnt/Passport/Kreatures-Engine`): `.gitmodules`, `Tooling/Starscape` gitlink, `docs/engine-roadmap.md`

**Interfaces:**
- Consumes: the finished `starscape` binary.
- Produces: a cooked example sky with numbers for the owner; the submodule and roadmap row.

- [ ] **Step 1: Full test run**

Run `ctest --test-dir /mnt/Passport/Worktrees/starscape/Build --output-on-failure`. Record the pass/fail counts. All must pass.

- [ ] **Step 2: Cook the example scene at final quality**

Run in the background (it may take tens of minutes) and keep the timestamped log:
```bash
nice -n 19 /mnt/Passport/Worktrees/starscape/Build/starscape cook /mnt/Passport/Worktrees/starscape/scenes/milky-way-from-earth.toml -o /mnt/Passport/Worktrees/starscape/Build/milky-way-from-earth.kre_starscape --preview /mnt/Passport/Worktrees/starscape/Build/milky-way-from-earth.hdr > /mnt/Passport/Worktrees/starscape/Build/milky-way-from-earth.log 2>&1
```
Then `starscape inspect` on the result. Report to the owner: wall time, F_res and its magnitude, star count, the points/diffuse split, the magnitude histogram, dropped stars, file size. The preview is **UNVERIFIED-VISUAL**: hand the owner the `.hdr` path; no gate judges it.

- [ ] **Step 3: Owner gate — publish**

Ask the owner before pushing: the repository is public. On approval:
```bash
git -C /mnt/Passport/Worktrees/starscape push origin starscape-v1:main
```
Verify with `git -C /mnt/Passport/Worktrees/starscape status -sb` that the branch is not `ahead`.

- [ ] **Step 4: Owner gate — register the submodule**

The Starscape clone already sits at `Tooling/Starscape` as an untracked directory. With the owner's approval and after the push (a gitlink must resolve on the remote):
```bash
git -C /mnt/Passport/Kreatures-Engine submodule add git@github.com:SpehleonLP/Starscape.git Tooling/Starscape
```
(`submodule add` adopts an existing clone of the same URL.) Then add a row to the Gate L table in `docs/engine-roadmap.md` (columns: Remaining outcome | Status | Scope / owner):

```
| Night sky from a cooked `.kre_starscape` | tool landed — engine side open | [Starscape](../Tooling/Starscape/docs/spec.md) cooks the diffuse cubemap + point stars. Remaining: a spec for the engine reader replacing the glTF path in `src/World/Atmosphere/starscape.cpp`, the Gaussian-PSF star shader and diffuse-cubemap draw, and deleting `Shaders/Starscape.glb`. Later: BC6H diffuse encoding. |
```

Commit in the superproject:
```bash
git -C /mnt/Passport/Kreatures-Engine add .gitmodules Tooling/Starscape docs/engine-roadmap.md && git -C /mnt/Passport/Kreatures-Engine commit -m "Add Tooling/Starscape: offline night-sky cook

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Pushing the superproject is a separate owner decision (`scripts/push-tree.sh --dry-run` first).
