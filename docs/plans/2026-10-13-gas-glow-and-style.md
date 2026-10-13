# Gas Glow and Stylized Sky Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a gas-glow emission term to disks, a scene-wide star colour saturation setting, and a stylized example scene that reproduces the owner-approved target render, with both example scenes in the galactic frame.

**Architecture:**
- **Gas** is computed in `EvaluateLocal` from the thin disk's young-star density, two noise fields and the dust's own filament and bubble fields. It is returned as an RGB emissivity on `GalaxySample`. `Walk` sums it onto each `Step` and `MarchRay` adds it to the diffuse radiance behind the step's transmittance.
- **Saturation** is applied in `Tabulate` to each bin's chromaticity, renormalised to luminance 1.
- **Scenes** are data: TOML files under `scenes/`.

**Tech Stack:** C++23 (GCC 13), CMake ≥ 3.25 + Ninja, vendored glm, toml++ and doctest.

**Spec:** `docs/specs/2026-10-13-gas-glow-and-style.md`. Read it first; it is the authority. The base model is `docs/spec.md`.

## Global Constraints

- **Repository:** `/mnt/Passport/Worktrees/starscape`, branch `starscape-v1`. Paths below are relative to it.
- **Start from the committed tree.** The working tree carries an uncommitted look-dev prototype (nebula spheres, env-var knobs). Task 1 stashes it. None of it lands:
  - `src/nebula.*`;
  - `BillowFbm`, `Chord`, `NebulaDepth`, `Step::dustDepth`;
  - every `NEB_*`, `GAS_GAIN` and `CHROMA` read.
- **Style:**
  - C++23, `-Wall -Wextra -Wpedantic` clean;
  - tabs and Allman braces, with **every** control-flow body braced;
  - comments say *why*, not *what*.
- **Parameters:**
  - Containers and structs larger than a couple of registers are passed `const&` or `&&`.
  - `glm` vectors and spans go by value.
- **Errors** are returned as values (`Expected<T>`). No exceptions, no asserts on input.
- **Units:**
  - parsecs;
  - L_V☉;
  - magnitude-0 flux;
  - `kFluxAt10pc` (population.h) is 10^(−0.4·4.83).
- **Gas noise constants, exact:**
  - clump wavelength 40 pc, wisp wavelength 25 pc, 5 octaves each;
  - salts `0x6A5` (clumps) and `0x6A6` (wisps), XORed with `galaxy.seed`;
  - clumps `LogNormal(1.2, Fbm/kFbmSigma)`;
  - wisps `6 · max(0, RidgedFbm)⁶`;
  - texture `(0.6·clumps + wisps) · sqrt(F/⟨F⟩)`;
  - colour mix `smoothstep(1, 3, clumps)`.

  These match the prototype whose ⟨G⟩ was measured (0.6045 at `filament` 3.5, seed 1), so `glow = 1.5` reproduces the target.
- **Build** (Starscape is outside the engine's jobserver wrapper; keep it low priority):
  - build: `nice -n 19 ionice -c 3 cmake --build /mnt/Passport/Worktrees/starscape/Build -j 4`
  - one test: `/mnt/Passport/Worktrees/starscape/Build/starscape-tests --test-case='<name>'`
  - all: `ctest --test-dir /mnt/Passport/Worktrees/starscape/Build --output-on-failure`
- **Mutations:** every test that names a mutation must be shown red under that mutation once. Then revert the mutation and `touch` the file so Ninja rebuilds it. Record "mutation verified: <what>" in the commit message.
- **Physics tests are never weakened to pass.** If an existing physics test goes red, stop and report the numbers.
- **Do not push.** Pushing Starscape is owner-gated.

## Review Focus

1. **A scene sets `glow` on a disk with no `arm_population`.** Expect a clear error naming `arm_population`, not a silently dark sky. Task 2 test "Gas keys are read and checked".
2. **A colour of all zeros or with a negative component.** Expect an error, never a division by zero in normalisation. Task 2, same test.
3. **The observer sits inside the dust bubble while young stars surround it.** Expect no gas glow from inside the cavity. Task 2 test "Gas shares the dust's bubble".
4. **Saturation large enough to push a channel below zero.** Expect a clamp and renormalisation, never a negative or non-unit-luminance chroma. Task 1 test "Saturation pushes star colours from white and keeps their light" at s = 4, which must find at least one clipped bin.
5. **Gas on a galaxy with dust disabled.** The filament factor must be 1, not a stale or zero `filamentMean`. Task 3 test "Gas glow is dimmed by the dust in front of it" runs its source with dust disabled.

---

### Task 1: Clean start and star colour saturation

**Files:**
- Modify: `src/scene.h`, `src/scene.cpp`, `src/population.h`, `src/population.cpp`, `src/march.cpp`
- Test: `tests/test_population.cpp`, `tests/test_scene.cpp`, `tests/test_march.cpp`

**Interfaces:**
- Produces:
  - `struct Style { double starSaturation = 1.0; };`, a member `Style style;` of `Scene`;
  - `PopulationTable Tabulate(Population const& population, double saturation = 1.0);`
  - `PrepareScene` passes `scene.style.starSaturation`.

- [ ] **Step 1: Stash the prototype and confirm a clean baseline**

```bash
cd /mnt/Passport/Worktrees/starscape
git stash push -u -m "nebula/gas look-dev prototype (2026-10-13): do not land"
git status --short          # expect: nothing
nice -n 19 ionice -c 3 cmake --build Build -j 4
ctest --test-dir Build --output-on-failure
```
Expected: the tree is clean and every test passes. Record the passing count in the ledger.

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_population.cpp`:

```cpp
TEST_CASE("Saturation pushes star colours from white and keeps their light")
{
	// Mutation: skipping the renormalisation in Saturate reddens the luminance check on the clipped bins.
	double const s = 4.0;
	int clipped = 0;
	for (Population const& population : PresetPopulations())
	{
		PopulationTable const physical = Tabulate(population);
		PopulationTable const vivid = Tabulate(population, s);
		for (int bin = 0; bin < kBins; ++bin)
		{
			size_t const i = size_t(bin);
			CHECK(Luminance(vivid.chroma[i]) == doctest::Approx(1.0).epsilon(1e-12));
			glm::dvec3 const pushed = glm::dvec3(1.0) + s * (physical.chroma[i] - glm::dvec3(1.0));
			if (pushed.x >= 0.0 && pushed.y >= 0.0 && pushed.z >= 0.0)
			{
				CHECK(glm::length(vivid.chroma[i] - glm::dvec3(1.0)) == doctest::Approx(s * glm::length(physical.chroma[i] - glm::dvec3(1.0))).epsilon(1e-9));
			}
			else
			{
				++clipped;
				CHECK(vivid.chroma[i].x >= 0.0);
				CHECK(vivid.chroma[i].y >= 0.0);
				CHECK(vivid.chroma[i].z >= 0.0);
			}
			CHECK(Luminance(vivid.lightBelow[i + 1]) == doctest::Approx(Luminance(physical.lightBelow[i + 1])).epsilon(1e-12));
		}
	}
	CHECK(clipped > 0);   // the renormalisation is exercised, not just the in-gamut identity
}
```

Add `#include <glm/geometric.hpp>` to its includes.

Append to `tests/test_scene.cpp`:

```cpp
TEST_CASE("Style is read and range-checked")
{
	// Mutation: a default star_saturation other than 1 reddens the first check.
	std::string const scene = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n";
	CHECK(ParseScene(scene, "test.toml")->style.starSaturation == 1.0);
	CHECK(ParseScene(scene + "[style]\nstar_saturation = 1.5\n", "test.toml")->style.starSaturation == 1.5);
	CHECK(ErrorOf(scene + "[style]\nstar_saturation = 5.0\n").find("style.star_saturation") != std::string::npos);
	CHECK(ErrorOf(scene + "[style]\nstar_saturation = -1.0\n").find("must be in") != std::string::npos);
	CHECK(ErrorOf(scene + "[style]\nsaturation = 1.5\n").find("unknown key") != std::string::npos);
}
```

Append to `tests/test_march.cpp`:

```cpp
TEST_CASE("Star saturation moves colour, not light or counts")
{
	// Mutation: PrepareScene not passing starSaturation to Tabulate reddens the spread check.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	scene.galaxies[0].dust.enabled = false;   // grey transmittance: only then is luminance exactly kept
	PreparedScene const physical = PrepareScene(scene);
	scene.style.starSaturation = 4.0;
	PreparedScene const vivid = PrepareScene(scene);
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir = glm::normalize(glm::dvec3(-1.0, 0.1, 0.02));
	RayResult a;
	RayResult b;
	MarchRay(physical, params, dir, 1e-4, 1.0, 1.0, a);
	MarchRay(vivid, params, dir, 1e-4, 1.0, 1.0, b);
	CHECK(Luminance(b.radiance) == doctest::Approx(Luminance(a.radiance)).epsilon(1e-12));
	CHECK(b.budget == a.budget);
	auto spread = [](glm::dvec3 c) { return (std::max({c.x, c.y, c.z}) - std::min({c.x, c.y, c.z})) / Luminance(c); };
	CHECK(spread(b.radiance) > 2.0 * spread(a.radiance));
}
```

Add `#include <algorithm>` and `#include <glm/geometric.hpp>` to `tests/test_march.cpp` if absent.

- [ ] **Step 3: Run the tests to see them fail**

Run: `nice -n 19 ionice -c 3 cmake --build Build -j 4`
Expected: compile errors. `Tabulate` takes one argument and `Scene` has no `style` member.

- [ ] **Step 4: Implement**

In `src/scene.h`, after `struct Observer`:

```cpp
// Scene-wide departures from the physical sky, for a stylised look.
struct Style
{
	double starSaturation = 1.0;   // 1: physical star colours; above 1 pushes them away from white
};
```

In `struct Scene`, after `Observer observer;`:

```cpp
	Style                   style;
```

In `src/scene.cpp`, after `ReadObserver`:

```cpp
Expected<void> ReadStyle(TableReader& root, Style& style)
{
	STARSCAPE_TRY_ASSIGN(table, root.Table("style"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(root, "style", *table));
	STARSCAPE_TRY(reader.Number("star_saturation", style.starSaturation, 0.0, 4.0));
	return reader.Finish();
}
```

In `ParseScene`, after `STARSCAPE_TRY(ReadObserver(reader, scene.observer));`:

```cpp
	STARSCAPE_TRY(ReadStyle(reader, scene.style));
```

In `src/population.h`, replace the `Tabulate` declaration:

```cpp
// `saturation` scales each bin's colour away from white (1: physical); every bin keeps its luminance.
PopulationTable Tabulate(Population const& population, double saturation = 1.0);
```

In `src/population.cpp`, add an anonymous-namespace helper above `Tabulate`:

```cpp
namespace
{

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
```

Then change `Tabulate`'s signature to `PopulationTable Tabulate(Population const& population, double saturation)`, and its chroma line to:

```cpp
		table.chroma[i]     = Saturate(BlackbodyChroma(population.temperature[i]), saturation);
```

Add `#include <glm/common.hpp>` to `population.cpp` if absent (for `glm::max`). If `population.cpp` already has an anonymous namespace, put `Saturate` inside it instead of opening a second one.

In `src/march.cpp` `PrepareScene`, change the populations loop body to:

```cpp
		prepared.populations.push_back(Tabulate(population, scene.style.starSaturation));
```

- [ ] **Step 5: Run the tests**

Run: `nice -n 19 ionice -c 3 cmake --build Build -j 4 && ctest --test-dir Build --output-on-failure`
Expected: every test passes, including the three new ones.

- [ ] **Step 6: Verify the mutations**

Apply each mutation, rebuild, run its test case and see it red. Then revert and `touch` the file:
- In `Saturate`, return `pushed` without dividing. Run `--test-case='Saturation pushes star colours from white and keeps their light'`; expect it red.
- In `PrepareScene`, call `Tabulate(population)`. Run `--test-case='Star saturation moves colour, not light or counts'`; expect it red.
- In `Style`, set the default to `1.5`. Run `--test-case='Style is read and range-checked'`; expect it red.

- [ ] **Step 7: Commit**

```bash
git add src/scene.h src/scene.cpp src/population.h src/population.cpp src/march.cpp tests/test_population.cpp tests/test_scene.cpp tests/test_march.cpp
git commit -m "Style: star colour saturation that keeps every bin's light

mutation verified: no renormalisation; PrepareScene ignoring saturation; default saturation 1.5"
```

---

### Task 2: The gas model in the galaxy, and its scene keys

**Files:**
- Modify: `src/galaxy.h`, `src/galaxy.cpp`, `src/scene.cpp`
- Test: `tests/test_galaxy.cpp`, `tests/test_scene.cpp`

**Interfaces:**
- Consumes: `Luminance` and `kFluxAt10pc` from `population.h`.
- Produces:
  - `struct Gas { double glow; glm::dvec3 colour; glm::dvec3 coreColour; double textureMean; double youngLight; };`, a member `Gas gas;` at the end of `Disk`;
  - `GalaxySample::gas`, a `glm::dvec3` emissivity in radiance per parsec before transmittance;
  - `PrepareGalaxy` normalises `gas.colour` and `gas.coreColour` to luminance 1 and measures `gas.textureMean`;
  - `gas.youngLight` is set by `PrepareScene` in Task 3. Tests set it directly.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_galaxy.cpp`, and add `#include "population.h"` to its includes:

```cpp
namespace
{

// A disk uniform over any test box, all young stars lighting gas; dust off unless a test turns it on.
Galaxy GasDisk(double glow)
{
	Galaxy galaxy;
	galaxy.disk = {.enabled = true, .density = 1.0, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .arms = {}, .population = 1, .armPopulation = 0,
		.youngFraction = 0.01};
	galaxy.disk.gas.glow = glow;
	galaxy.disk.gas.youngLight = 1.0;   // PrepareScene fills this from the population table
	return galaxy;
}

}

TEST_CASE("Gas glows at glow times the young stars' light, on average")
{
	// Mutation: dropping the division by textureMean leaves the ratio at ⟨G⟩ ≈ 0.88 of glow at filament 2.
	Galaxy galaxy = GasDisk(2.0);
	// Dust present only for its filaments; scale height 400 gives the Milky Way's 200 pc filament wavelength.
	galaxy.dust = {.enabled = true, .extinction = 0.001, .scaleLength = 1e12, .scaleHeight = 400.0, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 2.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	RandomStream rng(11, Purpose::Test, 0, 0, 0);
	double gas = 0.0;
	double young = 0.0;
	int const count = 200000;
	for (int i = 0; i < count; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		Evaluate(prepared, p, sample);
		REQUIRE(sample.count == 2);
		gas += Luminance(sample.gas);
		young += sample.parts[1].density * 100.0 * kFluxAt10pc;
	}
	CHECK(gas / young == doctest::Approx(2.0).epsilon(0.03));
}

TEST_CASE("Gas shares the dust's bubble")
{
	// Mutation: leaving the bubble factor off the gas lights the cavity.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.dust = {.enabled = true, .extinction = 0.002, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0,
		.bubbleCenter = {100.0, 0.0, 0.0}, .bubbleRadius = 300.0, .bubbleHeight = 800.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	GalaxySample inside;
	Evaluate(prepared, {100.0, 0.0, 0.0}, inside);
	REQUIRE(inside.count == 2);   // young stars are there to light it
	CHECK(inside.gas == glm::dvec3(0.0));
	GalaxySample outside;
	Evaluate(prepared, {100.0, 0.0, 900.0}, outside);
	CHECK(Luminance(outside.gas) > 0.0);
}

TEST_CASE("Gas colours are normalised and the densest clumps take the core colour")
{
	// Mutation: swapping the mix's endpoints makes pure-core samples the majority and reddens the ratio check.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.disk.gas.colour = {2.0, 0.0, 0.0};
	galaxy.disk.gas.coreColour = {0.0, 0.0, 3.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	CHECK(Luminance(prepared.galaxy.disk.gas.colour) == doctest::Approx(1.0));
	CHECK(Luminance(prepared.galaxy.disk.gas.coreColour) == doctest::Approx(1.0));
	CHECK(prepared.galaxy.disk.gas.colour.y == 0.0);
	CHECK(prepared.galaxy.disk.gas.coreColour.x == 0.0);
	RandomStream rng(13, Purpose::Test, 0, 0, 0);
	int red = 0;
	int blue = 0;
	for (int i = 0; i < 20000; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		Evaluate(prepared, p, sample);
		REQUIRE(sample.gas.y == 0.0);
		red += sample.gas.x > 0.0 && sample.gas.z == 0.0 ? 1 : 0;
		blue += sample.gas.z > 0.0 && sample.gas.x == 0.0 ? 1 : 0;
	}
	// Clumps are LN(1.2): below 1 for ~73% of points, above 3 for ~6.5%.
	CHECK(blue > 0);
	CHECK(red > 5 * blue);
}

TEST_CASE("Without glow there is no gas")
{
	GalaxySample sample;
	EvaluateLocal(GasDisk(0.0), {300.0, -200.0, 100.0}, sample);
	CHECK(sample.gas == glm::dvec3(0.0));
}
```

Append to `tests/test_scene.cpp`:

```cpp
TEST_CASE("Gas keys are read and checked")
{
	// Mutation: deleting the arm_population check in ReadDisk reddens the last check.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk.gas]\n";
	Expected<Scene> const scene = ParseScene(head + "glow = 1.5\ncolour = [1.0, 0.2, 0.3]\ncore_colour = [0.2, 1.0, 0.8]\n", "test.toml");
	REQUIRE(scene.has_value());
	Gas const& gas = scene->galaxies[0].disk.gas;
	CHECK(gas.glow == 1.5);
	CHECK(gas.colour == glm::dvec3(1.0, 0.2, 0.3));
	CHECK(gas.coreColour == glm::dvec3(0.2, 1.0, 0.8));
	CHECK(ParseScene(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n", "test.toml")->galaxies[0].disk.gas.glow == 0.0);
	CHECK(ErrorOf(head + "glow = -1.0\n").find("disk.gas.glow") != std::string::npos);
	CHECK(ErrorOf(head + "glow = 1001.0\n").find("must be in") != std::string::npos);
	CHECK(ErrorOf(head + "colour = [-0.1, 1.0, 1.0]\n").find("disk.gas.colour") != std::string::npos);
	CHECK(ErrorOf(head + "core_colour = [0.0, 0.0, 0.0]\n").find("positive luminance") != std::string::npos);
	CHECK(ErrorOf(head + "colour = \"red\"\n").find("three numbers") != std::string::npos);
	CHECK(ErrorOf(head + "glw = 1.0\n").find("gas.glw") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.thick_disk.gas]\nglow = 1.0\n").find("unknown key") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\ndensity = 0.1\ncutoff = 5000\npopulation = \"old-disk\"\n[galaxy.disk.gas]\nglow = 1.0\n")
		.find("arm_population") != std::string::npos);
}
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `nice -n 19 ionice -c 3 cmake --build Build -j 4`
Expected: compile errors. `Disk` has no `gas` member and `GalaxySample` has no `gas`.

- [ ] **Step 3: Implement the data**

In `src/galaxy.h`, before `struct Disk`:

```cpp
// Gas lit by a disk's young stars. Its light is `glow` times theirs, averaged over space: it gathers where
// they cluster and onto the dust's filaments, and the dust's bubble clears it too.
struct Gas
{
	double     glow = 0.0;
	glm::dvec3 colour{1.0, 0.18, 0.32};      // hydrogen red; PrepareGalaxy normalises it to luminance 1
	glm::dvec3 coreColour{0.25, 1.0, 0.85};  // the densest clumps; normalised likewise
	double     textureMean = 1.0;            // derived by PrepareGalaxy: ⟨G⟩, divided out like Dust::filamentMean
	double     youngLight  = 0.0;            // derived by PrepareScene: the arm population's mean L_V per star
};
```

At the end of `struct Disk`, after `associationMean`:

```cpp
	Gas    gas;                         // thin disks only; needs armPopulation
```

In `struct GalaxySample`, after `extinction`:

```cpp
	glm::dvec3                  gas{0.0};           // glowing gas: radiance per parsec, before transmittance
```

- [ ] **Step 4: Implement the model**

In `src/galaxy.cpp`, add `#include "population.h"`.

Next to the other salts and constants:

```cpp
// The gas texture's noise. Fixed rather than scene keys: the owner-approved look was tuned on exactly these.
constexpr uint64_t kGasClumpSalt = 0x6A5ull;
constexpr uint64_t kGasWispSalt = 0x6A6ull;
constexpr double kGasClumpWavelength = 40.0;   // parsecs: a bright nebula's core
constexpr double kGasWispWavelength = 25.0;
constexpr int kGasOctaves = 5;
```

Add these helpers to the anonymous namespace, after `AssociationMean`. They lift the dust's wavelength and bubble out of `EvaluateLocal` so the gas reads the same fields:

```cpp
double DustWavelength(Galaxy const& galaxy)
{
	return 0.5 * (galaxy.disk.arms.count > 0 ? galaxy.disk.arms.width : galaxy.dust.scaleHeight);
}

// Mean one over space: the dust's filaments with the extinction's normalisation.
double DustFilaments(Galaxy const& galaxy, glm::dvec3 p)
{
	return DustField(galaxy.dust, p / DustWavelength(galaxy), galaxy.seed) / galaxy.dust.filamentMean;
}

// Empty within 70% of the bubble's semi-axes, ramping smoothly to one at its edge; one without a bubble.
double BubbleFactor(Dust const& dust, glm::dvec3 p)
{
	if (dust.bubbleRadius <= 0.0)
	{
		return 1.0;
	}
	glm::dvec3 const d = p - dust.bubbleCenter;
	double const q = std::sqrt((d.x * d.x + d.y * d.y) / (dust.bubbleRadius * dust.bubbleRadius) + d.z * d.z / (dust.bubbleHeight * dust.bubbleHeight));
	double const t = std::clamp((q - 0.7) / 0.3, 0.0, 1.0);
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

void NormaliseColours(Gas& gas)
{
	gas.colour /= Luminance(gas.colour);
	gas.coreColour /= Luminance(gas.coreColour);
}
```

In `PrepareGalaxy`, after the `associationMean` lines:

```cpp
	prepared.galaxy.disk.gas.textureMean = GasTextureMean(prepared.galaxy, prepared.galaxy.disk);
	NormaliseColours(prepared.galaxy.disk.gas);
```

In `EvaluateLocal`, reset the gas next to `out.extinction = 0.0;`:

```cpp
	out.gas = glm::dvec3(0.0);
```

In `addDisk`, replace the final `push(midplane * young * associations * Sech2(p.z / (2.0 * youngHeight)), disk.armPopulation);` with:

```cpp
		double const youngDensity = midplane * young * associations * Sech2(p.z / (2.0 * youngHeight));
		push(youngDensity, disk.armPopulation);
		if (disk.gas.glow > 0.0 && youngDensity > 0.0)
		{
			double const clumps = GasClumps(galaxy, p);
			double const texture = GasTexture(galaxy, p, clumps) / disk.gas.textureMean;
			// Dense clumps turn from hydrogen red to the core colour; both have luminance 1, so the mix keeps the light.
			double const t = std::clamp((clumps - 1.0) / 2.0, 0.0, 1.0);
			glm::dvec3 const colour = glm::mix(disk.gas.colour, disk.gas.coreColour, t * t * (3.0 - 2.0 * t));
			double const light = disk.gas.glow * disk.gas.youngLight * youngDensity * texture * BubbleFactor(galaxy.dust, p);
			out.gas += (light * 100.0 * kFluxAt10pc) * colour;
		}
```

In the dust block of `EvaluateLocal`, replace the inline wavelength, filaments and bubble code with the helpers. This preserves behaviour; the existing dust and bubble tests guard it:

```cpp
		double const filaments = DustFilaments(galaxy, p);
		out.extinction = dust.extinction * std::exp(-radius / dust.scaleLength) * Sech2(p.z / (2.0 * dust.scaleHeight)) * laneModulation * filaments
			* BubbleFactor(dust, p);
```

Delete the now-unused `wavelength` local and the old `if (dust.bubbleRadius > 0.0) { ... }` block, but keep both comments explaining why. The log-normal comment goes above `DustFilaments`. The "Without the cavity…" comment goes above `BubbleFactor`.

- [ ] **Step 5: Implement the scene keys**

In `src/scene.cpp`, after `ReadArms`:

```cpp
Expected<void> ReadColour(TableReader& reader, std::string_view key, glm::dvec3& colour)
{
	STARSCAPE_TRY(reader.Vector(key, colour));
	// Normalisation divides by the luminance, so a black or negative colour must stop here.
	if (colour.x < 0.0 || colour.y < 0.0 || colour.z < 0.0 || !(Luminance(colour) > 0.0))
	{
		return reader.Fail(key, "needs non-negative components and a positive luminance");
	}
	return {};
}

Expected<void> ReadGas(TableReader& parent, Gas& gas)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("gas"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "gas", *table));
	STARSCAPE_TRY(reader.Number("glow", gas.glow, 0.0, 1000.0));
	STARSCAPE_TRY(ReadColour(reader, "colour", gas.colour));
	STARSCAPE_TRY(ReadColour(reader, "core_colour", gas.coreColour));
	return reader.Finish();
}
```

In `ReadDisk`, inside `if (thin)` after `STARSCAPE_TRY(ReadArms(reader, disk.arms));` and the existing `arm_young_fraction` check:

```cpp
			STARSCAPE_TRY(ReadGas(reader, disk.gas));
			if (disk.enabled && disk.gas.glow > 0.0 && disk.armPopulation < 0)
			{
				return reader.Fail("gas", "glow needs arm_population to be set: young stars light the gas");
			}
```

- [ ] **Step 6: Run the tests**

Run: `nice -n 19 ionice -c 3 cmake --build Build -j 4 && ctest --test-dir Build --output-on-failure`
Expected: every test passes. If the calibration test misses by more than 3%, do not widen the tolerance. Report the measured ratio and ⟨G⟩.

- [ ] **Step 7: Verify the mutations**

For each mutation: apply it, rebuild, run the named test and see it red, then revert and `touch`:
- In `EvaluateLocal`, drop `/ disk.gas.textureMean`. Run "Gas glows at glow times the young stars' light, on average".
- Drop `* BubbleFactor(galaxy.dust, p)` from the gas `light`. Run "Gas shares the dust's bubble".
- Swap `disk.gas.colour` and `disk.gas.coreColour` in the `glm::mix`. Run "Gas colours are normalised and the densest clumps take the core colour".
- Delete the `glow needs arm_population` check. Run "Gas keys are read and checked".

- [ ] **Step 8: Commit**

```bash
git add src/galaxy.h src/galaxy.cpp src/scene.cpp tests/test_galaxy.cpp tests/test_scene.cpp
git commit -m "Gas: young-star glow on the dust's filaments, outside its bubble

mutation verified: no textureMean division; no bubble on gas; mix endpoints swapped; arm_population check removed"
```

---

### Task 3: Gas in the march

**Files:**
- Modify: `src/march.cpp`
- Test: `tests/test_march.cpp`

**Interfaces:**
- Consumes: `GalaxySample::gas`, `Gas::youngLight` and `Tabulate(population, saturation)` from Tasks 1 and 2.
- Produces: `PrepareScene` sets `disk.gas.youngLight`. `Step` carries `glm::dvec3 gas`. `MarchRay` adds `weight · ds · T ⊙ gas` to `radiance`, and nothing to `budget` or `resolvedRadiance`. `HistogramRay` and `SampleStar` ignore it.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_march.cpp`:

```cpp
TEST_CASE("Gas glow is dimmed by the dust in front of it")
{
	// Mutation: adding the gas without the step's transmittance reddens this test.
	Scene scene = BaseScene();
	int const population = AddSingleBinPopulation(scene, 14);
	Galaxy source;   // dust disabled: the gas's filament factor must be 1
	source.position = {1.0e6, 0.0, 0.0};
	source.disk = {.enabled = true, .density = 0.05, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 20000.0, .arms = {},
		.population = population, .armPopulation = population, .youngFraction = 0.5};
	scene.galaxies.push_back(source);
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir(1.0, 0.0, 0.0);
	// The gas's own radiance: the same ray with and without glow, so the starlight cancels exactly.
	auto gasRadiance = [&](Scene const& base)
	{
		RayResult dark;
		MarchRay(PrepareScene(base), params, dir, HUGE_VAL, 1.0, 1.0, dark);
		Scene lit = base;
		lit.galaxies[0].disk.gas.glow = 1.0;
		RayResult glowing;
		MarchRay(PrepareScene(lit), params, dir, HUGE_VAL, 1.0, 1.0, glowing);
		return glowing.radiance - dark.radiance;
	};
	glm::dvec3 const clear = gasRadiance(scene);
	REQUIRE(Luminance(clear) > 0.0);
	double const alpha = 0.002;
	double const radius = 1000.0;
	scene.galaxies.push_back(UniformDust(alpha, radius));
	glm::dvec3 const dimmed = gasRadiance(scene);
	double const tau = alpha * (radius - scene.observer.sMin) / 1.0857;
	CHECK(dimmed.x / clear.x == doctest::Approx(std::exp(-0.75 * tau)).epsilon(0.005).scale(0.0));
	CHECK(dimmed.y / clear.y == doctest::Approx(std::exp(-1.00 * tau)).epsilon(0.005).scale(0.0));
	CHECK(dimmed.z / clear.z == doctest::Approx(std::exp(-1.32 * tau)).epsilon(0.005).scale(0.0));
}

TEST_CASE("Gas adds diffuse light and leaves the stars alone")
{
	// Mutation: adding the gas to resolvedRadiance instead of radiance reddens this test.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	PreparedScene const plain = PrepareScene(scene);
	scene.galaxies[0].disk.gas.glow = 1.5;
	PreparedScene const glowing = PrepareScene(scene);
	CHECK(glowing.galaxies[0].galaxy.disk.gas.youngLight == doctest::Approx(Luminance(glowing.populations[0].lightBelow[kBins])));
	MarchParams const params = ParamsFor(Quality::Quick);
	glm::dvec3 const dir = glm::normalize(glm::dvec3(-1.0, 0.1, 0.02));
	RayResult a;
	RayResult b;
	MarchRay(plain, params, dir, 1e-4, 1.0, 1.0, a);
	MarchRay(glowing, params, dir, 1e-4, 1.0, 1.0, b);
	CHECK(Luminance(b.radiance) > 1.01 * Luminance(a.radiance));
	CHECK(b.budget == a.budget);
	CHECK(b.resolvedRadiance == a.resolvedRadiance);
	FluxHistogram ha;
	FluxHistogram hb;
	HistogramRay(plain, params, dir, 1.0, ha);
	HistogramRay(glowing, params, dir, 1.0, hb);
	CHECK(ha.count == hb.count);
}
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `nice -n 19 ionice -c 3 cmake --build Build -j 4 && Build/starscape-tests --test-case='Gas*'`
Expected: both cases fail. `youngLight` is 0, so there's no gas, `clear` is zero and the `youngLight` check is red.

- [ ] **Step 3: Implement**

In `src/march.cpp`, add a member to `struct Step` after `parts`:

```cpp
	glm::dvec3                    gas;   // glowing gas: radiance per parsec at the midpoint, before transmittance
```

In `Walk`, declare `glm::dvec3 gas(0.0);` next to `double extinction = 0.0;`. Add `gas += sample.gas;` after `extinction += sample.extinction;`. Then pass it last in the `Step{...}`:

```cpp
			visit(Step{mid, ds, s, s + ds, TransmittanceForDepth(depth + 0.5 * extinction * ds), parts, gas});
```

In `MarchRay`, make this the first statement of the `Walk` lambda:

```cpp
		// Gas is diffuse light only: it has no stars to resolve, so it never touches the budget or the threshold.
		out.radiance += (weight * step.ds) * (step.transmittance * step.gas);
```

Replace `PrepareScene`'s two loops with populations first, then galaxies:

```cpp
	for (Population const& population : scene.populations)
	{
		prepared.populations.push_back(Tabulate(population, scene.style.starSaturation));
	}
	for (Galaxy const& galaxy : scene.galaxies)
	{
		PreparedGalaxy ready = PrepareGalaxy(galaxy);
		Disk& disk = ready.galaxy.disk;
		if (disk.armPopulation >= 0)
		{
			// Gas glows at a multiple of its young stars' light, and only their table knows how bright they are.
			disk.gas.youngLight = Luminance(prepared.populations[size_t(disk.armPopulation)].lightBelow[kBins]);
		}
		prepared.galaxies.push_back(std::move(ready));
	}
```

- [ ] **Step 4: Run the tests**

Run: `nice -n 19 ionice -c 3 cmake --build Build -j 4 && ctest --test-dir Build --output-on-failure`
Expected: every test passes.

- [ ] **Step 5: Verify the mutations**

Apply each, see it red, then revert and `touch`:
- In `MarchRay`, drop `step.transmittance *` from the gas line. Run "Gas glow is dimmed by the dust in front of it".
- Add the gas line to `out.resolvedRadiance` instead of `out.radiance`. Run "Gas adds diffuse light and leaves the stars alone".

- [ ] **Step 6: Commit**

```bash
git add src/march.cpp tests/test_march.cpp
git commit -m "March: gas glow joins the diffuse light behind the dust

mutation verified: gas without transmittance; gas added as resolved light"
```

---

### Task 4: Example scenes, docs and the owner's visual check

**Files:**
- Modify: `scenes/milky-way-from-earth.toml`, `CMakeLists.txt`, `docs/spec.md`, `README.md`
- Create: `scenes/stylized-milky-way.toml`
- Delete: `docs/specs/2026-10-13-gas-glow-and-style.md` (folded into `docs/spec.md`; this follows the repo's documentation lifecycle)

**Interfaces:**
- Consumes: the scene keys from Tasks 1 and 2.

- [ ] **Step 1: Put the Earth scene in the galactic frame**

In `scenes/milky-way-from-earth.toml`, replace the first two comment lines and the observer's axes:

```toml
# A Milky-Way-like galaxy seen from a Sun-like position: 8.2 kpc out, 20 pc above the plane.
# The sky frame is galactic: the band runs along the equator and the galactic centre is at
# longitude 0, the middle of the preview. A consumer orients the sky by rotating the skybox.
```

```toml
[observer]
position      = [8200.0, 0.0, 20.0]
spin_axis     = [0.0, 0.0, 1.0]
zero_meridian = [-1.0, 0.0, 0.0]
```

- [ ] **Step 2: Create the stylized scene**

`scenes/stylized-milky-way.toml`:

```toml
# The Earth scene pushed toward a fantastical look that keeps the galaxy's structure: darker, more
# filamentary dust, bolder arms, young-star gas glowing red with teal cores, and more vivid star colour.
# Target: docs/images/2026-10-13-stylized-target.png. Galactic frame, as in milky-way-from-earth.toml.
[output]
face_size    = 256
target_stars = 4000
quality      = "final"
preview      = { width = 8192, march = false }

[observer]
position      = [8200.0, 0.0, 20.0]
spin_axis     = [0.0, 0.0, 1.0]
zero_meridian = [-1.0, 0.0, 0.0]

[style]
star_saturation = 1.5

[[galaxy]]
preset = "milky-way"
seed   = 1

  [galaxy.disk.arms]
  contrast = 2.25

  # 1.5 reproduces the target: its prototype ran at 2.5 before the texture's mean (0.6045 here) was divided out.
  [galaxy.disk.gas]
  glow = 1.5

  [galaxy.dust]
  bubble   = { center = [8200.0, 0.0, 20.0], radius = 300.0, height = 800.0 }
  contrast = 6.0
  filament = 3.5

[[galaxy]]
preset   = "small-spiral"
position = [780000.0, 120000.0, -40000.0]
rotation = [77.0, 0.0, 35.0]
seed     = 2
```

- [ ] **Step 3: Cook the stylized scene in the CLI tests**

In `CMakeLists.txt`, after the `cli-draft` line:

```cmake
starscape_cli_test(cli-draft-stylized 0 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/scenes/stylized-milky-way.toml -o ${CLI_OUT}/stylized.kre_starscape --draft --preview ${CLI_OUT}/stylized.hdr)
```

Run: `nice -n 19 ionice -c 3 cmake --build Build -j 4 && ctest --test-dir Build --output-on-failure`
Expected: every test passes, including `cli-draft` and `cli-draft-stylized`.

- [ ] **Step 4: Render for the owner (visual acceptance, not a gate)**

```bash
mkdir -p Build-prof/accept
nice -n 19 Build/starscape cook scenes/stylized-milky-way.toml -o Build-prof/accept/stylized.kre_starscape --preview Build-prof/accept/stylized.hdr --quick --width 2048
```

Tone-map it to PNG at a fixed exposure: white at its 99.9th percentile of luminance, black at the 5th, gamma 2.2. Stack it above `docs/images/2026-10-13-stylized-target.png`, scaled to the same width, into `Build-prof/accept/compare.png`. The target used the same mapping, so a matching look means matching numbers.

Then print:
`UNVERIFIED-VISUAL: Build-prof/accept/compare.png — stylized scene vs owner-approved target; owner to judge`

A large mismatch in overall brightness or gas strength is a finding for the controller. Do not retune keys on your own.

- [ ] **Step 5: Fold the spec into `docs/spec.md` and delete it**

Edit `docs/spec.md`:
- **Non-goals:** replace "Emission nebulae (HII regions)" with "Placed or near nebulae: a 40° nebula looked right in a still frame but dominated the sky in a 3D viewer".
- **§1:** after the **Dust** subsection, add a **Gas** subsection. Copy the new spec's "Why" bullets on structure and nebulae (two sentences), and its "Model" and "In the march" text and formula verbatim.
- **§4 Scene (TOML):**
  - add `[style]` with `star_saturation` to the example, after `[observer]`;
  - add the `[galaxy.disk.gas]` block under the example's `[galaxy.disk]`;
  - add a bullet each for `[style]` (with the saturation formula and the luminance note) and for the gas keys' ranges and errors.
- **§4 example and observer prose:** say the example scenes use the galactic frame and why: consumers rotate the skybox, and the preview stays legible.
- **§5:** add the new spec's test rows to the table.
- **§6:** list `scenes/stylized-milky-way.toml` and `docs/images/`.

Then:

```bash
git rm docs/specs/2026-10-13-gas-glow-and-style.md
```

In `README.md`'s scene-format paragraph, add `[style]` (scene-wide stylisation such as `star_saturation`) to the list of tables. Name `scenes/stylized-milky-way.toml` next to the Earth example.

Check that nothing still links to the deleted file:

```bash
grep -rn "2026-10-13-gas-glow-and-style" --include=*.md . | grep -v "^./docs/plans/"
```
Expected: no output. The plan itself may keep its Spec line.

- [ ] **Step 6: Commit**

```bash
git add scenes/ CMakeLists.txt docs/spec.md README.md
git commit -m "Scenes: galactic frame, a stylized Milky Way; spec folds in gas and style"
```

- [ ] **Step 7: Drop the prototype stash only after the owner accepts the look**

Leave `stash@{0}` ("nebula/gas look-dev prototype") in place. The controller asks the owner after showing `compare.png`. Dropping it is the owner's call.
