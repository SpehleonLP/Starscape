# Nebula Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Four new emission terms (a gas palette, large-scale gas, a bubble-rim
glow and dust reflection) as scene keys that default to off, with the stylized
scene reproducing `docs/images/2026-10-15-nebula-target.png`.

**Architecture:**
- **Where each term lives.** Every term is computed in `EvaluateLocal`
  (`src/galaxy.cpp`) from fields the galaxy already has, and added to
  `GalaxySample::gas`. The march already adds that behind each step's
  transmittance, so `src/march.cpp` does not change.
- **Prepared values.** Colours are normalised, and the large-scale field's mean
  measured, in `PrepareGalaxy`.
- **Scene keys** are read in `src/scene.cpp` next to the existing gas, dust and
  bubble keys.

**Tech Stack:** C++23 (GCC 13), CMake ≥ 3.25 + Ninja, vendored glm, toml++ and doctest.

**Spec:** `docs/specs/2026-10-15-nebula-layer.md`. Read it first; it is the
authority. The base model is `docs/spec.md`, and the spike that the target came
from is branch `spike/nebula-layer` at 5d09f61 (reference only; never merge it).

## Global Constraints

- **Repository:**
  - Starscape lives at `/mnt/Passport/Kreatures-Engine/Tooling/Starscape`. It is
    a submodule of Kreatures-Engine, and its branch is `main`.
  - Work in a new worktree at `/mnt/Passport/Worktrees/starscape-nebula`, on
    branch `nebula-layer` from `main` (Task 1). Paths below are relative to it.
  - Never touch the Kreatures-Engine superproject, and never push.
- **Build and test:**
  - `cmake -S . -B Build -G Ninja && nice -n 19 cmake --build Build`.
  - Run everything with `ctest --test-dir Build --output-on-failure`.
  - Run one case with `Build/starscape-tests -tc="<name>"`.
- **Style:**
  - C++23, `-Wall -Wextra -Wpedantic` clean;
  - tabs and Allman braces, with **every** control-flow body braced;
  - comments say *why*, not *what*.
- **Parameters:** containers and structs are passed `const&`; glm vectors go by
  value.
- **Bounded collections** are `std::array` plus a count, never `std::vector`.
  The GPU port (the next spec) mirrors these structs.
- **Errors** are returned as `Expected<T>`. No exceptions, no asserts on input.
- **Units:**
  - parsecs, L_V☉, magnitude-0 flux;
  - an emissivity is `light · 100 · kFluxAt10pc` times a colour of luminance 1,
    exactly as the gas term does it.
- **Exact constants from the spike (needed to reproduce the target):**

  | Term | Salt | Field |
  |---|---|---|
  | Palette | `0x9A1E77Eull` | 3 octaves |
  | Rim mask | `0x81B1ull` | 150 pc, 3 octaves |
  | Rim wisps | `0x81B0ull` | 30 pc, 5 octaves |
  | Large-scale gas | `0xB16B16ull` | 4 octaves |

  Every salt is XORed with `galaxy.seed`.
- **Mutation rule:** every new test's comment names its mutation. Apply the
  mutation by hand, watch the test go red, and revert, *before* committing the
  task. Record "mutation checked" in the commit message.

## Review Focus

1. **A palette entry that is black or negative.** Normalisation divides by the
   luminance, so such an entry must be a scene error naming
   `palette[i]`, never a NaN in the sky. Owned by Task 2.
2. **A scene that sets `rim` but disables the dust** (`[galaxy.dust] enabled = false`
   with a bubble). Disabled dust is gone, bubble included, so the rim must be
   exactly zero. Owned by Task 4.
3. **A `palette_wavelength` or `large_wavelength` far larger than the galaxy**
   (up to the 100 000 pc limit). The terms must stay finite, and the large-scale
   mean must not depend on the wavelength (it is measured in noise space).
   Owned by Task 3.
4. **The thick disk.** It has no gas table, so the new gas keys under
   `[galaxy.thick_disk.gas]` must stay unknown keys, and the thick disk must
   never emit the new terms. Owned by Task 2.
5. **A palette of other than four entries.** The ramp is cyclic, so with three
   entries the last blends into the first. It must never read past
   `paletteCount` into the zeroed rest of the array, which would darken the gas.
   Owned by Task 2.

---

### Task 1: Worktree and the pre-change reference

**Files:** none in the repository; one scratch file holding two hashes.

The spec requires every existing scene to cook exactly as before. This task
records today's bytes, so Task 6 can prove it.

- [ ] **Step 1: Create the worktree**

```bash
cd /mnt/Passport/Kreatures-Engine/Tooling/Starscape
git worktree add /mnt/Passport/Worktrees/starscape-nebula -b nebula-layer main
cd /mnt/Passport/Worktrees/starscape-nebula
cmake -S . -B Build -G Ninja && nice -n 19 cmake --build Build
ctest --test-dir Build --output-on-failure
```

Expected: every test passes. Record the count; that is the baseline.

- [ ] **Step 2: Record draft-cook hashes of both example scenes**

```bash
mkdir -p Build/reference
for s in milky-way-from-earth stylized-milky-way; do
  Build/starscape cook scenes/$s.toml -o Build/reference/$s.kre_starscape --draft
done
sha256sum Build/reference/*.kre_starscape | tee Build/reference/pre-nebula.sha256
```

Expected: two hashes. Do not commit anything; `Build/` is ignored.

---

### Task 2: The palette

**Files:**
- Modify: `src/galaxy.h` (struct `Gas`)
- Modify: `src/galaxy.cpp` (new `PaletteColour`, `NormaliseColours`, the gas colour in `EvaluateLocal`)
- Modify: `src/scene.cpp` (`ReadGas`, new `ReadPalette`)
- Test: `tests/test_galaxy.cpp`, `tests/test_scene.cpp`

**Interfaces:**
- Produces:
  - `inline constexpr int kMaxPalette = 8;` in `galaxy.h`;
  - in `Gas`: `std::array<glm::dvec3, kMaxPalette> palette{}`, `int paletteCount = 0`,
    `double paletteWavelength = 1500.0`;
  - in `galaxy.cpp`, internal: `glm::dvec3 PaletteColour(Galaxy const& galaxy, Gas const& gas, glm::dvec3 p)`.
    It returns `gas.colour` when `paletteCount == 0`. Tasks 3 and 4 call it.

- [ ] **Step 1: Write the failing galaxy test**

Append to `tests/test_galaxy.cpp`, after `"Disabled dust leaves the gas texture alone"`:

```cpp
TEST_CASE("A palette colours the gas by a large-scale field without changing its light")
{
	// Mutation: dropping the per-entry normalisation in NormaliseColours leaves entries off Y = 1.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.disk.gas.palette = {glm::dvec3(2.0, 0.0, 0.0), glm::dvec3(0.0, 3.0, 0.0), glm::dvec3(0.0, 0.0, 4.0), glm::dvec3(1.0, 1.0, 0.0)};
	galaxy.disk.gas.paletteCount = 4;
	galaxy.disk.gas.paletteWavelength = 1500.0;
	// A core colour equal to nothing in the palette would hide which entry a point took; make it one.
	galaxy.disk.gas.coreColour = galaxy.disk.gas.palette[0];
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	for (int i = 0; i < 4; ++i)
	{
		CHECK(Luminance(prepared.galaxy.disk.gas.palette[size_t(i)]) == doctest::Approx(1.0).epsilon(1e-12));
	}
	Galaxy plain = galaxy;
	plain.disk.gas.paletteCount = 0;
	PreparedGalaxy const reference = PrepareGalaxy(plain);
	RandomStream rng(17, Purpose::Test, 0, 0, 0);
	std::array<int, 4> nearest{};
	int const count = 100000;
	for (int i = 0; i < count; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		GalaxySample without;
		Evaluate(prepared, p, sample);
		Evaluate(reference, p, without);
		REQUIRE(Luminance(without.gas) > 0.0);
		// Colour never changes the amount of light: the palette's gas is exactly as bright as the plain gas.
		REQUIRE(Luminance(sample.gas) == doctest::Approx(Luminance(without.gas)).epsilon(1e-12));
		glm::dvec3 const colour = sample.gas / Luminance(sample.gas);
		int best = 0;
		for (int e = 1; e < 4; ++e)
		{
			if (glm::distance(colour, prepared.galaxy.disk.gas.palette[size_t(e)]) < glm::distance(colour, prepared.galaxy.disk.gas.palette[size_t(best)]))
			{
				best = e;
			}
		}
		++nearest[size_t(best)];
	}
	// Every entry appears. The spec asked for 10% each, but with 0.9 entries per σ the entry two steps
	// from the field's mean is reached only in the tails (|z| > 1.1), so the floor is 1%. Task 7 corrects
	// the spec's row when it folds in.
	for (int e = 0; e < 4; ++e)
	{
		CHECK(nearest[size_t(e)] > count / 100);
	}
}

TEST_CASE("A three-colour palette wraps from its last entry to its first")
{
	// Mutation: wrapping with `% kMaxPalette` instead of `% paletteCount` blends into a zeroed entry, and
	// the gas's luminance drops below one.
	Galaxy galaxy = GasDisk(1.0);
	galaxy.disk.gas.palette = {glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(0.0, 1.0, 0.0), glm::dvec3(0.0, 0.0, 1.0)};
	galaxy.disk.gas.paletteCount = 3;
	// A green core: any core mixing then shows as green, so red with blue and no green can only be the wrap.
	galaxy.disk.gas.coreColour = {0.0, 1.0, 0.0};
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	Galaxy plain = galaxy;
	plain.disk.gas.paletteCount = 0;
	PreparedGalaxy const reference = PrepareGalaxy(plain);
	RandomStream rng(43, Purpose::Test, 0, 0, 0);
	bool blendedLastIntoFirst = false;
	for (int i = 0; i < 100000; ++i)
	{
		glm::dvec3 const p(rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5, rng.NextDouble() * 1e6 - 5e5);
		GalaxySample sample;
		GalaxySample without;
		Evaluate(prepared, p, sample);
		Evaluate(reference, p, without);
		REQUIRE(Luminance(without.gas) > 0.0);
		// Colour never changes the amount of light: the palette's gas is exactly as bright as the plain gas.
		REQUIRE(Luminance(sample.gas) == doctest::Approx(Luminance(without.gas)).epsilon(1e-12));
		blendedLastIntoFirst = blendedLastIntoFirst || (sample.gas.x > 0.0 && sample.gas.z > 0.0 && sample.gas.y == 0.0);
	}
	CHECK(blendedLastIntoFirst);
}

TEST_CASE("The thick disk never takes the thin disk's palette")
{
	// Mutation: passing galaxy.disk.gas instead of the thick disk's own gas into its gas block lights it;
	// a thick disk alone must stay dark.
	Galaxy galaxy;
	galaxy.thickDisk = {.enabled = true, .density = 1.0, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .arms = {}, .population = 1};
	galaxy.disk.gas.glow = 1.0;
	galaxy.disk.gas.youngLight = 1.0;
	galaxy.disk.gas.paletteCount = 2;
	galaxy.disk.gas.palette = {glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0)};
	GalaxySample sample;
	EvaluateLocal(PrepareGalaxy(galaxy).galaxy, {300.0, -200.0, 100.0}, sample);
	CHECK(sample.gas == glm::dvec3(0.0));
}
```

Add `#include <glm/geometric.hpp>` and `#include <array>` to the test file's
includes.

- [ ] **Step 2: Run it to verify it fails to compile**

Run: `nice -n 19 cmake --build Build`
Expected: errors that `Gas` has no member `palette`, `paletteCount` or `paletteWavelength`.

- [ ] **Step 3: Add the fields**

In `src/galaxy.h`, above `struct Gas`:

```cpp
// The most colours a gas palette takes: enough for any ramp, small enough to live in a fixed array.
inline constexpr int kMaxPalette = 8;
```

and in `struct Gas`, after `coreColour`:

```cpp
	// A cyclic ramp replacing `colour`, chosen by a field of `paletteWavelength`, so that neighbouring
	// clouds differ in hue. Empty (count 0): `colour` alone. PrepareGalaxy normalises each to luminance 1.
	std::array<glm::dvec3, kMaxPalette> palette{};
	int        paletteCount = 0;
	double     paletteWavelength = 1500.0;   // parsecs
```

- [ ] **Step 4: Normalise and evaluate the palette**

In `src/galaxy.cpp`, add the salt beside the gas constants:

```cpp
// The palette's hue field. Fixed rather than a scene key: the owner-approved nebula look was tuned on it.
constexpr uint64_t kPaletteSalt = 0x9A1E77Eull;
```

Replace `NormaliseColours` with:

```cpp
void NormaliseColours(Gas& gas)
{
	gas.colour /= Luminance(gas.colour);
	gas.coreColour /= Luminance(gas.coreColour);
	for (int i = 0; i < gas.paletteCount; ++i)
	{
		gas.palette[size_t(i)] /= Luminance(gas.palette[size_t(i)]);
	}
}
```

Add after `GasTextureMean`:

```cpp
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
```

In `EvaluateLocal`'s gas block, replace the colour line:

```cpp
			glm::dvec3 const colour = glm::mix(disk.gas.colour, disk.gas.coreColour, t * t * (3.0 - 2.0 * t));
```

with:

```cpp
			glm::dvec3 const colour = glm::mix(PaletteColour(galaxy, disk.gas, p), disk.gas.coreColour, t * t * (3.0 - 2.0 * t));
```

- [ ] **Step 5: Run the galaxy tests**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="*palette*,*Gas*,*gas*"`
Expected: PASS, including the existing gas tests. With no palette, `PaletteColour`
returns `colour` itself, so they are unchanged.

- [ ] **Step 6: Write the failing scene test**

Append to `tests/test_scene.cpp`:

```cpp
TEST_CASE("A gas palette is read and checked")
{
	// Mutation: deleting the luminance check in ReadPalette lets a black entry through to a division by zero.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk.gas]\nglow = 1.0\n";
	Expected<Scene> const scene = ParseScene(head + "palette = [[1.0, 0.2, 0.3], [0.2, 1.0, 0.8], [1.0, 0.7, 0.2]]\npalette_wavelength = 900.0\n", "test.toml");
	REQUIRE(scene.has_value());
	Gas const& gas = scene->galaxies[0].disk.gas;
	CHECK(gas.paletteCount == 3);
	CHECK(gas.palette[2] == glm::dvec3(1.0, 0.7, 0.2));
	CHECK(gas.paletteWavelength == 900.0);
	CHECK(ParseScene(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n", "test.toml")->galaxies[0].disk.gas.paletteCount == 0);
	CHECK(ParseScene(head + "palette = []\n", "test.toml")->galaxies[0].disk.gas.paletteCount == 0);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3]]\n").find("0 or 2 to 8") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0]]\n").find("0 or 2 to 8") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3], [0.0, 0.0, 0.0]]\n").find("palette[1]") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3], [-1.0, 1.0, 1.0]]\n").find("positive luminance") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3], [1.0, 1.0]]\n").find("three numbers") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [1.0, 0.2, 0.3]\n").find("three numbers") != std::string::npos);
	CHECK(ErrorOf(head + "palette = \"warm\"\n").find("expected an array") != std::string::npos);
	CHECK(ErrorOf(head + "palette_wavelength = 5.0\n").find("disk.gas.palette_wavelength") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.thick_disk.gas]\npalette = []\n").find("unknown key") != std::string::npos);
}
```

- [ ] **Step 7: Run it to verify it fails**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="A gas palette is read and checked"`
Expected: FAIL, `palette` reported as an unknown key.

- [ ] **Step 8: Read the palette**

In `src/scene.cpp`, add before `ReadGas`:

```cpp
// An array of colours, each checked like `colour`: normalisation divides by the luminance.
Expected<void> ReadPalette(TableReader& reader, Gas& gas)
{
	STARSCAPE_TRY_ASSIGN(array, reader.Array("palette"));
	if (array == nullptr)
	{
		return {};
	}
	if (array->size() == 1 || array->size() > size_t(kMaxPalette))
	{
		return reader.Fail(array->source(), "palette", std::format("needs 0 or 2 to {} colours", kMaxPalette));
	}
	for (size_t i = 0; i < array->size(); ++i)
	{
		toml::node const& element = (*array)[i];
		std::string const key = std::format("palette[{}]", i);
		toml::array const* rgb = element.as_array();
		if (rgb == nullptr || rgb->size() != 3)
		{
			return reader.Fail(element.source(), key, "expected an array of three numbers");
		}
		glm::dvec3 colour(0.0);
		for (size_t c = 0; c < 3; ++c)
		{
			std::optional<double> const value = (*rgb)[c].is_boolean() ? std::nullopt : (*rgb)[c].value<double>();
			if (!value || !std::isfinite(*value))
			{
				return reader.Fail((*rgb)[c].source(), key, "expected an array of three numbers");
			}
			colour[int(c)] = *value;
		}
		if (colour.x < 0.0 || colour.y < 0.0 || colour.z < 0.0 || !(Luminance(colour) > 0.0))
		{
			return reader.Fail(element.source(), key, "needs non-negative components and a positive luminance");
		}
		gas.palette[i] = colour;
	}
	gas.paletteCount = int(array->size());
	return {};
}
```

In `ReadGas`, before `return reader.Finish();`:

```cpp
	STARSCAPE_TRY(ReadPalette(reader, gas));
	STARSCAPE_TRY(reader.Number("palette_wavelength", gas.paletteWavelength, 10.0, 100000.0));
```

- [ ] **Step 9: Run all tests**

Run: `nice -n 19 cmake --build Build && ctest --test-dir Build --output-on-failure`
Expected: everything passes.

- [ ] **Step 10: Check both mutations, then commit**

The mutations are named in each new test's comment:
- drop the palette loop in `NormaliseColours`;
- delete the luminance check in `ReadPalette`.

Apply each, see red, and revert. Then:

```bash
git add src/galaxy.h src/galaxy.cpp src/scene.cpp tests/test_galaxy.cpp tests/test_scene.cpp
git commit -m "Gas: a palette chosen by a large-scale hue field (mutation checked)"
```

---

### Task 3: Large-scale gas

**Files:**
- Modify: `src/galaxy.h` (struct `Gas`)
- Modify: `src/galaxy.cpp` (`LargeGasMean`, the term in `EvaluateLocal`, `PrepareGalaxy`)
- Modify: `src/scene.cpp` (`ReadGas`, `ReadDisk`)
- Test: `tests/test_galaxy.cpp`, `tests/test_scene.cpp`

**Interfaces:**
- Consumes: `PaletteColour` (Task 2).
- Produces:
  - in `Gas`: `double largeGlow = 0.0`, `double largeWavelength = 100.0`;
  - `double largeMean = 1.0`, derived by `PrepareGalaxy`.

- [ ] **Step 1: Write the failing calibration test**

Append to `tests/test_galaxy.cpp`:

```cpp
TEST_CASE("Large-scale gas glows at large_glow times the young stars' light, on average")
{
	// Mutation: lighting the term by the disk's whole midplane density instead of its young density puts
	// the ratio at large_glow / young_fraction. The second wavelength spans few noise cells across the box,
	// which only a mean measured in noise space keeps calibrated.
	for (double const wavelength : {100.0, 20000.0})
	{
		Galaxy galaxy = GasDisk(0.0);
		galaxy.disk.gas.largeGlow = 2.0;
		galaxy.disk.gas.largeWavelength = wavelength;
		PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
		RandomStream rng(19, Purpose::Test, 0, 0, 0);
		double gas = 0.0;
		double young = 0.0;
		for (int i = 0; i < 200000; ++i)
		{
			glm::dvec3 const p(rng.NextDouble() * 1e8 - 5e7, rng.NextDouble() * 1e8 - 5e7, rng.NextDouble() * 1e8 - 5e7);
			GalaxySample sample;
			Evaluate(prepared, p, sample);
			REQUIRE(sample.count == 2);
			REQUIRE(std::isfinite(Luminance(sample.gas)));
			gas += Luminance(sample.gas);
			young += sample.parts[1].density * 100.0 * kFluxAt10pc;
		}
		CHECK(gas / young == doctest::Approx(2.0).epsilon(0.03));
	}
}

TEST_CASE("The large-scale field's mean is measured in noise space and divided out")
{
	// Mutations: replacing LargeGasMean's walk with a constant 1 reddens the first three checks; dropping
	// the division by largeMean in EvaluateLocal reddens the last. At σ 1 the mean is 1.007, too close to
	// one for the calibration test's 3% to see, so the division is pinned directly.
	Galaxy galaxy = GasDisk(0.0);
	galaxy.disk.gas.largeGlow = 1.0;
	PreparedGalaxy const a = PrepareGalaxy(galaxy);
	CHECK(a.galaxy.disk.gas.largeMean == doctest::Approx(1.006899).epsilon(1e-5));   // measured at seed 1, 2^18 points
	galaxy.seed = 2;
	CHECK(PrepareGalaxy(galaxy).galaxy.disk.gas.largeMean != a.galaxy.disk.gas.largeMean);
	galaxy.seed = 1;
	galaxy.disk.gas.largeWavelength = 5000.0;
	CHECK(PrepareGalaxy(galaxy).galaxy.disk.gas.largeMean == a.galaxy.disk.gas.largeMean);   // noise space: no wavelength
	PreparedGalaxy doubled = a;
	doubled.galaxy.disk.gas.largeMean *= 2.0;
	glm::dvec3 const p(300.0, -200.0, 100.0);
	GalaxySample once;
	GalaxySample halved;
	Evaluate(a, p, once);
	Evaluate(doubled, p, halved);
	REQUIRE(Luminance(once.gas) > 0.0);
	CHECK(Luminance(halved.gas) == doctest::Approx(0.5 * Luminance(once.gas)).epsilon(1e-12));
}
```

- [ ] **Step 2: Run it to verify it fails to compile**

Run: `nice -n 19 cmake --build Build`
Expected: errors that `Gas` has no `largeGlow`, `largeWavelength` or `largeMean`.

- [ ] **Step 3: Add the fields**

In `struct Gas`, after `paletteWavelength`:

```cpp
	// A second gas band at the size of giant cloud complexes. `largeGlow` is its light as a multiple of the
	// young stars' own, like `glow`.
	double     largeGlow = 0.0;
	double     largeWavelength = 100.0;  // parsecs
	double     largeMean = 1.0;          // derived by PrepareGalaxy: the field's mean in noise space, divided out
```

- [ ] **Step 4: Measure the mean and add the term**

In `src/galaxy.cpp`, beside the other salts:

```cpp
constexpr uint64_t kLargeGasSalt = 0xB16B16ull;
constexpr int kLargeGasOctaves = 4;
```

After `PaletteColour`:

```cpp
double LargeGasField(glm::dvec3 noiseSpace, uint64_t seed)
{
	return LogNormal(1.0, Fbm(noiseSpace, kLargeGasOctaves, seed ^ kLargeGasSalt) / kFbmSigma);
}

// Measured in noise space, so the mean holds at any wavelength: a box in parsecs would span a single
// noise cell at the largest wavelengths and measure nothing.
double LargeGasMean(Galaxy const& galaxy, Gas const& gas)
{
	if (gas.largeGlow <= 0.0)
	{
		return 1.0;
	}
	return SpatialMean([&](glm::dvec3 q) { return LargeGasField(q, galaxy.seed); });
}
```

In `PrepareGalaxy`, after the `textureMean` line:

```cpp
	prepared.galaxy.disk.gas.largeMean = LargeGasMean(prepared.galaxy, prepared.galaxy.disk.gas);
```

In `EvaluateLocal`'s `addDisk`, replace the whole gas block

```cpp
		if (disk.gas.glow > 0.0 && youngDensity > 0.0)
		{
			...
		}
```

with:

```cpp
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
```

The thick disk has no gas keys (Task 2's test pins it), so both strengths are
zero there and the block is skipped.

- [ ] **Step 5: Run the galaxy tests**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="*arge*,*Gas*,*gas*,*palette*"`
Expected: PASS.

- [ ] **Step 6: Write the failing scene test**

Append to `tests/test_scene.cpp`:

```cpp
TEST_CASE("Large-scale gas keys are read and checked")
{
	// Mutation: deleting the arm_population check for large_glow in ReadDisk reddens the last check.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk.gas]\n";
	Expected<Scene> const scene = ParseScene(head + "large_glow = 1.5\nlarge_wavelength = 120.0\n", "test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->galaxies[0].disk.gas.largeGlow == 1.5);
	CHECK(scene->galaxies[0].disk.gas.largeWavelength == 120.0);
	CHECK(ParseScene(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n", "test.toml")->galaxies[0].disk.gas.largeGlow == 0.0);
	CHECK(ErrorOf(head + "large_glow = -1.0\n").find("disk.gas.large_glow") != std::string::npos);
	CHECK(ErrorOf(head + "large_wavelength = 200000.0\n").find("must be in") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\ndensity = 0.1\ncutoff = 5000\npopulation = \"old-disk\"\n[galaxy.disk.gas]\nlarge_glow = 1.0\n")
		.find("large_glow needs arm_population") != std::string::npos);
}
```

- [ ] **Step 7: Run it to verify it fails**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="Large-scale gas keys are read and checked"`
Expected: FAIL, `large_glow` reported as an unknown key.

- [ ] **Step 8: Read the keys**

In `ReadGas`, after `palette_wavelength`:

```cpp
	STARSCAPE_TRY(reader.Number("large_glow", gas.largeGlow, 0.0, 1000.0));
	STARSCAPE_TRY(reader.Number("large_wavelength", gas.largeWavelength, 10.0, 100000.0));
```

In `ReadDisk`, after the existing `glow` check:

```cpp
		if (disk.enabled && disk.gas.largeGlow > 0.0 && disk.armPopulation < 0)
		{
			return reader.Fail("gas", "large_glow needs arm_population to be set: young stars light the gas");
		}
```

- [ ] **Step 9: Run all tests, check the mutations, commit**

Run: `nice -n 19 cmake --build Build && ctest --test-dir Build --output-on-failure`
Expected: everything passes.

Mutations:
- light the term by `midplane` instead of `youngDensity`;
- replace `LargeGasMean`'s walk with `return 1.0;`;
- drop the division by `largeMean` in `EvaluateLocal`;
- delete the `largeGlow` arm_population check.

```bash
git add src/galaxy.h src/galaxy.cpp src/scene.cpp tests/test_galaxy.cpp tests/test_scene.cpp
git commit -m "Gas: a large-scale band at giant-cloud scale (mutation checked)"
```

---

### Task 4: The bubble rim

**Files:**
- Modify: `src/galaxy.h` (struct `Dust`)
- Modify: `src/galaxy.cpp` (`BubbleRadius`, `BubbleFactor`, `RimGlow`, `EvaluateLocal`)
- Modify: `src/scene.cpp` (`ReadBubble`, `ReadDust`, `ReadGalaxy`)
- Test: `tests/test_galaxy.cpp`, `tests/test_scene.cpp`

**Interfaces:**
- Consumes: `PaletteColour` (Task 2).
- Produces:
  - in `Dust`: `double bubbleRim = 0.0`, appended after `filamentMean` so
    existing designated initialisers stay valid;
  - `ReadDust(TableReader&, Dust&, Disk const& thin)`. Task 5 extends it.

- [ ] **Step 1: Write the failing galaxy tests**

Append to `tests/test_galaxy.cpp`:

```cpp
namespace
{

// GasDisk with a bubble lit on its rim; the dust is real but faint, so the rim is the only gas.
Galaxy RimDisk(double rim)
{
	Galaxy galaxy = GasDisk(0.0);
	galaxy.dust = {.enabled = true, .extinction = 1e-6, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0,
		.bubbleCenter = {0.0, 0.0, 0.0}, .bubbleRadius = 300.0, .bubbleHeight = 800.0};
	galaxy.dust.bubbleRim = rim;
	return galaxy;
}

double RimLight(PreparedGalaxy const& prepared, glm::dvec3 p)
{
	GalaxySample sample;
	Evaluate(prepared, p, sample);
	return Luminance(sample.gas);
}

}

TEST_CASE("The rim glows on a shell just inside the bubble's edge, and nowhere else")
{
	// Mutation: centring the shell at q = 1 moves the peak out of [0.83, 0.87].
	PreparedGalaxy const prepared = PrepareGalaxy(RimDisk(4.0));
	RandomStream rng(23, Purpose::Test, 0, 0, 0);
	for (int i = 0; i < 20000; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		// q = 1 / |(d.x/300, d.y/300, d.z/800)| at unit scale; pick q away from the shell, then scale to it.
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		for (double const q : {0.2, 0.6, 0.64, 1.06, 1.5, 3.0})
		{
			CHECK(RimLight(prepared, d * (q / unit)) == 0.0);
		}
	}
	// Along radial lines through lit patches, the brightest point sits on the shell.
	int lines = 0;
	for (int i = 0; i < 2000 && lines < 50; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		double bestQ = 0.0;
		double best = 0.0;
		for (double q = 0.66; q <= 1.04; q += 0.002)
		{
			double const light = RimLight(prepared, d * (q / unit));
			if (light > best)
			{
				best = light;
				bestQ = q;
			}
		}
		if (best > 0.0)
		{
			++lines;
			CHECK(bestQ >= 0.80);
			CHECK(bestQ <= 0.90);
		}
	}
	CHECK(lines == 50);
}
```

The wisp texture W varies along a line, so a single line's peak can sit a little
off q = 0.85. The test therefore allows [0.80, 0.90] per line. Then add the
spec's tighter check on the shell profile S itself, by averaging over many lines:

```cpp
TEST_CASE("Averaged over the wall, the rim peaks at q = 0.85")
{
	// Mutation: centring the shell at q = 1 moves the averaged peak out of [0.83, 0.87].
	PreparedGalaxy const prepared = PrepareGalaxy(RimDisk(4.0));
	RandomStream rng(29, Purpose::Test, 0, 0, 0);
	std::array<double, 61> profile{};   // q from 0.70 to 1.00 in steps of 0.005
	for (int i = 0; i < 4000; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		for (size_t k = 0; k < profile.size(); ++k)
		{
			profile[k] += RimLight(prepared, d * ((0.70 + 0.005 * double(k)) / unit));
		}
	}
	size_t const peak = size_t(std::max_element(profile.begin(), profile.end()) - profile.begin());
	CHECK(0.70 + 0.005 * double(peak) >= 0.83);
	CHECK(0.70 + 0.005 * double(peak) <= 0.87);
}

TEST_CASE("The rim lights patches of the wall, not all of it")
{
	// Mutation: dropping the − 0.5 offset in the rim's mask lights most of the wall.
	PreparedGalaxy const prepared = PrepareGalaxy(RimDisk(4.0));
	RandomStream rng(31, Purpose::Test, 0, 0, 0);
	int dark = 0;
	int const count = 100000;
	for (int i = 0; i < count; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		dark += RimLight(prepared, d * (0.85 / unit)) == 0.0 ? 1 : 0;
	}
	CHECK(double(dark) / count > 0.40);
	CHECK(double(dark) / count < 0.90);
}

TEST_CASE("Without dust there is no rim")
{
	// Mutation: dropping the dust.enabled guard in RimGlow lights a bubble whose dust is gone.
	Galaxy galaxy = RimDisk(4.0);
	galaxy.dust.enabled = false;
	PreparedGalaxy const prepared = PrepareGalaxy(galaxy);
	RandomStream rng(37, Purpose::Test, 0, 0, 0);
	for (int i = 0; i < 2000; ++i)
	{
		glm::dvec3 const d = glm::normalize(glm::dvec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
		double const unit = std::sqrt((d.x * d.x + d.y * d.y) / (300.0 * 300.0) + d.z * d.z / (800.0 * 800.0));
		CHECK(RimLight(prepared, d * (0.85 / unit)) == 0.0);
	}
}
```

Add `#include <algorithm>` to the test file.

- [ ] **Step 2: Run them to verify they fail to compile**

Run: `nice -n 19 cmake --build Build`
Expected: errors that `Dust` has no member `bubbleRim`.

- [ ] **Step 3: Add the field**

In `struct Dust`, after `filamentMean`:

```cpp
	// Gas glowing on the bubble's wall, lit in patches: arcs across the sky from inside it. Its strength
	// is a multiple of the young stars' light at the wall's radius; it needs the thin disk's arm population.
	double bubbleRim = 0.0;
```

- [ ] **Step 4: Implement the rim**

In `src/galaxy.cpp`, beside the salts:

```cpp
// The rim's patches and wisps. Fixed rather than scene keys: the owner-approved nebula look was tuned on them.
constexpr uint64_t kRimMaskSalt = 0x81B1ull;
constexpr uint64_t kRimWispSalt = 0x81B0ull;
```

Split `BubbleFactor` around a shared radius:

```cpp
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
```

After `LargeGasMean`:

```cpp
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
```

In `EvaluateLocal`, after `addDisk(galaxy.thickDisk);`:

```cpp
	if (r <= galaxy.disk.cutoff)
	{
		out.gas += RimGlow(galaxy, radius, p);
	}
```

- [ ] **Step 5: Run the galaxy tests**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="*rim*,*Rim*,*bubble*"`
Expected: PASS, including `"A dust bubble empties its core and leaves the outside alone"`.

- [ ] **Step 6: Write the failing scene test**

Append to `tests/test_scene.cpp`:

```cpp
TEST_CASE("The bubble's rim is read and needs young stars")
{
	// Mutation: deleting the arm_population check in ReadBubble reddens the last check.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.dust]\n";
	Expected<Scene> const scene = ParseScene(head + "bubble = { center = [8200.0, 0.0, 20.0], radius = 300.0, height = 800.0, rim = 4.0 }\n", "test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->galaxies[0].dust.bubbleRim == 4.0);
	CHECK(ParseScene(head + "bubble = { center = [8200.0, 0.0, 20.0], radius = 300.0, height = 800.0 }\n", "test.toml")->galaxies[0].dust.bubbleRim == 0.0);
	CHECK(ErrorOf(head + "bubble = { center = [0.0, 0.0, 0.0], radius = 300.0, height = 800.0, rim = 1001.0 }\n").find("dust.bubble.rim") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\ndensity = 0.1\ncutoff = 5000\npopulation = \"old-disk\"\n[galaxy.dust]\nbubble = { center = [0.0, 0.0, 0.0], radius = 300.0, height = 800.0, rim = 1.0 }\n")
		.find("rim needs the thin disk's arm_population") != std::string::npos);
}
```

- [ ] **Step 7: Run it to verify it fails**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="The bubble's rim is read and needs young stars"`
Expected: FAIL, `bubble.rim` reported as an unknown key.

- [ ] **Step 8: Read the key**

Change `ReadBubble` and `ReadDust` to take the thin disk. The thin disk is read
before the dust, so its `armPopulation` is known.

```cpp
Expected<void> ReadBubble(TableReader& parent, Dust& dust, Disk const& thin)
{
	...existing body up to the semi-axis check...
	STARSCAPE_TRY(reader.Number("rim", dust.bubbleRim, 0.0, 1000.0));
	if (dust.bubbleRim > 0.0 && thin.armPopulation < 0)
	{
		return reader.Fail("rim", "rim needs the thin disk's arm_population: young stars light the wall");
	}
	return reader.Finish();
}
```

In `ReadDust`, change the signature to `Expected<void> ReadDust(TableReader& parent, Dust& dust, Disk const& thin)`
and the call inside it to `STARSCAPE_TRY(ReadBubble(reader, dust, thin));`. In
`ReadGalaxy`, change the call to `STARSCAPE_TRY(ReadDust(reader, galaxy.dust, galaxy.disk));`.

- [ ] **Step 9: Run all tests, check the mutations, commit**

Run: `nice -n 19 cmake --build Build && ctest --test-dir Build --output-on-failure`
Expected: everything passes.

Mutations:
- shell centred at 1.0;
- drop `- 0.5` in the mask;
- drop `!dust.enabled`;
- delete the rim's arm_population check.

```bash
git add src/galaxy.h src/galaxy.cpp src/scene.cpp tests/test_galaxy.cpp tests/test_scene.cpp
git commit -m "Dust: the bubble's rim glows in patches (mutation checked)"
```

---

### Task 5: Reflection

**Files:**
- Modify: `src/galaxy.h` (struct `Dust`)
- Modify: `src/galaxy.cpp` (`addDisk` returns the young density, the term, `PrepareGalaxy`)
- Modify: `src/scene.cpp` (`ReadDust`)
- Test: `tests/test_galaxy.cpp`, `tests/test_scene.cpp`

**Interfaces:**
- Consumes: `ReadDust(TableReader&, Dust&, Disk const& thin)` (Task 4).
- Produces: in `Dust`, after `bubbleRim`:
  - `double reflection = 0.0`;
  - `glm::dvec3 reflectionColour{0.45, 0.62, 1.0}`, normalised to Y = 1 by
    `PrepareGalaxy`.

- [ ] **Step 1: Write the failing galaxy tests**

Append to `tests/test_galaxy.cpp`:

```cpp
namespace
{

// GasDisk under uniform dust that reflects; no other gas.
Galaxy ReflectingDisk(double extinction)
{
	Galaxy galaxy = GasDisk(0.0);
	galaxy.dust = {.enabled = true, .extinction = extinction, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 1e9, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0};
	galaxy.dust.reflection = 2.0;
	galaxy.dust.reflectionColour = {0.9, 1.24, 2.0};   // twice the default: normalisation must bring it to Y = 1
	return galaxy;
}

}

TEST_CASE("Reflection tracks the dust and the young stars")
{
	// Mutation: using the dust's filament field without α₀ (out.extinction / dust.extinction) breaks the doubling.
	glm::dvec3 const p(300.0, -200.0, 100.0);
	GalaxySample once;
	GalaxySample twice;
	Evaluate(PrepareGalaxy(ReflectingDisk(0.001)), p, once);
	Evaluate(PrepareGalaxy(ReflectingDisk(0.002)), p, twice);
	REQUIRE(Luminance(once.gas) > 0.0);
	CHECK(Luminance(twice.gas) == doctest::Approx(2.0 * Luminance(once.gas)).epsilon(1e-12));
	// At 1 mag/kpc the reflected light is `reflection` times the young stars' own.
	CHECK(Luminance(once.gas) == doctest::Approx(2.0 * once.parts[1].density * 100.0 * kFluxAt10pc).epsilon(1e-12));
	Galaxy dark = ReflectingDisk(0.001);
	dark.dust.enabled = false;
	GalaxySample none;
	Evaluate(PrepareGalaxy(dark), p, none);
	CHECK(none.gas == glm::dvec3(0.0));
}

TEST_CASE("Reflected light takes its normalised colour")
{
	// Mutation: skipping the reflection colour's normalisation doubles the luminance.
	PreparedGalaxy const prepared = PrepareGalaxy(ReflectingDisk(0.001));
	CHECK(Luminance(prepared.galaxy.dust.reflectionColour) == doctest::Approx(1.0).epsilon(1e-12));
	GalaxySample sample;
	Evaluate(prepared, {300.0, -200.0, 100.0}, sample);
	glm::dvec3 const expected = prepared.galaxy.dust.reflectionColour * Luminance(sample.gas);
	CHECK(sample.gas.x == doctest::Approx(expected.x).epsilon(1e-12));
	CHECK(sample.gas.y == doctest::Approx(expected.y).epsilon(1e-12));
	CHECK(sample.gas.z == doctest::Approx(expected.z).epsilon(1e-12));
}
```

The second `CHECK` in the first test derives as follows. With uniform dust,
α_V = 0.001 mag/pc, so α_V · 1000 = 1. `GasDisk` has young density
`parts[1].density`, and young light 1. The emissivity's luminance is therefore
`reflection` · 1 · 1 · n_young · 100 · kFluxAt10pc.

- [ ] **Step 2: Run them to verify they fail to compile**

Run: `nice -n 19 cmake --build Build`
Expected: errors that `Dust` has no member `reflection` or `reflectionColour`.

- [ ] **Step 3: Add the fields**

In `struct Dust`, after `bubbleRim`:

```cpp
	// Dust lit by nearby young stars scatters their light. `reflection` is that light as a multiple of the
	// young stars' own where the dust runs at 1 mag/kpc. Added light, a stylisation: the stars lose nothing.
	double     reflection = 0.0;
	glm::dvec3 reflectionColour{0.45, 0.62, 1.0};   // PrepareGalaxy normalises it to luminance 1
```

- [ ] **Step 4: Implement it**

In `PrepareGalaxy`, after `NormaliseColours(prepared.galaxy.disk.gas);`:

```cpp
	prepared.galaxy.dust.reflectionColour /= Luminance(prepared.galaxy.dust.reflectionColour);
```

In `EvaluateLocal`, change `addDisk` to return the disk's young density:
- give the lambda a `-> double` return type;
- make the early return `return 0.0;`;
- end it with `return youngDensity;`.

Then replace the two calls with:

```cpp
	// The thin disk's young stars light the reflection below.
	double const youngDensity = addDisk(galaxy.disk);
	addDisk(galaxy.thickDisk);
```

At the end of the dust block, after `out.extinction = ...;`:

```cpp
		if (dust.reflection > 0.0 && youngDensity > 0.0)
		{
			// α_V in magnitudes per kiloparsec, so `reflection` is a multiple of the young light at 1 mag/kpc.
			double const light = dust.reflection * (out.extinction * 1000.0) * galaxy.disk.gas.youngLight * youngDensity;
			out.gas += (light * 100.0 * kFluxAt10pc) * dust.reflectionColour;
		}
```

- [ ] **Step 5: Run the galaxy tests**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="*eflect*"`
Expected: PASS.

- [ ] **Step 6: Write the failing scene test**

Append to `tests/test_scene.cpp`:

```cpp
TEST_CASE("Reflection keys are read and need young stars")
{
	// Mutation: deleting the arm_population check for reflection in ReadDust reddens the last check.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.dust]\n";
	Expected<Scene> const scene = ParseScene(head + "reflection = 2.0\nreflection_colour = [0.4, 0.6, 1.0]\n", "test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->galaxies[0].dust.reflection == 2.0);
	CHECK(scene->galaxies[0].dust.reflectionColour == glm::dvec3(0.4, 0.6, 1.0));
	CHECK(ParseScene(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n", "test.toml")->galaxies[0].dust.reflection == 0.0);
	CHECK(ErrorOf(head + "reflection = -1.0\n").find("dust.reflection") != std::string::npos);
	CHECK(ErrorOf(head + "reflection_colour = [0.0, 0.0, 0.0]\n").find("positive luminance") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\ndensity = 0.1\ncutoff = 5000\npopulation = \"old-disk\"\n[galaxy.dust]\nreflection = 1.0\n")
		.find("reflection needs the thin disk's arm_population") != std::string::npos);
}
```

- [ ] **Step 7: Run it to verify it fails**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="Reflection keys are read and need young stars"`
Expected: FAIL, `reflection` reported as an unknown key.

- [ ] **Step 8: Read the keys**

In `ReadDust`, after `filament`:

```cpp
	STARSCAPE_TRY(reader.Number("reflection", dust.reflection, 0.0, 1000.0));
	STARSCAPE_TRY(ReadColour(reader, "reflection_colour", dust.reflectionColour));
	if (dust.enabled && dust.reflection > 0.0 && thin.armPopulation < 0)
	{
		return reader.Fail("reflection", "reflection needs the thin disk's arm_population: young stars light the dust");
	}
```

- [ ] **Step 9: Run all tests, check the mutations, commit**

Run: `nice -n 19 cmake --build Build && ctest --test-dir Build --output-on-failure`
Expected: everything passes.

Mutations:
- `out.extinction / dust.extinction` in place of `out.extinction`;
- drop the colour normalisation;
- delete the reflection arm_population check.

```bash
git add src/galaxy.h src/galaxy.cpp src/scene.cpp tests/test_galaxy.cpp tests/test_scene.cpp
git commit -m "Dust: reflection of the young stars' light (mutation checked)"
```

---

### Task 6: The terms in the march, the stylized scene, and inert defaults

**Files:**
- Modify: `tests/test_march.cpp`
- Modify: `scenes/stylized-milky-way.toml`

The spec's "Defaults are inert" test is realised by two pins, not one new test:
- **Step 5's byte comparison** of the realistic scene against Task 1's reference;
- **the existing `"Gas colours are normalised and the densest clumps take the core colour"`**.
  It requires pure-red gas to have no green. With an empty palette, a palette
  path that took over would add colour and redden it.

Every other term multiplies by its strength, which defaults to zero, so no
mutation can make it non-inert at its default. A permanent test of that would be
a pin that cannot fail. Task 7 rewrites the spec's row to match.

**Interfaces:**
- Consumes: every field from Tasks 2–5.

- [ ] **Step 1: Extend the stars-alone test**

In `tests/test_march.cpp`, test `"Gas adds diffuse light and leaves the stars alone"`,
replace the opening lines, from `Scene scene = MilkyWayScene(...)` through
`PreparedScene const glowing = PrepareScene(scene);`, with:

```cpp
	// Mutation: adding any term (gas, large-scale, rim or reflection) to resolvedRadiance instead of
	// radiance reddens this test. Every term is on at once, around a bubble both scenes share, so only the
	// light terms differ between them.
	Scene scene = MilkyWayScene({8200.0, 0.0, 20.0});
	Galaxy& galaxy = scene.galaxies[0];
	galaxy.dust.bubbleCenter = {8200.0, 0.0, 20.0};
	galaxy.dust.bubbleRadius = 300.0;
	galaxy.dust.bubbleHeight = 800.0;
	PreparedScene const plain = PrepareScene(scene);
	galaxy.disk.gas.glow = 1.5;
	galaxy.disk.gas.palette = {glm::dvec3(1.0, 0.18, 0.32), glm::dvec3(0.25, 1.0, 0.85)};
	galaxy.disk.gas.paletteCount = 2;
	galaxy.disk.gas.largeGlow = 1.5;
	galaxy.dust.bubbleRim = 4.0;
	galaxy.dust.reflection = 2.0;
	PreparedScene const glowing = PrepareScene(scene);
```

Keep the rest of the test as it is. The ray starts at the bubble's centre, so it
crosses the shell.

- [ ] **Step 2: Add the attenuation test for the three light terms**

Append to `tests/test_march.cpp`:

```cpp
TEST_CASE("Large-scale gas, the rim and reflection are dimmed by the dust in front of them")
{
	// Mutation: adding any one term to the gas outside the step's transmittance (e.g. straight into
	// the radiance) leaves that term undimmed, and its section fails.
	struct Light
	{
		char const* name;
		void (*apply)(Galaxy&);
	};
	Light const lights[] = {
		{"large", [](Galaxy& g) { g.disk.gas.largeGlow = 1.0; }},
		{"rim", [](Galaxy& g) { g.dust.bubbleRim = 4.0; }},
		{"reflection", [](Galaxy& g) { g.dust.reflection = 2.0; }},
	};
	for (Light const& light : lights)
	{
		CAPTURE(light.name);
		Scene scene = BaseScene();
		int const population = AddSingleBinPopulation(scene, 14);
		Galaxy source;
		source.position = {1.0e6, 0.0, 0.0};
		source.disk = {.enabled = true, .density = 0.05, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 20000.0, .arms = {},
			.population = population, .armPopulation = population, .youngFraction = 0.5};
		// Faint dust with a bubble at the source's centre, so the rim and the reflection have something to light.
		source.dust = {.enabled = true, .extinction = 1e-6, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = 20000.0, .laneOffset = 0.0, .contrast = 1.0, .filament = 0.0,
			.bubbleCenter = {0.0, 0.0, 0.0}, .bubbleRadius = 300.0, .bubbleHeight = 800.0};
		scene.galaxies.push_back(source);
		MarchParams const params = ParamsFor(Quality::Quick);
		glm::dvec3 const dir(1.0, 0.0, 0.0);
		auto termRadiance = [&](Scene const& base)
		{
			RayResult dark;
			MarchRay(PrepareScene(base), params, dir, HUGE_VAL, 1.0, 1.0, dark);
			Scene lit = base;
			light.apply(lit.galaxies[0]);
			RayResult glowing;
			MarchRay(PrepareScene(lit), params, dir, HUGE_VAL, 1.0, 1.0, glowing);
			return glowing.radiance - dark.radiance;
		};
		glm::dvec3 const clear = termRadiance(scene);
		REQUIRE(Luminance(clear) > 0.0);
		double const alpha = 0.002;
		double const radius = 1000.0;
		scene.galaxies.push_back(UniformDust(alpha, radius));
		glm::dvec3 const dimmed = termRadiance(scene);
		double const tau = alpha * (radius - scene.observer.sMin) / 1.0857;
		CHECK(dimmed.x / clear.x == doctest::Approx(std::exp(-0.75 * tau)).epsilon(0.005).scale(0.0));
		CHECK(dimmed.y / clear.y == doctest::Approx(std::exp(-1.00 * tau)).epsilon(0.005).scale(0.0));
		CHECK(dimmed.z / clear.z == doctest::Approx(std::exp(-1.32 * tau)).epsilon(0.005).scale(0.0));
	}
}
```

The source galaxy's young density at its centre is `0.05 · 0.5 · sech²(0)`. The
rim's reference density uses `youngFraction + armYoungFraction = 0.5`, so it is
non-zero.

Before relying on the rim section, check that the ray crosses a lit patch:
`REQUIRE(Luminance(clear) > 0.0)` fails if it does not. If it fails for the rim,
move `bubbleCenter` along y in steps of 50 pc until the ray hits a patch, and
note in the test's comment why the offset is there.

- [ ] **Step 3: Run the march tests**

Run: `nice -n 19 cmake --build Build && Build/starscape-tests -tc="*dimmed*,*stars alone*"`
Expected: PASS.

- [ ] **Step 4: Update the stylized scene**

In `scenes/stylized-milky-way.toml`:

Under `[galaxy.disk.gas]`, after `glow = 1.5`:

```toml
  # The nebula layer (docs/images/2026-10-15-nebula-target.png): neighbouring clouds take different hues,
  # and a second band of gas at giant-cloud scale fills the arms. large_glow is the spike's 1.5 times
  # the field's mean (1.006899 at seed 1), which the cook now divides out.
  palette            = [[1.0, 0.18, 0.32], [0.25, 1.0, 0.85], [1.0, 0.70, 0.22], [0.62, 0.30, 1.0]]
  palette_wavelength = 1500.0
  large_glow         = 1.51
  large_wavelength   = 100.0
```

Under `[galaxy.dust]`:

```toml
  reflection = 2.0
```

and change the `bubble` line to:

```toml
  bubble   = { center = [8200.0, 0.0, 20.0], radius = 300.0, height = 800.0, rim = 4.0 }
```

Update the file's header comment to mention the nebula layer and its target image.

- [ ] **Step 5: Prove the existing scenes are untouched**

Run:
```bash
nice -n 19 cmake --build Build
mkdir -p Build/after
Build/starscape cook scenes/milky-way-from-earth.toml -o Build/after/milky-way-from-earth.kre_starscape --draft
sha256sum Build/after/milky-way-from-earth.kre_starscape
grep milky-way-from-earth Build/reference/pre-nebula.sha256
```
Expected: the two hashes are equal. The stylized scene's hash changes by design.
A different realistic hash is a defect: find the term that reached it, fix it,
and add a unit test that pins it.

- [ ] **Step 6: Run the whole suite**

Run: `ctest --test-dir Build --output-on-failure`
Expected: everything passes, including `cli-draft-stylized`.

- [ ] **Step 7: Commit**

```bash
git add tests/test_march.cpp scenes/stylized-milky-way.toml
git commit -m "Scenes: the stylized sky takes the nebula layer; new terms pinned in the march (mutation checked)"
```

---

### Task 7: Visual acceptance, the spec fold-in, and cleanup

**Files:**
- Modify: `docs/spec.md`
- Delete: `docs/specs/2026-10-15-nebula-layer.md`

- [ ] **Step 1: Render the stylized scene beside the target**

Run:

```bash
Build/starscape cook scenes/stylized-milky-way.toml -o Build/accept/stylized.kre_starscape --quick --width 2048 --preview Build/accept/stylized.hdr
```

Write `Build/accept/compare.py` with this content, then run
`python3 -I Build/accept/compare.py Build/accept/stylized.hdr docs/images/2026-10-15-nebula-target.png Build/accept/compare.png`:

```python
import sys, numpy as np
from PIL import Image, ImageDraw
def load(path):
    data = open(path, 'rb').read()
    pos = data.index(b'\n\n') + 2
    end = data.index(b'\n', pos)
    dims = data[pos:end].split()
    h, w = int(dims[1]), int(dims[3])
    pos = end + 1
    out = np.zeros((h, w, 4), np.uint8)
    for y in range(h):
        if data[pos] == 2 and data[pos + 1] == 2:
            pos += 4
            for c in range(4):
                x = 0
                while x < w:
                    n = data[pos]; pos += 1
                    if n > 128:
                        n -= 128; out[y, x:x + n, c] = data[pos]; pos += 1
                    else:
                        out[y, x:x + n, c] = np.frombuffer(data[pos:pos + n], np.uint8); pos += n
                    x += n
        else:
            out[y] = np.frombuffer(data[pos:pos + 4 * w], np.uint8).reshape(w, 4); pos += 4 * w
    e = out[..., 3].astype(np.float64)
    return out[..., :3] * np.where(e > 0, np.ldexp(1.0, (e - 136).astype(np.int32)), 0.0)[..., None]
hdr, target, dest = sys.argv[1:4]
rgb = load(hdr); lum = rgb @ [0.2126, 0.7152, 0.0722]
floor, ref = np.percentile(lum, 5), np.percentile(lum, 99.9)
new = Image.fromarray((np.clip((rgb - floor) / (ref - floor), 0, 1) ** (1 / 2.2) * 255 + 0.5).astype(np.uint8))
old = Image.open(target).convert('RGB').resize(new.size, Image.LANCZOS)
for im, label in ((old, 'target (spike)'), (new, 'stylized-milky-way.toml')):
    ImageDraw.Draw(im).text((8, 6), label, fill=(255, 255, 0))
canvas = Image.new('RGB', (new.width, new.height * 2)); canvas.paste(old, (0, 0)); canvas.paste(new, (0, new.height)); canvas.save(dest)
```

Print `UNVERIFIED-VISUAL: Build/accept/compare.png awaits the owner`
and hand the image and `Build/accept/stylized.hdr` to the owner, for two checks:
1. Does it match the target?
2. In Blender, with the HDR as the world, does the rim read as part of the sky?

The owner's answer decides Step 2. Do not continue without it.

- [ ] **Step 2: Apply the owner's verdict on the rim**

- **If the rim passed:** nothing to do.
- **If it failed:** set `rim = 0.0` in `scenes/stylized-milky-way.toml` (the term
  stays available), re-run Step 1's render for the record, and commit with the
  message `Scenes: the rim stays out of the stylized sky (owner, Blender check)`.

- [ ] **Step 3: Fold the spec into `docs/spec.md`**

Edit `docs/spec.md`:
- **§1 Gas:** after the existing *In the march* text, add the spec's §1–§4
  (shared pieces, rim, large-scale gas, reflection) as subsections *Palette*,
  *Large-scale gas*, *Bubble rim* and *Reflection*. Keep their formulas, and put
  the salts table from the spec's §6 under the palette.
- **§4 Scene (TOML):** add the new keys to the example block, and the spec's
  §5 rules to the bullet list. Extend the stylized-scene bullet with the new
  values and `docs/images/2026-10-15-nebula-target.png`.
- **§5 Testing:** add one row per test added in Tasks 2–6, in the existing
  table's format. Two rows differ from the draft spec:
  - **Palette:** each entry is nearest for at least 1% of points, not 10%
    (Task 2's comment gives the reason).
  - **"Defaults are inert":** this row becomes the realistic scene's byte
    comparison at landing, plus the existing gas-colour test (Task 6 explains
    why no permanent test is possible).
- **Non-goals:** change "Placed or near nebulae: …" to keep the ruling, and add
  that the bubble rim is near gas admitted as part of the medium, owner-checked
  in a 3D viewer on the date of Step 1.

Then delete the spec file: `git rm docs/specs/2026-10-15-nebula-layer.md`.
`docs/specs/` is then empty again; remove the directory if git leaves it.

- [ ] **Step 4: Run the whole suite again and commit**

Run: `ctest --test-dir Build --output-on-failure`
Expected: everything passes.

```bash
git add docs/spec.md
git commit -m "Spec: fold in the nebula layer"
```

- [ ] **Step 5: Hand back for the owner's branch review**

Report:
- the branch, `nebula-layer`, and its commit list;
- test counts against Task 1's baseline;
- the realistic-scene hash check from Task 6.

The owner reviews the branch before it fast-forwards into `main`. Merging is
not part of this plan.

After the owner has merged it, delete the spike:

```bash
git -C /mnt/Passport/Kreatures-Engine/Tooling/Starscape worktree remove --force /mnt/Passport/Worktrees/starscape-nebula-spike
git -C /mnt/Passport/Kreatures-Engine/Tooling/Starscape branch -D spike/nebula-layer
```

`--force` is needed only because the spike worktree holds an ignored `Build/`
and `Build-asan/`. Look at `git -C /mnt/Passport/Worktrees/starscape-nebula-spike status --short`
first: it must print nothing.
