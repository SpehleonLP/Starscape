# Gas glow and a stylized sky — design spec

Status: draft 2026-10-13, for owner review. Extends `docs/spec.md`; on landing,
its sections fold into that spec and this file is deleted.

## Why

The realistic sky is correct but muted: a brown-grey band. The owner wants a
more fantastical sky that still feels vast. Look-development on throwaway
prototypes settled what works:

- **Keep the galaxy's structure; exaggerate its parameters.** Space HDRIs and
  Mario Galaxy backdrops are a flat layer of noise clouds behind an even
  sprinkle of white dots. They read as painted backdrops, not as something you
  are small inside. The cues that read as scale are structural: dots to grain to
  glow, a band that says "inside a disk", a centre that is elsewhere, dust that
  occludes, and stars with a power-law brightness spread and real colour.
- **Nebulae belong to the band, not on top of it.** Separate nebula objects
  (spheres with their own shaders) read as pasted on, however they were tuned.
  Gas glow driven by the galaxy's own young-star and dust fields reads as part of
  the same sky.
- **Near nebulae are out.** A painted nebula 40° across looked right in a still
  frame, but in a 3D viewer it dominated the sky. Dropped from this spec.
- **The view is unrotated.** Example scenes use the galactic frame, with the band
  along the equator and the centre in the middle of the preview. Orientation is
  the consumer's job: rotating the skybox gives the same result and the preview
  stays legible.

The acceptance target is the owner-approved render
`docs/images/2026-10-13-stylized-target.png` (prototype panel `q`).

## Scope

1. A gas-glow emission term on disks with a young population.
2. A star colour saturation setting.
3. A stylized example scene matching the target, plus the galactic frame for both
   example scenes.
4. Deletion of the nebula prototype.

Non-goals: near or placed nebulae, absorption by gas (the dust field already
absorbs), reflection nebulae, tone mapping or a colour-grading LUT (parked),
constellations (their own spec).

## 1. Gas glow

Gas lit by young stars glows where they cluster. It shares the dust's filaments,
because gas and dust are one medium, and it shares the dust's bubble cavity.

### Model

For a thin disk with an `arm_population`, at galaxy-frame point x:

    j_gas(x) = glow · L̄ · n_young(x) · G(x) / ⟨G⟩ · B(x) · c(x) · 100 · 10^(−1.932)

- **n_young** is the young stars' density exactly as the disk already computes
  it: midplane × young share × association field × sech²(z / 2h_young).
- **L̄ = Σ f_i L_i** is the arm population's mean V luminosity per star. So
  `glow` is the gas's light as a multiple of the young stars' own diffuse light,
  averaged over space before the cavity.
- **G** is the gas texture:

      G = (0.6 · D + W) · sqrt(F / ⟨F⟩)
      D = LN(1.2, fbm(x / 40 pc, 5 octaves) / σ_fbm)        // soft clumps, mean one
      W = 6 · max(0, ridged(x / 25 pc, 5 octaves))⁶          // thin bright wisps

  F is the dust filament field at the dust's own wavelength and seed, so gas
  knots sit in the dust's filaments. If the galaxy has no dust, F/⟨F⟩ = 1. The
  noise uses `galaxy.seed` with its own salts.
- **⟨G⟩** is G's spatial mean, measured in `PrepareGalaxy` with the same 2¹⁸-point
  `SpatialMean` walk as the dust's ⟨LN⟩, and only when `glow` > 0. This keeps
  `glow` meaning the same thing at any texture setting.
- **B** is the dust bubble's cavity factor (zero within 70% of its semi-axes,
  smoothstep to one at the edge), and one without a bubble.
- **c** is the colour: mix(ĉ, ĉ_core, smoothstep(1, 3, D)). ĉ and ĉ_core are the
  scene's `colour` and `core_colour`, each normalised to Rec.709 Y = 1, so the
  colour never changes the light's amount. Dense clumps turn teal.

The emissivity is in the same units as a star bin's j_i (spec §1, Normalisation).

### In the march

`EvaluateLocal` returns the summed gas emissivity of a point as an RGB value on
`GalaxySample`. L̄ reaches the galaxy as a derived field that `PrepareScene` fills
from the population tables, so the galaxy code needs no population lookup.
`Walk` sums it over galaxies onto the `Step`. `MarchRay` adds

    radiance += weight · gas · T · ds

The gas is diffuse light only:

- the flux histogram, the threshold F_res and the star draws ignore it;
- the cook's diffuse-flux statistic includes it;
- `valueScale` already adapts to the brightest stored value.

Thick disks and other components without an arm population carry no gas.

### Scene keys

```toml
  [galaxy.disk.gas]
  glow        = 0.0                  # off by default in every preset
  colour      = [1.0, 0.18, 0.32]    # hydrogen red; normalised to Y = 1
  core_colour = [0.25, 1.0, 0.85]    # dense clumps; normalised to Y = 1
```

- The keys are under the thin disk's table (`[galaxy.disk.gas]`), like `arms`. The thick disk has no young
  population to light gas, so a `gas` table there is an unknown key.
- `glow` ranges over [0, 1000].
- Colour components must be ≥ 0, with Y > 0.
- `glow` > 0 on a disk without `arm_population` is an error, as
  `arm_young_fraction` already is.
- `glow` = 0 (every preset) makes the term exactly zero, so existing skies are unchanged.

## 2. Star colour saturation

```toml
[style]
star_saturation = 1.0    # 1: physical; >1 pushes colours away from white
```

- The setting applies in `Tabulate`, to every population's bin chromaticity:

      ĉ' = max(0, 1 + s · (ĉ − 1)),  then renormalised to Y = 1

  Chromaticities have Y = 1, so white is (1, 1, 1).
- Points and diffuse light both use the tables, so both change together.
- Each bin keeps its luminance and the resolution split reads only luminance and T_V, so star counts,
  F_res and the point budgets are unchanged. Luminance is exactly unchanged only where the transmittance
  is grey: dust dims blue more than red, so a bluer star loses a little more light behind dust.
- Range [0, 4]. 0 gives grey stars, which is harmless.
- `[style]` is a new top-level table for scene-wide stylisation, so that later
  knobs have a home.

## 3. Example scenes

- **`scenes/milky-way-from-earth.toml`** keeps its realism but drops `spin_axis`
  and `zero_meridian` in favour of the galactic frame:
  - `spin_axis = [0, 0, 1]`;
  - `zero_meridian` pointing at the galactic centre.

  The keys themselves stay in the format.
- **New `scenes/stylized-milky-way.toml`:** the Earth scene in the galactic frame
  plus these settings:

  | Key | Value |
  |---|---|
  | `[galaxy.dust] contrast` | 6.0 |
  | `[galaxy.dust] filament` | 3.5 |
  | `[galaxy.disk.arms] contrast` | 2.25 |
  | `[galaxy.disk.gas] glow` | re-tuned (see below) |
  | `[style] star_saturation` | 1.5 |

  The prototype ran at gas gain 2.5 without dividing by ⟨G⟩. Its ⟨G⟩ at
  `filament` 3.5 and seed 1 measures 0.6045, so `glow = 1.5` reproduces the
  target's brightness. The plan renders the scene beside the target for the owner.

## 4. Deleting the prototype

The working tree carries uncommitted look-dev code. None of it lands. The
implementation starts from the committed tree (5fa8ec4) and must not contain:

- `src/nebula.h` and `src/nebula.cpp`, and their CMake entry;
- `BillowFbm`;
- `Chord`, `BuildChords`, `NebulaDepth`, `tChords` and `Step::dustDepth`;
- the `depthAt` bookkeeping in `MarchRay`;
- every `NEB_*`, `GAS_GAIN` and `CHROMA` environment read.

`docs/spec.md`'s non-goal "Emission nebulae (HII regions)" becomes "Placed or near
nebulae" when this spec folds in.

## 5. Testing

The tests follow the repository's rule: each must be shown to fail under its
listed mutation, recorded in the test's comment.

| Test | Assertion | Mutation that must fail it |
|---|---|---|
| Gas calibration | Uniform young disk, dust filaments at `filament` 2 (⟨G⟩ = 0.88), 2 × 10⁵ random points: Σ gas luminance / Σ young-star diffuse luminance = `glow` within 3% | skip the ⟨G⟩ division |
| Gas cavity | Points inside 70% of the bubble's semi-axes, with young stars present: gas emissivity exactly 0 | leave the cavity factor off the gas |
| Gas colour | Prepared colours have Y = 1 and keep their direction; with pure red and pure blue colours, samples that are purely `colour` (D < 1, ~73%) outnumber those purely `core_colour` (D > 3, ~6.5%) by more than 5× and both occur | swap the mix's endpoints |
| Gas attenuation | Gas region behind a dust slab of optical depth τ: per-channel gas radiance ratio with/without the slab = exp(−k_c τ) within 0.5% | add gas without the step's transmittance |
| Gas leaves the stars alone | One Milky Way ray with `glow` 1.5 vs 0: radiance rises, while point budgets, resolved radiance and the flux histogram are identical | add the gas to `resolvedRadiance` instead of `radiance` |
| Saturation default | An absent `[style]` gives `star_saturation` 1, the physical sky | a default other than 1 |
| Saturation keeps luminance | s = 4 (so some bins clip at zero): every bin's chroma has Y = 1 to 10⁻¹²; unclipped bins sit 4× as far from white; a dust-free Milky Way ray keeps its luminance to 10⁻¹² and its point budgets exactly, while its channels spread apart | skip the renormalisation |
| Scene keys | Range and type errors for `glow`, `colour`, `core_colour` and `star_saturation`; `glow` > 0 without `arm_population` names the key | delete any one validation |
| Example scenes | The CLI draft-cooks both files under `scenes/` (the Earth scene already has `cli-draft`; the stylized one gets its own) | — |

**Visual acceptance:** the stylized scene's preview beside the target image, for
the owner to judge. This is not a gate, so the test log prints an
UNVERIFIED-VISUAL marker for it.
