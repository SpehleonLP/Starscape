# Nebula layer: rim, palette, large-scale gas and reflection — design spec

Status: draft 2026-10-15, for owner review. Extends `docs/spec.md`; on landing,
its sections fold into that spec and this file is deleted.

## Why

The stylized sky (`docs/images/2026-10-13-stylized-target.png`) looks good but
not fantastical. The owner asked for "another layer" and chose, from a spike, a
bigger and more colourful nebula layer over a second galaxy, a planet or a
runtime overlay:

- A second galaxy already works: scenes take several `[[galaxy]]` tables and the
  march handles them together. Distant galaxies look poor, which is a separate
  problem and out of scope here.
- A planet or moon would be animated, so it belongs to the runtime, not the cook.

**Nebulae stay part of the medium.** `docs/spec.md` §Gas rules out placed nebula
objects because they read as pasted on top of the band. Every term here is driven
by fields the galaxy already has (young stars, dust, the bubble). None of them
places an object.

The spike (`spike/nebula-layer` at 5d09f61, a local branch never meant to land)
cooked a matrix of the four terms. The owner judged them with the stylized
target's tone mapping (black at the 5th percentile, white at the 99.9th, gamma
2.2) and chose all four, with large-scale gas and reflection turned down:

| Term | What the owner saw | Chosen setting |
|---|---|---|
| Bubble rim | Coloured arcs above and below the band | strength 4 |
| Palette | Neighbouring clouds in red, teal, gold and violet | on, 1.5 kpc field |
| Large-scale gas | Bigger pink and magenta clouds through the band | 100 pc at strength 1.5 |
| Reflection | At 8 it bleached the band white; at 2 it adds a faint blue lift | strength 2 |

The acceptance target is `docs/images/2026-10-15-nebula-target.png`, a quick cook
at face 256 with exactly those settings.

**The near-nebula risk.** `docs/spec.md` rules out near nebulae: a painted
nebula 40° across looked right in a still frame but dominated the sky in a 3D
viewer. The rim is near gas, wrapped around the observer. It differs in being
faint, patchy and lit by the medium's own fields, but that does not prove it
survives a 3D viewer. Its acceptance therefore includes the owner looking at it
in Blender (§6).

## Scope

1. A bubble-rim emission term on the dust bubble.
2. A palette that replaces the gas's single `colour` with a cyclic ramp chosen by
   a large-scale field.
3. A large-scale gas term on the thin disk.
4. A reflection term on the dust.
5. The stylized scene gains all four at the chosen settings.

Every term defaults to off, so every existing scene cooks exactly as before.

Non-goals:
- Distant galaxies' look.
- The resolution question (face size, encodings).
- The GPU cook (`2026-10-15-gpu-march.md`). This spec lands on the CPU first, and
  the GPU cook ports whatever field model is current.
- Scattering that conserves energy: reflection here is added light, a stylisation.
- Tone mapping.

## 1. Shared pieces

### Reference young light

Terms 1, 3 and 4 scale with the young stars' light, so that their strengths mean
something at any galaxy setting, exactly like `glow`. With L̄ the arm
population's mean V luminosity per star (already derived for the gas) and n_young
the young density at x as the disk computes it, the young light density is
L̄ · n_young(x).

Every emissivity below is in the units of a star bin's j_i and enters the march
like the gas: added to `GalaxySample::gas`, then `radiance += weight · gas · T · ds`.
None of these terms touches the flux histogram, the threshold F_res, the point
budgets or the star draws.

### Palette colour

    P(x) = ramp(z(x)),   z = fbm(x / λ_p, 3 octaves) / σ_fbm

- **ramp** walks the palette's N entries cyclically. The index is
  u = (0.9 · z) mod N, with mod giving a result in [0, N). Between entry
  ⌊u⌋ and the next, the weight is the smoothstep of u's fractional part.
  The 0.9 means a one-σ change in the field moves about one palette entry.
- Every entry is normalised to Rec.709 Y = 1 when the scene is prepared, so a
  blend of two entries also has Y = 1: colour never changes the amount of light.
- λ_p is `palette_wavelength`, 1500 pc in the target. The noise uses
  `galaxy.seed` with its own salt.

Without a palette, P(x) is the gas's `colour`, so the terms that use P need no
special case.

## 2. Bubble rim

Gas lit on the bubble's wall: a shell just inside the dust's ramp, so most of
it sits in front of the wall's dust.

    j_rim(x) = rim · L̄ · n_ref(x) · S(q) · M(x) · W(x) · P(4x)

- **q** is the bubble's normalised radius:
  q = sqrt((d_x² + d_y²) / r² + d_z² / h²), where d = x − bubble centre and r, h
  are the bubble's semi-axes. The dust's cavity already uses this q, with its
  ramp from 0.7 to 1.
- **S(q) = exp(−((q − 0.85) / 0.04)²)** for |q − 0.85| ≤ 0.2, and exactly 0
  outside that range, so the term costs nothing away from the wall and the
  cavity test below can be exact.
- **M(x) = clamp(fbm(x / 150 pc, 3) / σ_fbm − 0.5, 0, 1)** lights only patches
  of the wall. The eye finds arcs where a patch's edge runs along the line of
  sight.
- **W(x) = 0.3 + 6 · max(0, ridged(x / 30 pc, 5))⁶** breaks the patches into
  wisps.
- **n_ref(x) = n_0 · exp(−R / R_d) · (young_fraction + arm_young_fraction)** is
  the young density the disk would have at that galactocentric radius R, without
  arms, associations or the sech². The rim is not confined to the young disk's
  thickness: the bubble opens into chimneys, so the wall wraps the whole sky.
- **P(4x)** samples the palette at four times its frequency, which gives the rim
  its own, finer colour changes. Without a palette it is the gas `colour`.

The term belongs to a galaxy with both a bubble and a thin disk with
`arm_population` (for L̄). M and W use `galaxy.seed` with their own salts.

## 3. Large-scale gas

A second gas band, at the size of giant cloud complexes rather than single
nebulae:

    j_large(x) = large_glow · L̄ · n_young(x) · LN(1.0, fbm(x / λ_l, 4) / σ_fbm) / ⟨LN⟩ · B(x) · P(x)

- λ_l is `large_wavelength`, 100 pc in the target.
- B is the bubble's cavity factor, as for the gas.
- ⟨LN⟩ is the field's spatial mean, measured in `PrepareGalaxy` with the
  `SpatialMean` walk, and only when `large_glow` > 0. So `large_glow`, like
  `glow`, is this term's light as a multiple of the young stars' own.
- The spike did not divide by ⟨LN⟩. The stylized scene's `large_glow` is
  therefore 1.5 · ⟨LN⟩ measured at λ_l = 100 pc and seed 1, which reproduces the
  target. The plan measures ⟨LN⟩ and records the value in the scene's comment,
  as `glow` = 1.5 was handled.

## 4. Reflection

Dust lit by nearby young stars scatters their light, bluish:

    j_refl(x) = reflection · (α_V(x) · 1000 pc) · L̄ · n_young(x) · ĉ_refl

- α_V(x) is the dust's local extinction rate in magnitudes per parsec, bubble
  factor included. α_V · 1000 pc is that rate in magnitudes per kiloparsec, so
  `reflection` is the scattered light as a multiple of the young stars' own
  wherever the dust runs at 1 mag/kpc.
- n_young is the thin disk's young density at x, so dust far from young stars
  stays dark.
- ĉ_refl is `reflection_colour`, default (0.45, 0.62, 1.0), normalised to Y = 1.
- The term is added light, not light taken from the stars. Energy conservation
  across scattering is a non-goal.

## 5. Scene keys

```toml
  [galaxy.disk.gas]
  glow               = 0.0
  colour             = [1.0, 0.18, 0.32]
  core_colour        = [0.25, 1.0, 0.85]
  palette            = []                 # list of RGB; empty: `colour` alone
  palette_wavelength = 1500.0             # parsecs
  large_glow         = 0.0                # off in every preset
  large_wavelength   = 100.0              # parsecs

  [galaxy.dust]
  reflection         = 0.0                # off in every preset
  reflection_colour  = [0.45, 0.62, 1.0]
  bubble             = { center = [...], radius = 300.0, height = 800.0, rim = 0.0 }
```

- **`palette`:** 0 entries, or 2 to 8. Each entry follows the colour rules
  (components ≥ 0, Y > 0). Entries are written as an ordinary TOML array of
  arrays. When non-empty, P replaces `colour` everywhere the gas used it, and
  the core mix still applies on top: mix(P, ĉ_core, smoothstep(1, 3, D)).
- **`palette_wavelength`** and **`large_wavelength`:** in [10, 100 000] pc.
- **`large_glow`**, **`reflection`** and **`rim`:** in [0, 1000].
- `large_glow` > 0 or `rim` > 0 without an `arm_population` on the thin disk is
  an error naming the key, as `glow` already is. So is `reflection` > 0, which
  needs the thin disk's young density. `rim` lives in the bubble's inline table
  because it is meaningless without one.
- At every default the terms are exactly zero, and an empty palette takes the
  existing `colour` path, so existing scenes are byte-identical.

The stylized scene adds:

| Key | Value |
|---|---|
| `[galaxy.disk.gas] palette` | red (1.0, 0.18, 0.32), teal (0.25, 1.0, 0.85), gold (1.0, 0.70, 0.22), violet (0.62, 0.30, 1.0) |
| `[galaxy.disk.gas] palette_wavelength` | 1500 |
| `[galaxy.disk.gas] large_glow` | 1.5 · ⟨LN⟩ (§3) |
| `[galaxy.disk.gas] large_wavelength` | 100 |
| `[galaxy.dust] reflection` | 2 |
| `[galaxy.dust.bubble] rim` | 4 |

## 6. Testing

Each test must be shown to fail under its listed mutation, recorded in the
test's comment.

| Test | Assertion | Mutation that must fail it |
|---|---|---|
| Defaults are inert | Stylized preset fields with every new key at its default, 10⁴ random points: `GalaxySample::gas` equals the existing gas term evaluated alone, exactly; and the CLI draft cook with the new keys absent equals one with each set explicitly to its default, byte for byte | evaluate any term at its default (e.g. `reflection` 0 still mixing `reflection_colour` into the gas) |
| Palette luminance | Every prepared palette entry has Y = 1 to 10⁻¹²; 10⁵ random points give P(x) with Y = 1 to 10⁻¹²; with four pure, distinct entries, each is the nearest colour for at least 10% of points | drop the per-entry normalisation |
| Rim placement | Points with \|q − 0.85\| > 0.2 have exactly zero rim emissivity; along a radial line through a lit patch, the emissivity peaks within q ∈ [0.83, 0.87] | centre the shell at q = 1 |
| Rim patchiness | Over 10⁵ points on the shell (q = 0.85), the fraction with M = 0 is between 40% and 90% | drop the − 0.5 offset in M |
| Large-scale calibration | Uniform young disk, 2 × 10⁵ random points: Σ large-scale gas luminance / Σ young-star diffuse luminance = `large_glow` within 3% | skip the ⟨LN⟩ division |
| Reflection tracks dust | At fixed young density, doubling `extinction` doubles the reflection emissivity to 10⁻¹²; with dust disabled it is exactly zero | use the dust's filament field without α₀ |
| Reflection colour | Reflection-only emissivity is proportional to normalised `reflection_colour` with Y = 1 | skip its normalisation |
| New terms leave the stars alone | One Milky Way ray with all four terms on vs off: radiance rises, while point budgets, resolved radiance and the flux histogram are identical | add any term to `resolvedRadiance` |
| New terms are attenuated | Each term's region behind a dust slab of optical depth τ: per-channel ratio with/without the slab = exp(−k_c τ) within 0.5% | add one term without the step's transmittance |
| Scene keys | Range and type errors for each new key; a palette of 1 or 9 entries; `large_glow`, `rim` or `reflection` > 0 without `arm_population`, each naming its key | delete any one validation |
| Example scenes | The CLI draft-cooks the updated stylized scene | — |

**Visual acceptance (not gates; the log prints UNVERIFIED-VISUAL):**
1. The stylized scene's preview beside `docs/images/2026-10-15-nebula-target.png`.
2. The rim in a 3D viewer: the owner loads the stylized scene's preview HDR as
   Blender's world and looks around. The rim must read as part of the sky, not
   as a layer over it. If it fails, the rim is dropped from the stylized scene
   while the term stays available.

**Salts.** Each term's noise uses the same salt as the spike (palette
`0x9A1E77E`, rim mask `0x81B1`, rim wisps `0x81B0`, large-scale `0xB16B16`).
Otherwise the target image cannot be reproduced.

## 7. Deleting the spike

None of `spike/nebula-layer` lands. The implementation starts from `main` and must
not contain any `SPIKE_*` environment read, `SpikeKnobs`, `SpikePalette` or
`SpikeBubbleQ`. After landing, the spike branch and its worktree
(`/mnt/Passport/Worktrees/starscape-nebula-spike`) are deleted.
