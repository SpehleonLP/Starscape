# Starscape — design spec

Status: approved 2026-10-09; refined during planning (valueScale, explicit cutoffs, arm young fraction, C++23, own .hdr writer).

Starscape is an offline command-line tool that cooks a night sky for a given
vantage point inside (or outside) one or more procedural galaxies. It writes a
`.kre_starscape` file — a diffuse HDR cubemap of the unresolved galactic light
plus a short list of resolved point stars — and a high-resolution equirectangular
Radiance `.hdr` preview.

It replaces the Kreatures engine's `Shaders/Starscape.glb`, a point cloud cooked
by an earlier experiment whose sky came out uniform. That uniformity was caused by
geometry bugs, recorded here so the new model does not repeat them:

1. The observer sat at 250 000 units from the centre of a galaxy 10 000 units
   across — about fifty galactic radii outside it.
2. The arm density ended in `+ distance / diameter`, which grows without bound,
   so far from the galaxy the density was a flat ~25 everywhere.
3. Only a ±4000 × ±500 box around the observer was sampled, so the core and arms
   were never sampled at all.
4. The hash `fract(sin(dot(p, k)) * 43758)` lost precision at those coordinates,
   corrupting the luminosity draw.

Even a correct model gives near-uniform *resolved* stars from inside a disk: the
stars bright enough to be points are local. The Milky Way's band, core and dust
lanes are unresolved light. Starscape therefore cooks both, and splits the light
between them so that nothing is counted twice or dropped.

Starscape is a public repository and depends on nothing private. Its consumers
include `include/kre_starscape.h` and read the file themselves.

## Goals

- Any vantage point: inside a disk (a band), above it (a face-on spiral), near a
  core, in a halo, or in intergalactic space.
- Scenes are additive: any number of galaxies, each with its own parameters and
  placement. Emission adds; dust attenuates whatever lies behind it, including
  light from other galaxies.
- Few point stars (a budget, default 4000); the diffuse layer carries the rest of
  the light.
- Physically checkable output: star counts, total flux and extinction match
  analytic expectations under machine-checked tests.
- Deterministic: the same scene produces byte-identical output on any thread count.
- Cook time is not a constraint. Correctness and quality are.

## Non-goals (this spec)

- The engine-side reader, shaders and removal of `Starscape.glb` (next spec).
- Block compression of the diffuse layer (BC6H). The format reserves an encoding
  value for it.
- Emission nebulae (HII regions), runtime grain, animated or time-varying skies.
- A GPU implementation. The CPU cook is the reference.

## 1. Galaxy model

### Units and frames

Distances are parsecs. Luminosities are **visual-band** luminosities in solar
units (L_V☉), which avoids bolometric corrections: a star's colour comes from its
temperature, its brightness from L_V.

- **Galaxy frame**: disk in the xy-plane, +z toward the galaxy's north pole.
- **Scene frame**: each galaxy is placed by a position and an Euler rotation
  (degrees, applied Z then Y then X).
- **Sky frame**: the observer's celestial frame. +z is the planet's spin axis
  (`observer.spin_axis`), +x is `observer.zero_meridian` projected orthogonal to
  +z, and +y = z × x. All output directions are in the sky frame; a consumer
  rotates the sky about +z by sidereal angle.

### Components

A galaxy is a sum of components. Each component has a number-density function
n(x) in stars/pc³ (galaxy frame) and names a stellar population (see *Stellar populations*).

| Component | Density n(x) | Default population |
|---|---|---|
| Thin disk | n₀ · exp(−R/R_d) · sech²(z / 2h_z) · A(R, φ) | `young-disk` in the arms, `old-disk` between them |
| Thick disk (optional) | n₀ₜ · exp(−R/R_dₜ) · sech²(z / 2h_zₜ) | `old-disk` |
| Bulge / bar | Hernquist: n_b · a / (r′ (r′ + a)³), r′ the radius after scaling the bar axes (floored at 1 pc) | `bulge` |
| Halo | n_h · (1 + r/r_h)^(−3.5) | `halo` |

R is the cylindrical radius, φ the azimuth and r the spherical radius in the
galaxy frame. Each density has a hard cutoff radius (default 6 R_d for disks,
10 a for the bulge, 30 r_h for the halo) that defines its bounding sphere; the
march uses those spheres to skip empty space (see *The march*).

**Arm modulation** A(R, φ). For m arms with pitch angle ψ, the arm phase at
radius R is θ_k(R) = φ₀ + 2πk/m + ln(R / R₀) / tan ψ. Let d be the distance from
(R, φ) to the nearest arm, measured along the circle of radius R (arc length,
parsecs). Then

    A = 1 + (contrast − 1) · exp(−d² / 2w²) · (1 + clumping · fbm(x))

clamped to ≥ 0, where w is the arm width and fbm is 4-octave gradient noise in
the galaxy frame (base wavelength = w), seeded by the galaxy's seed. The arm
fraction a = exp(−d² / 2w²) also blends the thin disk's population: weight
`arm_young_fraction` × a (default 0.05) goes to its `arm_population`
(`young-disk`), the rest to its `population` (`old-disk`). Young stars are a small
share by number but, at ~26× the mean luminosity, dominate an arm's light.

Every component has an explicit `cutoff` radius in parsecs (presets set them).
The bulge radius r′ is floored at 1 pc so an observer at the core stays finite.

### Dust

Dust is a separate field giving visual extinction per parsec, α_V(x) in
magnitudes/pc (the solar neighbourhood is about 0.001 mag/pc):

    α_V = α₀ · exp(−R/R_dust) · sech²(z / 2h_dust) · D(R, φ) · (1 + filament · fbm_dust(x))

clamped to ≥ 0. D is the arm modulation computed with arm phase offset inward by
`lane_offset · w` (lanes sit on the inner, concave edge of each arm) and its own
contrast. fbm_dust has 6 octaves with base wavelength w/2, so lanes break into
filaments.

Extinction is per channel. The optical depth of channel c along a path is

    τ_c = k_c · ∫ α_V ds / 1.0857

with k = (0.75, 1.00, 1.32) for linear Rec.709 R, G, B (approximating a Cardelli
R_V = 3.1 curve at 610, 550 and 465 nm). The transmittance is T_c = exp(−τ_c);
T_V uses k = 1.

### Stellar populations

A population is a table over **24 luminosity bins**, log-spaced in L_V from
10⁻⁴ to 10⁶ L_V☉ (edges at 10^(−4 + 10i/24)). Each bin i holds:

- f_i — fraction of the population's stars in the bin (Σ f_i = 1);
- T_i — the bin's representative effective temperature, kelvin.

A bin's representative luminosity L_i is the geometric centre of its edges.
Presets are derived in `presets.cpp` from a Kroupa (2001) IMF, a piecewise
mass–luminosity relation, the population's main-sequence turnoff mass (none for
`young-disk`; ~1 M☉ for `old-disk`; ~0.85 M☉ for `bulge` and `halo`) and a giant
branch (stars at the turnoff mass, spread over L_V 10–1000 at T 4000–4900 K, with
a number fraction set by the giant-branch lifetime). The derivation is a comment
beside each table; the tables are data, not computed at runtime.

Colour: the bin's chromaticity is the blackbody spectrum at T_i integrated
against the CIE 1931 2° observer, converted to linear Rec.709 RGB and normalised
so that Y = 1. A star's RGB flux is its V-band flux times this chromaticity.

### Normalisation

Flux unit: a magnitude-0 star. A source of luminosity L (L_V☉) at distance s
(pc) behind transmittance T_V has flux

    F = L · T_V · (10 / s)² · 10^(−0.4 · 4.83)

since the Sun's absolute visual magnitude is 4.83. Diffuse radiance is in
mag-0 flux per steradian. The emissivity of bin i of a component is
j_i(x) = n(x) · f_i · L_i · 100 · 10^(−1.932) per steradian per pc: a volume
element s² dΩ ds holds n dV stars whose fluxes each carry 1/s², so the s²
cancels and radiance accumulates as Σ j_i T ds. (The 4π is already inside the
solar normalisation.)

## 2. Splitting the light between points and diffuse

### Threshold

One apparent-flux threshold F_res. At distance s along a ray, bin i's stars are
**resolved** if L_i · T_V(s) · (10/s)² · 10^(−1.932) ≥ F_res and **diffuse**
otherwise. Distance and dust handle themselves: a supergiant is resolved across
kiloparsecs, a red dwarf only within a few parsecs, another galaxy's stars almost
never.

F_res comes from the scene: either `limiting_magnitude` m_lim directly
(F_res = 10^(−0.4 m_lim)), or `target_stars` N, in which case a **calibration
pass** finds it. The calibration pass is the march of *The march* at face size 64 with
one sample per texel; instead of classifying, it histograms each step's expected
star count by log₁₀ of the step's apparent flux (512 bins over 30 decades). The
cumulative count from the bright end gives N(>F); F_res is the smallest F with
N(>F) ≤ N, interpolated within the bin.

### The march

The cubemap has 6 · S² texels (S = `face_size`). For each texel the cook
traces `spp` stratified, jittered directions within the texel (spp = 16 at
`final`, 1 at `quick`). Along each direction:

1. Intersect the ray with every component's bounding sphere and merge the
   intervals. Gaps between intervals are empty: no emission, no dust.
2. Within each interval, step from s_start = max(s_min, interval start) with
   step ds = clamp(α · s, ds_min, ds_max) (α = 0.02, ds_min = 0.01 pc,
   ds_max = 5 pc at `final`; α = 0.1, ds_max = 50 pc at `quick`), evaluating
   at the midpoint of each step.
3. Maintain T_c(s) by integrating α_V with the same steps.
4. At each step, for each galaxy component and each bin i:
   - if resolved: add n(x) · f_i · s² · dΩ_sample · ds to the texel's **point
     budget** B[texel][i], where dΩ_sample is the texel's solid angle / spp;
   - otherwise: add j_i(x) · T_c(s) · chroma_i,c · ds / spp to the texel's
     **diffuse radiance**, channel c.

s_min defaults to 0.05 pc so the observer's own star system is excluded; the
engine draws its own sun. The march ends at the last interval's exit.

### Drawing the stars

For each texel and bin, the star count is Poisson(B[texel][i]). For the k-th
drawn star:

- its direction is uniform in the texel's solid angle (rejection-sampled on the
  cube face against the solid-angle Jacobian);
- its distance is sampled by **re-marching that direction** and inverting the
  cumulative resolved count along it for bin i (linear within a march step; a
  star the interpolation leaves a hair under F_res is raised to F_res) (only drawn stars re-march, so
  memory stays at one float per texel per bin, 37 MB at S = 256);
- its luminosity is the bin's L_i — the same value the classification used, so
  the star is resolved at its sampled distance by construction and the split is
  exact at the threshold (apparent magnitudes are still continuous, since
  distance and dust vary continuously);
- its flux is the formula of *Normalisation* with T_c at the sampled distance,
  giving a reddened RGB flux.

Points plus diffuse equal the full sky's light in expectation. Lowering
`target_stars` moves light from points to diffuse; it never loses light.

### Randomness

All random draws use a counter-based generator (Philox 4×32-10) keyed on
(seed, purpose, texel, bin, index). Noise fields are seeded per galaxy. No draw
depends on thread scheduling, so output is byte-identical on any thread count.
Work is distributed over texels with a shared atomic counter; each texel's
results are written to its own slot.

## 3. File format

### `.kre_starscape`

Little-endian. A fixed 64-byte header, then the diffuse payload, then the star
array. All records are written and read with one memcpy; padding is explicit
`uint8_t` arrays. Offsets are from the start of the file.

```cpp
struct kre_starscape_header {       // 64 bytes
    uint32_t magic;                 // 'KSTR' = 0x5254534B
    uint32_t version;               // 1
    uint64_t sceneHash;             // FNV-1a 64 of the scene file bytes
    uint32_t diffuseEncoding;       // 1 = RGB9E5; 2 reserved for BC6H
    uint32_t diffuseFaceSize;       // cubemap face edge in texels, power of two
    uint32_t diffuseMipCount;       // full chain down to 1×1
    uint32_t starCount;
    uint64_t diffuseOffset;
    uint64_t diffuseBytes;
    uint64_t starsOffset;
    float    fluxThreshold;         // F_res, mag-0 flux units (physical)
    float    valueScale;            // stored radiance and flux = physical × valueScale
};

struct kre_starscape_star {         // 16 bytes; uploads as a vertex buffer as-is
    float    dir[3];                // unit vector, sky frame
    uint32_t fluxRGB9E5;            // linear Rec.709 flux, shared-exponent
};
```

**Diffuse payload**: mip-major, then face (+X, −X, +Y, −Y, +Z, −Z), then rows
top to bottom, texels left to right, one `uint32_t` RGB9E5 per texel. Face
orientation follows the Vulkan cube map face selection table (Vulkan spec,
"Cube Map Face Selection"), with the sky-frame direction as the lookup vector.
Mips are 2×2 box filters weighted by texel solid angle.

RGB9E5 encoding and decoding follow the Vulkan specification's reference
algorithm for `VK_FORMAT_E5B9G9R9_UFLOAT_PACK32` (round to nearest, clamp
negative to zero, saturate at the format maximum). The same encoding packs a
star's RGB flux; the star's vertex shader unpacks it.

Stars are sorted by luminance (Y, the V-band proxy) flux, brightest first.

Sizes: S = 256 is about 2 MB with mips; S = 512 about 8 MB. The star array at
the default budget is 64 KB.

### Reader

`include/kre_starscape.h` is header-only and depends only on the C++ standard
library. It defines both structs with `static_assert`s on size and offset,
`EncodeRGB9E5` / `DecodeRGB9E5`, the face-mapping function, and

```cpp
kre_starscape::Result kre_starscape::Parse(std::span<const uint8_t> file, View& out);
```

`Parse` validates magic, version, encoding, face size, mip count against face
size, that every offset and size lies within the file and that the regions do
not overlap, and that `diffuseBytes` equals the size the header implies. It
returns an error value naming the failed check; it never asserts or throws. A
`View` holds spans into the caller's buffer and does not copy.

### Preview `.hdr`

Radiance RGBE, equirectangular, width W (default 8192) by W/2. Pixel (u, v)
centres map to longitude λ = 2π (u + ½)/W − π and latitude
β = π/2 − π (v + ½)/(W/2); the direction is (cos β cos λ, cos β sin λ, sin β)
in the sky frame.

The preview is rendered **from the cooked data**: the diffuse cubemap sampled
bilinearly with seamless face edges, plus every star splatted as a Gaussian of
FWHM 1.5 px. Each splat's weights are divided by their own sum and by each
pixel's solid angle, so a star's integrated flux is exact at every latitude.
`--exposure EV` scales by 2^EV. With `preview.march = true` the diffuse layer of
the preview is instead marched at the preview's own resolution, for use as a
standalone skybox elsewhere; the `.kre_starscape` is unchanged.

## 4. Scene file and CLI

### Scene (TOML)

```toml
[output]
face_size    = 256
target_stars = 4000           # or: limiting_magnitude = 6.5 (exactly one of the two)
quality      = "final"        # or "quick"
preview      = { width = 8192, march = false }

[observer]
position      = [8200.0, 0.0, 20.0]   # parsecs, scene frame
spin_axis     = [0.0, -0.46, 0.89]    # sky-frame +z, scene frame
zero_meridian = [1.0, 0.0, 0.0]       # sky-frame +x after orthogonalisation
s_min         = 0.05

[[galaxy]]
preset   = "milky-way"
position = [0.0, 0.0, 0.0]
rotation = [0.0, 0.0, 0.0]
seed     = 1

  [galaxy.disk]
  scale_length = 2600.0
  scale_height = 300.0
  arms = { count = 4, pitch_deg = 12.0, width = 400.0, contrast = 2.5, clumping = 0.6 }

  [galaxy.dust]
  scale_height = 100.0
  lane_offset  = 0.15
  filament     = 0.7

[[galaxy]]
preset   = "small-spiral"
position = [780000.0, 120000.0, -40000.0]
```

- Galaxy presets: `milky-way`, `small-spiral`, `elliptical`, `dwarf`. A preset
  supplies every field; the scene overrides any of them. A component table
  (`disk`, `thick_disk`, `bulge`, `halo`, `dust`) that is present enables that
  component unless it sets `enabled = false`.
- `[output] seed` (default 1) keys the jitter and star draws; galaxies keep
  their own `seed` for noise.
- Diffuse range: RGB9E5 holds full 9-bit relative precision over ~32 stops
  (1.5 × 10⁻⁵ to 65408); `valueScale` is the power of two that puts the sky's
  brightest stored value in (16384, 32768], giving the dark end every other stop.
- Population presets: `young-disk`, `old-disk`, `bulge`, `halo`; components name
  theirs with `population = "..."`.
- Unknown keys, wrong types, out-of-range values (negative scale lengths, pitch
  outside (0°, 90°), non-power-of-two face size, a zero spin axis, zero
  meridian parallel to the spin axis) and having both or neither of
  `target_stars` and `limiting_magnitude` are errors with file, line and column.

### CLI

```
starscape cook    <scene.toml> -o <out.kre_starscape> [--preview <out.hdr>] [--threads N] [--quick]
starscape preview <in.kre_starscape> -o <out.hdr> [--width N] [--exposure EV]
starscape inspect <in.kre_starscape>
```

- `cook` prints timestamped progress (phase, per-face percentage, ETA) and ends
  with F_res, the star count, and the share of total flux in points vs diffuse.
  `--quick` overrides `quality`. Output goes to a temporary file in the
  destination directory and is renamed into place on success.
- `preview` re-renders a preview from an existing file.
- `inspect` prints the header, a flux histogram of the stars by magnitude, the
  total diffuse flux (radiance integrated over solid angle) and the
  points/diffuse split.
- Exit codes: 0 success, 1 invalid input (scene or file), 2 I/O failure. Errors
  propagate as values to `main`, which prints each once.

## 5. Testing

Tests use vendored doctest and run under `ctest`. Each test below must be shown
to fail under its listed mutation before it is relied on; the mutation is
recorded in the test's comment.

| Test | Assertion | Mutation that must fail it |
|---|---|---|
| Euclidean counts | Infinite uniform dust-free single-bin medium: fitted slope of log N(>F) vs log F is −1.5 ± 0.05 over 3 decades | drop the s² in the point budget |
| Energy conservation | Milky Way scene: total sky flux (points + diffuse) at budgets 500, 4000, 20 000 agrees within 3σ Poisson noise of the points | classify a step as both resolved and diffuse, or as neither |
| Distant galaxy | Dust-free galaxy at 1 Mpc: total sky flux = Σ L · (10/d)² · 10^(−1.932) within 1 % | restore an unbounded density term |
| Dust slab | Uniform slab of optical depth τ in front of a point source: per-channel flux ratio = exp(−k_c τ) within 0.5 % | apply k_G to all channels |
| Inside a disk | Milky Way, observer at R = 8.2 kpc, z = 20 pc, 20 000 stars: mean diffuse radiance within 5° of the galactic plane ≥ 10× within 5° of the poles; point density per steradian within 10° of the plane vs above 60° latitude in [1.5, 4] (wider bands than the diffuse check so the counts are not noise-limited) | replace the sech² with 1 |
| Determinism | Quick cook with 1 thread and with 8 threads: byte-identical files | seed draws by thread index |
| Format round trip | Write then `Parse`: every field and payload byte equal | — |
| Format rejection | `Parse` rejects every truncation length, each offset/size pushed past EOF, overlapping regions, bad magic/version/encoding, mip count inconsistent with face size | delete any one validation check |
| RGB9E5 | Encode/decode matches the Vulkan reference algorithm over a logarithmic sweep of 10⁶ values plus the edge cases (0, max, negative, NaN → 0) | wrong exponent bias |
| Preview flux | One star of known flux at latitudes 0°, 60°, 85°: preview pixels summed × solid angle = flux within 0.5 % | omit the per-pixel solid angle |

## 6. Repository layout

```
Starscape/
  CMakeLists.txt                 # C++23; targets: starscape, starscape-tests
  include/kre_starscape.h        # public format + reader
  src/
    cli.cpp                      # argument parsing, top-level error printing
    scene.h/.cpp                 # TOML -> Scene, validation
    presets.h/.cpp               # galaxy and population presets
    galaxy.h/.cpp                # components, arm modulation, dust field, bounding spheres
    population.h/.cpp            # bins, chromaticity
    common.h                     # Error, Expected, TRY macros
    noise.h/.cpp                 # gradient noise, fbm
    rng.h                        # Philox 4x32-10
    march.h/.cpp                 # intervals, per-ray march, flux histogram, star re-march
    parallel.h                   # deterministic ParallelFor
    cook.h/.cpp                  # calibration, texel march, Poisson star drawing
    writer.h/.cpp                # mips, file assembly, atomic rename
    preview.h/.cpp               # equirect render, splats, RGBE write
  scenes/milky-way-from-earth.toml
  third_party/                   # vendored single headers, versions in third_party/VERSIONS
    glm/  toml.hpp  doctest.h    # the Radiance .hdr writer is ~30 lines in preview.cpp
  tests/
  docs/spec.md
  README.md                      # usage, units, citations
  LICENSE                        # MIT
```

Dependencies are vendored files rather than submodules so the repository clones
in one step and consumers need not recurse into it.

## 7. Integration with Kreatures (this spec)

- `Tooling/Starscape` is added to the Kreatures superproject as a submodule.
- `docs/engine-roadmap.md` gains a row for the sky replacement, pointing here and
  naming the engine-side spec as the remainder.

## References

- O. Elfwering, *Procedurally Generating an Artificial Galaxy*, Lund University,
  Department of Statistics, 2016 — the starting point of the earlier experiment.
- L. Hernquist, "An analytical model for spherical galaxies and bulges", ApJ 356
  (1990) — bulge profile.
- P. Kroupa, "On the variation of the initial mass function", MNRAS 322 (2001) —
  IMF for the population tables.
- J. A. Cardelli, G. C. Clayton, J. S. Mathis, "The relationship between infrared,
  optical, and ultraviolet extinction", ApJ 345 (1989) — channel extinction ratios.
- P. C. van der Kruit, L. Searle, "Surface photometry of edge-on spiral galaxies",
  A&A 95 (1981) — exponential/sech² disk.
- J. K. Salmon et al., "Parallel random numbers: as easy as 1, 2, 3", SC11 (2011)
  — Philox.
