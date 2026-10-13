# Starscape

Cooks a night sky for any vantage point in a scene of procedural galaxies: a diffuse HDR
cubemap of the light you cannot resolve into stars, plus the few thousand stars you can.
Output is a `.kre_starscape` file (format: [`docs/spec.md`](docs/spec.md); the reader is
[`include/kre_starscape.h`](include/kre_starscape.h), header-only, standard library only) and an
optional equirectangular Radiance `.hdr` preview. The Kreatures engine consumes `.kre_starscape`.

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
starscape cook scene.toml -o draft.kre_starscape --draft --preview draft.hdr   # seconds: quick, face <= 64, preview <= 1024 wide
starscape preview sky.kre_starscape -o bright.hdr --exposure 4
starscape inspect sky.kre_starscape
```

Options: `--threads N` (cook, preview; output is byte-identical for any thread count),
`--quick` and `--draft` (cook), `--width N` (cook and preview: preview width) and `--exposure EV` (preview). Exit codes: 0 success,
1 invalid input or usage, 2 I/O failure.

## Scene format

A scene is TOML: `[output]` (face size, `target_stars` or `limiting_magnitude`, quality,
preview), `[observer]` (position in parsecs, spin axis, zero meridian), an optional `[style]`
(scene-wide stylisation such as `star_saturation`) and one or more
`[[galaxy]]` tables starting from a preset (`milky-way`, `small-spiral`, `elliptical`,
`dwarf`) with any field overridden. Unknown keys are errors. See
[`scenes/milky-way-from-earth.toml`](scenes/milky-way-from-earth.toml) for the realistic sky,
[`scenes/stylized-milky-way.toml`](scenes/stylized-milky-way.toml) for a fantastical one, and
[`docs/spec.md`](docs/spec.md) for the model, formulas and units. The spec describes the main keys; every
accepted key, its range and its preset default is defined in [`src/scene.cpp`](src/scene.cpp), and
[`scenes/`](scenes/) shows them in use.

## Units

Parsecs; visual luminosity in solar units; star flux in units of a magnitude-0 star
(flux = 10^(-0.4 m)); diffuse radiance in magnitude-0 flux per steradian; linear Rec.709.
Stored values are multiplied by the header's `valueScale`; divide by it to read them.

## How the light is split

One apparent-flux threshold decides, at every step of every ray, whether a luminosity bin's
stars are bright enough to be points. Points and diffuse light add up to the whole sky, so
asking for fewer stars moves light into the texture instead of losing it.

## References

- O. Elfwering, *Procedurally Generating an Artificial Galaxy*, Lund University, 2016 - the starting point.
- L. Hernquist, ApJ 356 (1990) - bulge profile.
- P. Kroupa, MNRAS 322 (2001) - initial mass function.
- J. A. Cardelli, G. C. Clayton, J. S. Mathis, ApJ 345 (1989) - extinction curve.
- P. C. van der Kruit, L. Searle, A&A 95 (1981) - disk profile.
- C. Wyman, P.-P. Sloan, P. Shirley, JCGT 2(2) (2013) - analytic CIE colour matching functions.
- J. K. Salmon et al., SC11 (2011) - Philox.

## License

MIT.
