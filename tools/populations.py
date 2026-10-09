#!/usr/bin/env python3
"""Regenerates the stellar population tables in src/presets.cpp.

Each population is a number fraction and a light-weighted effective temperature
per luminosity bin (24 bins over log10 L_V = -4..6, matching src/population.h).

  python3 tools/populations.py            # print the C++ initialisers
  python3 tools/populations.py --report   # integrated colour temperature and mean L_V

Model:
- Kroupa (2001) IMF: xi(m) ~ m^-1.3 on 0.08-0.5 Msun, m^-2.3 above, to 150 Msun.
- Main sequence: L_bol = 0.23 m^2.3 (m < 0.43), m^4 (< 2), 1.4 m^3.5 (< 55), 32000 m;
  R = m^0.8 (m < 1) else m^0.57; T = 5772 (L / R^2)^0.25.
- L_V = L_bol x (blackbody flux under the CIE y-bar curve / sigma T^4), normalised so the
  Sun has L_V = 1.
- A turnoff mass removes the main sequence above it. Stars in the first `giant_span`
  (fractional) of mass above the turnoff are evolved: a share `clump_share` sits in the
  red clump (L_V 50, 4800 K); the rest spread over the giant branch from L_V 3 to 1000 with
  dN/dlog L ~ L^-0.75 (stars climb the branch ever faster), cooling from 5100 K to 3900 K
  in log L. Heavier stars are remnants and, like every star below the lowest bin edge,
  are counted in bin 0 (they add number but no light).
- Bins with no stars repeat the previous bin's temperature.
"""
import sys
import numpy as np

BINS, LOG_MIN, SPAN = 24, -4.0, 10.0
EDGES = LOG_MIN + SPAN * np.arange(BINS + 1) / BINS

h, c, k = 6.626e-34, 2.998e8, 1.381e-23
LAM = np.linspace(360e-9, 830e-9, 941)
NM = LAM * 1e9


def _g(x, mu, s1, s2):
    return np.exp(-0.5 * ((x - mu) / np.where(x < mu, s1, s2)) ** 2)


# CIE 1931 colour matching functions, Wyman, Sloan & Shirley (2013) multi-lobe fit.
XBAR = 1.056 * _g(NM, 599.8, 37.9, 31.0) + 0.362 * _g(NM, 442.0, 16.0, 26.7) - 0.065 * _g(NM, 501.1, 20.4, 26.2)
YBAR = 0.821 * _g(NM, 568.8, 46.9, 40.5) + 0.286 * _g(NM, 530.9, 16.3, 31.1)
ZBAR = 1.217 * _g(NM, 437.0, 11.8, 36.0) + 0.681 * _g(NM, 459.0, 26.0, 13.8)
XYZ_TO_709 = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])


def planck(T):
    return 1.0 / LAM ** 5 / np.expm1(h * c / (LAM * k * T))


def visual_fraction(T):
    """Flux under y-bar per unit bolometric flux, unnormalised."""
    dlam = LAM[1] - LAM[0]
    visual = 2 * h * c ** 2 * (planck(T) * YBAR).sum() * dlam
    return visual / (5.670e-8 * T ** 4 / np.pi)


SUN_V = visual_fraction(5772.0)


def lv(lbol, T):
    return lbol * visual_fraction(T) / SUN_V


def rgb(T):
    B = planck(T)
    v = XYZ_TO_709 @ np.array([(B * XBAR).sum(), (B * YBAR).sum(), (B * ZBAR).sum()])
    return v / v[1]


def colour_temperature(r_over_g):
    Ts = np.linspace(2500, 12000, 400)
    rg = np.array([rgb(T)[0] for T in Ts])
    return np.interp(-r_over_g, -rg, Ts)


def imf(m):
    return np.where(m < 0.5, (m / 0.5) ** -1.3, (m / 0.5) ** -2.3)


def main_sequence(m):
    lbol = np.where(m < 0.43, 0.23 * m ** 2.3, np.where(m < 2, m ** 4, np.where(m < 55, 1.4 * m ** 3.5, 32000 * m)))
    R = np.where(m < 1, m ** 0.8, m ** 0.57)
    T = 5772.0 * (lbol / R ** 2) ** 0.25
    return lbol, T


def population(turnoff=None, giant_span=0.08, clump_share=0.35):
    m = np.geomspace(0.08, 150.0, 200000)
    dm = np.gradient(m)
    w = imf(m) * dm
    count = np.zeros(BINS)
    light_T = np.zeros(BINS)
    light = np.zeros(BINS)

    def add(n, L, T):
        b = np.clip(np.floor((np.log10(L) - LOG_MIN) / SPAN * BINS).astype(int), -1, BINS)
        visible = (b >= 0) & (b < BINS)
        np.add.at(count, np.where(visible, b, 0), n)   # below the range: bin 0
        np.add.at(light, b[visible], n[visible] * L[visible])
        np.add.at(light_T, b[visible], n[visible] * L[visible] * T[visible])

    lbol, T = main_sequence(m)
    L = np.array([lv(lb, t) for lb, t in zip(lbol, T)])
    ms = np.ones_like(m, bool) if turnoff is None else m <= turnoff
    add(w[ms], L[ms], T[ms])
    if turnoff is not None:
        evolved = (m > turnoff) & (m <= turnoff * (1 + giant_span))
        giants = w[evolved].sum()
        add(np.array([giants * clump_share]), np.array([50.0]), np.array([4800.0]))
        g = np.geomspace(3.0, 1000.0, 400)
        p = g ** -0.75
        p *= giants * (1 - clump_share) / p.sum()
        Tg = np.interp(np.log10(g), [np.log10(3.0), 3.0], [5100.0, 3900.0])
        add(p, g, Tg)
        remnants = w[m > turnoff * (1 + giant_span)].sum()
        count[0] += remnants
    fraction = count / count.sum()
    temperature = np.zeros(BINS)
    last = 2680.0
    for i in range(BINS):
        if light[i] > 0:
            last = light_T[i] / light[i]
        temperature[i] = last
    return fraction, temperature


def report(name, fraction, temperature):
    Lc = 10 ** (EDGES[:-1] + SPAN / BINS / 2)
    light = fraction * Lc
    col = sum(f * rgb(T) for f, T in zip(light, temperature))
    print('%-11s mean L_V %.3f  colour %.0f K  giants carry %.0f%% of the light' % (
        name, light.sum(), colour_temperature(col[0] / col[1]), 100 * light[12:].sum() / light.sum()))


POPULATIONS = {
    'young-disk': dict(turnoff=None),
    # The thin disk formed stars for ~10 Gyr, so A/F stars up to ~1.6 Msun (~2.5 Gyr) are still on
    # the main sequence; a single old turnoff makes the disk's light far too red.
    'old-disk': dict(turnoff=1.6, giant_span=0.08),
    'bulge': dict(turnoff=0.85, giant_span=0.08),
}

if __name__ == '__main__':
    for name, args in POPULATIONS.items():
        f, T = population(**args)
        if '--report' in sys.argv:
            report(name, f, T)
        else:
            print('// %s' % name)
            print('{' + ', '.join('%.6e' % x for x in f) + '},')
            print('{' + ', '.join('%.0f' % x for x in T) + '}')
