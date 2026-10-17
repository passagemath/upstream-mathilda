#!/usr/bin/env python3
"""Generate the piecewise-polynomial tables behind the real machine Airy kernels.

Writes src/special_functions/sf_machine_airy_tab.c.  Run with the benchmark venv:

    benchmarks/.venv/bin/python -I tools/sf_accuracy/gen_airy_tables.py [--check]

Method (Chebyshev interpolation, e.g. Trefethen, "Approximation Theory and
Approximation Practice", ch. 3-4; the region split follows Cephes `airy`):

  * The real line is cut into three fitted regions plus two asymptotic tails:
        [-10, 0]   Ai, Ai', Bi, Bi' themselves (oscillatory, entire)
        [ 0,  2]   Ai, Ai', Bi, Bi' themselves (entire, no cancellation)
        [ 2, 10]   the SCALED functions Ai e^zeta, Ai' e^zeta, Bi e^-zeta,
                   Bi' e^-zeta  (zeta = 2/3 x^(3/2)), which vary slowly, so the
                   exponential is carried exactly by exp() and never fitted
        |x| > 10   the Poincare asymptotic expansions (in the C file).
  * Each region is split into equal-width intervals.  On each interval the
    function is interpolated at Chebyshev points of the first kind at 50 digits,
    the Chebyshev series is truncated where its tail falls below 2^-60 of the
    local magnitude, and converted EXACTLY (in mpmath) to monomial coefficients
    in the local variable t = (x - centre)/halfwidth in [-1, 1], then rounded to
    double.  The C side evaluates by Horner with fma; the --check pass simulates
    exactly that evaluation in Python (math.fma) and reports the worst error in
    ulps of the local magnitude over a dense grid.
"""
from __future__ import annotations

import argparse
import math
import os

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
WORKTREE = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(WORKTREE, "src", "special_functions", "sf_machine_airy_tab.c")
OUT_H = os.path.join(WORKTREE, "src", "special_functions", "sf_machine_airy_tab.h")
mp.mp.dps = 50


def zeta(x):
    return mp.mpf(2) / 3 * x * mp.sqrt(x)


FUNCS = {
    "ai":  lambda x: mp.airyai(x),
    "aip": lambda x: mp.airyai(x, 1),
    "bi":  lambda x: mp.airybi(x),
    "bip": lambda x: mp.airybi(x, 1),
}
SCALED = {
    "ai":  lambda x: mp.airyai(x) * mp.exp(zeta(x)),
    "aip": lambda x: mp.airyai(x, 1) * mp.exp(zeta(x)),
    "bi":  lambda x: mp.airybi(x) * mp.exp(-zeta(x)),
    "bip": lambda x: mp.airybi(x, 1) * mp.exp(-zeta(x)),
}

# (name, lo, hi, nintervals, scaled?)
REGIONS = [
    ("neg", -10.0, 0.0, 40, False),
    ("mid", 0.0, 2.0, 4, False),
    ("pos", 2.0, 10.0, 8, True),
]
TOL = mp.mpf(2) ** -60
# One degree per region (the max any function needs there), so the C side
# indexes every table of a region with the same stride.
REGION_DEG = {"neg": 12, "mid": 13, "pos": 16}


def cheb_coeffs(f, a, b, n):
    """Chebyshev coefficients c_0..c_n of the degree-n interpolant on [a, b]."""
    a, b = mp.mpf(a), mp.mpf(b)
    c, h = (a + b) / 2, (b - a) / 2
    nodes = [mp.cos(mp.pi * (j + mp.mpf(1) / 2) / (n + 1)) for j in range(n + 1)]
    fv = [f(c + h * t) for t in nodes]
    co = []
    for k in range(n + 1):
        s = mp.fsum(fv[j] * mp.cos(mp.pi * k * (j + mp.mpf(1) / 2) / (n + 1)) for j in range(n + 1))
        co.append(s * 2 / (n + 1))
    co[0] /= 2
    return co


def cheb_to_mono(co):
    """Exact conversion of sum c_k T_k(t) to monomial coefficients."""
    n = len(co) - 1
    T = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]]
    for k in range(2, n + 1):
        prev, prev2 = T[k - 1], T[k - 2]
        nxt = [mp.mpf(0)] * (k + 1)
        for i, v in enumerate(prev):
            nxt[i + 1] += 2 * v
        for i, v in enumerate(prev2):
            nxt[i] -= v
        T.append(nxt)
    mono = [mp.mpf(0)] * (n + 1)
    for k, ck in enumerate(co):
        for i, v in enumerate(T[k]):
            mono[i] += ck * v
    return mono


def local_mag(name, f, x, scaled):
    if scaled:
        return abs(f(x))
    v = f(x)
    d = mp.diff(f, x)
    return max(abs(v), abs(d) / mp.sqrt(max(1, abs(x))))


def fit_interval(name, f, a, b, scaled):
    # magnitude floor over the interval: smallest local magnitude on a grid
    grid = [a + (b - a) * i / 16 for i in range(17)]
    mag = min(local_mag(name, f, mp.mpf(x), scaled) for x in grid)
    n = 40
    co = cheb_coeffs(f, a, b, n)
    # truncate the tail below TOL * mag
    keep = n
    tail = mp.mpf(0)
    for k in range(n, -1, -1):
        tail += abs(co[k])
        if tail > TOL * mag:
            keep = k
            break
    deg = max(keep, 1)
    co = cheb_coeffs(f, a, b, deg)          # re-interpolate at the final degree
    mono = cheb_to_mono(co)
    return [float(m) for m in mono], deg, mag


def horner(c, t):
    s = c[-1]
    for v in reversed(c[:-1]):
        s = math.fma(s, t, v)
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    tables = {}
    for rname, lo, hi, nint, scaled in REGIONS:
        w = (hi - lo) / nint
        for fname in FUNCS:
            f = (SCALED if scaled else FUNCS)[fname]
            ivs = []
            for i in range(nint):
                ia, ib = lo + i * w, lo + (i + 1) * w
                mono, deg, mag = fit_interval(fname, f, ia, ib, scaled)
                ivs.append((ia, ib, mono, deg))
            maxdeg = REGION_DEG[rname]
            assert max(iv[3] for iv in ivs) <= maxdeg, (rname, fname, [iv[3] for iv in ivs])
            # pad every interval to the region's max degree (uniform stride)
            padded = []
            for ia, ib, mono, deg in ivs:
                f2 = (SCALED if scaled else FUNCS)[fname]
                co = cheb_coeffs(f2, ia, ib, maxdeg)
                padded.append([float(m) for m in cheb_to_mono(co)])
            tables[(rname, fname)] = (lo, hi, nint, maxdeg, padded)
            print(f"{rname:4s} {fname:4s} intervals={nint} degree={maxdeg} "
                  f"(per-interval {[iv[3] for iv in ivs]})", flush=True)
            if a.check:
                worst = 0.0
                for i, coeffs in enumerate(padded):
                    ia = lo + i * w
                    for j in range(41):
                        x = ia + w * j / 40
                        c = ia + w / 2
                        t = (x - c) / (w / 2)
                        got = horner(coeffs, t)
                        ref = f(mp.mpf(x))
                        mag = local_mag(fname, f, mp.mpf(x), scaled)
                        e = float(abs(got - ref) / mag) / 2.0 ** -52
                        worst = max(worst, e)
                print(f"        check: worst {worst:.3f} ulp(local magnitude)", flush=True)

    with open(OUT, "w") as fh:
        fh.write("/* GENERATED by tools/sf_accuracy/gen_airy_tables.py -- do not edit.\n"
                 " *\n"
                 " * Piecewise monomial coefficients (local variable t in [-1, 1]) for the real\n"
                 " * machine Airy kernels in sf_machine_bessel.c.  See the generator for the\n"
                 " * method (Chebyshev interpolation at 50 digits, exact conversion to the\n"
                 " * monomial basis, tail truncated below 2^-60 of the local magnitude). */\n"
                 "#include \"sf_machine_airy_tab.h\"\n\n")
        for (rname, fname), (lo, hi, nint, deg, padded) in tables.items():
            sym = f"sf_airy_tab_{rname}_{fname}"
            fh.write(f"/* {fname} on [{lo:g}, {hi:g}]: {nint} intervals, degree {deg} */\n")
            fh.write(f"const double {sym}[{nint * (deg + 1)}] = {{\n")
            for coeffs in padded:
                fh.write("    " + ", ".join(f"{v:.17e}" for v in coeffs) + ",\n")
            fh.write("};\n\n")
    with open(OUT_H, "w") as fh:
        fh.write("/* GENERATED by tools/sf_accuracy/gen_airy_tables.py -- do not edit.\n"
                 " * Geometry of the piecewise Airy tables in sf_machine_airy_tab.c. */\n"
                 "#ifndef MATHILDA_SF_MACHINE_AIRY_TAB_H\n#define MATHILDA_SF_MACHINE_AIRY_TAB_H\n\n")
        for rname, lo, hi, nint, scaled in REGIONS:
            R = rname.upper()
            fh.write(f"#define SF_AIRY_{R}_LO  ({lo!r})\n#define SF_AIRY_{R}_HI  ({hi!r})\n"
                     f"#define SF_AIRY_{R}_N   {nint}\n#define SF_AIRY_{R}_DEG {REGION_DEG[rname]}\n")
            for fname in FUNCS:
                n = nint * (REGION_DEG[rname] + 1)
                fh.write(f"extern const double sf_airy_tab_{rname}_{fname}[{n}];\n")
            fh.write("\n")
        fh.write("#endif\n")
    print("wrote", OUT, OUT_H)


if __name__ == "__main__":
    main()
