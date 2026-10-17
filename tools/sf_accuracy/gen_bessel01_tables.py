#!/usr/bin/env python3
"""Generate the tables behind the real machine I_0, I_1, K_0, K_1 kernels.

Writes src/special_functions/sf_machine_bessel_tab.{c,h}.  Run with the
benchmark venv:

    benchmarks/.venv/bin/python -I tools/sf_accuracy/gen_bessel01_tables.py [--check]

Regions (x > 0; the Cephes i0/i1/k0/k1 split, with the fits regenerated here):

  I_0, I_1   x <= 2       ascending series, t = x^2/4:  I_0 = sum t^k/k!^2,
                          I_1 = (x/2) sum t^k/(k!(k+1)!)  -- positive terms, so
                          Horner is accurate to ~1 ulp (DLMF 10.25.2)
             2 < x < 32   I e^-x, piecewise Chebyshev (entire, slowly varying)
             x >= 32      Hankel: I = e^x/sqrt(x) sum (-1)^k c a_k(nu)/x^k,
                          c = 1/sqrt(2 pi)  (DLMF 10.40.1)
  K_0, K_1   x <= 1/2     the log series (DLMF 10.31.2 / 10.31.1):
                          K_0 = -ln(x/2) I_0 + sum psi(k+1) t^k/k!^2
                          K_1 = 1/x + ln(x/2) I_1
                                - (x/4) sum (psi(k+1)+psi(k+2)) t^k/(k!(k+1)!)
             1/2 < x < 32 K e^x, piecewise Chebyshev
             x >= 32      Hankel: K = e^-x/sqrt(x) sum c a_k(nu)/x^k,
                          c = sqrt(pi/2)  (DLMF 10.40.2)

Piecewise fits: Chebyshev interpolation at 50 digits on equal power-of-two
width intervals, tail truncated below 2^-60, converted exactly to the monomial
basis in t = (x - centre)/halfwidth.  --check simulates the C evaluation in
Python (math.fma) over each whole function and reports ulps vs mpmath.
"""
from __future__ import annotations

import argparse
import math
import os

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
WORKTREE = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(WORKTREE, "src", "special_functions", "sf_machine_bessel_tab.c")
OUT_H = os.path.join(WORKTREE, "src", "special_functions", "sf_machine_bessel_tab.h")
mp.mp.dps = 50
TOL = mp.mpf(2) ** -60

SCALED = {
    "i0": lambda x: mp.besseli(0, x) * mp.exp(-x),
    "i1": lambda x: mp.besseli(1, x) * mp.exp(-x),
    "k0": lambda x: mp.besselk(0, x) * mp.exp(x),
    "k1": lambda x: mp.besselk(1, x) * mp.exp(x),
}
# (region, function family, lo, hi, n intervals)
REGIONS = [
    ("ia", ("i0", "i1"), 2.0, 8.0, 6),
    ("ib", ("i0", "i1"), 8.0, 32.0, 6),
    ("ka", ("k0", "k1"), 0.5, 1.0, 2),
    ("kb", ("k0", "k1"), 1.0, 2.0, 2),
    ("kc", ("k0", "k1"), 2.0, 8.0, 6),
    ("kd", ("k0", "k1"), 8.0, 32.0, 6),
]
NASYM = 22       # Hankel terms kept for x >= 32
NSER = 15        # ascending-series terms for x <= 2


def cheb_coeffs(f, a, b, n):
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


def needed_degree(f, a, b):
    mag = min(abs(f(mp.mpf(a) + (mp.mpf(b) - a) * i / 16)) for i in range(17))
    co = cheb_coeffs(f, a, b, 40)
    tail = mp.mpf(0)
    for k in range(40, -1, -1):
        tail += abs(co[k])
        if tail > TOL * mag:
            return max(k, 1)
    return 1


def horner(c, t):
    s = c[-1]
    for v in reversed(c[:-1]):
        s = math.fma(s, t, v)
    return s


def hankel_a(nu, n):
    mu = 4 * mp.mpf(nu) ** 2
    a = [mp.mpf(1)]
    for k in range(1, n):
        a.append(a[-1] * (mu - (2 * k - 1) ** 2) / (8 * k))
    return a


def series_tables():
    fac = [mp.factorial(k) for k in range(NSER + 2)]
    psi = [mp.digamma(k + 1) for k in range(NSER + 2)]
    i0 = [1 / fac[k] ** 2 for k in range(NSER)]
    i1 = [1 / (fac[k] * fac[k + 1]) for k in range(NSER)]
    k0 = [psi[k] / fac[k] ** 2 for k in range(NSER)]
    k1 = [(psi[k] + psi[k + 1]) / (fac[k] * fac[k + 1]) for k in range(NSER)]
    return {"ser_i0": i0, "ser_i1": i1, "ser_k0": k0, "ser_k1": k1}


def asym_tables():
    cI = 1 / mp.sqrt(2 * mp.pi)
    cK = mp.sqrt(mp.pi / 2)
    out = {}
    for nu in (0, 1):
        a = hankel_a(nu, NASYM)
        out[f"asy_i{nu}"] = [cI * (-1) ** k * a[k] for k in range(NASYM)]
        out[f"asy_k{nu}"] = [cK * a[k] for k in range(NASYM)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    fits = {}
    region_deg = {}
    for rname, fams, lo, hi, nint in REGIONS:
        w = (hi - lo) / nint
        deg = 0
        for fn in fams:
            for i in range(nint):
                deg = max(deg, needed_degree(SCALED[fn], lo + i * w, lo + (i + 1) * w))
        region_deg[rname] = deg
        for fn in fams:
            tabs = []
            for i in range(nint):
                co = cheb_coeffs(SCALED[fn], lo + i * w, lo + (i + 1) * w, deg)
                tabs.append([float(m) for m in cheb_to_mono(co)])
            fits[(rname, fn)] = tabs
        print(f"{rname}: [{lo:g},{hi:g}] n={nint} degree={deg}", flush=True)
    ser = series_tables()
    asy = asym_tables()

    with open(OUT, "w") as fh:
        fh.write("/* GENERATED by tools/sf_accuracy/gen_bessel01_tables.py -- do not edit.\n"
                 " * Tables for the real machine I_0, I_1, K_0, K_1 kernels in\n"
                 " * sf_machine_bessel.c; see the generator for the method. */\n"
                 "#include \"sf_machine_bessel_tab.h\"\n\n")
        for (rname, fn), tabs in fits.items():
            n = len(tabs) * len(tabs[0])
            fh.write(f"const double sf_bes_tab_{rname}_{fn}[{n}] = {{\n")
            for t in tabs:
                fh.write("    " + ", ".join(f"{v:.17e}" for v in t) + ",\n")
            fh.write("};\n\n")
        for name, vals in list(ser.items()) + list(asy.items()):
            fh.write(f"const double sf_bes_{name}[{len(vals)}] = {{\n    ")
            fh.write(",\n    ".join(f"{float(v):.17e}" for v in vals))
            fh.write("\n};\n\n")
    with open(OUT_H, "w") as fh:
        fh.write("/* GENERATED by tools/sf_accuracy/gen_bessel01_tables.py -- do not edit. */\n"
                 "#ifndef MATHILDA_SF_MACHINE_BESSEL_TAB_H\n#define MATHILDA_SF_MACHINE_BESSEL_TAB_H\n\n")
        for rname, fams, lo, hi, nint in REGIONS:
            R = rname.upper()
            fh.write(f"#define SF_BES_{R}_LO  ({lo!r})\n#define SF_BES_{R}_HI  ({hi!r})\n"
                     f"#define SF_BES_{R}_N   {nint}\n#define SF_BES_{R}_DEG {region_deg[rname]}\n")
            for fn in fams:
                fh.write(f"extern const double sf_bes_tab_{rname}_{fn}[{nint * (region_deg[rname] + 1)}];\n")
            fh.write("\n")
        fh.write(f"#define SF_BES_NSER  {NSER}\n#define SF_BES_NASYM {NASYM}\n")
        for name, vals in list(ser.items()) + list(asy.items()):
            fh.write(f"extern const double sf_bes_{name}[{len(vals)}];\n")
        fh.write("\n#endif\n")
    print("wrote", OUT, OUT_H)

    if a.check:
        check(fits, region_deg, ser, asy)


# ---- a Python model of the C evaluation, for --check -------------------------
def model(fn, x, fits, region_deg, ser, asy):
    S = {k: [float(v) for v in vals] for k, vals in ser.items()}
    A = {k: [float(v) for v in vals] for k, vals in asy.items()}
    isI = fn[0] == "i"
    nu = int(fn[1])
    if x >= 32.0:
        s = horner(A[f"asy_{fn}"], 1.0 / x)
        return (math.exp(x) if isI else math.exp(-x)) * (s / math.sqrt(x))
    if isI and x <= 2.0:
        t = x * x * 0.25
        s = horner(S[f"ser_{fn}"], t)
        return s if nu == 0 else 0.5 * x * s
    if not isI and x <= 0.5:
        t = x * x * 0.25
        l = math.log(0.5 * x)
        i = horner(S[f"ser_i{nu}"], t)
        if nu == 0:
            return math.fma(-l, i, horner(S["ser_k0"], t))
        i *= 0.5 * x
        return 1.0 / x + math.fma(l, i, -0.25 * x * horner(S["ser_k1"], t))
    for rname, fams, lo, hi, nint in REGIONS:
        if fn in fams and lo <= x < hi:
            w = (hi - lo) / nint
            i = min(int((x - lo) / w), nint - 1)
            c = lo + (i + 0.5) * w
            t = (x - c) * (2.0 / w)
            v = horner(fits[(rname, fn)][i], t)
            return v * (math.exp(x) if isI else math.exp(-x))
    raise ValueError(x)


def check(fits, region_deg, ser, asy):
    import random
    rng = random.Random(7)
    for fn in ("i0", "i1", "k0", "k1"):
        ref = {"i0": lambda x: mp.besseli(0, x), "i1": lambda x: mp.besseli(1, x),
               "k0": lambda x: mp.besselk(0, x), "k1": lambda x: mp.besselk(1, x)}[fn]
        worst, wx, errs = 0.0, None, []
        for _ in range(3000):
            x = math.exp(rng.uniform(math.log(1e-3), math.log(700.0)))
            v = model(fn, x, fits, region_deg, ser, asy)
            r = ref(mp.mpf(x))
            if not math.isfinite(v) or abs(r) > 1e300 or abs(r) < 1e-300:
                continue
            u = 2.0 ** (math.frexp(float(abs(r)))[1] - 53)
            e = float(abs(v - r)) / u
            errs.append(e)
            if e > worst:
                worst, wx = e, x
        errs.sort()
        print(f"  check {fn}: median {errs[len(errs)//2]:.2f} max {worst:.2f} ulp (at x={wx!r})")


if __name__ == "__main__":
    main()
