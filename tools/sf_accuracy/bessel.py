#!/usr/bin/env python3
"""Accuracy gate for the "bessel" family of machine special-function kernels.

Heads: BesselJ[2, x], BesselY[1, x], BesselI[0, x], BesselK[0, x] (the benchmark
orders), the same four at a random REAL order nu (the *nu cases), and AiryAi,
AiryBi, AiryAiPrime, AiryBiPrime.

Drives `<worktree>/Mathilda -file` on a generated .m that prints each value with
ToString[InputForm[SetPrecision[r, 17]]] (17 digits round-trip a double), and
compares against mpmath at 30 digits.

Two error measures per point, both in ulps:
  rel     |f~ - f| / ulp(|f|)                 the plain relative error
  scaled  |f~ - f| / ulp(M),  M = max(|f|, h |f'|)
          the error against the LOCAL MAGNITUDE of the function.  h is the local
          length scale (1 for Bessel, 1/sqrt|x| on the oscillatory Airy side), so
          M is the envelope of an oscillating function: away from a zero M = |f|
          and scaled == rel; at a zero the relative error is meaningless and the
          brief asks for "absolute error scaled by the local magnitude".

Modes:
  scalar   one statement per point:  F[x`]          (the scalar builtin path)
  packed   F[NDArray[{x1`, ...}]]                    (the ND kernel path)

Usage:
  benchmarks/.venv/bin/python -I tools/sf_accuracy/bessel.py \
      [--heads BesselJ,AiryAi] [--regimes MR,C,MRW,CW] [--mode scalar|packed]
      [--n 2000] [--bin <worktree>/Mathilda] [--worst 3]

Regimes: MR / C are the benchmark domains (benchmarks/97-special-functions/
registry.py); MRW / CW a wider band, which may be looser if documented.
"""
from __future__ import annotations

import argparse
import math
import os
import random
import re
import statistics
import subprocess
import sys
import tempfile
import zlib

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
WORKTREE = os.path.dirname(os.path.dirname(HERE))
mp.mp.dps = 30


def wl_real(x: float) -> str:
    s = repr(float(x))
    if "e" in s or "E" in s:
        m, _, e = s.lower().partition("e")
        return f"{m}`*^{int(e)}"
    return s + "`"


def wl_complex(z: complex) -> str:
    return f"Complex[{wl_real(z.real)}, {wl_real(z.imag)}]"


def h_bessel(x):
    return 1.0


def h_airy(x):
    r = abs(complex(x))
    return 1.0 / math.sqrt(max(1.0, r))


AI = lambda z: mp.airyai(z)
BI = lambda z: mp.airybi(z)
AIP = lambda z: mp.airyai(z, 1)
BIP = lambda z: mp.airybi(z, 1)

# name -> dict(template, fn(nu, z), order (fixed float | (lo, hi) | None),
#              MR, C, MRW, CW domains, h).
# Real domains (lo, hi); complex (re_lo, re_hi, im_lo, im_hi).
HEADS = {
    "BesselJ": dict(t="BesselJ[{1}, {0}]", f=mp.besselj, nu=2, h=h_bessel,
                    MR=(0.5, 15.0), C=(0.5, 8.0, -4.0, 4.0),
                    MRW=(1e-3, 80.0), CW=(-30.0, 30.0, -20.0, 20.0)),
    "BesselY": dict(t="BesselY[{1}, {0}]", f=mp.bessely, nu=1, h=h_bessel,
                    MR=(0.5, 15.0), C=(0.5, 8.0, 0.1, 4.0),
                    MRW=(1e-3, 80.0), CW=(-30.0, 30.0, -20.0, 20.0)),
    "BesselI": dict(t="BesselI[{1}, {0}]", f=mp.besseli, nu=0, h=h_bessel,
                    MR=(0.1, 8.0), C=(0.5, 4.0, -3.0, 3.0),
                    MRW=(-60.0, 600.0), CW=(-30.0, 30.0, -30.0, 30.0)),
    "BesselK": dict(t="BesselK[{1}, {0}]", f=mp.besselk, nu=0, h=h_bessel,
                    MR=(0.1, 8.0), C=(0.5, 4.0, 0.1, 3.0),
                    MRW=(1e-3, 600.0), CW=(-30.0, 30.0, -30.0, 30.0)),
    "BesselJn": dict(t="BesselJ[{1}, {0}]", f=mp.besselj, nu=(-12, 12), h=h_bessel,
                     MR=(0.5, 15.0), C=(0.5, 8.0, -4.0, 4.0),
                     MRW=(-80.0, 80.0), CW=(-30.0, 30.0, -20.0, 20.0), intorder=True),
    "BesselYn": dict(t="BesselY[{1}, {0}]", f=mp.bessely, nu=(-12, 12), h=h_bessel,
                     MR=(0.5, 15.0), C=(0.5, 8.0, 0.1, 4.0),
                     MRW=(1e-2, 80.0), CW=(-30.0, 30.0, -20.0, 20.0), intorder=True),
    "BesselJnu": dict(t="BesselJ[{1}, {0}]", f=mp.besselj, nu=(-4.7, 9.3), h=h_bessel,
                      MR=(0.5, 15.0), C=(0.5, 8.0, -4.0, 4.0),
                      MRW=(1e-2, 80.0), CW=(-30.0, 30.0, -20.0, 20.0)),
    "BesselYnu": dict(t="BesselY[{1}, {0}]", f=mp.bessely, nu=(-4.7, 9.3), h=h_bessel,
                      MR=(0.5, 15.0), C=(0.5, 8.0, 0.1, 4.0),
                      MRW=(1e-2, 80.0), CW=(-30.0, 30.0, -20.0, 20.0)),
    "BesselInu": dict(t="BesselI[{1}, {0}]", f=mp.besseli, nu=(-4.7, 9.3), h=h_bessel,
                      MR=(0.1, 8.0), C=(0.5, 4.0, -3.0, 3.0),
                      MRW=(1e-2, 400.0), CW=(-30.0, 30.0, -30.0, 30.0)),
    "BesselKnu": dict(t="BesselK[{1}, {0}]", f=mp.besselk, nu=(-4.7, 9.3), h=h_bessel,
                      MR=(0.1, 8.0), C=(0.5, 4.0, 0.1, 3.0),
                      MRW=(1e-2, 400.0), CW=(-30.0, 30.0, -30.0, 30.0)),
    "AiryAi": dict(t="AiryAi[{0}]", f=AI, nu=None, h=h_airy,
                   MR=(-5.0, 5.0), C=(-3.0, 3.0, -3.0, 3.0),
                   MRW=(-60.0, 100.0), CW=(-15.0, 15.0, -15.0, 15.0)),
    "AiryBi": dict(t="AiryBi[{0}]", f=BI, nu=None, h=h_airy,
                   MR=(-5.0, 4.0), C=(-3.0, 2.0, -3.0, 3.0),
                   MRW=(-60.0, 100.0), CW=(-15.0, 15.0, -15.0, 15.0)),
    "AiryAiPrime": dict(t="AiryAiPrime[{0}]", f=AIP, nu=None, h=h_airy,
                        MR=(-5.0, 5.0), C=(-3.0, 3.0, -3.0, 3.0),
                        MRW=(-60.0, 100.0), CW=(-15.0, 15.0, -15.0, 15.0)),
    "AiryBiPrime": dict(t="AiryBiPrime[{0}]", f=BIP, nu=None, h=h_airy,
                        MR=(-5.0, 4.0), C=(-3.0, 2.0, -3.0, 3.0),
                        MRW=(-60.0, 100.0), CW=(-15.0, 15.0, -15.0, 15.0)),
}


def ulp(x: float) -> float:
    x = abs(float(x))
    if x == 0.0 or not math.isfinite(x):
        return 2.0 ** -1074
    if x < 2.0 ** -1022:
        return 2.0 ** -1074
    return 2.0 ** (math.frexp(x)[1] - 53)


NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def parse_value(s: str):
    s = s.strip().replace("*^", "e").replace("`", "")
    if s in ("", "Indeterminate", "ComplexInfinity") or "[" in s:
        return None
    m = re.fullmatch(rf"({NUM})\s*([-+])\s*({NUM})\*I", s)
    if m:
        return mp.mpc(mp.mpf(m.group(1)), mp.mpf(m.group(3)) * (1 if m.group(2) == "+" else -1))
    m = re.fullmatch(rf"({NUM})\*I", s)
    if m:
        return mp.mpc(0, mp.mpf(m.group(1)))
    m = re.fullmatch(rf"({NUM})", s)
    if m:
        return mp.mpf(m.group(1))
    return None


def sample(spec, regime, n, rng):
    dom = spec[regime]
    pts = []
    for _ in range(n):
        if regime in ("MR", "MRW"):
            lo, hi = dom
            if lo > 0 and hi / lo > 1e3:
                x = math.exp(rng.uniform(math.log(lo), math.log(hi)))
            else:
                x = rng.uniform(lo, hi)
        else:
            a, b, c, d = dom
            x = complex(rng.uniform(a, b), rng.uniform(c, d))
        nu = spec["nu"]
        if isinstance(nu, tuple):
            if spec.get("intorder"):
                nu = rng.randint(nu[0], nu[1])
            else:
                nu = rng.uniform(*nu)
        pts.append((x, nu))
    return pts


def lit_order(nu):
    if nu is None:
        return ""
    if isinstance(nu, int):
        return str(nu)
    if float(nu).is_integer() and abs(nu) < 1e6:
        return str(int(nu))
    return wl_real(nu)


def run_mathilda(binary, spec, pts, mode):
    tmpl = spec["t"]
    lines = []
    if mode == "scalar":
        for x, nu in pts:
            lit = wl_complex(x) if isinstance(x, complex) else wl_real(x)
            lines.append(f"Print[ToString[InputForm[SetPrecision[{tmpl.format(lit, lit_order(nu))}, 17]]]]")
    else:
        # one NDArray per order value (the benchmark orders are fixed)
        groups = {}
        for i, (x, nu) in enumerate(pts):
            groups.setdefault(nu, []).append(i)
        order = []
        for nu, idx in groups.items():
            lits = ", ".join(wl_complex(pts[i][0]) if isinstance(pts[i][0], complex)
                             else wl_real(pts[i][0]) for i in idx)
            lines.append(f"pts$ = NDArray[{{{lits}}}];")
            lines.append(f"res$ = {tmpl.format('pts$', lit_order(nu))};")
            lines.append('If[Head[res$] =!= NDArray, Print["#DECLINED"]];')
            lines.append("Scan[Print[ToString[InputForm[SetPrecision[#, 17]]]] &, Normal[res$]]")
            order.extend(idx)
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        f.write("\n".join(lines) + "\n")
        path = f.name
    try:
        out = subprocess.run([binary, "-file", path], capture_output=True, text=True,
                             timeout=3600).stdout
    finally:
        os.unlink(path)
    vals = [v for v in out.strip().split("\n")]
    declined = sum(1 for v in vals if v == "#DECLINED")
    vals = [v for v in vals if v != "#DECLINED"]
    if declined:
        print(f"      [{declined} NDArray group(s) declined; values are the List path]")
    if len(vals) != len(pts):
        raise RuntimeError(f"got {len(vals)} lines for {len(pts)} points:\n{out[:2000]}")
    if mode != "scalar":
        reordered = [None] * len(pts)
        for k, i in enumerate(order):
            reordered[i] = vals[k]
        vals = reordered
    return vals


def reference(spec, x, nu):
    f = spec["f"]
    z = mp.mpc(x.real, x.imag) if isinstance(x, complex) else mp.mpf(x)
    if nu is None:
        g = f
    else:
        n = mp.mpf(nu)
        g = lambda t: f(n, t)
    v = g(z)
    try:
        d = mp.diff(g, z)
        scale = max(abs(v), abs(d) * spec["h"](x))
    except Exception:
        scale = abs(v)
    return v, scale


def evaluate(binary, head, regime, n, mode, seed, worst):
    spec = HEADS[head]
    rng = random.Random(seed ^ zlib.crc32((head + regime).encode()))
    pts = sample(spec, regime, n, rng)
    raw = run_mathilda(binary, spec, pts, mode)
    rels, scs, bad = [], [], []
    for (x, nu), s in zip(pts, raw):
        v, scale = reference(spec, x, nu)
        if not mp.isfinite(abs(v)) or abs(v) > 1e300 or (0 < abs(v) < 1e-300):
            continue                      # outside the double range: not a machine answer
        got = parse_value(s)
        if got is None:
            bad.append((x, nu, s))
            continue
        err = abs(mp.mpc(got) - mp.mpc(v))
        mag = abs(v)
        rel = float(err) / ulp(float(mag)) if mag != 0 else float("inf")
        sc = float(err) / ulp(float(max(mag, scale)))
        rels.append((rel, x, nu, s, v, sc))
        scs.append(sc)
    if not rels:
        return {"head": head, "regime": regime, "n": n, "unevaluated": len(bad), "bad": bad[:3]}
    rr = sorted(r[0] for r in rels)
    res = {
        "head": head, "regime": regime, "n": len(rels),
        "rel_median": statistics.median(rr), "rel_max": rr[-1],
        "scaled_median": statistics.median(scs), "scaled_max": max(scs),
        "unevaluated": len(bad),
    }
    res["worst"] = [(r[1], r[2], r[3], mp.nstr(r[4], 20), round(r[5], 2))
                    for r in sorted(rels, key=lambda t: -t[5])[:worst]]
    res["bad"] = bad[:3]
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--heads", default="BesselJ,BesselY,BesselI,BesselK,AiryAi,AiryBi,"
                                       "AiryAiPrime,AiryBiPrime")
    ap.add_argument("--regimes", default="MR,C,MRW,CW")
    ap.add_argument("--mode", default="scalar", choices=("scalar", "packed"))
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--bin", default=os.path.join(WORKTREE, "Mathilda"))
    ap.add_argument("--worst", type=int, default=3)
    a = ap.parse_args()
    print(f"{'head':12s} {'reg':4s} {'n':>5s} {'rel med':>8s} {'rel max':>10s} "
          f"{'scl med':>8s} {'scl max':>8s} {'uneval':>6s}")
    for h in a.heads.split(","):
        for r in a.regimes.split(","):
            res = evaluate(a.bin, h, r, a.n, a.mode, a.seed, a.worst)
            if "rel_median" not in res:
                print(f"{h:12s} {r:4s} all {res['unevaluated']} unevaluated")
                for b in res["bad"]:
                    print(f"      uneval: x={b[0]!r} nu={b[1]!r} -> {b[2][:80]}")
                continue
            print(f"{h:12s} {r:4s} {res['n']:5d} {res['rel_median']:8.2f} {res['rel_max']:10.2f} "
                  f"{res['scaled_median']:8.2f} {res['scaled_max']:8.2f} {res['unevaluated']:6d}")
            for w in res["worst"]:
                print(f"      worst: x={w[0]!r} nu={w[1]!r} got={w[2]} ref={w[3]} ({w[4]} scaled ulp)")
            for b in res["bad"]:
                print(f"      uneval: x={b[0]!r} nu={b[1]!r} -> {b[2][:80]}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
