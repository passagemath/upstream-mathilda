#!/usr/bin/env python3
"""Accuracy gate for the "expint" family of machine special-function kernels.

Heads: SinIntegral, CosIntegral, SinhIntegral, CoshIntegral, ExpIntegralEi,
ExpIntegralE[2, x], LogIntegral, EllipticK, EllipticE (complete), ProductLog,
LegendreP[3, x], LegendreQ[2, x].

Drives `<worktree>/Mathilda -file` on a generated .m that prints each value with
ToString[InputForm[SetPrecision[r, 17]]] (17 digits round-trip a double), and
compares against mpmath at 30 digits.

Two error measures per point, both in units of ulp(reference):
  rel     |f~ - f| / ulp(|f|)                       the plain relative error
  scaled  |f~ - f| / ulp(max(|f|, |x f'(x)|))       the error relative to the
          change a 1-ulp perturbation of the ARGUMENT makes.  Equal to `rel`
          except near a zero of f, where the relative error is meaningless and
          the brief asks for "absolute error scaled by the local magnitude".

Modes:
  scalar   one statement per point:  F[x`]        (the scalar builtin path)
  packed   F[{x1`, x2`, ...}] on a packed buffer   (the ND kernel / helper path)

Usage:
  benchmarks/.venv/bin/python -I tools/sf_accuracy/expint.py \
      [--heads SinIntegral,CosIntegral] [--regimes MR,C,MRW,CW] [--mode scalar|packed]
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


# name -> (template, mp function, MR domain, C domain, MR wide, C wide)
# Domains: real (lo, hi); complex (re_lo, re_hi, im_lo, im_hi).
HEADS = {
    "SinIntegral":  ("SinIntegral[{0}]", mp.si,  (0.1, 15.0), (0.5, 6.0, -3.0, 3.0),
                     (-60.0, 60.0), (-30.0, 30.0, -12.0, 12.0)),
    "CosIntegral":  ("CosIntegral[{0}]", mp.ci,  (0.1, 15.0), (0.5, 6.0, 0.1, 3.0),
                     (1e-3, 80.0), (-30.0, 30.0, -12.0, 12.0)),
    "SinhIntegral": ("SinhIntegral[{0}]", mp.shi, (0.1, 8.0), (0.5, 4.0, -2.0, 2.0),
                     (-60.0, 60.0), (-12.0, 12.0, -30.0, 30.0)),
    "CoshIntegral": ("CoshIntegral[{0}]", mp.chi, (0.1, 8.0), (0.5, 4.0, 0.1, 2.0),
                     (1e-3, 80.0), (-12.0, 12.0, -30.0, 30.0)),
    "ExpIntegralEi": ("ExpIntegralEi[{0}]", mp.ei, (0.1, 8.0), (0.5, 4.0, -3.0, 3.0),
                      (-60.0, 80.0), (-30.0, 30.0, -30.0, 30.0)),
    "ExpIntegralE": ("ExpIntegralE[2, {0}]", (lambda z: mp.expint(2, z)), (0.2, 8.0),
                     (0.5, 4.0, -3.0, 3.0), (1e-3, 60.0), (-10.0, 30.0, -20.0, 20.0)),
    "LogIntegral":  ("LogIntegral[{0}]", mp.li, (2.0, 20.0), (2.0, 10.0, -3.0, 3.0),
                     (1e-3, 1e6), (-30.0, 30.0, -30.0, 30.0)),
    "EllipticK":    ("EllipticK[{0}]", mp.ellipk, (0.0, 0.95), (-0.9, 0.9, -0.9, 0.9),
                     (-50.0, 0.9999), (-20.0, 20.0, -20.0, 20.0)),
    "EllipticE":    ("EllipticE[{0}]", mp.ellipe, (0.0, 0.95), (-0.9, 0.9, -0.9, 0.9),
                     (-50.0, 1.0), (-20.0, 20.0, -20.0, 20.0)),
    "ProductLog":   ("ProductLog[{0}]", mp.lambertw, (-0.3, 10.0), (-0.3, 5.0, -2.0, 2.0),
                     (-0.36787944117144233, 1e6), (-30.0, 30.0, -30.0, 30.0)),
    "LegendreP":    ("LegendreP[3, {0}]", (lambda z: mp.legendre(3, z)), (-0.95, 0.95),
                     (-0.9, 0.9, -0.9, 0.9), (-20.0, 20.0), (-10.0, 10.0, -10.0, 10.0)),
    "LegendreQ":    ("LegendreQ[2, {0}]", (lambda z: mp.legenq(2, 0, z, type=2)), (-0.9, 0.9),
                     (-0.8, 0.8, -0.8, 0.8), (-0.9999, 0.9999), (-5.0, 5.0, -5.0, 5.0)),
}

# Complex-domain references whose principal branch differs from mpmath's default.
# LegendreQ type 1 (Mathematica's LegendreQ[n, z]) off the cut is
# P_n(z) (log(1+z) - log(1-z))/2 - W_{n-1}(z); mpmath's legenq(type=2) agrees with
# it inside |x| < 1 on the real line only, so the complex reference is built here.
def _legq2_type1(z):
    z = mp.mpmathify(z)
    p2 = (3 * z * z - 1) / 2
    return p2 * (mp.log(1 + z) - mp.log(1 - z)) / 2 - 3 * z / 2


COMPLEX_REF = {"LegendreQ": _legq2_type1}


def ulp(x: float) -> float:
    x = abs(float(x))
    if x == 0.0 or not math.isfinite(x):
        return 2.0 ** -1074
    return 2.0 ** (math.frexp(x)[1] - 53)


NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def parse_value(s: str):
    s = s.strip().replace("*^", "e")
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


def sample(regime, dom, n, rng):
    pts = []
    if regime in ("MR", "MRW"):
        lo, hi = dom
        if regime == "MRW" and lo > 0 and hi / lo > 1e3:          # log-spaced band
            for _ in range(n):
                pts.append(math.exp(rng.uniform(math.log(lo), math.log(hi))))
        else:
            for _ in range(n):
                pts.append(rng.uniform(lo, hi))
    else:
        a, b, c, d = dom
        for _ in range(n):
            pts.append(complex(rng.uniform(a, b), rng.uniform(c, d)))
    return pts


def run_mathilda(binary, head, pts, mode):
    tmpl = HEADS[head][0]
    lines = []
    if mode == "scalar":
        for p in pts:
            lit = wl_complex(p) if isinstance(p, complex) else wl_real(p)
            lines.append(f"Print[ToString[InputForm[SetPrecision[{tmpl.format(lit)}, 17]]]]")
    else:
        lits = ", ".join(wl_complex(p) if isinstance(p, complex) else wl_real(p) for p in pts)
        # A visible NDArray forces the element-wise kernel; if ANY element
        # declines, the whole array falls back to the List path (reported).
        lines.append(f"pts$ = NDArray[{{{lits}}}];")
        lines.append(f"res$ = {tmpl.format('pts$')};")
        lines.append('If[Head[res$] =!= NDArray, Print["#DECLINED"]];')
        lines.append("Scan[Print[ToString[InputForm[SetPrecision[#, 17]]]] &, Normal[res$]]")
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        f.write("\n".join(lines) + "\n")
        path = f.name
    try:
        out = subprocess.run([binary, "-file", path], capture_output=True, text=True,
                             timeout=1800).stdout
    finally:
        os.unlink(path)
    vals = out.strip().split("\n")
    if vals and vals[0] == "#DECLINED":
        print(f"      [{head}: the NDArray kernel declined an element; values are the List path]")
        vals = vals[1:]
    if len(vals) != len(pts):
        raise RuntimeError(f"{head}: got {len(vals)} lines for {len(pts)} points:\n{out[:2000]}")
    return vals


def reference(head, p):
    fn = HEADS[head][1]
    if isinstance(p, complex) and head in COMPLEX_REF:
        fn = COMPLEX_REF[head]
    x = mp.mpc(p.real, p.imag) if isinstance(p, complex) else mp.mpf(p)
    v = fn(x)
    try:
        d = mp.diff(fn, x)
        scale = max(abs(v), abs(x * d))
    except Exception:
        scale = abs(v)
    return v, scale


def evaluate(binary, head, regime, n, mode, seed, worst):
    tmpl, fn, dmr, dc, dmrw, dcw = HEADS[head]
    dom = {"MR": dmr, "C": dc, "MRW": dmrw, "CW": dcw}[regime]
    rng = random.Random(seed ^ zlib.crc32((head + regime).encode()))
    pts = sample(regime, dom, n, rng)
    raw = run_mathilda(binary, head, pts, mode)
    rels, scs, bad = [], [], []
    for p, s in zip(pts, raw):
        v, scale = reference(head, p)
        got = parse_value(s)
        if got is None:
            bad.append((p, s))
            continue
        if regime in ("MR", "MRW") and isinstance(v, mp.mpc) and abs(v.imag) > 0:
            # off the real axis: compare as complex
            pass
        err = abs(mp.mpc(got) - mp.mpc(v))
        mag = abs(v)
        rel = float(err) / ulp(float(mag)) if mag != 0 else float("inf")
        sc = float(err) / ulp(float(max(mag, scale)))
        rels.append((rel, p, s, v))
        scs.append((sc, p, s, v))
    if not rels:
        return {"head": head, "regime": regime, "n": n, "unevaluated": len(bad)}
    rr = sorted(r[0] for r in rels)
    res = {
        "head": head, "regime": regime, "n": len(rels),
        "rel_median": statistics.median(rr), "rel_max": rr[-1],
        "scaled_median": statistics.median(t[0] for t in scs), "scaled_max": max(t[0] for t in scs),
        "unevaluated": len(bad),
    }
    res["worst"] = [(r[1], r[2], mp.nstr(r[3], 20), round(r[0], 2))
                    for r in sorted(rels, key=lambda t: -t[0])[:worst]]
    ws = max(scs, key=lambda t: t[0])
    res["worst_scaled"] = (ws[1], ws[2], mp.nstr(ws[3], 20), round(ws[0], 2))
    res["bad"] = bad[:3]
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--heads", default=",".join(HEADS))
    ap.add_argument("--regimes", default="MR,C,MRW,CW")
    ap.add_argument("--mode", default="scalar", choices=("scalar", "packed"))
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--bin", default=os.path.join(WORKTREE, "Mathilda"))
    ap.add_argument("--worst", type=int, default=3)
    a = ap.parse_args()
    print(f"{'head':14s} {'reg':4s} {'n':>5s} {'rel med':>8s} {'rel max':>10s} "
          f"{'scl med':>8s} {'scl max':>8s} {'uneval':>6s}")
    for h in a.heads.split(","):
        for r in a.regimes.split(","):
            res = evaluate(a.bin, h, r, a.n, a.mode, a.seed, a.worst)
            if "rel_median" not in res:
                print(f"{h:14s} {r:4s} all {res['unevaluated']} unevaluated")
                continue
            print(f"{h:14s} {r:4s} {res['n']:5d} {res['rel_median']:8.2f} {res['rel_max']:10.2f} "
                  f"{res['scaled_median']:8.2f} {res['scaled_max']:8.2f} {res['unevaluated']:6d}")
            for w in res["worst"]:
                print(f"      worst: x={w[0]!r} got={w[1]} ref={w[2]} ({w[3]} ulp)")
            w = res["worst_scaled"]
            print(f"      worst scaled: x={w[0]!r} got={w[1]} ref={w[2]} ({w[3]} ulp)")
            for b in res["bad"]:
                print(f"      uneval: x={b[0]!r} -> {b[1][:80]}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
