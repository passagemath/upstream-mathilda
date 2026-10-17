#!/usr/bin/env python3
"""Accuracy gate for the hypergeometric-family machine kernels.

Drives `Mathilda -file` on a generated script that evaluates each head at
machine-precision points (machine literals with a trailing backtick, complex as
Complex[a`, b`]) and prints ToString[InputForm[SetPrecision[r, 17]]], then
compares against mpmath at 30 digits.

Two error measures per point, in units of the double ulp of the reference:
  rel     |m - ref| / ulp(|ref|)
  scaled  |m - ref| / ulp(local magnitude), local magnitude = max |f| over a
          +-30% neighbourhood of the point (13 samples). This is the brief's
          "absolute error scaled by the local magnitude" near a zero of an
          oscillatory function; it equals `rel` where f is monotone.

Usage (benchmark venv, isolated):
  benchmarks/.venv/bin/python -I tools/sf_accuracy/hyper.py [--bin PATH]
        [--n N] [--only SUBSTR] [--json OUT]
"""
import argparse
import json
import math
import os
import random
import subprocess
import sys
import tempfile

import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def mlit(x):
    """Machine literal for a double: shortest repr plus a backtick."""
    r = repr(float(x))
    if "e" in r or "E" in r:
        m, e = r.lower().split("e")
        return f"{m}`*^{int(e)}"
    return r + "`"


def marg(v):
    if isinstance(v, complex):
        return f"Complex[{mlit(v.real)}, {mlit(v.imag)}]"
    if isinstance(v, str):
        return v                       # exact parameter, verbatim
    return mlit(v)


def parse_val(s):
    s = s.strip()
    if s.startswith("$FAIL") or "[" in s or not s:
        return None
    s = s.replace("*^", "e").replace("`", "")
    # forms: "a", "a + b*I", "a - b*I", "b*I"
    if "I" in s:
        t = s.replace(" ", "")
        if t.endswith("*I"):
            t = t[:-2]
        elif t.endswith("I"):
            t = t[:-1] or "1"
        # split at the last +/- that is not an exponent sign
        k = len(t) - 1
        while k > 0:
            if t[k] in "+-" and t[k - 1] not in "eE":
                break
            k -= 1
        if k <= 0:
            return complex(0.0, float(t))
        return complex(float(t[:k]), float(t[k:]))
    return float(s)


def ulp(x):
    x = abs(float(x))
    if x == 0.0 or not math.isfinite(x):
        return 5e-324
    return math.ulp(x)


# ---------------------------------------------------------------- the cases
# Each case: label, Mathilda template with {z}, list of exact/mpmath params,
# mpmath function f(z), sampler -> list of points (float or complex).

def U(lo, hi):
    return lambda rng: rng.uniform(lo, hi)


def B(xl, xh, yl, yh):
    return lambda rng: complex(rng.uniform(xl, xh), rng.uniform(yl, yh))


def F0(b):
    return lambda z: mp.hyp0f1(b, z)


def F1(a, b):
    return lambda z: mp.hyp1f1(a, b, z)


def F2(a, b, c):
    return lambda z: mp.hyp2f1(a, b, c, z)


half = mp.mpf(1) / 2


def D(num, den=1):
    """A parameter as a machine engine sees it: the exact rational num/den
    rounded to the nearest double (3/10 is evaluated as 0.299999...989). The
    reference is the function AT that double, so parameter rounding -- which no
    machine kernel can undo -- is not charged as kernel error."""
    return mp.mpf(float(mp.mpf(num) / den))


CASES = [
    # ---- benchmark domains (the gate)
    ("0F1(;2;x) bench MR", "Hypergeometric0F1[2, {z}]", F0(2), U(-5, 5), True),
    ("0F1(;2;z) bench C", "Hypergeometric0F1[2, {z}]", F0(2), B(-3, 3, -3, 3), True),
    ("1F1(1;2;x) bench MR", "Hypergeometric1F1[1, 2, {z}]", F1(1, 2), U(-5, 5), True),
    ("1F1(1;2;z) bench C", "Hypergeometric1F1[1, 2, {z}]", F1(1, 2), B(-3, 3, -3, 3), True),
    ("2F1(1,2;3;x) bench MR", "Hypergeometric2F1[1, 2, 3, {z}]", F2(1, 2, 3), U(-0.9, 0.9), True),
    ("2F1(1,2;3;z) bench C", "Hypergeometric2F1[1, 2, 3, {z}]", F2(1, 2, 3), B(-0.8, 0.8, -0.8, 0.8), True),
    ("2F2(1,2;3,4;x) bench MR", "HypergeometricPFQ[{{1, 2}}, {{3, 4}}, {z}]",
     lambda z: mp.hyper([1, 2], [3, 4], z), U(-3, 3), True),
    ("2F2(1,2;3,4;z) bench C", "HypergeometricPFQ[{{1, 2}}, {{3, 4}}, {z}]",
     lambda z: mp.hyper([1, 2], [3, 4], z), B(-2, 2, -2, 2), True),
    ("U(1,2,x) bench MR", "HypergeometricU[1, 2, {z}]", lambda z: mp.hyperu(1, 2, z), U(0.2, 8), True),
    ("U(1,2,z) bench C", "HypergeometricU[1, 2, {z}]", lambda z: mp.hyperu(1, 2, z), B(0.3, 4, -2, 2), True),
    ("M(1,1/2,x) bench MR", "WhittakerM[1, 1/2, {z}]", lambda z: mp.whitm(1, half, z), U(0.5, 8), True),
    ("M(1,1/2,z) bench C", "WhittakerM[1, 1/2, {z}]", lambda z: mp.whitm(1, half, z), B(0.5, 4, -2, 2), True),
    ("W(1,1/2,x) bench MR", "WhittakerW[1, 1/2, {z}]", lambda z: mp.whitw(1, half, z), U(0.5, 8), True),
    ("W(1,1/2,z) bench C", "WhittakerW[1, 1/2, {z}]", lambda z: mp.whitw(1, half, z), B(0.5, 4, -2, 2), True),
    # ---- wider bands and general parameters (documented, looser allowed)
    ("0F1(;2;x) wide", "Hypergeometric0F1[2, {z}]", F0(2), U(-400, 2000), False),
    ("0F1(;0.3;x)", "Hypergeometric0F1[3/10, {z}]", F0(D(3, 10)), U(-60, 60), False),
    ("0F1(;7.5;x)", "Hypergeometric0F1[15/2, {z}]", F0(D(15, 2)), U(-200, 200), False),
    ("0F1(;-1.5;x)", "Hypergeometric0F1[-3/2, {z}]", F0(D(-3, 2)), U(-30, 30), False),
    ("0F1(;2;z) wide C", "Hypergeometric0F1[2, {z}]", F0(2), B(-15, 15, -15, 15), False),
    ("1F1(1;2;x) wide", "Hypergeometric1F1[1, 2, {z}]", F1(1, 2), U(-300, 300), False),
    ("1F1(0.3;1.7;x)", "Hypergeometric1F1[3/10, 17/10, {z}]", F1(D(3, 10), D(17, 10)), U(-60, 60), False),
    ("1F1(-2.5;3.2;x)", "Hypergeometric1F1[-5/2, 16/5, {z}]", F1(D(-5, 2), D(16, 5)), U(-60, 60), False),
    ("1F1(2.5;0.7;x)", "Hypergeometric1F1[5/2, 7/10, {z}]", F1(D(5, 2), D(7, 10)), U(-60, 60), False),
    ("1F1(-3;2;x) poly", "Hypergeometric1F1[-3, 2, {z}]", F1(-3, 2), U(-30, 30), False),
    ("1F1(0.3;1.7;z) C", "Hypergeometric1F1[3/10, 17/10, {z}]", F1(D(3, 10), D(17, 10)), B(-12, 12, -12, 12), False),
    ("2F1(1,2;3;x) wide", "Hypergeometric2F1[1, 2, 3, {z}]", F2(1, 2, 3), U(-30, 0.999), False),
    ("2F1(.3,.7;1.9;x)", "Hypergeometric2F1[3/10, 7/10, 19/10, {z}]",
     F2(D(3, 10), D(7, 10), D(19, 10)), U(-30, 0.999), False),
    ("2F1(1.5,2.25;3.1;x)", "Hypergeometric2F1[3/2, 9/4, 31/10, {z}]",
     F2(D(3, 2), D(9, 4), D(31, 10)), U(-30, 0.999), False),
    ("2F1(.5,.5;1;x) m=0", "Hypergeometric2F1[1/2, 1/2, 1, {z}]", F2(half, half, 1), U(-30, 0.999), False),
    ("2F1(.5,1;2.5;x) m=1", "Hypergeometric2F1[1/2, 1, 5/2, {z}]", F2(half, 1, D(5, 2)), U(-30, 0.999), False),
    ("2F1(1.5,2;2.5;x) m=-1", "Hypergeometric2F1[3/2, 2, 5/2, {z}]", F2(D(3, 2), 2, D(5, 2)), U(-30, 0.999), False),
    ("2F1(1,2;3;x>1) cut", "Hypergeometric2F1[1, 2, 3, {z}]", F2(1, 2, 3), U(1.001, 20), False),
    ("2F1(1,2;3;z) wide C", "Hypergeometric2F1[1, 2, 3, {z}]", F2(1, 2, 3), B(-4, 4, -4, 4), False),
    ("2F1(.3,.7;1.9;z) C", "Hypergeometric2F1[3/10, 7/10, 19/10, {z}]",
     F2(D(3, 10), D(7, 10), D(19, 10)), B(-3, 3, -3, 3), False),
    ("2F2(1,2;3,4;x) wide", "HypergeometricPFQ[{{1, 2}}, {{3, 4}}, {z}]",
     lambda z: mp.hyper([1, 2], [3, 4], z), U(-60, 60), False),
    ("3F2(.5,1,1.5;2,2.5;x)", "HypergeometricPFQ[{{1/2, 1, 3/2}}, {{2, 5/2}}, {z}]",
     lambda z: mp.hyper([half, 1, D(3, 2)], [2, D(5, 2)], z), U(-0.95, 0.95), False),
    ("U(0.7,1.3,x)", "HypergeometricU[7/10, 13/10, {z}]", lambda z: mp.hyperu(D(7, 10), D(13, 10), z), U(0.05, 40), False),
    ("U(1.5,3,x)", "HypergeometricU[3/2, 3, {z}]", lambda z: mp.hyperu(D(3, 2), 3, z), U(0.05, 40), False),
    ("U(0.7,1.3,z) C", "HypergeometricU[7/10, 13/10, {z}]", lambda z: mp.hyperu(D(7, 10), D(13, 10), z), B(0.1, 10, -8, 8), False),
    ("M(0.3,0.7,x)", "WhittakerM[3/10, 7/10, {z}]", lambda z: mp.whitm(D(3, 10), D(7, 10), z), U(0.05, 40), False),
    ("M(0.3,0.7,z) C", "WhittakerM[3/10, 7/10, {z}]", lambda z: mp.whitm(D(3, 10), D(7, 10), z), B(-8, 8, -8, 8), False),
    ("W(0.3,0.7,x)", "WhittakerW[3/10, 7/10, {z}]", lambda z: mp.whitw(D(3, 10), D(7, 10), z), U(0.05, 40), False),
    ("W(0.3,0.7,z) C", "WhittakerW[3/10, 7/10, {z}]", lambda z: mp.whitw(D(3, 10), D(7, 10), z), B(0.1, 10, -8, 8), False),
]


def ref_val(f, z):
    try:
        v = f(mp.mpc(z.real, z.imag) if isinstance(z, complex) else mp.mpf(z))
    except Exception:
        return None
    if isinstance(z, complex) or (isinstance(v, mp.mpc) and abs(v.imag) > 0):
        v = mp.mpc(v)
        return complex(float(v.real), float(v.imag)), v
    return float(v), mp.mpf(v)


def local_mag(f, z, refmag):
    mags = [refmag]
    for k in range(-6, 7):
        if k == 0:
            continue
        zz = z * (1 + 0.05 * k)
        try:
            v = f(mp.mpc(zz.real, zz.imag) if isinstance(zz, complex) else mp.mpf(zz))
            mags.append(float(abs(v)))
        except Exception:
            pass
    return max(m for m in mags if math.isfinite(m))


def run_case(binpath, case, n, seed):
    label, tmpl, f, sampler, _gate = case
    rng = random.Random(seed)
    pts = [sampler(rng) for _ in range(n)]
    lines = []
    for z in pts:
        expr = tmpl.format(z=marg(z))
        lines.append(f'Print[ToString[InputForm[SetPrecision[{expr}, 17]]]];')
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write("\n".join(lines) + "\n")
        path = fh.name
    try:
        out = subprocess.run([binpath, "-file", path], capture_output=True, text=True,
                             timeout=3600).stdout.splitlines()
    finally:
        os.unlink(path)
    rel, sca, fails, undef = [], [], 0, 0
    worst = (0.0, None, None, None)
    for z, s in zip(pts, out):
        r = ref_val(f, z)
        if r is None:
            continue
        rv, rhi = r
        m = parse_val(s)
        if m is None:
            undef += 1
            continue
        if isinstance(rv, complex) or isinstance(m, complex):
            mm = complex(m)
            err = abs(mp.mpc(mm.real, mm.imag) - mp.mpc(rhi))
        else:
            err = abs(mp.mpf(m) - rhi)
        mag = float(abs(rhi))
        if not math.isfinite(mag) or mag == 0.0:
            continue
        e_rel = float(err) / ulp(mag)
        lm = local_mag(f, z, mag) if e_rel > 2.0 else mag
        e_sca = float(err) / ulp(lm)
        rel.append(e_rel)
        sca.append(e_sca)
        if e_sca > worst[0]:
            worst = (e_sca, z, s, rv)
    if len(out) < len(pts):
        fails = len(pts) - len(out)

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else float("nan")
    return {
        "label": label, "n": len(pts), "evaluated": len(rel), "unevaluated": undef,
        "missing": fails,
        "rel_median": med(rel), "rel_max": max(rel) if rel else float("nan"),
        "scaled_median": med(sca), "scaled_max": max(sca) if sca else float("nan"),
        "worst": [worst[0], repr(worst[1]), worst[2], repr(worst[3])],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=os.path.join(ROOT, "Mathilda"))
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--only", default="")
    ap.add_argument("--json", default="")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    res = []
    print(f"{'case':28s} {'n':>5s} {'unev':>5s} {'rel.med':>8s} {'rel.max':>10s} "
          f"{'sc.med':>8s} {'sc.max':>10s}  worst")
    for case in CASES:
        if args.only and args.only not in case[0]:
            continue
        r = run_case(args.bin, case, args.n if case[4] else max(args.n // 2, 200),
                     args.seed)
        res.append(r)
        print(f"{r['label']:28s} {r['evaluated']:5d} {r['unevaluated']:5d} "
              f"{r['rel_median']:8.2f} {r['rel_max']:10.2f} {r['scaled_median']:8.2f} "
              f"{r['scaled_max']:10.2f}  z={r['worst'][1]} got={r['worst'][2]}", flush=True)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
