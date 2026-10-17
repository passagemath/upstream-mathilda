#!/usr/bin/env python3
"""Accuracy gate for the zeta-family machine kernels (Zeta, HurwitzZeta,
LerchPhi, PolyLog), measured through the REAL scalar builtins.

For every (head, regime) it draws N machine-double argument tuples, writes a .m
file that prints each value as ToString[InputForm[SetPrecision[F[...], 17]]],
runs `<Mathilda> -file`, and compares against mpmath at 90 digits for Zeta /
HurwitzZeta and 50 for LerchPhi / PolyLog. (30 or even 50 digits are not
enough for the zeta functions: mpmath's Hurwitz zeta works to an ABSOLUTE
tolerance, so a tiny result -- 3.7e-45 at s = 27, a = 46 -- is wrong in the
15th digit at dps = 50, which first showed up here as false 1e-10 "errors".)

Errors are reported in ulps:
  * real results: |got - exact| / ulp(exact)
  * complex results: |got - exact| / (|exact| 2^-52)
Next to the plain figure, `cond` divides each error by the function's
condition number kappa = max(1, |x f'(x) / f(x)| summed over the arguments):
the error a backward-stable evaluation would be entitled to. That is the
"absolute error scaled by the local magnitude" view the brief asks for near
zeros / cancellation (a complex Hurwitz sum whose value is small against its
terms), and is shown for transparency, not used to excuse a plain-ulp failure.

Usage:
  benchmarks/.venv/bin/python -I tools/sf_accuracy/zeta.py [--bin PATH] [--n N]
                                                          [--only Zeta,PolyLog]
"""
import argparse
import math
import os
import random
import re
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BIN = os.path.normpath(os.path.join(HERE, "..", "..", "Mathilda"))
DPS = 50


def lit(x):
    """A machine-double literal Mathilda reads back bit-exactly."""
    return repr(float(x)) + "`"


def clit(re_, im_):
    return "Complex[%s, %s]" % (lit(re_), lit(im_))


NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def parse(line):
    """Parse one InputForm line into a complex (None if not numeric)."""
    s = line.strip().replace("`", "")
    m = re.fullmatch(r"(%s)" % NUM, s)
    if m:
        return complex(float(m.group(1)), 0.0), False
    m = re.fullmatch(r"(%s) ([-+]) (%s)\*I" % (NUM, NUM), s)
    if m:
        im_ = float(m.group(3)) * (-1 if m.group(2) == "-" else 1)
        return complex(float(m.group(1)), im_), True
    m = re.fullmatch(r"(%s)\*I" % NUM, s)
    if m:
        return complex(0.0, float(m.group(1))), True
    return None, False


# ---------------------------------------------------------------- the cases
# Each case: (label, head, n_args_template, generator, reference, cond)
#   generator(rnd) -> tuple of argument values (float or complex)
#   reference(*mp_args) -> exact value
#   cond(*mp_args, value) -> condition number (>= 1)

def U(r, lo, hi):
    return r.uniform(lo, hi)


def C(r, a, b, c, d):
    return complex(r.uniform(a, b), r.uniform(c, d))


def k_zeta(s, v):
    return abs(s * mp.zeta(s, 1, 1) / v)


def k_hurwitz(s, a, v):
    return abs(s * mp.zeta(s, a, 1) / v) + abs(a * s * mp.zeta(s + 1, a) / v)


def k_polylog(s, z, v):
    return abs(mp.polylog(s - 1, z) / v)


CASES = [
    # --- Zeta
    ("Zeta MR bench s in [1.5, 9]", "Zeta[%s]",
     lambda r: (U(r, 1.5, 9),), lambda s: mp.zeta(s), k_zeta),
    ("Zeta R wide s in [0.5, 60]", "Zeta[%s]",
     lambda r: (U(r, 0.5, 60),), lambda s: mp.zeta(s), k_zeta),
    ("Zeta R wide s in [-40, 0.5]", "Zeta[%s]",
     lambda r: (U(r, -40, 0.5),), lambda s: mp.zeta(s), k_zeta),
    ("Zeta C bench Re[1.5,5] Im[-3,3]", "Zeta[%s]",
     lambda r: (C(r, 1.5, 5, -3, 3),), lambda s: mp.zeta(s), k_zeta),
    ("Zeta C wide Re[0.5,10] Im[-40,40]", "Zeta[%s]",
     lambda r: (C(r, 0.5, 10, -40, 40),), lambda s: mp.zeta(s), k_zeta),
    ("Zeta C wide Re[-10,0.5] Im[-20,20]", "Zeta[%s]",
     lambda r: (C(r, -10, 0.5, -20, 20),), lambda s: mp.zeta(s), k_zeta),
    # --- HurwitzZeta
    ("HurwitzZeta MR bench s[1.5,6] a[0.5,4]", "HurwitzZeta[%s, %s]",
     lambda r: (U(r, 1.5, 6), U(r, 0.5, 4)), lambda s, a: mp.zeta(s, a), k_hurwitz),
    ("HurwitzZeta R wide s[0.5,40] a[0.01,50]", "HurwitzZeta[%s, %s]",
     lambda r: (U(r, 0.5, 40), math.exp(U(r, math.log(0.01), math.log(50)))),
     lambda s, a: mp.zeta(s, a), k_hurwitz),
    ("HurwitzZeta C bench", "HurwitzZeta[%s, %s]",
     lambda r: (C(r, 1.5, 4, -2, 2), C(r, 0.5, 3, -1, 1)),
     lambda s, a: mp.zeta(s, a), k_hurwitz),
    ("HurwitzZeta C wide", "HurwitzZeta[%s, %s]",
     lambda r: (C(r, 0.5, 10, -20, 20), C(r, 0.05, 10, -5, 5)),
     lambda s, a: mp.zeta(s, a), k_hurwitz),
    # --- LerchPhi
    ("LerchPhi MR bench z[0,.9] s[2,4] a[.5,2]", "LerchPhi[%s, %s, %s]",
     lambda r: (U(r, 0, .9), U(r, 2, 4), U(r, .5, 2)),
     lambda z, s, a: mp.lerchphi(z, s, a), None),
    ("LerchPhi R wide z[-.99,.99] s[-3,8] a[.1,10]", "LerchPhi[%s, %s, %s]",
     lambda r: (U(r, -.99, .99), U(r, -3, 8), U(r, .1, 10)),
     lambda z, s, a: mp.lerchphi(z, s, a), None),
    ("LerchPhi C bench z box .6, s[2,4] a[.5,2]", "LerchPhi[%s, %s, %s]",
     lambda r: (C(r, -.6, .6, -.6, .6), U(r, 2, 4), U(r, .5, 2)),
     lambda z, s, a: mp.lerchphi(z, s, a), None),
    ("LerchPhi C wide (complex z, s, a)", "LerchPhi[%s, %s, %s]",
     lambda r: (C(r, -.7, .7, -.7, .7), C(r, -2, 6, -3, 3), C(r, .2, 5, -2, 2)),
     lambda z, s, a: mp.lerchphi(z, s, a), None),
    # --- PolyLog
    ("PolyLog[2,x] MR bench x[-0.9,0.9]", "PolyLog[2, %s]",
     lambda r: (U(r, -.9, .9),), lambda x: mp.polylog(2, x),
     lambda x, v: k_polylog(2, x, v)),
    ("PolyLog[2,x] R wide x[-50,1]", "PolyLog[2, %s]",
     lambda r: (U(r, -50, 1),), lambda x: mp.polylog(2, x),
     lambda x, v: k_polylog(2, x, v)),
    ("PolyLog[2,x] R wide x[1,50] (cut)", "PolyLog[2, %s]",
     lambda r: (U(r, 1, 50),), lambda x: mp.polylog(2, mp.mpc(x, 0) - mp.mpc(0, mp.mpf(10) ** -40)),
     None),
    ("PolyLog[2,z] C bench box .8", "PolyLog[2, %s]",
     lambda r: (C(r, -.8, .8, -.8, .8),), lambda z: mp.polylog(2, z),
     lambda z, v: k_polylog(2, z, v)),
    ("PolyLog[2,z] C wide box 20", "PolyLog[2, %s]",
     lambda r: (C(r, -20, 20, -20, 20),), lambda z: mp.polylog(2, z),
     lambda z, v: k_polylog(2, z, v)),
    ("PolyLog[3,x] R x[-1,1]", "PolyLog[3, %s]",
     lambda r: (U(r, -1, 1),), lambda x: mp.polylog(3, x), None),
    ("PolyLog[s,x] R s[-3,8] x[-1,1]", "PolyLog[%s, %s]",
     lambda r: (U(r, -3, 8), U(r, -1, 1)), lambda s, x: mp.polylog(s, x), None),
    ("PolyLog[3,z] C wide box 20", "PolyLog[3, %s]",
     lambda r: (C(r, -20, 20, -20, 20),), lambda z: mp.polylog(3, z), None),
    ("PolyLog[s,z] C s[-2,6] box 2", "PolyLog[%s, %s]",
     lambda r: (U(r, -2, 6), C(r, -2, 2, -2, 2)), lambda s, z: mp.polylog(s, z), None),
]


def fmt_arg(x):
    return clit(x.real, x.imag) if isinstance(x, complex) else lit(x)


def to_mp(x):
    return mp.mpc(x.real, x.imag) if isinstance(x, complex) else mp.mpf(x)


def ulps(got, exact, cplx):
    if cplx or abs(exact.imag) > 0:
        ae = abs(exact)
        if ae == 0:
            return 0.0 if got == 0 else math.inf
        return float(abs(mp.mpc(got.real, got.imag) - exact) / (ae * mp.mpf(2) ** -52))
    e = mp.re(exact)
    u = math.ulp(float(e)) if e != 0 else 5e-324
    return float(abs(mp.mpf(got.real) - e) / u)


def reference(case_idx, pts):
    label, _, _, ref, cond = CASES[case_idx]
    mp.mp.dps = 90 if label.startswith(("Zeta", "HurwitzZeta")) else DPS
    out = []
    for p in pts:
        a = [to_mp(x) for x in p]
        try:
            v = ref(*a)
            k = 1.0
            if cond is not None:
                k = max(1.0, float(cond(*a, v)))
        except Exception:  # noqa: BLE001 -- a reference failure skips the point
            v, k = None, 1.0
        out.append((v, k))
    return out


def run_case(args):
    idx, n, binpath, seed = args
    label, tmpl, gen, _, _ = CASES[idx]
    rnd = random.Random(seed + idx)
    pts = [gen(rnd) for _ in range(n)]
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        for p in pts:
            call = tmpl % tuple(fmt_arg(x) for x in p)
            f.write("Print[ToString[InputForm[SetPrecision[%s, 17]]]];\n" % call)
        mfile = f.name
    try:
        out = subprocess.run([binpath, "-file", mfile], capture_output=True,
                             text=True, timeout=3600).stdout.splitlines()
    finally:
        os.unlink(mfile)
    refs = reference(idx, pts)
    errs, conds, unevaluated = [], [], 0
    worst = (0.0, None)
    for p, line, (v, k) in zip(pts, out, refs):
        got, cplx = parse(line)
        if got is None or v is None:
            unevaluated += 1
            continue
        u = ulps(got, v, cplx or any(isinstance(x, complex) for x in p))
        errs.append(u)
        conds.append(u / k)
        if u > worst[0]:
            worst = (u, p)
    return label, errs, conds, unevaluated, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=os.environ.get("MATHILDA_BIN", DEFAULT_BIN))
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--only", default="")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    heads = [h for h in a.only.split(",") if h]
    sel = [i for i, c in enumerate(CASES)
           if not heads or any(c[0].startswith(h + " ") or c[0].startswith(h + "[") for h in heads)]
    print("%-46s %6s %5s %7s %7s %8s %8s" % ("case", "n", "unev", "median", "p99", "max", "cond.max"))
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for label, errs, conds, unev, worst in ex.map(
                run_case, [(i, a.n, a.bin, a.seed) for i in sel]):
            if not errs:
                print("%-46s all unevaluated" % label)
                continue
            s = sorted(errs)
            print("%-46s %6d %5d %7.2f %7.2f %8.2f %8.2f   worst=%s" % (
                label, len(s), unev, statistics.median(s), s[int(0.99 * (len(s) - 1))],
                s[-1], max(conds), worst[1]))
            sys.stdout.flush()


if __name__ == "__main__":
    main()
