#!/usr/bin/env python3
"""Accuracy gate for the Gamma-family machine kernels (src/special_functions/sf_machine_gamma.c).

Drives the Mathilda binary (`-file`) over >= 2000 random machine-precision points
per head per regime, and compares against mpmath at 30 digits.

    benchmarks/.venv/bin/python -I tools/sf_accuracy/gamma.py [--bin ./Mathilda] [--n 2000]

Reported per (head, regime):
  rel   median / p99 / max of |f - ref| / |ref|, in units of 2^-53 ("ulp").
  scl   max of |f - ref| / max(|ref|, |z f'(z)|) in the same units: near a zero of
        f this is the absolute error scaled by the local magnitude (the size of
        the input-rounding error, |z f'| eps), and elsewhere it equals rel.
  dec   points where Mathilda did not return a number (should be 0 in-domain).

The benchmark domains are those of benchmarks/97-special-functions/registry.py;
the wider bands probe the regime boundaries of the kernels.
"""
import argparse
import os
import random
import statistics
import subprocess
import sys
import tempfile

import mpmath as mp

mp.mp.dps = 30
U = mp.mpf(2) ** -53


def lit(x):
    return repr(float(x)) + "`"


def clit(z):
    return "Complex[%s, %s]" % (lit(z.real), lit(z.imag))


def arg_lit(a):
    if isinstance(a, complex):
        return clit(a)
    if isinstance(a, int):
        return str(a)
    return lit(a)


# head -> (WL template, mpmath value, mpmath derivative wrt the varying first argument or None)
HEADS = {
    "Gamma":      ("Gamma[{0}]", lambda z: mp.gamma(z), lambda z: mp.gamma(z) * mp.digamma(z)),
    "LogGamma":   ("LogGamma[{0}]", lambda z: mp.loggamma(z), lambda z: mp.digamma(z)),
    "PolyGamma0": ("PolyGamma[0, {0}]", lambda z: mp.digamma(z), lambda z: mp.polygamma(1, z)),
    "PolyGamma1arg": ("PolyGamma[{0}]", lambda z: mp.digamma(z), lambda z: mp.polygamma(1, z)),
    "PolyGamma2": ("PolyGamma[2, {0}]", lambda z: mp.polygamma(2, z), lambda z: mp.polygamma(3, z)),
    "PolyGamma1": ("PolyGamma[1, {0}]", lambda z: mp.polygamma(1, z), lambda z: mp.polygamma(2, z)),
    "Beta":       ("Beta[{0}, {1}]", lambda a, b: mp.beta(a, b), None),
    "Pochhammer3": ("Pochhammer[{0}, 3]", lambda a: mp.rf(a, 3), None),
    "Pochhammer": ("Pochhammer[{0}, {1}]", lambda a, n: mp.rf(a, n), None),
}

R = random.uniform


def rr(lo, hi):
    return lambda: R(lo, hi)


def cc(xlo, xhi, ylo, yhi):
    return lambda: complex(R(xlo, xhi), R(ylo, yhi))


def ii(lo, hi):
    return lambda: float(random.randint(lo, hi))


# (head, regime label, [generator per slot], benchmark-domain?)
REGIMES = [
    ("Gamma", "MR [0.3,6] (libm tgamma)", [rr(0.3, 6)], True),
    ("Gamma", "C box", [cc(0.5, 4, -3, 3)], True),
    ("Gamma", "C wide", [cc(-10, 30, -20, 20)], False),
    ("LogGamma", "MR [0.3,8] (libm lgamma)", [rr(0.3, 8)], True),
    ("LogGamma", "C box", [cc(0.5, 4, -3, 3)], True),
    ("LogGamma", "C wide", [cc(-25, 30, -20, 20)], False),
    ("PolyGamma0", "MR [0.3,8]", [rr(0.3, 8)], True),
    ("PolyGamma0", "R wide [1e-3,40]", [rr(1e-3, 40)], False),
    ("PolyGamma0", "R negative [-20,0]", [rr(-20, 0)], False),
    ("PolyGamma0", "C box", [cc(0.5, 4, -3, 3)], True),
    ("PolyGamma0", "C wide", [cc(-10, 30, -20, 20)], False),
    ("PolyGamma1arg", "MR [0.3,8]", [rr(0.3, 8)], True),
    ("PolyGamma1arg", "C box", [cc(0.5, 4, -3, 3)], True),
    ("PolyGamma2", "MR [0.3,8]", [rr(0.3, 8)], True),
    ("PolyGamma2", "R wide [1e-3,50]", [rr(1e-3, 50)], False),
    ("PolyGamma2", "C box", [cc(0.5, 4, -3, 3)], True),
    ("PolyGamma2", "C wide (Re >= 0)", [cc(0, 30, -20, 20)], False),
    ("PolyGamma1", "R wide [1e-3,50]", [rr(1e-3, 50)], False),
    ("Beta", "MR [0.3,5]^2", [rr(0.3, 5), rr(0.3, 5)], True),
    ("Beta", "R wide [-10,30]^2", [rr(-10, 30), rr(-10, 30)], False),
    ("Beta", "C box", [cc(0.5, 3, -2, 2), cc(0.5, 3, -2, 2)], True),
    ("Beta", "C wide", [cc(-5, 10, -5, 5), cc(-5, 10, -5, 5)], False),
    ("Pochhammer3", "MR [0.3,5]", [rr(0.3, 5)], True),
    ("Pochhammer3", "C box", [cc(0.5, 3, -2, 2)], True),
    ("Pochhammer", "R wide, real n", [rr(-20, 20), rr(-20, 20)], False),
    ("Pochhammer", "R wide, integer n", [rr(-20, 20), ii(-30, 30)], False),
    ("Pochhammer", "C wide, complex n", [cc(-5, 10, -5, 5), cc(-5, 5, -3, 3)], False),
]


def run_mathilda(binary, exprs):
    lines = []
    for e in exprs:
        lines.append('r = Quiet[%s]; Print[If[NumberQ[r], ToString[InputForm[SetPrecision[Re[r], 17]]] '
                     '<> " " <> ToString[InputForm[SetPrecision[Im[r], 17]]], "X"]];' % e)
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        f.write("\n".join(lines) + "\n")
        path = f.name
    try:
        out = subprocess.run([binary, "-file", path], capture_output=True, text=True, timeout=3600)
    finally:
        os.unlink(path)
    res = out.stdout.strip().split("\n")
    if len(res) != len(exprs):
        sys.exit("mathilda produced %d lines for %d inputs\n%s" % (len(res), len(exprs), out.stderr[-2000:]))
    return res


def parse_num(s):
    return float(s.replace("*^", "e"))


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--bin", default=os.path.join(here, "..", "..", "Mathilda"))
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    random.seed(args.seed)
    only = set(args.only.split(",")) if args.only else None
    print("%-14s %-26s %5s %4s  %8s %8s %8s   %8s" % ("head", "regime", "n", "dec", "rel-med", "rel-p99", "rel-max", "scl-max"))
    worst_bench = 0.0
    for head, label, gens, bench in REGIMES:
        if only and head not in only:
            continue
        tmpl, fval, fder = HEADS[head]
        pts = [[g() for g in gens] for _ in range(args.n)]
        exprs = [tmpl.format(*[arg_lit(a) for a in p]) for p in pts]
        outs = run_mathilda(args.bin, exprs)
        rel, scl, dec = [], [], 0
        worst = None
        for p, o in zip(pts, outs):
            if o == "X":
                dec += 1
                continue
            re_s, im_s = o.split()
            v = mp.mpc(parse_num(re_s), parse_num(im_s))
            mpargs = [mp.mpc(a) if isinstance(a, complex) else mp.mpf(a) for a in p]
            try:
                t = fval(*mpargs)
            except Exception:
                dec += 1
                continue
            err = abs(v - t)
            e = float(err / abs(t) / U) if t != 0 else (0.0 if err == 0 else float("inf"))
            rel.append(e)
            if fder is not None:
                loc = max(abs(t), abs(mpargs[0] * fder(mpargs[0])))
                scl.append(float(err / loc / U))
            else:
                scl.append(e)
            if worst is None or e > worst[0]:
                worst = (e, p)
        rel.sort()
        if not rel:
            print("%-14s %-26s %5d %4d  (no values)" % (head, label, len(pts), dec))
            continue
        line = "%-14s %-26s %5d %4d  %8.2f %8.2f %8.2f   %8.2f" % (
            head, label, len(rel), dec, statistics.median(rel), rel[int(0.99 * (len(rel) - 1))], rel[-1], max(scl))
        print(line + ("" if not bench else "   [benchmark domain]"))
        if bench:
            worst_bench = max(worst_bench, rel[-1])
        sys.stdout.flush()
    print("\nworst max rel error on the benchmark domains: %.2f ulp" % worst_bench)


if __name__ == "__main__":
    main()
