#!/usr/bin/env python3.11
"""Randomised correctness sweep for every FactorInteger method.

For each method and size class, factor a batch of random inputs in Mathilda
(one process per (method, class), each input independently checked) and
verify every answer: the product of the returned prime powers must equal n,
and every base must pass sympy.isprime unless Mathilda flagged it with
FactorInteger::nofac (a bounded method is allowed to give up, never to lie).

    python3.11 verify_methods.py [--count 200] [--methods ECM,Fermat]
"""
import argparse
import os
import random
import re
import subprocess
import tempfile
from math import prod

from sympy import isprime, randprime

HERE = os.path.dirname(os.path.abspath(__file__))
MATHILDA = os.environ.get("MATHILDA", os.path.join(HERE, "..", "..", "Mathilda"))

# method -> size classes (bits of n) it is expected to handle quickly.
METHODS = {
    "Automatic":         [8, 20, 32, 48, 63, 64, 65, 80, 100, 128, 160, 200],
    "TrialDivision":     [8, 20, 32, 40],
    "PollardRho":        [20, 32, 48, 63, 64, 80, 100],
    "Fermat":            [20, 32, 40],
    "Lehman":            [20, 32, 40, 42],
    "ShanksSquareForms": [20, 32, 48, 62, 80],
    "CFRAC":             [40, 64, 80, 100],
    "Dixon":             [40, 50, 64],
    "QuadraticSieve":    [40, 64, 80, 100, 120],
    "PollardPMinus1":    [20, 40, 64, 100],
    "WilliamsPPlus1":    [20, 40, 64, 100],
    "ECM":               [40, 64, 100, 128],
    "BlakeRationalBaseDescent": [16, 24],
}


def gen(bits, rng):
    """Mix of shapes: random, semiprime, prime, prime power, smooth."""
    kind = rng.randrange(5)
    if kind == 0 or bits < 12:
        return rng.randrange(1 << (bits - 1), 1 << bits)
    if kind == 1:
        a = bits // 2
        return randprime(1 << (a - 1), 1 << a) * randprime(1 << (bits - a - 1), 1 << (bits - a))
    if kind == 2:
        return randprime(1 << (bits - 1), 1 << bits)
    if kind == 3:
        e = rng.choice([2, 3, 5])
        b = max(2, bits // e)
        return randprime(1 << (b - 1), 1 << b) ** e
    m = 1
    while m.bit_length() < bits:
        m *= rng.choice([2, 3, 5, 7, 11, 13, 101, 997, 65537])
    return m


def run(method, ns):
    opt = "" if method == "Automatic" else f', Method -> "{method}"'
    lst = ", ".join(str(n) for n in ns)
    script = (f"ns = {{{lst}}};\n"
              f"Do[r = FactorInteger[n{opt}]; Print[\"R \", StringRiffle[ToString /@ Flatten[r], \" \"]], {{n, ns}}];\n")
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False, dir=os.environ.get("TMPDIR")) as f:
        f.write(script)
        path = f.name
    try:
        p = subprocess.run([MATHILDA, "-file", path], capture_output=True, text=True, timeout=1800)
    finally:
        os.unlink(path)
    rows = [l[2:] for l in p.stdout.splitlines() if l.startswith("R ")]
    nofac = set(re.findall(r"FactorInteger::nofac: (\d+) is composite", p.stderr))
    return rows, nofac, p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", type=int, default=100)
    ap.add_argument("--methods", default="")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    want = a.methods.split(",") if a.methods else list(METHODS)
    bad_total = 0
    for m in want:
        for bits in METHODS[m]:
            ns = [gen(bits, rng) for _ in range(a.count)]
            rows, nofac, rc = run(m, ns)
            bad, gaveup = 0, 0
            if len(rows) != len(ns):
                print(f"{m:26s} {bits:4d}b  CRASH/SHORT rc={rc} rows={len(rows)}/{len(ns)}")
                bad_total += 1
                continue
            for n, row in zip(ns, rows):
                t = [int(x) for x in row.split()]
                pairs = list(zip(t[0::2], t[1::2]))
                if prod(p ** e for p, e in pairs) != n:
                    bad += 1
                    print(f"  WRONG PRODUCT {m} n={n} -> {pairs}")
                    continue
                for p, _ in pairs:
                    if p > 1 and not isprime(p):
                        if str(p) in nofac:
                            gaveup += 1
                        else:
                            bad += 1
                            print(f"  SILENT COMPOSITE {m} n={n} p={p}")
            bad_total += bad
            print(f"{m:26s} {bits:4d}b  n={len(ns)} wrong={bad} gave_up={gaveup}", flush=True)
    print("TOTAL WRONG:", bad_total)
    raise SystemExit(1 if bad_total else 0)


if __name__ == "__main__":
    main()
