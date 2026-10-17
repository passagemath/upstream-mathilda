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
    "BlakeRationalBaseDescent": [16, 24, 40],
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


def planted(rng, qbits, pbits, lattice=False):
    """A semiprime N = P q with P ~ c (a/b)^k: [RBD] Theorem 1's hypotheses
    hold (|Delta| < (a/b)^k / q), or, with lattice=True, the miss at the
    aligned depth is ~2^(qbits/2), far outside the gcd window but inside the
    Howgrave-Graham tolerance.  Returns (N, P, q, a, b)."""
    from fractions import Fraction
    from math import gcd as g
    while True:
        b = rng.randrange(1, 12)
        a = rng.randrange(b + 1, 40)
        if g(a, b) == 1:
            break
    c = rng.randrange(1, 1 << 16)
    q = randprime(1 << (qbits - 1), 1 << qbits)
    r = Fraction(a, b)
    k = 1
    while (c * r ** k).numerator // (c * r ** k).denominator < (1 << pbits):
        k += 1
    base = (c * a ** k) // (b ** k)
    tol = (a ** k) // (b ** k * q)              # (a/b)^k / q
    while True:
        if lattice:
            miss = rng.randrange(1 << (qbits // 3), 1 << (qbits // 3 + 1))
            delta = (a ** k * miss) // (b ** k * q)    # miss = Delta (b/a)^k q
        else:
            delta = rng.randrange(-(tol // 2), tol // 2 + 1) if tol > 2 else 0
        P = base + delta
        from sympy import nextprime
        P = nextprime(P)
        if lattice or abs(P - base) < tol:
            return P * q, P, q, a, b


def run_rbd_planted(count, rng):
    """Each variant must return exactly {q, P} on planted instances."""
    cases = []
    for _ in range(count):
        qb = rng.choice([40, 64, 100])
        # P must exceed c q by a margin, or the theorem's tolerance
        # (a/b)^k / q is too small to hold a prime (c < 2^16).
        cases.append(("exact", planted(rng, qb, qb + 16 + rng.choice([40, 100, 200]))))
    for _ in range(max(1, count // 4)):
        cases.append(("lattice", planted(rng, 98, 99, lattice=True)))
    calls = []
    for kind, (N, P, q, a, b) in cases:
        if kind == "exact":
            calls += [f'FactorInteger`BlakeRationalBaseDescent[{N}, "Base" -> {a}/{b}]',
                      f'FactorInteger`BlakeRationalBaseDescentClassic[{N}, "Base" -> {a}/{b}]',
                      f'FactorInteger`BlakeRationalBaseDescent[{N}, "Height" -> 40]']
        else:
            calls += [f'FactorInteger`BlakeRationalBaseDescentLattice[{N}, "Base" -> {a}/{b}]']
    script = "".join(f'r = Quiet[{c}]; Print["R ", StringRiffle[ToString /@ Flatten[r], " "]];\n' for c in calls)
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False, dir=os.environ.get("TMPDIR")) as f:
        f.write(script)
        path = f.name
    try:
        p = subprocess.run([MATHILDA, "-file", path], capture_output=True, text=True, timeout=3600)
    finally:
        os.unlink(path)
    rows = [l[2:] for l in p.stdout.splitlines() if l.startswith("R ")]
    bad, i = 0, 0
    for kind, (N, P, q, a, b) in cases:
        for _ in range(3 if kind == "exact" else 1):
            row = rows[i] if i < len(rows) else ""
            got = sorted(int(x) for x in row.split()[0::2]) if row else []
            if got != sorted([P, q]):
                bad += 1
                print(f"  RBD {kind} call {calls[i][:60]}... N={N} base={a}/{b}: got {got}")
            i += 1
    print(f"BlakeRationalBaseDescent planted  cases={len(cases)} calls={len(calls)} wrong={bad}", flush=True)
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", type=int, default=100)
    ap.add_argument("--methods", default="")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rbd-planted", type=int, default=0,
                    help="also run N planted rational-base instances through every RBD variant")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    want = a.methods.split(",") if a.methods else list(METHODS)
    bad_total = 0
    if a.rbd_planted:
        bad_total += run_rbd_planted(a.rbd_planted, rng)
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
