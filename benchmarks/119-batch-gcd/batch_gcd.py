#!/usr/bin/env python3
"""Experiment 119 -- batch GCD over RSA moduli (group F).  See batch_gcd.m.

Baseline: the "app:" cases are plain Python ints (math.gcd, built-in big-int
multiply / divmod), as a pure-Python batch-GCD scan is written.  The "prim:"
cases use gmpy2 (GMP): the M-fold product as functools.reduce(gmpy2.mul) over
mpz, and gmpy2.gcd for the final GCD vector.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import prime_of, CHECK_MOD

import gmpy2
from functools import reduce
from math import gcd

require(["gmpy2:mul", "gmpy2:gcd", "math:gcd", "math:prod"])

M = 1024
bgP = [prime_of(5000 + 2 * i, 512) for i in range(1, M + 1)]
bgQ = [bgP[i - 18] if i % 37 == 0 else prime_of(5001 + 2 * i, 512) for i in range(1, M + 1)]
moduli = [p * q for p, q in zip(bgP, bgQ)]
moduli256 = moduli[:256]

# ---- the algorithm, written in the language ---------------------------------


def product_tree(ns):
    tree = [ns]
    while len(ns) > 1:
        ns = [ns[i] * ns[i + 1] if i + 1 < len(ns) else ns[i] for i in range(0, len(ns), 2)]
        tree.append(ns)
    return tree


def remainder_tree(tree):
    r = tree[-1]
    for level in reversed(tree[:-1]):
        r = [r[j // 2] % (level[j] * level[j]) for j in range(len(level))]
    return r


def batch_gcd(ns):
    z = remainder_tree(product_tree(ns))
    return [gcd(n, zi // n) for n, zi in zip(ns, z)]


def weak_digest(g):
    w = [x for x in g if x > 1]
    return len(w) * 10**10 + int(sum(w)) % CHECK_MOD


res = {}

# ---- prim: the big-integer kernels ------------------------------------------

mpz_moduli = [gmpy2.mpz(n) for n in moduli]
bench("prim: Times@@ 1024 moduli (1M-bit product) x4",
      lambda: [res.__setitem__("prod", reduce(gmpy2.mul, mpz_moduli)) for _ in range(4)])
check("prim: Times@@ 1024 moduli (1M-bit product) x4", int(res["prod"] % CHECK_MOD))

tree1024 = product_tree(moduli)
z1024 = remainder_tree(tree1024)
zq1024 = [gmpy2.mpz(z // n) for z, n in zip(z1024, moduli)]
bench("prim: GCD 1024-bit pairs x1024x8",
      lambda: [res.__setitem__("gf", [gmpy2.gcd(n, z) for n, z in zip(mpz_moduli, zq1024)]) for _ in range(8)])
check("prim: GCD 1024-bit pairs x1024x8", weak_digest(res["gf"]))

# ---- app: the two tree phases, then the whole scan ---------------------------

bench("app: product tree M=1024 x4", lambda: [res.__setitem__("tr", product_tree(moduli)) for _ in range(4)])
check("app: product tree M=1024 x4", res["tr"][-1][0] % CHECK_MOD + 10**10 * len(res["tr"]))

bench("app: remainder tree mod N^2 M=1024", lambda: res.__setitem__("zz", remainder_tree(tree1024)))
check("app: remainder tree mod N^2 M=1024", sum(res["zz"]) % CHECK_MOD)

bench("app: batch GCD M=256 x4", lambda: [res.__setitem__("g256", batch_gcd(moduli256)) for _ in range(4)])
check("app: batch GCD M=256 x4", weak_digest(res["g256"]))

bench("app: batch GCD M=1024", lambda: res.__setitem__("g1024", batch_gcd(moduli)))
check("app: batch GCD M=1024", weak_digest(res["g1024"]))


def naive_gcd(ns):
    m = len(ns)
    g = [1] * m
    for i in range(m - 1):
        for j in range(i + 1, m):
            d = gcd(ns[i], ns[j])
            if d > 1:
                g[i] *= d
                g[j] *= d
    return g


bench_once("app: naive pairwise GCD M=256 (32640 GCDs)", lambda: res.__setitem__("gn", naive_gcd(moduli256)))
check("app: naive pairwise GCD M=256 (32640 GCDs)", weak_digest(res["gn"]))
