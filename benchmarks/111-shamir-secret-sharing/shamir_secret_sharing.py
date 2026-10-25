#!/usr/bin/env python3
"""Experiment 111 -- Shamir secret sharing (group F).  See shamir_secret_sharing.m.

Baseline: plain Python ints + built-in pow(d, -1, p) for the "app:" protocol
(idiomatic pure-Python secret sharing).  The "prim:" interpolation is
sympy's GF(p) linear algebra: the Vandermonde system solved by
DomainMatrix.lu_solve over GF(p), whose first unknown is f(0).  sympy has no
GF(p) interpolation routine, and its Q-route (sympy.polys.polyfuncs.interpolate,
the analogue of the WL case) is unusable at this size: 1.7 s at 30 points,
minutes at 50 -- see the README.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, prime_of, CHECK_MOD

from math import prod
from sympy import GF
from sympy.polys.matrices import DomainMatrix

require(["builtins:pow", "sympy.polys.matrices:DomainMatrix"])

P = prime_of(1111, 256)
T = 50
N = 200


def poly_of(k):
    return [big_of(100000 * k + j, 256) % P for j in range(T)]


def horner_mod(cs, x):
    acc = 0
    for c in reversed(cs):
        acc = (acc * x + c) % P
    return acc


def split(cs):
    return [(i, horner_mod(cs, i)) for i in range(1, N + 1)]


def lagrange0(sh):
    xs = [x for x, _ in sh]
    total = 0
    for i in range(len(xs)):
        num, den = 1, 1
        for j in range(len(xs)):
            if j != i:
                num = (num * xs[j]) % P
                den = (den * (xs[j] - xs[i])) % P
        total += (sh[i][1] * num * pow(den, -1, P)) % P
    return total % P


def pick(sh):
    return [sh[k - 1] for k in range(3, 3 + 4 * (T - 1) + 1, 4)]


polys = [poly_of(k) for k in range(1, 21)]
res = {}

# ---- app ----------------------------------------------------------------

bench("app: split 256-bit n=200 t=50 x20",
      lambda: res.__setitem__("shares", [split(c) for c in polys]))
shares = res["shares"]
check("app: split 256-bit n=200 t=50 x20", sum(y for sh in shares for _, y in sh) % CHECK_MOD)

bench("app: reconstruct Lagrange t=50 x20",
      lambda: res.__setitem__("recs", [lagrange0(pick(sh)) for sh in shares]))
check("app: reconstruct Lagrange t=50 x20", sum(res["recs"]) % CHECK_MOD)

batch_shares = [[(x, (y + k) % P) for x, y in pick(shares[k % 20])] for k in range(1, 2001)]


def batch():
    xs0 = [x for x, _ in batch_shares[0]]
    w = [(prod(xs0[j] if j != i else 1 for j in range(T))
          * pow(prod(xs0[j] - xs0[i] if j != i else 1 for j in range(T)), -1, P)) % P
         for i in range(T)]
    res["brec"] = [sum(y * wi for (_, y), wi in zip(sh, w)) % P for sh in batch_shares]


bench("app: batch reconstruct t=50 x2000", batch)
check("app: batch reconstruct t=50 x2000", sum(res["brec"]) % CHECK_MOD)

# ---- prim ---------------------------------------------------------------


K = GF(P)


def interp0(sh):
    """f(0) of the interpolant through sh, by a Vandermonde solve over GF(p)."""
    n = len(sh)
    V = DomainMatrix([[K(pow(x, j, P)) for j in range(n)] for x, _ in sh], (n, n), K)
    b = DomainMatrix([[K(y)] for _, y in sh], (n, 1), K)
    return int(V.lu_solve(b).to_Matrix()[0, 0]) % P


bench("prim: InterpolatingPolynomial t=50 mod p x20",
      lambda: res.__setitem__("irecs", [interp0(pick(sh)) for sh in shares]))
check("prim: InterpolatingPolynomial t=50 mod p x20", sum(res["irecs"]) % CHECK_MOD)

bench_once("prim: InterpolatingPolynomial n=200 mod p x4",
           lambda: [res.__setitem__("i200", interp0(shares[r])) for r in range(4)])
check("prim: InterpolatingPolynomial n=200 mod p x4", res["i200"] % CHECK_MOD)
