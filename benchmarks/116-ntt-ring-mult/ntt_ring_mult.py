#!/usr/bin/env python3
"""Experiment 116 -- NTT ring multiplication (group F).  See ntt_ring_mult.m.

"app:" a pure-Python iterative radix-2 Cooley-Tukey NTT on lists of ints
(bit-reversed input, precomputed twiddles, psi pre-/post-multiplication for
the negacyclic wrap) -- the same butterflies as the .m half's vectorised
stages.  "prim:" numpy.convolve on int64 (n q^2 < 2^63 for both rings) folded
by x^n = -1, and sympy integer polynomials: (A*B).rem(x^n + 1) then mod q.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import small_list, CHECK_MOD

import numpy as np
from sympy import Poly, symbols
from sympy.ntheory import primitive_root

require(["numpy:convolve", "sympy:Poly", "sympy.ntheory:primitive_root"])

X = symbols("x")


def ntt_setup(q, n):
    g = primitive_root(q)
    psi = pow(g, (q - 1) // (2 * n), q)
    psi_inv = pow(psi, -1, q)
    w = psi * psi % q
    w_inv = pow(w, -1, q)
    k = n.bit_length() - 1
    perm = [int(format(i, "0%db" % k)[::-1], 2) for i in range(n)]
    n_inv = pow(n, -1, q)
    return {"q": q, "n": n, "perm": perm,
            "tw": [pow(w, j, q) for j in range(n // 2)],
            "twInv": [pow(w_inv, j, q) for j in range(n // 2)],
            "psiPow": [pow(psi, j, q) for j in range(n)],
            "postPow": [n_inv * pow(psi_inv, j, q) % q for j in range(n)]}


def ntt(a, tw, ctx):
    q, n = ctx["q"], ctx["n"]
    x = [a[i] for i in ctx["perm"]]
    h = 1
    while h < n:
        step = n // (2 * h)
        for s in range(0, n, 2 * h):
            for j in range(h):
                u = x[s + j]
                v = x[s + j + h] * tw[j * step] % q
                x[s + j] = (u + v) % q
                x[s + j + h] = (u - v) % q
        h *= 2
    return x


def ntt_mul(a, b, ctx):
    q = ctx["q"]
    fa = ntt([x * p % q for x, p in zip(a, ctx["psiPow"])], ctx["tw"], ctx)
    fb = ntt([x * p % q for x, p in zip(b, ctx["psiPow"])], ctx["tw"], ctx)
    c = ntt([x * y % q for x, y in zip(fa, fb)], ctx["twInv"], ctx)
    return [x * p % q for x, p in zip(c, ctx["postPow"])]


def conv_mul(a, b, q, n):
    full = np.convolve(a, b)
    lo = full[:n].copy()
    lo[:n - 1] -= full[n:]
    return (lo % q).tolist()


def poly_mul(a, b, q, n):
    pa = Poly(list(reversed(a)), X)
    pb = Poly(list(reversed(b)), X)
    r = (pa * pb).rem(Poly(X ** n + 1, X))
    cs = [int(c) % q for c in reversed(r.all_coeffs())]
    return cs + [0] * (n - len(cs))


def digest(cs):
    return sum(sum((i + 1) * c for i, c in enumerate(v)) for v in cs) % CHECK_MOD


ctxA = ntt_setup(12289, 1024)
ctxB = ntt_setup(7681, 256)
pairsA = [(small_list(1600 + 2 * k, 1024, 12289), small_list(1601 + 2 * k, 1024, 12289)) for k in range(1, 21)]
pairsB = [(small_list(1700 + 2 * k, 256, 7681), small_list(1701 + 2 * k, 256, 7681)) for k in range(1, 101)]
npA = [(np.array(a, dtype=np.int64), np.array(b, dtype=np.int64)) for a, b in pairsA]
npB = [(np.array(a, dtype=np.int64), np.array(b, dtype=np.int64)) for a, b in pairsB]
res = {}

bench("app: NTT negacyclic mult q=12289 n=1024 x20",
      lambda: res.__setitem__("A1", [ntt_mul(a, b, ctxA) for a, b in pairsA]))
check("app: NTT negacyclic mult q=12289 n=1024 x20", digest(res["A1"]))

bench("prim: ListConvolve negacyclic q=12289 n=1024 x4",
      lambda: res.__setitem__("A2", [conv_mul(a, b, 12289, 1024) for a, b in npA[:4]]))
check("prim: ListConvolve negacyclic q=12289 n=1024 x4", digest(res["A2"]))
assert res["A2"] == res["A1"][:4]

bench("app: NTT negacyclic mult q=7681 n=256 x100",
      lambda: res.__setitem__("B1", [ntt_mul(a, b, ctxB) for a, b in pairsB]))
check("app: NTT negacyclic mult q=7681 n=256 x100", digest(res["B1"]))

bench("prim: ListConvolve negacyclic q=7681 n=256 x100",
      lambda: res.__setitem__("B2", [conv_mul(a, b, 7681, 256) for a, b in npB]))
check("prim: ListConvolve negacyclic q=7681 n=256 x100", digest(res["B2"]))

pairsC = [(small_list(1900 + 2 * k, 128, 7681), small_list(1901 + 2 * k, 128, 7681)) for k in range(1, 3)]
bench_once("prim: PolynomialRemainder negacyclic q=7681 n=128 x2",
           lambda: res.__setitem__("C3", [poly_mul(a, b, 7681, 128) for a, b in pairsC]))
check("prim: PolynomialRemainder negacyclic q=7681 n=128 x2", digest(res["C3"]))
assert res["C3"] == [conv_mul(np.array(a), np.array(b), 7681, 128) for a, b in pairsC]
