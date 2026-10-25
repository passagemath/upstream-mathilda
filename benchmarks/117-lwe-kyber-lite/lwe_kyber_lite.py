#!/usr/bin/env python3
"""Experiment 117 -- Module-LWE / Kyber-like encryption (group F).  See lwe_kyber_lite.m.

Baseline: numpy int64 arrays, which is how lattice crypto is prototyped in
Python (pure-Python lists would time the interpreter, not the arithmetic).
The negacyclic product is np.convolve folded back with a sign flip -- the same
linear-convolution algorithm as the .m half's ListConvolve.  The modular
matrix inverse / solve use sympy's DomainMatrix over GF(q) (dense Gauss-Jordan
in pure Python; python-flint is not installed, so sympy has no C backend here).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import small_list, CHECK_MOD

import numpy as np
from sympy import GF
from sympy.polys.matrices import DomainMatrix

require(["numpy:convolve", "numpy:dot", "sympy.polys.matrices:DomainMatrix"])

q = 3329
N = 256
HALF = 1665
CBD_W = np.array([1, 1, -1, -1], dtype=np.int64)


def cbd(seed, length):
    return np.array(small_list(seed, 4 * length, 2), dtype=np.int64).reshape(length, 4) @ CBD_W


def uniform(seed, length):
    return np.array(small_list(seed, length, q), dtype=np.int64)


def pmul(a, b):
    c = np.convolve(a, b)
    r = c[:N].copy()
    r[:N - 1] -= c[N:]
    return r % q


def to_bits(w):
    return ((4 * (w % q) + q) // (2 * q)) % 2


def fold_hash(arrs):
    f = np.concatenate([np.asarray(a, dtype=np.int64).ravel() for a in arrs])
    return int(np.dot(np.arange(1, len(f) + 1, dtype=np.int64), f)) % CHECK_MOD


res = {}

# ---- prim: the ring product -----------------------------------------------

pa = uniform(11, N); pb = cbd(12, N)


def polymul_loop():
    for _ in range(200):
        res["pc"] = pmul(pa, pb)


bench("prim: negacyclic polymul n=256 x200", polymul_loop)
check("prim: negacyclic polymul n=256 x200", fold_hash([res["pc"]]))

# ---- prim: plain LWE on a 512 x 512 matrix ---------------------------------

lweA = uniform(21, 512 * 512).reshape(512, 512)
lweS = cbd(22, 512); lweE = cbd(23, 512)


def matvec_loop():
    for _ in range(200):
        res["b"] = (lweA @ lweS + lweE) % q


bench("prim: LWE Mod[A.s+e,q] 512x512 x200", matvec_loop)
check("prim: LWE Mod[A.s+e,q] 512x512 x200", fold_hash([res["b"]]))

lweR = np.array(small_list(24, 256 * 512, 2), dtype=np.int64).reshape(256, 512)
bench("prim: LWE batch Mod[R.A,q] 256x512x512",
      lambda: res.__setitem__("U", (lweR @ lweA) % q))
check("prim: LWE batch Mod[R.A,q] 256x512x512", fold_hash([res["U"]]))

# ---- prim: modular matrix inverse / solve ----------------------------------

K = GF(q)
invRows = np.array(small_list(31, 128 * 128, q)).reshape(128, 128).tolist()
invB = small_list(32, 128, q)
dM = DomainMatrix([[K(x) for x in row] for row in invRows], (128, 128), K)
dB = DomainMatrix([[K(x)] for x in invB], (128, 1), K)


def inv_loop():
    for _ in range(4):
        res["Mi"] = dM.inv()


def solve_loop():
    for _ in range(16):
        res["X"] = dM.lu_solve(dB)


# sympy's GF(q) arithmetic is pure Python (~5 s per 128 x 128 inverse), so
# these two run once, without warm-up.
bench_once("prim: Inverse[M,Modulus->q] 128x128 x4", inv_loop)
check("prim: Inverse[M,Modulus->q] 128x128 x4",
      fold_hash([[int(x) % q for row in res["Mi"].to_list() for x in row]]))
bench_once("prim: LinearSolve[M,b,Modulus->q] 128x128 x16", solve_loop)
check("prim: LinearSolve[M,b,Modulus->q] 128x128 x16",
      fold_hash([[int(x) % q for row in res["X"].to_list() for x in row]]))

# ---- app: Kyber-lite --------------------------------------------------------


def kyber_data(k, n_msg):
    b = 1000 * k
    return {
        "A": [[uniform(b + 10 * i + j, N) for j in range(1, k + 1)] for i in range(1, k + 1)],
        "s": [cbd(b + 100 + i, N) for i in range(1, k + 1)],
        "e": [cbd(b + 150 + i, N) for i in range(1, k + 1)],
        "r": [[cbd(b + 200 + 10 * t + i, N) for i in range(1, k + 1)] for t in range(1, n_msg + 1)],
        "e1": [[cbd(b + 500 + 10 * t + i, N) for i in range(1, k + 1)] for t in range(1, n_msg + 1)],
        "e2": [cbd(b + 800 + t, N) for t in range(1, n_msg + 1)],
        "m": [np.array(small_list(b + 900 + t, N, 2), dtype=np.int64) for t in range(1, n_msg + 1)],
    }


def kyber_keygen(a, s, e):
    k = len(s)
    return [(sum(pmul(a[i][j], s[j]) for j in range(k)) + e[i]) % q for i in range(k)]


def kyber_encrypt(a, t, r, e1, e2, m):
    k = len(r)
    u = [(sum(pmul(a[j][i], r[j]) for j in range(k)) + e1[i]) % q for i in range(k)]
    v = (sum(pmul(t[j], r[j]) for j in range(k)) + e2 + HALF * m) % q
    return u + [v]


def kyber_decrypt(s, ct):
    k = len(s)
    return to_bits(ct[k] - sum(pmul(s[j], ct[j]) for j in range(k)))


def kyber_run(d):
    t = kyber_keygen(d["A"], d["s"], d["e"])
    cts = [kyber_encrypt(d["A"], t, d["r"][j], d["e1"][j], d["e2"][j], d["m"][j])
           for j in range(len(d["m"]))]
    msg = [kyber_decrypt(d["s"], ct) for ct in cts]
    correct = sum(int(np.sum(mm == m0)) for mm, m0 in zip(msg, d["m"]))
    return correct * 10**10 + fold_hash([c for ct in cts for c in ct])


ky2 = kyber_data(2, 8)
bench("app: Kyber-lite k=2 keygen+enc+dec x8", lambda: res.__setitem__("k2", kyber_run(ky2)))
check("app: Kyber-lite k=2 keygen+enc+dec x8", res["k2"])

ky3 = kyber_data(3, 8)
bench("app: Kyber-lite k=3 keygen+enc+dec x8", lambda: res.__setitem__("k3", kyber_run(ky3)))
check("app: Kyber-lite k=3 keygen+enc+dec x8", res["k3"])

# Python has one spelling of the ring-vector sum; this row is the same program,
# timed under the label of the .m half's MapThread variant.
bench("app: Kyber-lite k=2 MapThread spelling x8", lambda: res.__setitem__("k2m", kyber_run(ky2)))
check("app: Kyber-lite k=2 MapThread spelling x8", res["k2m"])

# ---- app: Regev LWE, one bit per 512-vector ---------------------------------

lweMsg = np.array(small_list(25, 256, 2), dtype=np.int64)


def lwe_run():
    b = (lweA @ lweS + lweE) % q
    cts = []
    for t in range(256):
        r = lweR[t]
        cts.append(((r @ lweA) % q, int((r @ b + HALF * lweMsg[t]) % q)))
    dec = np.array([int(to_bits(v - int(lweS @ u))) for u, v in cts])
    correct = int(np.sum(dec == lweMsg))
    flat = []
    for u, v in cts:
        flat.append(u)
        flat.append(np.array([v], dtype=np.int64))
    return correct * 10**10 + fold_hash(flat)


bench("app: Regev LWE 512x512 enc+dec 256 bits", lambda: res.__setitem__("l", lwe_run()))
check("app: Regev LWE 512x512 enc+dec 256 bits", res["l"])
