#!/usr/bin/env python3
"""Experiment 118 -- quadratic-residuosity cryptography (group F).  See qr_prngs.m.

Baseline: the "app:" cases are pure Python ints with built-in pow() (how a
Python BBS / Goldwasser-Micali is written); the "prim:" symbol cases use
gmpy2.jacobi / gmpy2.kronecker (GMP's mpz_jacobi / mpz_kronecker).  The
Blum primes are found with gmpy2.is_prime(c, 25), the same deterministic
search as the .m half.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, small_list, CHECK_MOD

import gmpy2

require(["gmpy2:jacobi", "gmpy2:kronecker", "gmpy2:is_prime"])


def blum_prime(start):
    c = start - start % 4 + 3
    while not gmpy2.is_prime(c, 25):
        c += 4
    return c


P = blum_prime(big_of(1181, 512))
Q = blum_prime(big_of(1182, 512))
N = P * Q
X0 = big_of(1183, 1000) ** 2 % N

res = {}

# ---- app: Blum-Blum-Shub ---------------------------------------------------


def bbs_parity(x, n, steps):
    h = 0
    for _ in range(steps):
        x = x * x % n
        h = (2 * h + x % 2) % CHECK_MOD
    return h


def bbs_low10(x, n, steps):
    h = 0
    for _ in range(steps):
        x = x * x % n
        h = (1024 * h + x % 1024) % CHECK_MOD
    return h


bench("app: BBS-1024 parity x20000 bits", lambda: res.__setitem__("b1", bbs_parity(X0, N, 20000)))
check("app: BBS-1024 parity x20000 bits", res["b1"])

bench("app: BBS-1024 10 bits/step x20000 steps", lambda: res.__setitem__("b2", bbs_low10(X0, N, 20000)))
check("app: BBS-1024 10 bits/step x20000 steps", res["b2"])

# ---- app / prim: Goldwasser-Micali -----------------------------------------

Y = N - 1
msg = small_list(1184, 2000, 2)
rnd = [big_of(2000 + i, 1000) for i in range(1, 2001)]


def fold_bits(bits):
    return int("".join(map(str, bits)), 2) % CHECK_MOD


bench("app: GM-1024 encrypt x2000 bits x10",
      lambda: [res.__setitem__("c", [x * x * Y ** b % N for x, b in zip(rnd, msg)]) for _ in range(10)])
check("app: GM-1024 encrypt x2000 bits x10", sum(res["c"]) % CHECK_MOD)

gmC = res["c"]
bench("prim: JacobiSymbol GM-1024 decrypt x2000 bits x4",
      lambda: [res.__setitem__("d1", [(1 - gmpy2.jacobi(c, P)) // 2 for c in gmC]) for _ in range(4)])
check("prim: JacobiSymbol GM-1024 decrypt x2000 bits x4",
      fold_bits(res["d1"]) + (10**10 if res["d1"] == msg else 0))

E = (P - 1) // 2
bench("app: GM-1024 decrypt Euler criterion x2000 bits",
      lambda: res.__setitem__("d2", [0 if pow(c, E, P) == 1 else 1 for c in gmC]))
check("app: GM-1024 decrypt Euler criterion x2000 bits",
      fold_bits(res["d2"]) + (10**10 if res["d2"] == msg else 0))

# ---- prim: the symbol kernels on 1024-bit operands -------------------------

jacA = [big_of(3000 + i, 1024) for i in range(1, 5001)]
jacM = [big_of(9000 + i, 1024) | 1 for i in range(1, 5001)]


def jac_loop():
    for _ in range(4):
        res["j"] = [gmpy2.jacobi(a, m) for a, m in zip(jacA, jacM)]


bench("prim: JacobiSymbol 1024-bit x5000x4", jac_loop)
check("prim: JacobiSymbol 1024-bit x5000x4", sum(res["j"]) % CHECK_MOD)

kroA = [(-1) ** i * jacA[i - 1] for i in range(1, 5001)]
kroM = [jacM[i - 1] - i % 2 for i in range(1, 5001)]


def kro_loop():
    for _ in range(4):
        res["k"] = [gmpy2.kronecker(a, m) for a, m in zip(kroA, kroM)]


bench("prim: KroneckerSymbol signed/even 1024-bit x5000x4", kro_loop)
check("prim: KroneckerSymbol signed/even 1024-bit x5000x4", sum(res["k"]) % CHECK_MOD)
