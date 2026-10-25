#!/usr/bin/env python3
"""Experiment 103 -- finite-field Diffie-Hellman (group F).  See diffie_hellman.m.

Baseline: "prim:" cases use sympy.primitive_root (sympy has the only Python
primitive-root routine), gmpy2.powmod and gmpy2.jacobi.  "app:" cases are
plain Python ints and built-in pow(); the safe-prime search is driven by
gmpy2.is_prime(n, 25), the primality test a Python implementation uses.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, CHECK_MOD, MODP2048

import gmpy2
from sympy import primitive_root

require(["gmpy2:is_prime", "sympy:primitive_root", "gmpy2:powmod", "gmpy2:jacobi"])

P = MODP2048
Q = (P - 1) // 2
res = {}

# ---- app: safe-prime generation -------------------------------------------


def safe_prime(start):
    q = start | 1
    while not (gmpy2.is_prime(q, 25) and gmpy2.is_prime(2 * q + 1, 25)):
        q += 2
    return 2 * q + 1


bench("app: safe-prime search 256-bit x4",
      lambda: res.__setitem__("sp256", [safe_prime(big_of(s, 255)) for s in range(1, 5)]))
check("app: safe-prime search 256-bit x4", sum(res["sp256"]) % CHECK_MOD)

bench_once("app: safe-prime search 512-bit x2",
           lambda: res.__setitem__("sp512", [safe_prime(big_of(s, 511)) for s in range(9, 11)]))
check("app: safe-prime search 512-bit x2", sum(res["sp512"]) % CHECK_MOD)

# ---- primitive roots of safe primes ----------------------------------------

safes = [safe_prime(big_of(s, 255)) for s in range(1, 9)]

bench("prim: PrimitiveRoot 256-bit safe prime x800",
      lambda: [res.__setitem__("roots", [int(primitive_root(p)) for p in safes])
               for _ in range(100)])
check("prim: PrimitiveRoot 256-bit safe prime x800", sum(res["roots"]))


def safe_root(p):
    q, g = (p - 1) // 2, 2
    while g * g % p == 1 or pow(g, q, p) == 1:
        g += 1
    return g


bench("app: primitive root of safe prime 256-bit x800",
      lambda: [res.__setitem__("roots2", [safe_root(p) for p in safes]) for _ in range(100)])
check("app: primitive root of safe prime 256-bit x800", sum(res["roots2"]))

# ---- MODP-2048 key exchange -------------------------------------------------

privA = [big_of(s, 256) for s in range(6001, 6051)]
privB = [big_of(s, 256) for s in range(7001, 7051)]
gP = gmpy2.mpz(P)

bench("prim: PowerMod 2^a mod MODP-2048 (256-bit a) x100",
      lambda: res.__setitem__("pubs", [gmpy2.powmod(2, a, gP) for a in privA + privB]))
check("prim: PowerMod 2^a mod MODP-2048 (256-bit a) x100", int(sum(res["pubs"])) % CHECK_MOD)


def dh_exchange(a, b):
    pa = pow(2, a, P)
    pb = pow(2, b, P)
    sa = pow(pb, a, P)
    sb = pow(pa, b, P)
    return sa if sa == sb else 0


bench("app: MODP-2048 DH exchange x50",
      lambda: res.__setitem__("sec", [dh_exchange(a, b) for a, b in zip(privA, privB)]))
check("app: MODP-2048 DH exchange x50", sum(res["sec"]) % CHECK_MOD)

# ---- peer-key validation ------------------------------------------------------

cands = [pow(2, a, P) for a in privA] + [big_of(s, 2048) % P for s in range(8001, 8051)]


def valid(y):
    return 1 < y < P - 1 and pow(y, Q, P) == 1


bench("app: subgroup validation y^q mod p MODP-2048 x100",
      lambda: res.__setitem__("nv", sum(1 for y in cands if valid(y))))
check("app: subgroup validation y^q mod p MODP-2048 x100", res["nv"])

gcands = [gmpy2.mpz(y) for y in cands]
bench("prim: JacobiSymbol validation MODP-2048 x10000",
      lambda: [res.__setitem__("nj", sum(1 for y in gcands if 1 < y < P - 1 and gmpy2.jacobi(y, gP) == 1))
               for _ in range(100)])
check("prim: JacobiSymbol validation MODP-2048 x10000", res["nj"])
