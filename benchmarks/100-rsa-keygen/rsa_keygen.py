#!/usr/bin/env python3
"""Experiment 100 -- RSA key generation (group F).  See rsa_keygen.m.

Baseline: gmpy2 (GMP), the library a Python RSA implementation actually uses.
gmpy2.is_prime(n, 25) is 25 Miller-Rabin rounds; next_prime is GMP's
mpz_nextprime.  The "app:" keygen is the same incremental search as the .m
half, driven by gmpy2.is_prime.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, prime_of, CHECK_MOD

import gmpy2
from math import gcd, lcm

require(["gmpy2:next_prime", "gmpy2:is_prime", "gmpy2:powmod", "gmpy2:invert"])

E = 65537

# ---- prim ---------------------------------------------------------------

starts1024 = [big_of(s, 1024) for s in range(101, 105)]
res = {}
bench("prim: NextPrime 1024-bit x4",
      lambda: res.__setitem__("nps", [gmpy2.next_prime(x) for x in starts1024]))
check("prim: NextPrime 1024-bit x4", int(sum(res["nps"])) % CHECK_MOD)

primes2048 = [prime_of(s, 2048) for s in range(201, 209)]
bench("prim: PrimeQ 2048-bit prime x8",
      lambda: res.__setitem__("c", sum(1 for p in primes2048 if gmpy2.is_prime(p, 25))))
check("prim: PrimeQ 2048-bit prime x8", res["c"])

odds2048 = [big_of(s, 2048) | 1 for s in range(301, 701)]
bench("prim: PrimeQ 2048-bit odd x400",
      lambda: res.__setitem__("c", sum(1 for p in odds2048 if gmpy2.is_prime(p, 25))))
check("prim: PrimeQ 2048-bit odd x400", res["c"])

lams = [x for x in (big_of(s, 2048) for s in range(1001, 3101)) if gcd(x, E) == 1][:2000]
bench("prim: PowerMod[e,-1,lambda] 2048-bit x20000",
      lambda: [res.__setitem__("ds", [pow(E, -1, l) for l in lams]) for _ in range(10)])
check("prim: PowerMod[e,-1,lambda] 2048-bit x20000", sum(res["ds"]) % CHECK_MOD)

bench("prim: ModularInverse 2048-bit x20000",
      lambda: [res.__setitem__("ds2", [gmpy2.invert(E, l) for l in lams]) for _ in range(10)])
check("prim: ModularInverse 2048-bit x20000", int(sum(res["ds2"])) % CHECK_MOD)

# ---- app ----------------------------------------------------------------


def gen_prime(start):
    c = start | 1
    while not (gmpy2.is_prime(c, 25) and gcd(c - 1, E) == 1):
        c += 2
    return c


def rsa_key(bits, s):
    p = gen_prime(big_of(s, bits // 2))
    q = gen_prime(big_of(s + 1, bits // 2))
    lam = lcm(p - 1, q - 1)
    d = pow(E, -1, lam)
    return [p * q, d, d % (p - 1), d % (q - 1), pow(q, -1, p)]


def key_digest(keys):
    return sum(sum(k) for k in keys) % CHECK_MOD


bench("app: RSA-1024 keygen x8",
      lambda: res.__setitem__("k1", [rsa_key(1024, 2 * s) for s in range(1, 9)]))
check("app: RSA-1024 keygen x8", key_digest(res["k1"]))

bench_once("app: RSA-2048 keygen x2",
           lambda: res.__setitem__("k2", [rsa_key(2048, 2 * s + 50) for s in range(1, 3)]))
check("app: RSA-2048 keygen x2", key_digest(res["k2"]))
