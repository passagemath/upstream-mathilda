#!/usr/bin/env python3
"""Experiment 108 -- Paillier additively homomorphic encryption (group F).
See paillier.m.

"app:" cases: the same Paillier (g = n + 1) on plain Python ints with the
built-in pow() (idiomatic pure-Python crypto), math.lcm, and //.

"prim:" baseline: gmpy2.powmod (GMP mpz_powm) modulo the 4096-bit n^2.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, prime_of, CHECK_MOD

import gmpy2
from math import lcm

require(["gmpy2:powmod"])

# ---- keys ----------------------------------------------------------------

pqs = [(prime_of(2 * s + 3001, 1024), prime_of(2 * s + 3002, 1024)) for s in range(1, 9)]


def pl_L(u, n):
    return (u - 1) // n


def pl_keygen(pq):
    p, q = pq
    n = p * q
    n2 = n * n
    lam = lcm(p - 1, q - 1)
    mu = pow(pl_L(pow(n + 1, lam, n2), n), -1, n)
    return (n, lam, mu)


res = {}
bench("app: Paillier keygen x8 (given p, q)",
      lambda: res.__setitem__("keys", [pl_keygen(pq) for pq in pqs]))
check("app: Paillier keygen x8 (given p, q)", sum(sum(k) for k in res["keys"]) % CHECK_MOD)

nn, lam, mu = res["keys"][0]
nn2 = nn * nn

# ---- prim ----------------------------------------------------------------

pm_bases = [big_of(s, 4096) % nn2 for s in range(3101, 3151)]
bench("prim: PowerMod 4096-bit modulus x50",
      lambda: res.__setitem__("pms", [gmpy2.powmod(b, nn, nn2) for b in pm_bases]))
check("prim: PowerMod 4096-bit modulus x50", int(sum(res["pms"])) % CHECK_MOD)

# ---- app: encrypt / decrypt ----------------------------------------------


def pl_encrypt(m, r):
    return pow(nn + 1, m, nn2) * pow(r, nn, nn2) % nn2


def pl_encrypt_fast(m, r):
    return (1 + m * nn) * pow(r, nn, nn2) % nn2


def pl_decrypt(c):
    return pl_L(pow(c, lam, nn2), nn) * mu % nn


msgs = [big_of(s, 2048) % nn for s in range(3201, 3221)]
rs = [big_of(s, 2048) % nn for s in range(3301, 3321)]

bench("app: Paillier encrypt x20",
      lambda: res.__setitem__("cts", [pl_encrypt(m, r) for m, r in zip(msgs, rs)]))
cts = res["cts"]
check("app: Paillier encrypt x20", sum(cts) % CHECK_MOD)

bench("app: Paillier encrypt x20 (g=n+1: 1+mn)",
      lambda: res.__setitem__("ctsf", [pl_encrypt_fast(m, r) for m, r in zip(msgs, rs)]))
check("app: Paillier encrypt x20 (g=n+1: 1+mn)", sum(res["ctsf"]) % CHECK_MOD)

bench("app: Paillier decrypt x20",
      lambda: res.__setitem__("dec", [pl_decrypt(c) for c in cts]))
check("app: Paillier decrypt x20", sum(res["dec"]) % CHECK_MOD if res["dec"] == msgs else -1)

# ---- app: homomorphic sum of 1000 ciphertexts ----------------------------

rns = [pow(r, nn, nn2) for r in rs]
sum_msgs = [big_of(s, 2000) for s in range(4001, 5001)]
sum_cts = [(1 + sum_msgs[i - 1] * nn) * rns[i % 20] % nn2 for i in range(1, 1001)]


def run_sum():
    for _ in range(10):
        c = 1
        for x in sum_cts:
            c = c * x % nn2
        res["msum"] = pl_decrypt(c)


bench("app: homomorphic sum of 1000 ciphertexts + decrypt x10", run_sum)
check("app: homomorphic sum of 1000 ciphertexts + decrypt x10",
      res["msum"] % CHECK_MOD if res["msum"] == sum(sum_msgs) % nn else -1)

# ---- app: scalar multiplication ------------------------------------------

scalars = [big_of(s, 2048) % nn for s in range(3401, 3421)]
bench("app: Paillier scalar mult c^k x20 (k < n)",
      lambda: res.__setitem__("scs", [pow(c, k, nn2) for c, k in zip(cts, scalars)]))
check("app: Paillier scalar mult c^k x20 (k < n)",
      sum(res["scs"]) % CHECK_MOD
      if pl_decrypt(res["scs"][0]) == scalars[0] * msgs[0] % nn else -1)
