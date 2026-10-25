#!/usr/bin/env python3
"""Experiment 101 -- RSA-2048 encrypt / decrypt, plain vs CRT (group F).  See rsa_crt_decrypt.m.

Baseline: "prim:" cases use gmpy2 (GMP) -- gmpy2.powmod for one exponentiation,
gmpy2.powmod_base_list for the whole-list encryption (gmpy2's batch spelling of
the Listable PowerMod[list, e, n]), and sympy.ntheory.modular.crt for the
ChineseRemainder recombination (gmpy2 has no CRT).  "app:" cases are plain
Python ints and built-in pow(), as an idiomatic pure-Python RSA is written.
The key is generated untimed, exactly as in 100-rsa-keygen.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, CHECK_MOD

import gmpy2
from math import gcd, lcm
from sympy.ntheory.modular import crt

require(["gmpy2:powmod", "gmpy2:powmod_base_list", "sympy.ntheory.modular:crt",
         "gmpy2:is_prime"])

E = 65537

# ---- the key (untimed) ----------------------------------------------------


def gen_prime(start):
    c = start | 1
    while not (gmpy2.is_prime(c, 25) and gcd(c - 1, E) == 1):
        c += 2
    return c


P = gen_prime(big_of(1101, 1024))
Q = gen_prime(big_of(1102, 1024))
N = P * Q
D = pow(E, -1, lcm(P - 1, Q - 1))
DP = D % (P - 1)
DQ = D % (Q - 1)
QINV = pow(Q, -1, P)

msgs = [big_of(s, 2000) for s in range(2001, 6001)]
msgs200 = msgs[:200]
cts = [pow(m, E, N) for m in msgs]
cts50 = cts[:50]
cts200 = cts[:200]
gN = gmpy2.mpz(N)
res = {}

# ---- prim: encryption ------------------------------------------------------

bench("prim: PowerMod encrypt e=65537 2048-bit x4000",
      lambda: res.__setitem__("enc", [gmpy2.powmod(m, E, gN) for m in msgs]))
check("prim: PowerMod encrypt e=65537 2048-bit x4000", int(sum(res["enc"])) % CHECK_MOD)

bench("prim: PowerMod[list] encrypt e=65537 2048-bit x4000",
      lambda: res.__setitem__("encL", gmpy2.powmod_base_list(msgs, E, gN)))
check("prim: PowerMod[list] encrypt e=65537 2048-bit x4000", int(sum(res["encL"])) % CHECK_MOD)

# ---- prim: plain decryption ------------------------------------------------

bench("prim: PowerMod decrypt c^d mod n 2048-bit x50",
      lambda: res.__setitem__("dec", [gmpy2.powmod(c, D, gN) for c in cts50]))
check("prim: PowerMod decrypt c^d mod n 2048-bit x50", int(sum(res["dec"])) % CHECK_MOD)

# ---- app: CRT (Garner) decryption ------------------------------------------


def crt_decrypt(c):
    m1 = pow(c, DP, P)
    m2 = pow(c, DQ, Q)
    h = (QINV * (m1 - m2)) % P
    return m2 + h * Q


bench("app: RSA-2048 CRT decrypt x200",
      lambda: res.__setitem__("decC", [crt_decrypt(c) for c in cts200]))
check("app: RSA-2048 CRT decrypt x200", sum(res["decC"]) % CHECK_MOD)

# ---- recombination only ------------------------------------------------------

m1s = [pow(c, DP, P) for c in cts]
m2s = [pow(c, DQ, Q) for c in cts]
pairs = list(zip(m1s, m2s))

bench("prim: ChineseRemainder recombine 2048-bit x20000",
      lambda: [res.__setitem__("recC", [int(crt([P, Q], [a, b])[0]) for a, b in pairs])
               for _ in range(5)])
check("prim: ChineseRemainder recombine 2048-bit x20000", sum(res["recC"]) % CHECK_MOD)

bench("app: Garner recombine 2048-bit x20000",
      lambda: [res.__setitem__("recG", [b + ((QINV * (a - b)) % P) * Q for a, b in pairs])
               for _ in range(5)])
check("app: Garner recombine 2048-bit x20000", sum(res["recG"]) % CHECK_MOD)

# ---- app: full round trip -----------------------------------------------------

bench("app: RSA-2048 encrypt+CRT decrypt round trip x200",
      lambda: res.__setitem__("rt", sum(1 for m in msgs200 if crt_decrypt(pow(m, E, N)) == m)))
check("app: RSA-2048 encrypt+CRT decrypt round trip x200", res["rt"])
