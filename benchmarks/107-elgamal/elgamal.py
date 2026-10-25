#!/usr/bin/env python3
"""Experiment 107 -- ElGamal over the RFC 3526 MODP-2048 group (group F).
See elgamal.m.

"app:" cases: the same textbook ElGamal on plain Python ints with the
built-in pow() (idiomatic pure-Python crypto; pow(s, -1, p) is the inverse).
Products of all ciphertexts are formed whole (math.prod) and reduced once,
as the .m half's Times @@ does.

"prim:" baseline: gmpy2.powmod (GMP mpz_powm), the 2048-bit exponentiation
every ElGamal step reduces to.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, CHECK_MOD, MODP2048

import gmpy2
from math import prod

require(["gmpy2:powmod"])

p = MODP2048
g = 2

# ---- prim ----------------------------------------------------------------

pm_bases = [big_of(s, 2048) % p for s in range(2001, 2201)]
pm_exps = [big_of(s, 2048) for s in range(2201, 2401)]
res = {}
bench("prim: PowerMod 2048-bit x200",
      lambda: res.__setitem__("pms", [gmpy2.powmod(b, e, p) for b, e in zip(pm_bases, pm_exps)]))
check("prim: PowerMod 2048-bit x200", int(sum(res["pms"])) % CHECK_MOD)

# ---- app -----------------------------------------------------------------


def eg_keygen(seed):
    x = big_of(seed, 2048) % (p - 2) + 1
    return (x, pow(g, x, p))


def eg_encrypt(m, h, y):
    return (pow(g, y, p), m * pow(h, y, p) % p)


def eg_decrypt_inv(c, x):
    return c[1] * pow(pow(c[0], x, p), -1, p) % p


def eg_decrypt_fermat(c, x):
    return c[1] * pow(pow(c[0], x, p), p - 2, p) % p


bench("app: ElGamal keygen x50",
      lambda: res.__setitem__("keys", [eg_keygen(s) for s in range(901, 951)]))
check("app: ElGamal keygen x50", sum(a + b for a, b in res["keys"]) % CHECK_MOD)

priv_x, pub_h = res["keys"][0]
msgs = [big_of(s, 2048) % (p - 1) + 1 for s in range(1001, 1051)]
ephs = [big_of(s, 2048) % (p - 2) + 1 for s in range(1101, 1151)]

bench("app: ElGamal encrypt x50",
      lambda: res.__setitem__("cts", [eg_encrypt(m, pub_h, y) for m, y in zip(msgs, ephs)]))
cts = res["cts"]
check("app: ElGamal encrypt x50", sum(a + b for a, b in cts) % CHECK_MOD)

bench("app: ElGamal decrypt x50 (PowerMod[s,-1,p])",
      lambda: res.__setitem__("d1", [eg_decrypt_inv(c, priv_x) for c in cts]))
check("app: ElGamal decrypt x50 (PowerMod[s,-1,p])",
      sum(res["d1"]) % CHECK_MOD if res["d1"] == msgs else -1)

bench("app: ElGamal decrypt x50 (s^(p-2))",
      lambda: res.__setitem__("d2", [eg_decrypt_fermat(c, priv_x) for c in cts]))
check("app: ElGamal decrypt x50 (s^(p-2))",
      sum(res["d2"]) % CHECK_MOD if res["d2"] == msgs else -1)


def run_homo():
    for _ in range(20):
        prod_c = (prod(c[0] for c in cts) % p, prod(c[1] for c in cts) % p)
        res["pm"] = eg_decrypt_inv(prod_c, priv_x)


bench("app: homomorphic product of 50 ciphertexts + decrypt x20", run_homo)
check("app: homomorphic product of 50 ciphertexts + decrypt x20",
      res["pm"] % CHECK_MOD if res["pm"] == prod(msgs) % p else -1)

rands = [big_of(s, 2048) % (p - 2) + 1 for s in range(1201, 1251)]
bench("app: ElGamal re-randomise x50",
      lambda: res.__setitem__("rcts", [(c[0] * pow(g, r, p) % p, c[1] * pow(pub_h, r, p) % p)
                                       for c, r in zip(cts, rands)]))
check("app: ElGamal re-randomise x50",
      sum(a + b for a, b in res["rcts"]) % CHECK_MOD
      if [eg_decrypt_inv(c, priv_x) for c in res["rcts"]] == msgs else -1)
