#!/usr/bin/env python3
"""Experiment 109 -- RSA attacks (group F).  See rsa_attacks.m.

Baseline: plain Python ints + pow() for the "app:" attacks (idiomatic
pure-Python crypto).  The "prim:" continued-fraction cases use sympy:
continued_fraction for ContinuedFraction and continued_fraction_convergents
for Convergents (the head Mathilda lacks).  The Hastad cube root uses
gmpy2.iroot, Python's exact integer root primitive.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, CHECK_MOD

import gmpy2
from math import gcd, isqrt
from sympy import Rational
from sympy.ntheory.continued_fraction import (
    continued_fraction, continued_fraction_convergents)

require(["sympy.ntheory.continued_fraction:continued_fraction",
         "gmpy2:iroot"])

E1, E2 = 65537, 65539


def next_prime(x):
    return int(gmpy2.next_prime(x))


# ---- (a) Wiener ----------------------------------------------------------


def wiener_key(s):
    p, q = next_prime(big_of(s, 512)), next_prime(big_of(s + 1, 512))
    n, phi = p * q, (p - 1) * (q - 1)
    d = big_of(s + 2, 120) | 1
    while gcd(d, phi) != 1:
        d += 2
    e = pow(d, -1, phi)
    return (e, n, d)


def wiener(e, n):
    cf = continued_fraction(Rational(e, n))
    hm1, hm2, km1, km2, found = 1, 0, 0, 1, 0
    for a in cf:
        h, k = a * hm1 + hm2, a * km1 + km2
        if h > 0 and k > 0 and (e * k - 1) % h == 0:
            phi = (e * k - 1) // h
            b = n - phi + 1
            disc = b * b - 4 * n
            if disc >= 0 and isqrt(disc) ** 2 == disc:
                found = k
        hm2, hm1, km2, km1 = hm1, h, km1, k
    return found


wkeys = [wiener_key(1700 + 3 * s) for s in range(20)]
res = {}

bench("prim: ContinuedFraction e/N 1024-bit x600",
      lambda: [res.__setitem__("cfs", [list(continued_fraction(Rational(e, n))) for e, n, _ in wkeys]) for _ in range(30)])
check("prim: ContinuedFraction e/N 1024-bit x600",
      (sum(sum(cf) for cf in res["cfs"]) + sum(len(cf) for cf in res["cfs"])) % CHECK_MOD)


def wiener_from_conv(cv, e, n):
    found = 0
    for c in cv:
        c = Rational(c)
        h, k = int(c.p), int(c.q)
        if h > 0 and k > 0 and (e * k - 1) % h == 0:
            phi = (e * k - 1) // h
            b = n - phi + 1
            disc = b * b - 4 * n
            if disc >= 0 and isqrt(disc) ** 2 == disc:
                found = k
    return found


def convergents():
    for _ in range(2):
        res["cvrec"] = [wiener_from_conv(
            list(continued_fraction_convergents(continued_fraction(Rational(e, n)))), e, n)
            for e, n, _ in wkeys]


bench("prim: Convergents e/N 1024-bit x40", convergents)
check("prim: Convergents e/N 1024-bit x40", sum(res["cvrec"]) % CHECK_MOD)

bench("app: Wiener attack 1024-bit x20",
      lambda: res.__setitem__("wrec", [wiener(e, n) for e, n, _ in wkeys]))
check("app: Wiener attack 1024-bit x20", sum(res["wrec"]) % CHECK_MOD + CHECK_MOD * (sum(1 for r, (_, _, d) in zip(res["wrec"], wkeys) if r == d)))

# ---- (b) common modulus --------------------------------------------------

cmN = next_prime(big_of(1801, 512)) * next_prime(big_of(1802, 512))


def common_modulus(m):
    c1, c2 = pow(m, E1, cmN), pow(m, E2, cmN)
    g = gcd(E1, E2)
    # extended gcd
    a, b, old_a, old_b, r0, r1 = 0, 1, 1, 0, E1, E2
    while r1:
        qt = r0 // r1
        r0, r1 = r1, r0 - qt * r1
        old_a, a = a, old_a - qt * a
        old_b, b = b, old_b - qt * b
    return (pow(c1, old_a, cmN) * pow(c2, old_b, cmN)) % cmN


cm_msgs = [big_of(s, 900) % cmN for s in range(1810, 1910)]
bench("app: common-modulus attack 1024-bit x600",
      lambda: [res.__setitem__("cm", [common_modulus(m) for m in cm_msgs]) for _ in range(6)])
check("app: common-modulus attack 1024-bit x600", sum(res["cm"]) % CHECK_MOD + CHECK_MOD * (sum(1 for r, m in zip(res["cm"], cm_msgs) if r == m)))

# ---- (c) Hastad broadcast (e = 3) ---------------------------------------

h_mods = [next_prime(big_of(1900, 512)) * next_prime(big_of(1901, 512)),
          next_prime(big_of(1902, 512)) * next_prime(big_of(1903, 512)),
          next_prime(big_of(1904, 512)) * next_prime(big_of(1905, 512))]


def crt(rs, ms):
    M = 1
    for m in ms:
        M *= m
    x = 0
    for r, m in zip(rs, ms):
        Mi = M // m
        x += r * Mi * pow(Mi, -1, m)
    return x % M


def hastad(m):
    cs = [pow(m, 3, n) for n in h_mods]
    return int(gmpy2.iroot(crt(cs, h_mods), 3)[0])


h_msgs = [big_of(s, 300) for s in range(1950, 2000)]
bench("app: Hastad broadcast e=3 x1500",
      lambda: [res.__setitem__("h", [hastad(m) for m in h_msgs]) for _ in range(30)])
check("app: Hastad broadcast e=3 x1500", sum(res["h"]) % CHECK_MOD + CHECK_MOD * (sum(1 for r, m in zip(res["h"], h_msgs) if r == m)))
