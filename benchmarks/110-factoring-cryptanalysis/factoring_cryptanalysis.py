#!/usr/bin/env python3
"""Experiment 110 -- Factoring cryptanalysis (group F).  See factoring_cryptanalysis.m.

Baseline: Pollard rho (Brent), Pollard p-1 and Fermat as plain Python ("app:",
identical algorithm/order to the .m/.gp halves); the "prim:" case is
sympy.factorint, Python's native integer factoriser (gmpy2 is not a factoriser).
Each check is the sum of the recovered smaller factors mod CHECK_MOD, plus a
recovered-count guard.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, CHECK_MOD

import gmpy2
from math import gcd, isqrt, log
from sympy import factorint, primerange, prime, primepi

require(["sympy:factorint"])


def next_prime(x):
    return int(gmpy2.next_prime(x))


def smaller(n, f):
    return min(f, n // f)


# ---- app: Pollard rho (Brent) -------------------------------------------


def rho_brent(n):
    if n % 2 == 0:
        return 2
    y, c, m, g, r, q = 2, 1, 128, 1, 1, 1
    x = ys = 0
    while g == 1:
        x = y
        for _ in range(r):
            y = (y * y + c) % n
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = (y * y + c) % n
                q = (q * abs(x - y)) % n
            g = gcd(q, n)
            k += m
        r *= 2
    if g == n:
        g = 1
        while g == 1:
            ys = (ys * ys + c) % n
            g = gcd(abs(x - ys), n)
    return g


def weak_rho(s, small, bits):
    return next_prime(big_of(s, small)) * next_prime(big_of(s + 1, bits - small))


rho_ns = [weak_rho(2501, 28, 64), weak_rho(2511, 28, 80),
          weak_rho(2521, 28, 96), weak_rho(2531, 28, 112)]
res = {}

bench("app: Pollard rho (Brent) 28-bit factor x4",
      lambda: res.__setitem__("rhoF", [rho_brent(n) for n in rho_ns]))
check("app: Pollard rho (Brent) 28-bit factor x4",
      sum(smaller(n, f) for n, f in zip(rho_ns, res["rhoF"])) % CHECK_MOD + CHECK_MOD * (sum(1 for n, f in zip(rho_ns, res["rhoF"]) if n % f == 0)))

# ---- app: Pollard p-1 (stage 1) -----------------------------------------

PRIMES = list(primerange(2, 20000))


def pm1(n, B1):
    a = 2
    for pr in PRIMES:
        if pr > B1:
            break
        e = int(log(B1) / log(pr))
        a = pow(a, pr ** e, n)
    return gcd(a - 1, n)


def smooth_prime(b, t0):
    base = 2
    for p in primerange(2, b + 1):
        base *= p
    t = t0
    p = base * t + 1
    while not gmpy2.is_prime(p):
        t += 1
        p = base * t + 1
    return int(p)


pm1_keys = [(smooth_prime(40, 1) * next_prime(big_of(2601, 45)), 60),
            (smooth_prime(50, 3) * next_prime(big_of(2602, 55)), 80),
            (smooth_prime(60, 2) * next_prime(big_of(2603, 60)), 100),
            (smooth_prime(70, 5) * next_prime(big_of(2604, 64)), 120)]

bench("app: Pollard p-1 stage 1 x120",
      lambda: [res.__setitem__("pm1F", [pm1(n, b) for n, b in pm1_keys]) for _ in range(30)])
check("app: Pollard p-1 stage 1 x120",
      sum(smaller(n, f) for (n, _), f in zip(pm1_keys, res["pm1F"])) % CHECK_MOD + CHECK_MOD * (sum(1 for (n, _), f in zip(pm1_keys, res["pm1F"]) if 1 < f < n and n % f == 0)))

# ---- app: Fermat (close primes) -----------------------------------------


def fermat(n):
    a = isqrt(n)
    if a * a < n:
        a += 1
    b2 = a * a - n
    while isqrt(b2) ** 2 != b2:
        a += 1
        b2 = a * a - n
    return a - isqrt(b2)


def close_semiprime(s, pb, gap):
    p = next_prime(big_of(s, pb))
    return p * next_prime(p + gap)


fermat_ns = [close_semiprime(2701, 40, 2 ** 25), close_semiprime(2702, 40, 2 ** 25),
             close_semiprime(2703, 40, 2 ** 26), close_semiprime(2704, 40, 2 ** 24)]

bench("app: Fermat (close primes) x4",
      lambda: res.__setitem__("ferF", [fermat(n) for n in fermat_ns]))
check("app: Fermat (close primes) x4",
      sum(smaller(n, f) for n, f in zip(fermat_ns, res["ferF"])) % CHECK_MOD + CHECK_MOD * (sum(1 for n, f in zip(fermat_ns, res["ferF"]) if 1 < f < n and n % f == 0)))

# ---- prim: factorint on balanced semiprimes -----------------------------


def bal(s, bits):
    return next_prime(big_of(s, bits // 2)) * next_prime(big_of(s + 1, bits // 2))


def smallest_factor(n):
    return min(factorint(n))


bal80 = [bal(2800 + 2 * k, 80) for k in range(25)]
bench_once("prim: FactorInteger balanced 80-bit x25",
           lambda: res.__setitem__("sf80", [smallest_factor(n) for n in bal80]))
check("prim: FactorInteger balanced 80-bit x25", sum(res["sf80"]) % CHECK_MOD)

bal100 = [bal(2850 + 2 * k, 100) for k in range(10)]
bench_once("prim: FactorInteger balanced 100-bit x10",
           lambda: res.__setitem__("sf100", [smallest_factor(n) for n in bal100]))
check("prim: FactorInteger balanced 100-bit x10", sum(res["sf100"]) % CHECK_MOD)

bal120 = [bal(2880 + 2 * k, 120) for k in range(3)]
bench_once("prim: FactorInteger balanced 120-bit x3",
           lambda: res.__setitem__("sf120", [smallest_factor(n) for n in bal120]))
check("prim: FactorInteger balanced 120-bit x3", sum(res["sf120"]) % CHECK_MOD)

n160 = bal(2900, 160)
bench_once("prim: FactorInteger balanced 160-bit x1",
           lambda: res.__setitem__("sf160", smallest_factor(n160)))
check("prim: FactorInteger balanced 160-bit x1", res["sf160"] % CHECK_MOD)
