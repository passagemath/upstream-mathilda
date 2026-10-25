#!/usr/bin/env python3
"""Experiment 104 -- discrete logarithms (group F).  See discrete_log.m.

Baseline: "prim:" cases use sympy.ntheory.discrete_log (order via factoring,
Pohlig-Hellman, then trial / BSGS / Pollard rho per prime factor) -- the only
discrete-log routine in the standard Python stack.  "app:" cases are plain
Python ints, built-in pow() and dict, the same algorithms as the .m half.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, prime_of, CHECK_MOD

import gmpy2
from math import isqrt, prod
from sympy import factorint
from sympy.ntheory import discrete_log
from sympy.ntheory.modular import crt

require(["sympy.ntheory:discrete_log", "sympy.ntheory.modular:crt", "sympy:factorint"])

# ---- instances (untimed) ----------------------------------------------------


def subgroup(bits, s):
    """[p, r, g]: r = prime_of(s, bits), p = 2kr + 1 the first prime of that
    shape, g = 2^((p-1)/r) (next base if that is 1) of prime order r."""
    r = prime_of(s, bits)
    k = 1
    while not gmpy2.is_prime(2 * k * r + 1, 25):
        k += 1
    p = 2 * k * r + 1
    b = 2
    while pow(b, (p - 1) // r, p) == 1:
        b += 1
    return [p, r, pow(b, (p - 1) // r, p)]


def instance(bits, s, t):
    p, r, g = subgroup(bits, s)
    x = big_of(t, 64) % r
    return [p, r, g, pow(g, x, p), x]


inst24 = [instance(24, 9000 + i, 9100 + i) for i in range(4)]
inst32 = [instance(32, 9200 + i, 9300 + i) for i in range(2)]

# smooth p: p - 1 = 2k * prod of 8 primes just above 40000, 45000, ..., 75000
qs = [int(gmpy2.next_prime(40000 + 5000 * i)) for i in range(8)]
L = prod(qs)
k = 1
while not gmpy2.is_prime(2 * k * L + 1, 25):
    k += 1
PS = 2 * k * L + 1
FS = sorted(factorint(PS - 1).items())


def least_prim_root(p, fs):
    g = 2
    while any(pow(g, (p - 1) // q, p) == 1 for q, _ in fs):
        g += 1
    return g


GS = least_prim_root(PS, FS)
xs_smooth = [big_of(s, 128) % (PS - 1) for s in range(9401, 9411)]
hs_smooth = [pow(GS, x, PS) for x in xs_smooth]
res = {}

# ---- app: baby-step giant-step --------------------------------------------------


def bsgs_inc(g, h, p, r):
    """Baby steps inserted one at a time into a hash table."""
    m = isqrt(r) + 1
    tab = {}
    y = 1
    for i in range(m):
        tab[y] = i
        y = y * g % p
    c = pow(g, -m, p)
    y, j = h, 0
    while y not in tab:
        y = y * c % p
        j += 1
    return j * m + tab[y]


def bsgs(g, h, p, r):
    """Baby-step table built in bulk from the list of powers."""
    m = isqrt(r) + 1
    pw = [1] * m
    for i in range(1, m):
        pw[i] = pw[i - 1] * g % p
    tab = dict(zip(pw, range(m)))
    c = pow(g, -m, p)
    y, j = h, 0
    while y not in tab:
        y = y * c % p
        j += 1
    return j * m + tab[y]


bench("app: BSGS 24-bit subgroup (incremental table) x4",
      lambda: res.__setitem__("b1", [bsgs_inc(g, h, p, r) for p, r, g, h, _ in inst24]))
check("app: BSGS 24-bit subgroup (incremental table) x4", sum(res["b1"]) % CHECK_MOD)

bench("app: BSGS 24-bit subgroup (bulk table) x4",
      lambda: res.__setitem__("b2", [bsgs(g, h, p, r) for p, r, g, h, _ in inst24]))
check("app: BSGS 24-bit subgroup (bulk table) x4", sum(res["b2"]) % CHECK_MOD)

bench("app: BSGS 32-bit subgroup (bulk table) x2",
      lambda: res.__setitem__("b3", [bsgs(g, h, p, r) for p, r, g, h, _ in inst32]))
check("app: BSGS 32-bit subgroup (bulk table) x2", sum(res["b3"]) % CHECK_MOD)

# ---- app: Pollard rho for logarithms ---------------------------------------------


def rho_step(x, a, b, g, h, p, r):
    s = x % 3
    if s == 0:
        return x * x % p, 2 * a % r, 2 * b % r
    if s == 1:
        return x * g % p, (a + 1) % r, b
    return x * h % p, a, (b + 1) % r


def pollard_rho(g, h, p, r):
    """Floyd cycle finding on the 3-partition walk; restart from g^t on a
    degenerate collision."""
    t = 0
    while True:
        x, a, b = pow(g, t, p), t, 0
        X, A, B = x, a, b
        while True:
            x, a, b = rho_step(x, a, b, g, h, p, r)
            X, A, B = rho_step(X, A, B, g, h, p, r)
            X, A, B = rho_step(X, A, B, g, h, p, r)
            if x == X:
                break
        den = (b - B) % r
        if den != 0:
            return (A - a) * pow(den, -1, r) % r
        t += 1


bench_once("app: Pollard rho 32-bit subgroup x2",
           lambda: res.__setitem__("r1", [pollard_rho(g, h, p, r) for p, r, g, h, _ in inst32]))
check("app: Pollard rho 32-bit subgroup x2", sum(res["r1"]) % CHECK_MOD)

# ---- app: Pohlig-Hellman over a smooth p -------------------------------------------


def pohlig_hellman(g, h, p, fs):
    n = p - 1
    rs, ms = [], []
    for q, e in fs:
        gq = pow(g, n // q, p)
        xk = 0
        for j in range(e):
            hk = pow(pow(g, -xk, p) * h % p, n // q ** (j + 1), p)
            xk += bsgs(gq, hk, p, q) * q ** j
        rs.append(xk)
        ms.append(q ** e)
    return int(crt(ms, rs)[0])


bench("app: Pohlig-Hellman 131-bit smooth p x10",
      lambda: res.__setitem__("ph", [pohlig_hellman(GS, h, PS, FS) for h in hs_smooth]))
check("app: Pohlig-Hellman 131-bit smooth p x10", sum(res["ph"]) % CHECK_MOD)

# ---- prim ---------------------------------------------------------------------------

bench("prim: discrete log 24-bit subgroup x40",
      lambda: [res.__setitem__("d1", [discrete_log(p, h, g) for p, r, g, h, _ in inst24])
               for _ in range(10)])
check("prim: discrete log 24-bit subgroup x40", sum(res["d1"]) % CHECK_MOD)

bench("prim: discrete log 32-bit subgroup x10",
      lambda: [res.__setitem__("d2", [discrete_log(p, h, g) for p, r, g, h, _ in inst32])
               for _ in range(5)])
check("prim: discrete log 32-bit subgroup x10", sum(res["d2"]) % CHECK_MOD)

bench("prim: discrete log 131-bit smooth p x40",
      lambda: [res.__setitem__("d3", [discrete_log(PS, h, GS) for h in hs_smooth])
               for _ in range(4)])
check("prim: discrete log 131-bit smooth p x40", sum(res["d3"]) % CHECK_MOD)

