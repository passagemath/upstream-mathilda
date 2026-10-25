#!/usr/bin/env python3
"""Experiment 102 -- primality testing (group F).  See primality_testing.m.

Baseline: "prim:" cases use gmpy2.is_prime(n, 25) (GMP: trial division + 25
Miller-Rabin rounds) for PrimeQ and sympy.isprime for the "proven" row --
NOTE sympy.isprime is BPSW, not a proof: Python has no ECPP/APRCL in the
standard stack, so that row is the nearest available, not a like-for-like.
"app:" cases are the same Miller-Rabin / Baillie-PSW code as the .m half,
plain Python ints and built-in pow().
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, prime_of, CHECK_MOD

import gmpy2
import sympy
from math import gcd

require(["gmpy2:is_prime", "sympy:isprime", "sympy:jacobi_symbol"])

# ---- inputs (untimed) -------------------------------------------------------

primes = ([prime_of(s, 512) for s in range(401, 411)] +
          [prime_of(s, 1024) for s in range(411, 421)] +
          [prime_of(s, 2048) for s in range(421, 431)])
odds1024 = [big_of(s, 1024) | 1 for s in range(5001, 7001)]
odds500 = odds1024[:500]


def chernick(start):
    """First k >= start with 6k+1, 12k+1, 18k+1 all prime; their product."""
    k = start
    while not (gmpy2.is_prime(6 * k + 1, 25) and gmpy2.is_prime(12 * k + 1, 25)
               and gmpy2.is_prime(18 * k + 1, 25)):
        k += 1
    return (6 * k + 1) * (12 * k + 1) * (18 * k + 1)


carmichaels = [chernick(big_of(s, 64)) for s in range(3001, 3009)]
spsps = [3215031751, 2152302898747, 3474749660383, 341550071728321,
         3825123056546413051, 318665857834031151167461,
         3317044064679887385961981]
mersenne_c = [2**p - 1 for p in sympy.primerange(1000, 1200)]
fermat_c = [2**(2**k) + 1 for k in range(8, 12)]
pseudos = carmichaels + spsps + mersenne_c + fermat_c

# ---- app: the tests written in the language -------------------------------

MR_BASES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def mr_round(n, a, d, s):
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def miller_rabin(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    return all(a % n == 0 or mr_round(n, a, d, s) for a in MR_BASES)


def jacobi(a, n):
    return sympy.jacobi_symbol(a, n)


def half(x, n):
    return (x + n) // 2 if x % 2 else x // 2


def strong_lucas(n):
    """Strong Lucas probable-prime test with Selfridge's method A parameters."""
    D = 5
    while True:
        j = jacobi(D % n, n)
        if j == -1:
            break
        if j == 0 and abs(D) != n:
            return False
        D = -D - 2 if D > 0 else -D + 2
    P, Q = 1, (1 - D) // 4
    d, s = n + 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    bits = bin(d)[3:]
    U, V, Qk = 1, P, Q % n
    for b in bits:
        U, V = U * V % n, (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if b == "1":
            U, V = half((P * U + V) % n, n), half((D * U + P * V) % n, n)
            Qk = Qk * Q % n
    if U == 0 or V == 0:
        return True
    for _ in range(s - 1):
        V = (V * V - 2 * Qk) % n
        if V == 0:
            return True
        Qk = Qk * Qk % n
    return False


def is_square(n):
    r = gmpy2.isqrt(n)
    return r * r == n


def bpsw(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    if not mr_round(n, 2, d, s):
        return False
    if is_square(n):
        return False
    return strong_lucas(n)


# The same tests with early exits as loop conditions (the .m "no-Return"
# variants; see primality_testing.m).  Same algorithm, same verdicts.


def mr_round2(n, a, d, s):
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    r = 1
    while r < s and x != n - 1:
        x = x * x % n
        r += 1
    return x == n - 1


def miller_rabin2(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    return all(a % n == 0 or mr_round2(n, a, d, s) for a in MR_BASES)


def strong_lucas2(n):
    D = 5
    while jacobi(D % n, n) == 1:
        D = -D - 2 if D > 0 else -D + 2
    if jacobi(D % n, n) == 0 and abs(D) != n:
        return False
    P, Q = 1, (1 - D) // 4
    d, s = n + 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    U, V, Qk = 1, P, Q % n
    for b in bin(d)[3:]:
        U, V = U * V % n, (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if b == "1":
            U, V = half((P * U + V) % n, n), half((D * U + P * V) % n, n)
            Qk = Qk * Q % n
    if U == 0 or V == 0:
        return True
    r = 1
    while r < s and V != 0:
        V = (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        r += 1
    return V == 0


def bpsw2(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    if not mr_round2(n, 2, d, s):
        return False
    if is_square(n):
        return False
    return strong_lucas2(n)


res = {}

# ---- prim ---------------------------------------------------------------------

bench("prim: PrimeQ 512/1024/2048-bit primes x30",
      lambda: res.__setitem__("c", sum(1 for n in primes if gmpy2.is_prime(n, 25))))
check("prim: PrimeQ 512/1024/2048-bit primes x30", res["c"])

bench("prim: PrimeQ 1024-bit odd x2000",
      lambda: res.__setitem__("c", sum(1 for n in odds1024 if gmpy2.is_prime(n, 25))))
check("prim: PrimeQ 1024-bit odd x2000", res["c"])

bench("prim: PrimeQ pseudoprimes (Carmichael/spsp/Mersenne/Fermat) x47",
      lambda: res.__setitem__("c", sum(1 for n in pseudos if gmpy2.is_prime(n, 25))))
check("prim: PrimeQ pseudoprimes (Carmichael/spsp/Mersenne/Fermat) x47", res["c"])

M4423 = 2**4423 - 1
bench("prim: PrimeQ Mersenne 2^4423-1",
      lambda: res.__setitem__("c", 1 if gmpy2.is_prime(M4423, 25) else 0))
check("prim: PrimeQ Mersenne 2^4423-1", res["c"])

proof_primes = [prime_of(s, 256) for s in range(431, 435)]
bench_once("prim: proven primality 256-bit prime x4",
           lambda: res.__setitem__("c", sum(1 for n in proof_primes if sympy.isprime(n))))
check("prim: proven primality 256-bit prime x4", res["c"])

# ---- app ----------------------------------------------------------------------

bench("app: Miller-Rabin 12 bases primes x30",
      lambda: res.__setitem__("c", sum(1 for n in primes if miller_rabin(n))))
check("app: Miller-Rabin 12 bases primes x30", res["c"])

bench("app: Miller-Rabin 12 bases pseudoprimes x47",
      lambda: res.__setitem__("c", sum(1 for n in pseudos if miller_rabin(n))))
check("app: Miller-Rabin 12 bases pseudoprimes x47", res["c"])

bench("app: Baillie-PSW primes x30",
      lambda: res.__setitem__("c", sum(1 for n in primes if bpsw(n))))
check("app: Baillie-PSW primes x30", res["c"])

bench("app: Baillie-PSW pseudoprimes x47",
      lambda: res.__setitem__("c", sum(1 for n in pseudos if bpsw(n))))
check("app: Baillie-PSW pseudoprimes x47", res["c"])

bench("app: Baillie-PSW 1024-bit odd x500",
      lambda: res.__setitem__("c", sum(1 for n in odds500 if bpsw(n))))
check("app: Baillie-PSW 1024-bit odd x500", res["c"])

bench("app: Miller-Rabin 12 bases (no-Return) primes x30",
      lambda: res.__setitem__("c", sum(1 for n in primes if miller_rabin2(n))))
check("app: Miller-Rabin 12 bases (no-Return) primes x30", res["c"])

bench("app: Miller-Rabin 12 bases (no-Return) pseudoprimes x47",
      lambda: res.__setitem__("c", sum(1 for n in pseudos if miller_rabin2(n))))
check("app: Miller-Rabin 12 bases (no-Return) pseudoprimes x47", res["c"])

bench("app: Baillie-PSW (no-Return) primes x30",
      lambda: res.__setitem__("c", sum(1 for n in primes if bpsw2(n))))
check("app: Baillie-PSW (no-Return) primes x30", res["c"])

bench("app: Baillie-PSW (no-Return) 1024-bit odd x500",
      lambda: res.__setitem__("c", sum(1 for n in odds500 if bpsw2(n))))
check("app: Baillie-PSW (no-Return) 1024-bit odd x500", res["c"])
