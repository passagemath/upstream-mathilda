#!/usr/bin/env python3.11
"""Generate the FactorInteger shootout corpus (deterministic, seeded).

Writes corpus.json: a list of groups.  A *batch* group is timed as a whole
(many small inputs; per-call latency dominated), an *single* group is one hard
input per process (so a timeout or crash is a recorded outcome for that case
alone, never the whole run).

The harness verifies every engine's answer independently: the product of the
returned prime powers must equal n and every base must pass sympy.isprime.

Run with the Python that has sympy (python3.11 on the reference machine):
    /usr/local/opt/python@3.11/bin/python3.11 gen_corpus.py
"""
import json
import random
import math
from sympy import isprime, nextprime, primerange

SEED = 20261017
rng = random.Random(SEED)


def rp(digits=None, bits=None):
    """Random prime with an exact digit or bit count (seeded through rng)."""
    if bits is not None:
        lo, hi = 1 << (bits - 1), 1 << bits
    else:
        lo, hi = 10 ** (digits - 1), 10 ** digits
    while True:
        c = rng.randrange(lo, hi) | 1
        if isprime(c):
            return c


def smooth_prime(digits, B, plus=False):
    """Prime p with p-1 (or p+1 when plus) B-smooth, ~digits digits."""
    small = list(primerange(2, B))
    target = 10 ** (digits - 1)
    while True:
        m = 2
        while m < target:
            m *= rng.choice(small)
        p = m - 1 if plus else m + 1
        if isprime(p):
            return p


groups = []


def batch(name, desc, ns):
    groups.append({"name": name, "kind": "batch", "desc": desc,
                   "cases": [{"id": f"{name}-{i}", "n": str(n)} for i, n in enumerate(ns)]})


def single(name, desc, items):
    groups.append({"name": name, "kind": "single", "desc": desc,
                   "cases": [{"id": cid, "n": str(n)} for cid, n in items]})


# --- batch: small inputs, latency bound ------------------------------------
batch("u32", "2000 random integers in [2^31, 2^32)",
      [rng.randrange(1 << 31, 1 << 32) for _ in range(2000)])
batch("u63", "2000 random integers in [2^62, 2^63)",
      [rng.randrange(1 << 62, 1 << 63) for _ in range(2000)])
batch("semi63", "500 balanced semiprimes p*q, p,q ~ 2^31",
      [rp(bits=31) * rp(bits=32) for _ in range(500)])
batch("semi96", "20 balanced semiprimes p*q, p,q ~ 2^48",
      [rp(bits=48) * rp(bits=48) for _ in range(20)])
batch("rand40", "100 random 40-digit integers",
      [rng.randrange(10 ** 39, 10 ** 40) for _ in range(100)])


def smooth_num():
    small = list(primerange(2, 1 << 20))
    m = 1
    while m < 10 ** 59:
        m *= rng.choice(small)
    return m


batch("smooth", "200 ~60-digit products of random primes < 2^20", [smooth_num() for _ in range(200)])

# --- single: one hard input per process ------------------------------------
qs = []
for d in (30, 40, 50, 60, 70, 80, 90):
    for s in range(2):
        qs.append((f"qs{d}-{s}", rp(digits=d // 2) * rp(digits=d - d // 2)))
single("qs", "balanced semiprimes, 30-90 digits (quadratic-sieve territory)", qs)

ecm = []
for k in (20, 25, 30, 35):
    for s in range(2):
        ecm.append((f"ecm-p{k}-{s}", rp(digits=k) * rp(digits=150 - k)))
single("ecm", "150-digit composites hiding one p20..p35 factor (ECM territory)", ecm)

single("pm1", "p40 with p-1 10^5-smooth, times a p80 (Pollard p-1 territory)",
       [(f"pm1-{s}", smooth_prime(40, 10 ** 5) * rp(digits=80)) for s in range(3)])
single("pp1", "p40 with p+1 10^5-smooth, times a p80 (Williams p+1 territory)",
       [(f"pp1-{s}", smooth_prime(40, 10 ** 5, plus=True) * rp(digits=80)) for s in range(3)])


def close_pair(d, gap):
    p = rp(digits=d)
    return p * nextprime(p + gap)


single("fermat", "p*q with |p-q| tiny relative to sqrt(n) (Fermat/Lehman territory)",
       [("fermat-100d", close_pair(50, 10 ** 12)),
        ("fermat-80d", close_pair(40, 10 ** 10)),
        ("fermat-60d", close_pair(30, 10 ** 8))])

p30 = rp(digits=30)
m61 = (1 << 61) - 1
single("power", "perfect and near-perfect powers",
       [("pow-p30^7", p30 ** 7),
        ("pow-m61^3*q^2", m61 ** 3 * nextprime(10 ** 25) ** 2),
        ("pow-2^4096", 1 << 4096),
        ("pow-smooth^5", (2 ** 10 * 3 ** 7 * 1000003 * 998244353) ** 5)])

single("special", "classic named composites",
       [("F7=2^128+1", 2 ** 128 + 1),
        ("F8=2^256+1", 2 ** 256 + 1),
        ("M137=2^137-1", 2 ** 137 - 1),
        ("M149=2^149-1", 2 ** 149 - 1),
        ("30!+1", math.factorial(30) + 1),
        ("10^60+1", 10 ** 60 + 1)])

# --- method-specific groups (appended: earlier groups keep their values) ---
batch("rho96", "30 numbers p32 * q64 (Pollard rho territory)",
      [rp(bits=32) * rp(bits=64) for _ in range(30)])
batch("squfof60", "100 balanced semiprimes ~2^60 (SQUFOF territory)",
      [rp(bits=30) * rp(bits=30) for _ in range(100)])
batch("lehman40", "100 balanced semiprimes ~2^40 (Lehman territory)",
      [rp(bits=20) * rp(bits=20) for _ in range(100)])


def smooth_1e5():
    small = list(primerange(2, 10 ** 5))
    m = 1
    while m < 10 ** 39:
        m *= rng.choice(small)
    return m


batch("trial40", "50 ~40-digit products of primes < 10^5 (trial-division territory)",
      [smooth_1e5() for _ in range(50)])

with open("corpus.json", "w") as f:
    json.dump({"seed": SEED, "groups": groups}, f, indent=1)

print(sum(len(g["cases"]) for g in groups), "cases in", len(groups), "groups")
