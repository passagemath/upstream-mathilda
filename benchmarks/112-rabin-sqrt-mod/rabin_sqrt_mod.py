#!/usr/bin/env python3
"""Experiment 112 -- Rabin cryptosystem & modular square roots (group F).

See rabin_sqrt_mod.m.  Baseline: Tonelli-Shanks and Cipolla as plain Python
("app:"); the native primitive ("prim:") is sympy.ntheory.sqrt_mod -- gmpy2
has no modular square root.  The "prim: PowerMod[a,1/2,p]" and
"prim: SqrtMod" labels both exercise sympy.sqrt_mod here (PowerMod[a,1/2,p] is
the WL spelling of exactly this primitive; SqrtMod is a name WL lacks).  Every
root is canonicalised to min(r, p-r) so all systems agree on one integer.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, prime_of, CHECK_MOD

import gmpy2
from sympy.ntheory import sqrt_mod

require(["sympy.ntheory:sqrt_mod"])


def prime_mod(seed, bits, r, m):
    p = prime_of(seed, bits)
    while p % m != r:
        p = int(gmpy2.next_prime(p))
    return p


def canon(r, p):
    return min(r, p - r)


def tonelli(n, p):
    Q, S = p - 1, 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    if S == 1:
        return pow(n, (p + 1) // 4, p)
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    M, c, t, R = S, pow(z, Q, p), pow(n, Q, p), pow(n, (Q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = (tt * tt) % p
            i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M, c = i, (b * b) % p
        t, R = (t * c) % p, (R * b) % p
    return R


def qr(x, p):
    return x == 0 or pow(x, (p - 1) // 2, p) == 1


def fmul(a, b, w, p):
    x1, y1 = a
    x2, y2 = b
    return ((x1 * x2 + y1 * y2 * w) % p, (x1 * y2 + y1 * x2) % p)


def fpow(base, e, w, p):
    result, b = (1, 0), base
    while e > 0:
        if e & 1:
            result = fmul(result, b, w, p)
        b = fmul(b, b, w, p)
        e >>= 1
    return result


def cipolla(n, p):
    a0 = 1
    while qr((a0 * a0 - n) % p, p):
        a0 += 1
    w = (a0 * a0 - n) % p
    return fpow((a0, 1), (p + 1) // 2, w, p)[0]


p256 = prime_mod(1201, 256, 1, 8)
as256 = [pow(big_of(s, 256), 2, p256) for s in range(10, 110)]
p1024 = prime_mod(1251, 1024, 1, 8)
as1024 = [pow(big_of(s, 1024), 2, p1024) for s in range(10, 20)]
res = {}

bench("app: Tonelli-Shanks 256-bit p=1mod8 x100",
      lambda: res.__setitem__("ts256", [canon(tonelli(a, p256), p256) for a in as256]))
check("app: Tonelli-Shanks 256-bit p=1mod8 x100", sum(res["ts256"]) % CHECK_MOD)

bench_once("app: Tonelli-Shanks 1024-bit p=1mod8 x10",
           lambda: res.__setitem__("ts1024", [canon(tonelli(a, p1024), p1024) for a in as1024]))
check("app: Tonelli-Shanks 1024-bit p=1mod8 x10", sum(res["ts1024"]) % CHECK_MOD)

bench("app: Cipolla 256-bit p=1mod8 x100",
      lambda: res.__setitem__("cp256", [canon(cipolla(a, p256), p256) for a in as256]))
check("app: Cipolla 256-bit p=1mod8 x100", sum(res["cp256"]) % CHECK_MOD)

bench_once("app: Cipolla 1024-bit p=1mod8 x10",
           lambda: res.__setitem__("cp1024", [canon(cipolla(a, p1024), p1024) for a in as1024]))
check("app: Cipolla 1024-bit p=1mod8 x10", sum(res["cp1024"]) % CHECK_MOD)

bench("prim: PowerMod[a,1/2,p] 256-bit x500",
      lambda: res.__setitem__("pm256", [canon(int(sqrt_mod(a, p256)), p256) for _ in range(5) for a in as256][-100:]))
check("prim: PowerMod[a,1/2,p] 256-bit x500", sum(res["pm256"]) % CHECK_MOD)

bench("prim: PowerMod[a,1/2,p] 1024-bit x40",
      lambda: res.__setitem__("pm1024", [canon(int(sqrt_mod(a, p1024)), p1024) for _ in range(4) for a in as1024][-10:]))
check("prim: PowerMod[a,1/2,p] 1024-bit x40", sum(res["pm1024"]) % CHECK_MOD)

bench("prim: SqrtMod 256-bit x500",
      lambda: res.__setitem__("sm256", [canon(int(sqrt_mod(a, p256)), p256) for _ in range(5) for a in as256][-100:]))
check("prim: SqrtMod 256-bit x500", sum(res["sm256"]) % CHECK_MOD)

pr = prime_mod(1301, 256, 3, 4)
qr_ = prime_mod(1302, 256, 3, 4)
nr = pr * qr_
ypr = pow(pr, -1, qr_)
yqr = pow(qr_, -1, pr)


def rabin_roots(c):
    mp = pow(c, (pr + 1) // 4, pr)
    mq = pow(c, (qr_ + 1) // 4, qr_)
    u = (mp * qr_ * yqr + mq * pr * ypr) % nr
    v = (mp * qr_ * yqr - mq * pr * ypr) % nr
    return (canon(u, nr) + canon(v, nr)) % nr


msgs = [big_of(s, 256) % nr for s in range(500, 1000)]
cts = [pow(m, 2, nr) for m in msgs]
bench("app: Rabin decrypt 512-bit n x500",
      lambda: res.__setitem__("rr", [rabin_roots(c) for c in cts]))
check("app: Rabin decrypt 512-bit n x500", sum(res["rr"]) % CHECK_MOD)
