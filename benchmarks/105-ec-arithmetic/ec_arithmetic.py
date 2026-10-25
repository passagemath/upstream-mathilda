#!/usr/bin/env python3
"""Experiment 105 -- elliptic-curve point arithmetic on secp256k1 (group F).
See ec_arithmetic.m.

The "app:" cases are the same affine / Jacobian / ladder code as the .m half
on plain Python ints, inverse via the built-in pow(x, -1, p).

"prim:" baselines: Python has no EC group in the standard library.
  * scalar mult k*G -- the `cryptography` package (OpenSSL's secp256k1):
    ec.derive_private_key(k, SECP256K1()) computes k*G natively; timed when
    installed, else SKIP.
  * point add -- no installed library exposes a bare point addition
    (coincurve's PublicKey.combine_keys would; it is not installed), so the
    case is bench_if'd on coincurve and SKIPs.  The in-language code is NOT
    substituted under a "prim:" label.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, bench_if, check, check_if, require, have
from crypto_data import big_of, CHECK_MOD, SECP256K1_P, SECP256K1_N, SECP256K1_GX, SECP256K1_GY

require(["cryptography.hazmat.primitives.asymmetric.ec:SECP256K1", "coincurve:PublicKey"])

P = SECP256K1_P
N = SECP256K1_N
G = (SECP256K1_GX, SECP256K1_GY)
INFJ = (1, 1, 0)

# ---- affine coordinates (None is the point at infinity) ------------------


def dbl_a(pt):
    if pt is None or pt[1] == 0:
        return None
    x, y = pt
    l = 3 * x * x * pow(2 * y, -1, P) % P
    x3 = (l * l - 2 * x) % P
    return (x3, (l * (x - x3) - y) % P)


def add_a(p1, p2):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2:
        return None if (y1 + y2) % P == 0 else dbl_a(p1)
    l = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (l * l - x1 - x2) % P
    return (x3, (l * (x1 - x3) - y1) % P)


def mul_a(k, pt):
    r = None
    for b in bin(k)[2:]:
        r = dbl_a(r)
        if b == "1":
            r = add_a(r, pt)
    return r

# ---- Jacobian coordinates ------------------------------------------------


def to_j(pt):
    return INFJ if pt is None else (pt[0], pt[1], 1)


def from_j(pt):
    x, y, z = pt
    if z == 0:
        return None
    zi = pow(z, -1, P)
    return (x * zi * zi % P, y * zi * zi * zi % P)


def dbl_j(pt):
    x, y, z = pt
    if z == 0 or y == 0:
        return INFJ
    a = x * x % P
    b = y * y % P
    c = b * b % P
    d = 2 * ((x + b) ** 2 - a - c) % P
    e = 3 * a
    x3 = (e * e - 2 * d) % P
    return (x3, (e * (d - x3) - 8 * c) % P, 2 * y * z % P)


def add_j(p1, p2):
    if p1[2] == 0:
        return p2
    if p2[2] == 0:
        return p1
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    z1z1 = z1 * z1 % P
    z2z2 = z2 * z2 % P
    u1 = x1 * z2z2 % P
    u2 = x2 * z1z1 % P
    s1 = y1 * z2 * z2z2 % P
    s2 = y2 * z1 * z1z1 % P
    h = (u2 - u1) % P
    r = (s2 - s1) % P
    if h == 0:
        return dbl_j(p1) if r == 0 else INFJ
    h2 = h * h % P
    h3 = h * h2 % P
    u1h2 = u1 * h2 % P
    x3 = (r * r - h3 - 2 * u1h2) % P
    return (x3, (r * (u1h2 - x3) - s1 * h3) % P, z1 * z2 * h % P)


def mul_j(k, pt):
    r, pj = INFJ, to_j(pt)
    for b in bin(k)[2:]:
        r = dbl_j(r)
        if b == "1":
            r = add_j(r, pj)
    return from_j(r)


def ladder_j(k, pt):
    r0, r1 = INFJ, to_j(pt)
    for b in bin(k)[2:]:
        if b == "0":
            r1 = add_j(r0, r1)
            r0 = dbl_j(r0)
        else:
            r0 = add_j(r0, r1)
            r1 = dbl_j(r1)
    return from_j(r0)


def x_digest(pts):
    return sum(p[0] for p in pts) % CHECK_MOD


ks = [big_of(s, 256) % N for s in range(501, 565)]
res = {}

# ---- app: single operations, batched -------------------------------------

g2 = dbl_a(G)


def run_add_a():
    q = g2
    for _ in range(10000):
        q = add_a(q, G)
    res["qa"] = q


bench("app: affine point add x10000", run_add_a)
check("app: affine point add x10000", res["qa"][0] % CHECK_MOD)

gJ = to_j(G)


def run_add_j():
    q = to_j(g2)
    for _ in range(10000):
        q = add_j(q, gJ)
    res["qj"] = q


bench("app: Jacobian point add x10000", run_add_j)
check("app: Jacobian point add x10000", from_j(res["qj"])[0] % CHECK_MOD)


def run_dbl_a():
    q = G
    for _ in range(10000):
        q = dbl_a(q)
    res["da"] = q


bench("app: affine point double x10000", run_dbl_a)
check("app: affine point double x10000", res["da"][0] % CHECK_MOD)


def run_dbl_j():
    q = gJ
    for _ in range(10000):
        q = dbl_j(q)
    res["dj"] = q


bench("app: Jacobian point double x10000", run_dbl_j)
check("app: Jacobian point double x10000", from_j(res["dj"])[0] % CHECK_MOD)

# ---- app: scalar multiplication ------------------------------------------

bench("app: scalar mult affine double-and-add x64",
      lambda: res.__setitem__("ma", [mul_a(k, G) for k in ks]))
check("app: scalar mult affine double-and-add x64", x_digest(res["ma"]))

bench("app: scalar mult Jacobian double-and-add x64",
      lambda: res.__setitem__("mj", [mul_j(k, G) for k in ks]))
check("app: scalar mult Jacobian double-and-add x64", x_digest(res["mj"]))

bench("app: scalar mult Jacobian Montgomery ladder x64",
      lambda: res.__setitem__("ml", [ladder_j(k, G) for k in ks]))
check("app: scalar mult Jacobian Montgomery ladder x64", x_digest(res["ml"]))

# ---- prim: a native EC library -------------------------------------------


def run_prim_add():
    import coincurve  # not installed here: bench_if SKIPs before this runs
    raise NotImplementedError


bench_if("prim: point add x10000", "coincurve:PublicKey", run_prim_add)

EC_SPEC = "cryptography.hazmat.primitives.asymmetric.ec:SECP256K1"


def run_prim_mul():
    from cryptography.hazmat.primitives.asymmetric import ec
    curve = ec.SECP256K1()
    res["me"] = [ec.derive_private_key(k, curve).public_key().public_numbers().x
                 for k in ks]


bench_if("prim: scalar mult k*G x64", EC_SPEC, run_prim_mul)
check_if("prim: scalar mult k*G x64", EC_SPEC, lambda: sum(res["me"]) % CHECK_MOD)
