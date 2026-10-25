#!/usr/bin/env python3
"""Experiment 106 -- ECDSA sign / verify over secp256k1 (group F).  See ecdsa.m.

"app:" cases: the same Jacobian arithmetic, Shamir's-trick verifier, and
deterministic (NOT RFC 6979) LCG nonces as the .m half, on plain Python ints
with the built-in pow(x, -1, m).

"prim:" baselines:
  * SHA-256 of 1000 messages -- hashlib.sha256 over the UTF-8 bytes.
  * ECDSA verify -- the `cryptography` package (OpenSSL secp256k1), verifying
    the app-generated signatures with a prehashed 32-byte z.
  * ECDSA sign -- OpenSSL draws its own nonce, so it cannot reproduce the
    deterministic signatures; coincurve (custom nonce) is not installed, so the
    case is bench_if'd on coincurve and SKIPs.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, bench_if, check, check_if, require, have
from crypto_data import big_of, CHECK_MOD, SECP256K1_P, SECP256K1_N, SECP256K1_GX, SECP256K1_GY

import hashlib

EC_SPEC = "cryptography.hazmat.primitives.asymmetric.ec:SECP256K1"
require(["hashlib:sha256", EC_SPEC, "coincurve:PrivateKey"])

P = SECP256K1_P
N = SECP256K1_N
G = (SECP256K1_GX, SECP256K1_GY)
INFJ = (1, 1, 0)

# ---- Jacobian point arithmetic (as in 105-ec-arithmetic) -----------------


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


def shamir_j(u1, u2, g, q):
    """u1*g + u2*q in one joint double-and-add, g + q precomputed."""
    gj, qj = to_j(g), to_j(q)
    gq = add_j(gj, qj)
    n = max(u1.bit_length(), u2.bit_length())
    r = INFJ
    for i in range(n - 1, -1, -1):
        r = dbl_j(r)
        b1, b2 = (u1 >> i) & 1, (u2 >> i) & 1
        if b1 and b2:
            r = add_j(r, gq)
        elif b1:
            r = add_j(r, gj)
        elif b2:
            r = add_j(r, qj)
    return from_j(r)

# ---- ECDSA ---------------------------------------------------------------


def ecdsa_sign(z, d, k):
    r = mul_j(k, G)[0] % N
    s = pow(k, -1, N) * (z + r * d) % N
    return (r, s)


def ecdsa_verify(z, sig, q):
    r, s = sig
    if not (1 <= r < N and 1 <= s < N):
        return 0
    w = pow(s, -1, N)
    x = shamir_j(z * w % N, r * w % N, G, q)
    if x is None:
        return 0
    return 1 if x[0] % N == r else 0


priv_d = big_of(701, 256) % N
pub_q = mul_j(priv_d, G)

zs = [big_of(s, 256) for s in range(801, 849)]
nonces = [big_of(s, 256) % N for s in range(1801, 1849)]
res = {}

bench("app: ECDSA sign x48",
      lambda: res.__setitem__("sigs", [ecdsa_sign(z, priv_d, k) for z, k in zip(zs, nonces)]))
check("app: ECDSA sign x48", sum(r + s for r, s in res["sigs"]) % CHECK_MOD)
sigs = res["sigs"]

bench("app: ECDSA verify x48 (Shamir)",
      lambda: res.__setitem__("nv", sum(ecdsa_verify(z, sg, pub_q) for z, sg in zip(zs, sigs))))
check("app: ECDSA verify x48 (Shamir)", res["nv"])

bench("app: ECDSA verify tampered x48 (Shamir)",
      lambda: res.__setitem__("nb", sum(ecdsa_verify(z + 1, sg, pub_q) for z, sg in zip(zs, sigs))))
check("app: ECDSA verify tampered x48 (Shamir)", res["nb"])

# ---- prim: SHA-256 -------------------------------------------------------

msgs = ["message %d" % i for i in range(1, 1001)]
bench("prim: SHA-256 of 1000 messages",
      lambda: res.__setitem__("hs", [int.from_bytes(hashlib.sha256(m.encode()).digest(), "big")
                                     for m in msgs]))
check("prim: SHA-256 of 1000 messages", sum(res["hs"]) % CHECK_MOD)

# ---- prim: native ECDSA --------------------------------------------------


def run_prim_sign():
    import coincurve  # not installed here: bench_if SKIPs before this runs
    raise NotImplementedError


bench_if("prim: ECDSA sign x48", "coincurve:PrivateKey", run_prim_sign)


def run_prim_verify():
    from cryptography.hazmat.primitives.asymmetric import ec, utils
    from cryptography.hazmat.primitives import hashes
    from cryptography.exceptions import InvalidSignature
    pk = ec.EllipticCurvePublicNumbers(pub_q[0], pub_q[1], ec.SECP256K1()).public_key()
    alg = ec.ECDSA(utils.Prehashed(hashes.SHA256()))
    n = 0
    for z, (r, s) in zip(zs, sigs):
        try:
            pk.verify(utils.encode_dss_signature(r, s), z.to_bytes(32, "big"), alg)
            n += 1
        except InvalidSignature:
            pass
    res["pv"] = n


bench_if("prim: ECDSA verify x48", EC_SPEC, run_prim_verify)
check_if("prim: ECDSA verify x48", EC_SPEC, lambda: res["pv"])
