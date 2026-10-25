#!/usr/bin/env python3
"""Experiment 114 -- AES-128 and GF(2^8) (group F).  See aes_gf256.m.

The same AES-128 written twice in plain Python: idiomatically with the
^ & | << >> operators, and as the "(no-bitwise)" variant that XORs bytes by
lookup in a 256x256 table built (untimed) from base-2 digit arithmetic --
the identical algorithm to the .m / .gp halves.  The GF(2)[x] route is a
pure-Python polynomial extended Euclid over GF(2) (polynomials as int bit
masks, the idiomatic Python representation).  The library baseline is
pycryptodome's AES (C, AES-NI where available) in ECB mode.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import small_list, CHECK_MOD

from Crypto.Cipher import AES

require(["Crypto.Cipher.AES:new"])


def fold(bs):
    v = 0
    for b in bs:
        v = 256 * v + b
    return v % CHECK_MOD


fips_key = list(range(16))
fips_pt = list(bytes.fromhex("00112233445566778899aabbccddeeff"))
shift_idx = [r + 4 * ((c + r) % 4) for c in range(4) for r in range(4)]
ecb_key = small_list(1401, 16, 256)
ecb_flat = small_list(1402, 128 * 16, 256)
ecb_blocks = [ecb_flat[16 * i:16 * i + 16] for i in range(128)]

# =================== idiomatic (^ & | << >>) =================================


def xtime(b):
    return ((b << 1) & 255) ^ (27 if b >= 128 else 0)


def gmul(a, b):
    p = 0
    while b > 0:
        if b & 1:
            p ^= a
        a = xtime(a)
        b >>= 1
    return p


def ginv(a):
    r, s = 1, a
    for _ in range(7):
        s = gmul(s, s)
        r = gmul(r, s)
    return r


def rotl8(b, k):
    return ((b << k) | (b >> (8 - k))) & 255


def make_sbox():
    out = []
    for a in range(256):
        b = ginv(a)
        out.append(b ^ rotl8(b, 1) ^ rotl8(b, 2) ^ rotl8(b, 3) ^ rotl8(b, 4) ^ 99)
    return out


def expand_key(key, sb):
    w = [key[4 * i:4 * i + 4] for i in range(4)]
    rc = 1
    for i in range(4, 44):
        t = w[i - 1]
        if i % 4 == 0:
            t = [sb[b] for b in t[1:] + t[:1]]
            t[0] ^= rc
            rc = xtime(rc)
        w.append([x ^ y for x, y in zip(w[i - 4], t)])
    return [w[4 * r] + w[4 * r + 1] + w[4 * r + 2] + w[4 * r + 3] for r in range(11)]


def mix_columns(s):
    out = []
    for c in range(4):
        col = s[4 * c:4 * c + 4]
        t = col[0] ^ col[1] ^ col[2] ^ col[3]
        out += [col[r] ^ t ^ xtime(col[r] ^ col[(r + 1) % 4]) for r in range(4)]
    return out


def encrypt_block(pt, rk, sb):
    s = [x ^ y for x, y in zip(pt, rk[0])]
    for r in range(1, 10):
        u = [sb[b] for b in s]
        s = [x ^ y for x, y in zip(mix_columns([u[i] for i in shift_idx]), rk[r])]
    u = [sb[b] for b in s]
    return [x ^ y for x, y in zip([u[i] for i in shift_idx], rk[10])]

# =================== (no-bitwise): arithmetic + table lookup =================


def bits8(a):
    d = []
    for _ in range(8):
        a, r = divmod(a, 2)
        d.append(r)
    return d[::-1]


def from_bits(d):
    v = 0
    for b in d:
        v = 2 * v + b
    return v


XOR_TAB = [[from_bits([(x + y) % 2 for x, y in zip(bits8(a), bits8(b))])
            for b in range(256)] for a in range(256)]


def x8(a, b):
    return XOR_TAB[a][b]


def x8v(u, v):
    return [XOR_TAB[a][b] for a, b in zip(u, v)]


def xtime_n(b):
    return x8(2 * b - 256, 27) if b >= 128 else 2 * b


def gmul_n(a, b):
    p = 0
    while b > 0:
        if b % 2 == 1:
            p = x8(p, a)
        a = xtime_n(a)
        b //= 2
    return p


def ginv_n(a):
    r, s = 1, a
    for _ in range(7):
        s = gmul_n(s, s)
        r = gmul_n(r, s)
    return r


def rotl8_n(b, k):
    return (b * 2 ** k) % 256 + b // 2 ** (8 - k)


def make_sbox_n():
    out = []
    for a in range(256):
        b = ginv_n(a)
        v = b
        for y in (rotl8_n(b, 1), rotl8_n(b, 2), rotl8_n(b, 3), rotl8_n(b, 4), 99):
            v = x8(v, y)
        out.append(v)
    return out


def expand_key_n(key, sb):
    w = [key[4 * i:4 * i + 4] for i in range(4)]
    rc = 1
    for i in range(4, 44):
        t = w[i - 1]
        if i % 4 == 0:
            t = [sb[b] for b in t[1:] + t[:1]]
            t[0] = x8(t[0], rc)
            rc = xtime_n(rc)
        w.append(x8v(w[i - 4], t))
    return [w[4 * r] + w[4 * r + 1] + w[4 * r + 2] + w[4 * r + 3] for r in range(11)]


def mix_columns_n(s):
    out = []
    for c in range(4):
        col = s[4 * c:4 * c + 4]
        t = x8(x8(x8(col[0], col[1]), col[2]), col[3])
        out += [x8(x8(col[r], t), xtime_n(x8(col[r], col[(r + 1) % 4]))) for r in range(4)]
    return out


def encrypt_block_n(pt, rk, sb):
    s = x8v(pt, rk[0])
    for r in range(1, 10):
        u = [sb[b] for b in s]
        s = x8v(mix_columns_n([u[i] for i in shift_idx]), rk[r])
    u = [sb[b] for b in s]
    return x8v([u[i] for i in shift_idx], rk[10])

# =================== cases ===================================================


res = {}


def rep(n, f):
    for _ in range(n):
        f()


bench("app: AES S-box generation x4", lambda: rep(4, lambda: res.__setitem__("sb", make_sbox())))
check("app: AES S-box generation x4", fold(res["sb"]))

bench("app: AES S-box generation x4 (no-bitwise)",
      lambda: rep(4, lambda: res.__setitem__("sbN", make_sbox_n())))
check("app: AES S-box generation x4 (no-bitwise)", fold(res["sbN"]))
sbox, sboxN = res["sb"], res["sbN"]

bench("app: AES-128 FIPS-197 C.1 x100",
      lambda: rep(100, lambda: res.__setitem__("f", encrypt_block(fips_pt, expand_key(fips_key, sbox), sbox))))
check("app: AES-128 FIPS-197 C.1 x100", fold(res["f"]))
assert bytes(res["f"]).hex() == "69c4e0d86a7b0430d8cdb78070b4c55a"

bench("app: AES-128 FIPS-197 C.1 x100 (no-bitwise)",
      lambda: rep(100, lambda: res.__setitem__("fN", encrypt_block_n(fips_pt, expand_key_n(fips_key, sboxN), sboxN))))
check("app: AES-128 FIPS-197 C.1 x100 (no-bitwise)", fold(res["fN"]))


def ecb(blocks):
    rk = expand_key(ecb_key, sbox)
    return [encrypt_block(b, rk, sbox) for b in blocks]


def ecb_n(blocks):
    rk = expand_key_n(ecb_key, sboxN)
    return [encrypt_block_n(b, rk, sboxN) for b in blocks]


bench("app: AES-128 ECB encrypt 128 blocks", lambda: res.__setitem__("e", ecb(ecb_blocks)))
check("app: AES-128 ECB encrypt 128 blocks", sum(fold(c) for c in res["e"]) % CHECK_MOD)

bench("app: AES-128 ECB encrypt 128 blocks (no-bitwise)", lambda: res.__setitem__("eN", ecb_n(ecb_blocks)))
check("app: AES-128 ECB encrypt 128 blocks (no-bitwise)", sum(fold(c) for c in res["eN"]) % CHECK_MOD)

# ---- GF(2)[x] route: extended Euclid on bit-mask polynomials ---------------

AES_POLY = 0x11b


def pdivmod(a, b):
    q, db = 0, b.bit_length()
    while a.bit_length() >= db:
        s = a.bit_length() - db
        q ^= 1 << s
        a ^= b << s
    return q, a


def pmul(a, b):
    p = 0
    while b:
        if b & 1:
            p ^= a
        a <<= 1
        b >>= 1
    return p


def poly_inv(a):
    """s with s*a + t*m = 1 over GF(2)[x]; returned as the integer s(2)."""
    r0, r1, s0, s1 = a, AES_POLY, 1, 0
    while r1:
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, s0 ^ pmul(q, s1)
    return s0


bench("prim: GF(2^8) inverse via PolynomialExtendedGCD Modulus->2 x7650",
      lambda: rep(30, lambda: res.__setitem__("inv", [poly_inv(a) for a in range(1, 256)])))
check("prim: GF(2^8) inverse via PolynomialExtendedGCD Modulus->2 x7650",
      sum(a * v for a, v in zip(range(1, 256), res["inv"])) % CHECK_MOD)

# ---- the library ------------------------------------------------------------

lib_data = bytes(ecb_flat) * 512
bench("prim: AES-128 ECB library 65536 blocks",
      lambda: res.__setitem__("lib", AES.new(bytes(ecb_key), AES.MODE_ECB).encrypt(lib_data)))
check("prim: AES-128 ECB library 65536 blocks",
      sum(fold(res["lib"][16 * i:16 * i + 16]) for i in range(65536)) % CHECK_MOD)
