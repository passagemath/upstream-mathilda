#!/usr/bin/env python3
"""Experiment 115 -- ARX hashing and stream ciphers (group F).  See arx_hash_stream.m.

SHA-256, HMAC-SHA256 and the ChaCha20 block function in plain Python, twice:
idiomatically with ^ & | << >> ~, and as the "(no-bitwise)" variant (32-bit
XOR/AND as four 8-bit lookups in 256x256 tables built untimed from base-2
digit arithmetic; rotations/shifts by // and %) -- the identical algorithm
to the .m / .gp halves.  The native primitive is hashlib.sha256 (OpenSSL).
Every answer is also asserted against hashlib / hmac / the RFC 7539 vector.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import small_list, CHECK_MOD

import hashlib, hmac

require(["hashlib:sha256", "hmac:new"])

M32 = 2 ** 32
SHA_K = [int(x, 16) for x in "428a2f98 71374491 b5c0fbcf e9b5dba5 3956c25b 59f111f1 923f82a4 ab1c5ed5 d807aa98 12835b01 243185be 550c7dc3 72be5d74 80deb1fe 9bdc06a7 c19bf174 e49b69c1 efbe4786 0fc19dc6 240ca1cc 2de92c6f 4a7484aa 5cb0a9dc 76f988da 983e5152 a831c66d b00327c8 bf597fc7 c6e00bf3 d5a79147 06ca6351 14292967 27b70a85 2e1b2138 4d2c6dfc 53380d13 650a7354 766a0abb 81c2c92e 92722c85 a2bfe8a1 a81a664b c24b8b70 c76c51a3 d192e819 d6990624 f40e3585 106aa070 19a4c116 1e376c08 2748774c 34b0bcb5 391c0cb3 4ed8aa4a 5b9cca4f 682e6ff3 748f82ee 78a5636f 84c87814 8cc70208 90befffa a4506ceb bef9a3f7 c67178f2".split()]
SHA_H0 = [int(x, 16) for x in "6a09e667 bb67ae85 3c6ef372 a54ff53a 510e527f 9b05688c 1f83d9ab 5be0cd19".split()]


def from_digits(ds, base):
    v = 0
    for d in ds:
        v = v * base + d
    return v


def to_digits(v, base, n):
    out = []
    for _ in range(n):
        v, r = divmod(v, base)
        out.append(r)
    return out[::-1]


def sha_pad(msg):
    return msg + [128] + [0] * ((55 - len(msg)) % 64) + to_digits(8 * len(msg), 256, 8)


msgs = [[c + 32 for c in small_list(1500 + j, 256, 94)] for j in range(1, 9)]
abc = list(b"abc")
hmac_key = list(b"Jefe")
hmac_msg = list(b"what do ya want for nothing?")
chacha_key = list(range(32))
chacha_nonce = list(bytes.fromhex("000000090000004a00000000"))

# =================== idiomatic ==============================================


def rotr(x, n):
    return (x >> n) | ((x << (32 - n)) & (M32 - 1))


def rotl(x, n):
    return rotr(x, 32 - n)


def sha_compress(h, blk):
    w = [from_digits(blk[4 * i:4 * i + 4], 256) for i in range(16)] + [0] * 48
    for t in range(16, 64):
        s0 = rotr(w[t - 15], 7) ^ rotr(w[t - 15], 18) ^ (w[t - 15] >> 3)
        s1 = rotr(w[t - 2], 17) ^ rotr(w[t - 2], 19) ^ (w[t - 2] >> 10)
        w[t] = (w[t - 16] + s0 + w[t - 7] + s1) % M32
    a, b, c, d, e, f, g, hh = h
    for t in range(64):
        t1 = (hh + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) +
              ((e & f) ^ ((e ^ (M32 - 1)) & g)) + SHA_K[t] + w[t]) % M32
        t2 = ((rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c))) % M32
        hh, g, f, e, d, c, b, a = g, f, e, (d + t1) % M32, c, b, a, (t1 + t2) % M32
    return [(x + y) % M32 for x, y in zip(h, [a, b, c, d, e, f, g, hh])]


def sha256(msg):
    p, h = sha_pad(msg), SHA_H0
    for i in range(0, len(p), 64):
        h = sha_compress(h, p[i:i + 64])
    return from_digits(h, M32)


def hmac_sha256(key, msg):
    k = key + [0] * (64 - len(key))
    inner = sha256([b ^ 54 for b in k] + msg)
    return sha256([b ^ 92 for b in k] + to_digits(inner, 256, 32))


def le32(bs):
    return [from_digits(bs[4 * i:4 * i + 4][::-1], 256) for i in range(len(bs) // 4)]


def chacha_init(key, ctr, nonce):
    return [1634760805, 857760878, 2036477234, 1797285236] + le32(key) + [ctr] + le32(nonce)


CHACHA_ROUNDS = [(0, 4, 8, 12), (1, 5, 9, 13), (2, 6, 10, 14), (3, 7, 11, 15),
                 (0, 5, 10, 15), (1, 6, 11, 12), (2, 7, 8, 13), (3, 4, 9, 14)]


def quarter(s, idx):
    i, j, k, l = idx
    a, b, c, d = s[i], s[j], s[k], s[l]
    a = (a + b) % M32; d = rotl(d ^ a, 16)
    c = (c + d) % M32; b = rotl(b ^ c, 12)
    a = (a + b) % M32; d = rotl(d ^ a, 8)
    c = (c + d) % M32; b = rotl(b ^ c, 7)
    s[i], s[j], s[k], s[l] = a, b, c, d


def chacha_block(key, ctr, nonce):
    s0 = chacha_init(key, ctr, nonce)
    s = s0[:]
    for _ in range(10):
        for idx in CHACHA_ROUNDS:
            quarter(s, idx)
    out = []
    for x, y in zip(s, s0):
        out += to_digits((x + y) % M32, 256, 4)[::-1]
    return out

# =================== (no-bitwise) ===========================================


def bits8(a):
    return to_digits(a, 2, 8)


XOR_TAB = [[from_digits([(x + y) % 2 for x, y in zip(bits8(a), bits8(b))], 2)
            for b in range(256)] for a in range(256)]
AND_TAB = [[from_digits([x * y for x, y in zip(bits8(a), bits8(b))], 2)
            for b in range(256)] for a in range(256)]
BYTE_POW = [1, 256, 65536, 16777216]


def x32(a, b):
    return sum(p * XOR_TAB[a // p % 256][b // p % 256] for p in BYTE_POW)


def a32(a, b):
    return sum(p * AND_TAB[a // p % 256][b // p % 256] for p in BYTE_POW)


def n32(a):
    return M32 - 1 - a


def rotr_n(x, n):
    return x // 2 ** n + (x % 2 ** n) * 2 ** (32 - n)


def rotl_n(x, n):
    return rotr_n(x, 32 - n)


def sha_compress_n(h, blk):
    w = [from_digits(blk[4 * i:4 * i + 4], 256) for i in range(16)] + [0] * 48
    for t in range(16, 64):
        s0 = x32(x32(rotr_n(w[t - 15], 7), rotr_n(w[t - 15], 18)), w[t - 15] // 8)
        s1 = x32(x32(rotr_n(w[t - 2], 17), rotr_n(w[t - 2], 19)), w[t - 2] // 1024)
        w[t] = (w[t - 16] + s0 + w[t - 7] + s1) % M32
    a, b, c, d, e, f, g, hh = h
    for t in range(64):
        t1 = (hh + x32(x32(rotr_n(e, 6), rotr_n(e, 11)), rotr_n(e, 25)) +
              x32(a32(e, f), a32(n32(e), g)) + SHA_K[t] + w[t]) % M32
        t2 = (x32(x32(rotr_n(a, 2), rotr_n(a, 13)), rotr_n(a, 22)) +
              x32(x32(a32(a, b), a32(a, c)), a32(b, c))) % M32
        hh, g, f, e, d, c, b, a = g, f, e, (d + t1) % M32, c, b, a, (t1 + t2) % M32
    return [(x + y) % M32 for x, y in zip(h, [a, b, c, d, e, f, g, hh])]


def sha256_n(msg):
    p, h = sha_pad(msg), SHA_H0
    for i in range(0, len(p), 64):
        h = sha_compress_n(h, p[i:i + 64])
    return from_digits(h, M32)


def hmac_sha256_n(key, msg):
    k = key + [0] * (64 - len(key))
    inner = sha256_n([XOR_TAB[b][54] for b in k] + msg)
    return sha256_n([XOR_TAB[b][92] for b in k] + to_digits(inner, 256, 32))


def quarter_n(s, idx):
    i, j, k, l = idx
    a, b, c, d = s[i], s[j], s[k], s[l]
    a = (a + b) % M32; d = rotl_n(x32(d, a), 16)
    c = (c + d) % M32; b = rotl_n(x32(b, c), 12)
    a = (a + b) % M32; d = rotl_n(x32(d, a), 8)
    c = (c + d) % M32; b = rotl_n(x32(b, c), 7)
    s[i], s[j], s[k], s[l] = a, b, c, d


def chacha_block_n(key, ctr, nonce):
    s0 = chacha_init(key, ctr, nonce)
    s = s0[:]
    for _ in range(10):
        for idx in CHACHA_ROUNDS:
            quarter_n(s, idx)
    out = []
    for x, y in zip(s, s0):
        out += to_digits((x + y) % M32, 256, 4)[::-1]
    return out

# =================== cases ==================================================


res = {}


def rep(n, f):
    for _ in range(n):
        f()


bench("app: SHA-256 abc x40", lambda: rep(40, lambda: res.__setitem__("abc", sha256(abc))))
check("app: SHA-256 abc x40", res["abc"] % CHECK_MOD)
bench("app: SHA-256 abc x40 (no-bitwise)", lambda: rep(40, lambda: res.__setitem__("abcN", sha256_n(abc))))
check("app: SHA-256 abc x40 (no-bitwise)", res["abcN"] % CHECK_MOD)
assert res["abc"] == res["abcN"] == int(hashlib.sha256(b"abc").hexdigest(), 16)

bench("app: SHA-256 8x256-byte messages", lambda: res.__setitem__("m", [sha256(m) for m in msgs]))
check("app: SHA-256 8x256-byte messages", sum(res["m"]) % CHECK_MOD)
bench("app: SHA-256 8x256-byte messages (no-bitwise)", lambda: res.__setitem__("mN", [sha256_n(m) for m in msgs]))
check("app: SHA-256 8x256-byte messages (no-bitwise)", sum(res["mN"]) % CHECK_MOD)
assert res["m"] == res["mN"] == [int(hashlib.sha256(bytes(m)).hexdigest(), 16) for m in msgs]

bench("app: HMAC-SHA256 RFC4231 #2 x12", lambda: rep(12, lambda: res.__setitem__("t", hmac_sha256(hmac_key, hmac_msg))))
check("app: HMAC-SHA256 RFC4231 #2 x12", res["t"] % CHECK_MOD)
bench("app: HMAC-SHA256 RFC4231 #2 x12 (no-bitwise)", lambda: rep(12, lambda: res.__setitem__("tN", hmac_sha256_n(hmac_key, hmac_msg))))
check("app: HMAC-SHA256 RFC4231 #2 x12 (no-bitwise)", res["tN"] % CHECK_MOD)
assert res["t"] == res["tN"] == int("5bdcc146bf60754e6a042426089575c75a003f089d2739839dec58b964ec3843", 16)

bench("app: ChaCha20 keystream 48 blocks",
      lambda: res.__setitem__("ks", [chacha_block(chacha_key, c, chacha_nonce) for c in range(1, 49)]))
check("app: ChaCha20 keystream 48 blocks", sum(from_digits(b, 256) for b in res["ks"]) % CHECK_MOD)
bench("app: ChaCha20 keystream 48 blocks (no-bitwise)",
      lambda: res.__setitem__("ksN", [chacha_block_n(chacha_key, c, chacha_nonce) for c in range(1, 49)]))
check("app: ChaCha20 keystream 48 blocks (no-bitwise)", sum(from_digits(b, 256) for b in res["ksN"]) % CHECK_MOD)
assert res["ks"] == res["ksN"]
assert bytes(res["ks"][0]).hex().startswith("10f1e7e4d13b5915500fdd1fa32071c4c7d1f4c733c068030422aa9ac3d46c4e")

msg_bytes = [bytes(m) for m in msgs]
bench("prim: Hash SHA256 8x256-byte messages x100",
      lambda: rep(100, lambda: res.__setitem__("p", [int(hashlib.sha256(m).hexdigest(), 16) for m in msg_bytes])))
check("prim: Hash SHA256 8x256-byte messages x100", sum(res["p"]) % CHECK_MOD)
