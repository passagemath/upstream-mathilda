#!/usr/bin/env python3
"""Experiment 113 -- lattice attacks (group F).  See lattice_attacks.m.

Baseline: sympy's DomainMatrix.lll() (pure-Python exact LLL over ZZ,
delta = 3/4, the same parameter as Mathilda's LatticeReduce) -- the lattice
reduction a Python cryptanalysis script reaches for without fpylll.  Basis
vectors are ROWS, as in LatticeReduce.  Everything around the reduction
(key generation, CJLOSS basis, recovery scan, Merkle-Hellman decryption) is
plain Python ints + built-in pow().
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require
from crypto_data import big_of, prime_of, small_list, lcg_states, CHECK_MOD

from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

require(["sympy.polys.matrices:DomainMatrix"])


def lll(rows):
    """LLL-reduce the row basis `rows` (list of lists of ints)."""
    m = DomainMatrix([[ZZ(v) for v in r] for r in rows], (len(rows), len(rows[0])), ZZ)
    return [[int(v) for v in r] for r in m.lll().to_list()]


def dot(x, a):
    return sum(i * j for i, j in zip(x, a))


def from_bits(x):
    v = 0
    for b in x:
        v = 2 * v + b
    return v


# ---- input vectors: NOT raw LCG output ------------------------------------
# (see lattice_attacks.m: LCG output is linear and LLL sees it, so each LCG
# word goes through g^e mod P before truncation to a bits-bit integer)

def mix_bigs(seed, cnt, bits):
    p = prime_of(seed, bits + 64)
    return [pow(3, e, p) % (1 << (bits - 1)) + (1 << (bits - 1))
            for e in lcg_states(seed, cnt)]


# ---- (a) Merkle-Hellman knapsack ------------------------------------------

def mh_key(seed, n, k):
    r = mix_bigs(seed, n, k)
    w, acc = [], 0
    for ri in r:
        wi = acc + ri
        acc += wi
        w.append(wi)
    q = prime_of(seed, k + n + 1)
    m = big_of(seed + 1000, k + n - 1) % q
    return w, q, m, [m * wi % q for wi in w]


def mh_message(seed, n):
    return small_list(seed, n, 2)


def mh_decrypt(key, c):
    w, q, m, _ = key
    t = c * pow(m, -1, q) % q
    x = [0] * len(w)
    for i in range(len(w) - 1, -1, -1):
        if t >= w[i]:
            x[i] = 1
            t -= w[i]
    return from_bits(x)


def cjloss_basis(a, s):
    n = len(a)
    rows = [[2 if j == i else 0 for j in range(n)] + [n * a[i]] for i in range(n)]
    rows.append([1] * n + [n * s])
    return rows


def cjloss_recover(a, s, red):
    for v in red:
        if v[-1] == 0 and set(abs(e) for e in v[:-1]) == {1}:
            x = [(1 + e) // 2 for e in v[:-1]]
            if dot(x, a) == s:
                return from_bits(x)
            x = [(1 - e) // 2 for e in v[:-1]]
            if dot(x, a) == s:
                return from_bits(x)
    return -1


def mh_break(a, s):
    return cjloss_recover(a, s, lll(cjloss_basis(a, s)))


key24 = mh_key(11, 24, 56)
key32 = mh_key(13, 32, 40)
cts24 = [dot(mh_message(2000 + j, 24), key24[3]) for j in range(1, 21)]
cts32 = [dot(mh_message(3000 + j, 32), key32[3]) for j in range(1, 13)]
res = {}

bench_once("prim: LatticeReduce CJLOSS knapsack n=24 x20",
      lambda: res.__setitem__("b24", [mh_break(key24[3], s) for s in cts24]))
check("prim: LatticeReduce CJLOSS knapsack n=24 x20", sum(res["b24"]) % CHECK_MOD)

bench_once("prim: LatticeReduce CJLOSS knapsack n=32 x12",
      lambda: res.__setitem__("b32", [mh_break(key32[3], s) for s in cts32]))
check("prim: LatticeReduce CJLOSS knapsack n=32 x12", sum(res["b32"]) % CHECK_MOD)

ctsDec = [dot(mh_message(5000 + j, 32), key32[3]) for j in range(1, 2001)]
bench("app: Merkle-Hellman decrypt n=32 x2000",
      lambda: res.__setitem__("d", [mh_decrypt(key32, c) for c in ctsDec]))
check("app: Merkle-Hellman decrypt n=32 x2000", sum(res["d"]) % CHECK_MOD)

# ---- (b) approximate-GCD / hidden-number secret recovery ------------------

def agcd_samples(seed, t):
    p = prime_of(seed, 128)
    return p, [p * a + b for a, b in zip(mix_bigs(seed + 10, t + 1, 128),
                                         mix_bigs(seed + 500, t + 1, 16))]


def agcd_recover(xs):
    t, x0 = len(xs) - 1, xs[0]
    basis = [[2 ** 17] + xs[1:]]
    for i in range(1, t + 1):
        basis.append([-x0 if j == i else 0 for j in range(t + 1)])
    red = lll(basis)
    v = abs(red[0][0])
    if v % 2 ** 17 == 0 and v > 0:
        return x0 // (v // 2 ** 17)
    return -1


agcd = [agcd_samples(7000 + 1000 * j, 8) for j in range(1, 41)]
bench_once("app: approximate-GCD secret recovery dim=9 x40",
      lambda: res.__setitem__("ps", [agcd_recover(xs) for _, xs in agcd]))
check("app: approximate-GCD secret recovery dim=9 x40",
      sum(res["ps"]) % CHECK_MOD + (0 if res["ps"] == [p for p, _ in agcd] else 1))

# ---- (c) plain LLL on a square knapsack-type basis ------------------------

def qary_basis(seed, n, bits):
    q = prime_of(seed, bits)
    rows = [[q] + [0] * (n - 1)]
    for i, a in enumerate(mix_bigs(seed + 1, n - 1, bits), 1):
        rows.append([a % q] + [1 if j == i else 0 for j in range(1, n)])
    return rows


def absdet(rows):
    return abs(int(DomainMatrix([[ZZ(v) for v in r] for r in rows],
                                (len(rows), len(rows)), ZZ).det()))


qa30 = qary_basis(200000, 30, 300)
qa50 = qary_basis(300000, 50, 200)

bench_once("prim: LatticeReduce q-ary basis n=30 300-bit", lambda: res.__setitem__("r30", lll(qa30)))
check("prim: LatticeReduce q-ary basis n=30 300-bit", absdet(res["r30"]) % CHECK_MOD)

bench_once("prim: LatticeReduce q-ary basis n=50 200-bit", lambda: res.__setitem__("r50", lll(qa50)))
check("prim: LatticeReduce q-ary basis n=50 200-bit", absdet(res["r50"]) % CHECK_MOD)
