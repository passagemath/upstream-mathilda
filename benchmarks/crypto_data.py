"""benchmarks/crypto_data.py -- shared deterministic inputs for group F.

The Python twin of crypto_data.m / crypto_data.gp; see crypto_data.m for the
rationale.  Every input comes from one 64-bit LCG so all four systems compute
on bit-identical integers.
"""

MASK64 = (1 << 64) - 1
CHECK_MOD = 1000000007


def lcg_next(x):
    return (6364136223846793005 * x + 1442695040888963407) & MASK64


def lcg_states(seed, n):
    """n successive states after seed (the seed itself excluded)."""
    out, x = [], seed
    for _ in range(n):
        x = lcg_next(x)
        out.append(x)
    return out


def small_list(seed, n, m):
    """n values in [0, m) from the high 32 bits of each state."""
    return [(s >> 32) % m for s in lcg_states(seed, n)]


def big_of(seed, bits):
    """A bits-bit integer with its top bit set; least significant word first."""
    s = lcg_states(seed, -(-bits // 64))
    v = 0
    for i, w in enumerate(s):
        v |= w << (64 * i)
    return v % (1 << (bits - 1)) + (1 << (bits - 1))


def prime_of(seed, bits):
    """The first prime strictly above big_of(seed, bits)."""
    import gmpy2
    return int(gmpy2.next_prime(big_of(seed, bits)))


MODP2048 = int("ffffffffffffffffc90fdaa22168c234c4c6628b80dc1cd129024e088a67cc74020bbea63b139b22514a08798e3404ddef9519b3cd3a431b302b0a6df25f14374fe1356d6d51c245e485b576625e7ec6f44c42e9a637ed6b0bff5cb6f406b7edee386bfb5a899fa5ae9f24117c4b1fe649286651ece45b3dc2007cb8a163bf0598da48361c55d39a69163fa8fd24cf5f83655d23dca3ad961c62f356208552bb9ed529077096966d670c354e4abc9804f1746c08ca18217c32905e462e36ce3be39e772c180e86039b2783a2ec07a28fb5c55df06f4c52c9de2bcbf6955817183995497cea956ae515d2261898fa051015728e5a8aacaa68ffffffffffffffff", 16)

SECP256K1_P = 2**256 - 2**32 - 977
SECP256K1_N = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
SECP256K1_GX = int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16)
SECP256K1_GY = int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)
