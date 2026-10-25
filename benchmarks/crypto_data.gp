\\ benchmarks/crypto_data.gp -- shared deterministic inputs for group F.
\\ The PARI/GP twin of crypto_data.m / crypto_data.py; see crypto_data.m.
\\ Loaded after the harness:  read("../harness.gp"); read("../crypto_data.gp");

lcgNext(x) = (6364136223846793005 * x + 1442695040888963407) % 2^64;

\\ n successive states after seed (the seed itself excluded).
lcgStates(seed, n) = my(v = vector(n), x = seed); for (i = 1, n, x = lcgNext(x); v[i] = x); v;

\\ n values in [0, m) from the high 32 bits of each state.
smallList(seed, n, m) = apply(s -> (s \ 2^32) % m, lcgStates(seed, n));

\\ A bits-bit integer with its top bit set; least significant word first.
bigOf(seed, bits) =
{
  my(s = lcgStates(seed, ceil(bits / 64)));
  fromdigits(Vecrev(s), 2^64) % 2^(bits - 1) + 2^(bits - 1);
}

\\ The first prime strictly above bigOf(seed, bits) (nextprime is >=, hence +1).
primeOf(seed, bits) = nextprime(bigOf(seed, bits) + 1);

CHECKMOD = 1000000007;

MODP2048 = 0xffffffffffffffffc90fdaa22168c234c4c6628b80dc1cd129024e088a67cc74020bbea63b139b22514a08798e3404ddef9519b3cd3a431b302b0a6df25f14374fe1356d6d51c245e485b576625e7ec6f44c42e9a637ed6b0bff5cb6f406b7edee386bfb5a899fa5ae9f24117c4b1fe649286651ece45b3dc2007cb8a163bf0598da48361c55d39a69163fa8fd24cf5f83655d23dca3ad961c62f356208552bb9ed529077096966d670c354e4abc9804f1746c08ca18217c32905e462e36ce3be39e772c180e86039b2783a2ec07a28fb5c55df06f4c52c9de2bcbf6955817183995497cea956ae515d2261898fa051015728e5a8aacaa68ffffffffffffffff;

SECP256K1_P  = 2^256 - 2^32 - 977;
SECP256K1_N  = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141;
SECP256K1_GX = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798;
SECP256K1_GY = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8;
