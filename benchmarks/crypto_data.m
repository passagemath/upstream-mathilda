(* ==========================================================================
   benchmarks/crypto_data.m -- shared deterministic inputs for group F
   (experiments 100-119, cryptographic number theory).
   ==========================================================================

   Loaded after the harness:   Get["../harness.m"]; Get["../crypto_data.m"];

   The four systems (Mathilda, Mathematica, Python, PARI/GP) cannot draw the
   same random stream, so EVERY group-F input is derived from one 64-bit LCG
   (Knuth's MMIX constants), spelled identically in crypto_data.{m,py,gp}:

       x -> (6364136223846793005 x + 1442695040888963407) mod 2^64

   Only bignum arithmetic is used here -- no bitwise heads -- so the file loads
   in every system whatever its bitwise coverage.  Hex constants go through
   FromDigits[..., 16] because Mathilda does not parse the 16^^ff base-literal
   form.
   ========================================================================== *)

lcgNext[x_] := Mod[6364136223846793005 x + 1442695040888963407, 2^64];

(* n successive states after seed (the seed itself excluded). *)
lcgStates[seed_, n_] := Rest[NestList[lcgNext, seed, n]];

(* n small values in [0, m), taken from the HIGH 32 bits of each state (the
   low bits of a power-of-two LCG have short periods). *)
smallList[seed_, n_, m_] := Mod[Quotient[lcgStates[seed, n], 2^32], m];

(* A bits-bit integer with its top bit set: ceil(bits/64) states, least
   significant word first. *)
bigOf[seed_, bits_] := Module[{s = lcgStates[seed, Ceiling[bits/64]]},
  Mod[FromDigits[Reverse[s], 2^64], 2^(bits - 1)] + 2^(bits - 1)];

(* The first prime strictly above bigOf[seed, bits]. *)
primeOf[seed_, bits_] := NextPrime[bigOf[seed, bits]];

(* The check reducer: every big-integer check is folded to this residue so the
   value prints identically everywhere. *)
$CheckMod = 1000000007;

(* ---- fixed group / curve constants ---- *)

(* RFC 3526 group 14: the 2048-bit MODP safe prime, generator 2. *)
$MODP2048 = FromDigits["ffffffffffffffffc90fdaa22168c234c4c6628b80dc1cd129024e088a67cc74020bbea63b139b22514a08798e3404ddef9519b3cd3a431b302b0a6df25f14374fe1356d6d51c245e485b576625e7ec6f44c42e9a637ed6b0bff5cb6f406b7edee386bfb5a899fa5ae9f24117c4b1fe649286651ece45b3dc2007cb8a163bf0598da48361c55d39a69163fa8fd24cf5f83655d23dca3ad961c62f356208552bb9ed529077096966d670c354e4abc9804f1746c08ca18217c32905e462e36ce3be39e772c180e86039b2783a2ec07a28fb5c55df06f4c52c9de2bcbf6955817183995497cea956ae515d2261898fa051015728e5a8aacaa68ffffffffffffffff", 16];

(* secp256k1: y^2 = x^3 + 7 over F_p, base point G of prime order n. *)
$Secp256k1P  = 2^256 - 2^32 - 977;
$Secp256k1N  = FromDigits["fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16];
$Secp256k1Gx = FromDigits["79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16];
$Secp256k1Gy = FromDigits["483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16];
