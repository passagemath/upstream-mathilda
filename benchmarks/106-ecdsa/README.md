# 106 — ECDSA sign / verify over secp256k1 (group F)

ECDSA written in each language (`app:`) on the in-language Jacobian point
arithmetic of 105 (copied in, not imported): sign R = k·G, r = x(R) mod n,
s = k⁻¹(z + r·d) mod n; verify X = u1·G + u2·Q via Shamir's trick (one joint
double-and-add with G + Q precomputed). The nonce k is **not** RFC 6979 — it is a
deterministic `bigOf[s, 256] mod n` from `crypto_data`, so every system produces the
same signatures; the message hash z is an integer `bigOf[s, 256]` (portable, no
hashing). The `prim:` cases time the native pieces: SHA-256 of the message, and
sign / verify on a native curve group.

| case | what |
|---|---|
| `app: ECDSA sign x48` | 48 signatures; check Σ(r + s) mod 10⁹+7 |
| `app: ECDSA verify x48 (Shamir)` | verify the 48 signatures; check = number valid (48) |
| `app: ECDSA verify tampered x48 (Shamir)` | verify against z + 1; check = number valid (0) |
| `prim: SHA-256 of 1000 messages` | `Hash[s, "SHA256"]` / `hashlib.sha256` over `"message i"`; **probe-SKIP** in Mathilda (`Hash[s, "SHA256"]` returns unevaluated); PARI has no SHA (`skip`) |
| `prim: ECDSA sign x48` | PARI `ellmul` + `Mod^-1`; Python SKIP (OpenSSL picks its own nonce; coincurve, which accepts one, is not installed); **ABSENT** in Mathilda (no EC head) |
| `prim: ECDSA verify x48` | PARI `elladd(ellmul, ellmul)`; Python `cryptography` (OpenSSL) prehashed verify of the app signatures; **ABSENT** in Mathilda |

Baselines: Python = plain ints + built-in `pow`, `hashlib`, `cryptography`;
PARI = plain integers with `%`, `Mod(x, m)^-1`, `binary`, `ellmul` / `elladd`.

## Findings

Measured (ms; Mathilda / Python / PARI / Mathematica):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| app: ECDSA sign x48 | 1005 | 92 | 105 | 620 |
| app: ECDSA verify x48 (Shamir) | 1668 | 122 | 129 | 1061 |
| app: ECDSA verify tampered x48 (Shamir) | 3386 | 127 | 129 | 774 |
| prim: SHA-256 of 1000 messages | SKIP | 1.1 | SKIP | 485 |
| prim: ECDSA sign x48 | ABSENT | SKIP | 27 | ABSENT |
| prim: ECDSA verify x48 | ABSENT | 38 | 55 | ABSENT |

All checks agree in all four systems; the app signatures verify under OpenSSL.

1. **`Hash[s, "SHA256"]` returns unevaluated** (as do `"MD5"`, `"SHA"`, and the
   3-argument `"HexString"` form); only the default `Hash[s]` evaluates. Reproducer:
   `Hash["abc", "SHA256"]` → Mathilda `Hash[abc, SHA256]`; Mathematica
   `84342368487090800366523834928142263660104883695016514377462985829716817089965` (= `int(hashlib.sha256(b"abc").hexdigest(), 16)`).
2. **ABSENT: no elliptic-curve group head** (WL has none either).
3. **Progressive slowdown from unreclaimed `Module` temporaries** (see 105,
   finding 2): the tampered-verify case does exactly the work of the valid-verify
   case, but runs later in the same process and costs 2× (3386 vs 1668 ms); in
   Mathematica the two are equal. App ECDSA is 11–27× Python overall.
