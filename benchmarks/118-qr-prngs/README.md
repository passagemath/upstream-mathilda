# 118 — Quadratic-residuosity cryptography: BBS and Goldwasser–Micali (group F)

Two constructions resting on quadratic residuosity mod a 1024-bit Blum integer
n = pq (p ≡ q ≡ 3 mod 4, found by the same deterministic search from
`bigOf[1181|1182, 512]` stepping by 4 under each system's own primality test):
Blum–Blum–Shub written in the language (one parity bit per step, and the
log₂log₂n = 10-bits-per-step variant), and Goldwasser–Micali (encrypt
c = x² yᵇ mod n with pseudo-square y = n − 1; decrypt by the Legendre symbol of
c mod p, through the native symbol and through Euler's criterion). Plus the raw
symbol kernels on 1024-bit operands. Bit streams are checked as a binary fold
mod 10⁹+7 (decryptions add 10¹⁰ when they equal the message).

| case | what |
|---|---|
| `app: BBS-1024 parity x20000 bits` | x → x² mod n, emit x mod 2, folded |
| `app: BBS-1024 10 bits/step x20000 steps` | same, emit x mod 1024 |
| `app: GM-1024 encrypt x2000 bits x10` | `Mod[x^2 y^b, n]` over 2000 bits, ×10 |
| `prim: JacobiSymbol GM-1024 decrypt x2000 bits x4` | `JacobiSymbol[c, p]` / `gmpy2.jacobi` / `kronecker` |
| `app: GM-1024 decrypt Euler criterion x2000 bits` | `PowerMod[c, (p-1)/2, p] == 1` / `pow` / `Mod(c,p)^e` |
| `prim: JacobiSymbol 1024-bit x5000x4` | Jacobi on 5000 random (a, odd m) pairs |
| `prim: KroneckerSymbol signed/even 1024-bit x5000x4` | signed a, half the moduli even — **ABSENT** in Mathilda |

Baselines: Python app = plain ints + built-in `pow`; Python prim = `gmpy2.jacobi`
/ `gmpy2.kronecker`. PARI = plain integers `x^2 % n`, `Mod(c, p)^e`, and
`kronecker` (PARI's only symbol builtin). No portable twin is timed for
`KroneckerSymbol`: Jacobi is undefined for even moduli.

## Findings

Measured (ms): Mathilda / Python / PARI / Mathematica.

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| app: BBS-1024 parity x20000 bits | 46.9 | 76.3 | 23 | 60.3 |
| app: BBS-1024 10 bits/step x20000 steps | 47.2 | 74.1 | 23 | 60.7 |
| app: GM-1024 encrypt x2000 bits x10 | 58.4 | 112.6 | 23 | 56.3 |
| prim: JacobiSymbol GM-1024 decrypt x2000 bits x4 | 41.7 | 36.2 | 138 | 488 |
| app: GM-1024 decrypt Euler criterion x2000 bits | 108.3 | 1324 | 271 | 634 |
| prim: JacobiSymbol 1024-bit x5000x4 | 173.0 | 168.1 | 766 | 2587 |
| prim: KroneckerSymbol signed/even 1024-bit x5000x4 | ABSENT | 89.1 | 375 | 1405 |

1. **Absent head — `KroneckerSymbol`** (`Names["KroneckerSymbol"] === {}`).
   GMP has `mpz_kronecker`, which `JacobiSymbol` presumably already sits next to.
2. No wrong answers: all checks agree in all four systems. Mathilda is at
   parity with gmpy2 on the symbol kernels and ahead of PARI / Mathematica;
   the in-language BBS loop is ~2× PARI.
