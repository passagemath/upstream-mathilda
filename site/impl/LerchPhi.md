---
source: src/special_functions/lerchphi.c
references:
  - "DLMF §25.14 — the Lerch transcendent."
  - "A. Erdelyi et al., Higher Transcendental Functions, Vol. I, §1.11(8) — the large-z continuation."
---
**Algorithm.** `builtin_lerchphi` evaluates `LerchPhi[z, s, a] = Sum_{k>=0}
z^k/(k+a)^s`, the common generalization of `Zeta`, `HurwitzZeta` and `PolyLog`.
Fewer than three arguments emit `LerchPhi::argrx`. (1) **Machine-scalar gate**:
a three-argument machine-precision call is answered by the shared `double`
kernels of `sf_machine_zeta.c` — `sf_machine_lerchphi` (= `sf_lerch_r`) for real
arguments, `sf_lerch_c` for complex `z`, `s` or `a`. Both sum the defining series
for `|z| < 1` (with `zᵏ` carried in double-double, each `k + a` formed by TwoSum
with a first-order power correction, Neumaier-compensated, declining when the
terms exceed the sum by a cancellation factor), hand `z = 1` to the Hurwitz
kernel, and decline `|z| ≥ 1` or `Re a ≤ 0`. (2) Options `DoublyInfinite` and
`IncludeSingularTerm` are parsed (a non-option emits `LerchPhi::nonopt`). (3) A
packed / visible-`NDArray` argument in one slot is threaded over
`sf_machine_lerchphi` by `lerchphi_ndarray`; on any decline the call is delisted
and re-evaluated element-wise (never left inert). (4) **Exact reductions**
(`lerchphi_single`): `IncludeSingularTerm` at non-positive integer `a →
ComplexInfinity`; non-positive integer `a` drops the singular term and shifts
onto `PolyLog`; `z = 0 → a^-s`; `s = 0 → 1/(1-z)`; `z = 1 → Zeta[s, a]`;
positive integer `a →` the `PolyLog` shift; `z = -1 → 2^-s (Zeta[s, a/2] −
Zeta[s, (a+1)/2])`, with a digamma difference at `s = 1` and a Dirichlet-beta
closed form at half-integer `a`; negative integer `s = -n → (z d/dz + a)^n
[1/(1-z)]` via the Euler operator and `Together`. `DoublyInfinite -> True` is
`Φ(z,s,a) + z⁻¹ Φ(1/z,s,1-a)`. (5) **Numeric** (`lp_numeric`, at least one
inexact operand): for `|z| ≤ 1`, arb's `acb_dirichlet_lerch_phi` above machine
precision (default convention), else the complex-MPFR power series `lp_series_cx`
(symmetric power `((k+a)²)^(−s/2)` for `Re a < 0`); for `|z| > 1`, off the cut
`[1, ∞)`, with `s` not a positive integer, `a` not a non-positive integer and
`|Log z| < 2π`, the Erdélyi large-`z` continuation `lp_erdelyi_large_z` through
`HurwitzZeta` and `Gamma`; otherwise symbolic.

**Data structures.** `Expr`; the `lcx` (`mpfr_t` re/im) complex-MPFR toolkit;
the exact and continuation paths reuse the `PolyLog` / `Zeta` / `HurwitzZeta` /
`Gamma` builtins. ND: n-ary kernel `NDKN_LerchPhi` (arity 3, `REG_N`) plus the
in-builtin `lerchphi_ndarray` threader; `LerchPhi` is on `pack.c`'s `AWARE` list.
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Machine values ≈ 1–2 µs; the MPFR series term cap is
`wp*64 + 100000`, the Erdélyi continuation truncates at `NMAX = 2000` with an
optimal-truncation safeguard. `|z| > 1` stays symbolic for integer `s`, on the
cut `[1, ∞)`, or outside `|Log z| < 2π`. Known gap: the Erdélyi path is
currently inaccurate on the negative real axis `z < -1` (e.g. `LerchPhi[-3.0,
1.5, 0.5]`). `Compile[]` lowers at scalar shape, real and complex (complex
arguments use the n-ary `KERNN_C` opcode), but **not** at rank-1 array shape
(`Compiled -> False`); the ND kernel still serves packed / `NDArray` input at the
REPL.
