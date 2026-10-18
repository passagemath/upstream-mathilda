# LerchPhi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LerchPhi[z, s, a]`**

is the Lerch transcendent Phi(z, s, a) = Sum\_{k\>=0} z^k/(k + a)^s.

**`Zeta[s, a]`**

and z LerchPhi\[z, s, 1\] is PolyLog\[s, z\]. Exact reductions

<details>
<summary>Notes</summary>

It generalizes Zeta, HurwitzZeta and PolyLog: LerchPhi\[1, s, a\] is cover z = 0 (a^-s), s = 0 (1/(1-z)), z = +-1, positive integer a (a PolyLog form) and negative integer s (a rational function of z). The options DoublyInfinite -\> True (sum k from -Infinity to Infinity) and IncludeSingularTerm -\> True (keep the k + a = 0 term) are supported. Inexact arguments with |z| \< 1 evaluate numerically at machine or arbitrary (MPFR) precision; |z| \> 1 stays symbolic. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= LerchPhi[z, s, 1]
Out[1]= PolyLog[s, z]/z

In[2]:= LerchPhi[0.5, 3, 2.5]
Out[2]= 0.0794983

In[3]:= LerchPhi[-1, 2, 1/2]
Out[3]= 4 Catalan
```

### Applications (10)

At z = 1 the Lerch transcendent is a Hurwitz zeta

```mathematica
In[4]:= LerchPhi[1, 2, 1]
Out[4]= 1/6 Pi^2
```

A negative integer s gives a rational function

```mathematica
In[5]:= LerchPhi[z, -1, 1]
Out[5]= 1/(1 - z)^2
```

The alternating half-integer sum is 4 times Catalan's constant

```mathematica
In[6]:= LerchPhi[-1, 2, 1/2]
Out[6]= 4 Catalan
```

A = 0 drops the singular term and reduces to PolyLog[2, 1/2]

```mathematica
In[7]:= LerchPhi[1/2, 2, 0]
Out[7]= -1/2 Log[2]^2 + 1/12 Pi^2
```

A positive integer a shifts onto PolyLog

```mathematica
In[8]:= LerchPhi[2, 3, 1]
Out[8]= 1/2 PolyLog[3, 2]
```

Machine reals on the double-double series kernel

```mathematica
In[9]:= LerchPhi[{0.1, 0.5, 0.9}, 2.0, 1.5]
Out[9]= {0.461314, 0.554291, 0.763093}
```

Complex z on the machine kernel

```mathematica
In[10]:= LerchPhi[0.5 I, 2.0, 1.5]
Out[10]= 0.425862 + 0.0744738*I
```

Arbitrary precision via arb

```mathematica
In[11]:= N[LerchPhi[1/3, 3/2, 2/5], 40]
Out[11]= 4.1916478522421557664542573036590861911793
```

|z| > 1 off the cut: the Erdelyi continuation

```mathematica
In[12]:= LerchPhi[2.0 + 1.0 I, 1.5, 0.5]
Out[12]= 2.69391 + 1.45599*I
```

On the cut with integer s there is no continuation, so it stays unevaluated

```mathematica
In[13]:= LerchPhi[2.0, 2.0, 0.5]
Out[13]= LerchPhi[2.0, 2.0, 0.5]
```

## Algorithm

Mathilda -- the Lerch transcendent LerchPhi.

```text
  LerchPhi[z, s, a]   Phi(z, s, a) = Sum_{k>=0} z^k / (k + a)^s
                      (|z| < 1; analytic continuation elsewhere, branch cut
                       z in [1, Infinity))
```

For Re(a) < 0 the principal value uses the symmetric power ((k+a)^2)^(-s/2),

```text
and any term with k + a = 0 is excluded.  LerchPhi is the common
```

generalization of Zeta, HurwitzZeta and PolyLog:

```text
  Phi(1, s, a) = Zeta[s, a],      z Phi(z, s, 1) = PolyLog[s, z].
```

The evaluator routes each kind of argument to the cheapest exact or fastest numeric path (mirroring src/special_functions/{zeta,hurwitzzeta,polylog}.c):

```text
  exact reductions (any z, s, a):
      z = 0                  ->  a^-s                        (the k = 0 term)
      s = 0                  ->  1/(1 - z)                   (geometric sum)
      z = 1                  ->  Zeta[s, a]
      z = -1                 ->  2^-s (Zeta[s,a/2] - Zeta[s,(a+1)/2])
      a positive integer m   ->  z^-m (PolyLog[s,z] - Sum_{j<m} z^j j^-s)
      s negative integer -n  ->  (z d/dz + a)^n [1/(1-z)]    (rational in z)
  options:
      IncludeSingularTerm->True at a non-positive integer a  ->  ComplexInfinity
      DoublyInfinite->True   ->  Phi(z,s,a) + z^-1 Phi(1/z, s, 1-a)
  numeric (>= 1 inexact operand):
      |z| < 1                ->  complex-MPFR power series
      |z| > 1                ->  diverges, stays symbolic (no continuation)
  everything else            ->  stays symbolic (return NULL)
```

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Zeta](../../special-functions/Zeta/), [HurwitzZeta](../../special-functions/HurwitzZeta/), [PolyLog](../../special-functions/PolyLog/), [NDArray](../../linear-algebra/NDArray/)

- DLMF §25.14 — the Lerch transcendent.
- A. Erdelyi et al., Higher Transcendental Functions, Vol. I, §1.11(8) — the large-z continuation.
- Source: [`src/special_functions/lerchphi.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/lerchphi.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_integrate_intrep.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_intrep.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)
- Tests: [`tests/test_lerchphi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_lerchphi.c)

## Notes & additional examples

### Notes

- `LerchPhi[z, s, a] = Σ_{k≥0} zᵏ/(k+a)ˢ` generalizes `Zeta`, `HurwitzZeta`
  and `PolyLog`: `LerchPhi[1, s, a] = Zeta[s, a]` and
  `z LerchPhi[z, s, 1] = PolyLog[s, z]`. For `Re a < 0` it uses the symmetric
  power `((k+a)²)^(−s/2)` and drops any `k + a = 0` term, as Mathematica does.
- Many shapes reduce exactly: `z = 0`, `s = 0`, `z = ±1`, integer `a`, and
  negative integer `s`.
- Options: `DoublyInfinite -> True` gives the two-sided sum.
  `IncludeSingularTerm -> True` keeps the `k + a = 0` term, which makes the
  result `ComplexInfinity` at a non-positive integer `a`.
- Numeric surfaces: for `|z| < 1`, machine-precision real and complex arguments
  run on one `double` kernel, shared by scalar, packed/`NDArray` and scalar
  `Compile[]`. Rank-1 `Compile[]` does not lower. Higher precision with
  `|z| ≤ 1` uses arb's `acb_dirichlet_lerch_phi`. For `|z| > 1`, the Erdélyi
  continuation applies only within its domain.
- Known gap: values on the negative real axis with `z < -1` (e.g.
  `LerchPhi[-3.0, 1.5, 0.5]`) are currently inaccurate, so do not rely on them.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `Zeta`, `HurwitzZeta`, `PolyLog`.
