# Zeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Zeta[s]`**

is the Riemann zeta function zeta(s) = Sum\_{k\>=1} k^-s.

**`Zeta[s, a]`**

is the Hurwitz zeta function zeta(s, a) = Sum\_{k\>=0} (k+a)^-s.

<details>
<summary>Notes</summary>

Even positive integers give rational multiples of Pi^(2n), negative integers give rationals, Zeta\[0\] is -1/2, and Zeta\[1\] is ComplexInfinity; odd positive integers stay symbolic. Hurwitz zeta at a positive integer a reduces to Zeta\[s\] minus a finite power sum, and Zeta\[s, 1/2\] is (2^s - 1) Zeta\[s\]. Real, complex, machine and arbitrary-precision (MPFR) numeric arguments evaluate numerically via mpfr\_zeta (real Riemann) or an Euler-Maclaurin kernel. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Zeta[2]
Out[1]= 1/6 Pi^2

In[2]:= Series[Zeta[x], {x, 1, 2}] // Normal
Out[2]= EulerGamma + 1/(-1 + x) - StieltjesGamma[1] (-1 + x) + 1/2 StieltjesGamma[2] (-1 + x)^2

In[3]:= N[Zeta[2 + 3 I], 30]
Out[3]= 0.7980219851462757206222945007251 - 0.1137443080529385002159133658573*I

In[4]:= Zeta[-2.5, -0.5]
Out[4]= 0.169765
```

### Applications (9)

Even positive integers: a rational multiple of a power of Pi

```mathematica
In[5]:= Zeta[6]
Out[5]= 1/945 Pi^6
```

Negative integers give -BernoulliB[n+1]/(n+1); the even ones are the trivial zeros

```mathematica
In[6]:= Table[Zeta[-n], {n, 1, 6}]
Out[6]= {-1/12, 0, 1/120, 0, -1/252, 0}
```

Odd values stay symbolic; N reaches them through arb

```mathematica
In[7]:= N[Zeta[5], 50]
Out[7]= 1.03692775514336992633136548645703416805708091950191
```

The first nontrivial zero, to machine precision

```mathematica
In[8]:= Zeta[0.5 + 14.134725141734693 I]
Out[8]= 1.16675e-16 - 7.32888e-16*I
```

A complex argument at 30 digits

```mathematica
In[9]:= N[Zeta[2 + 3 I], 30]
Out[9]= 0.7980219851462757206222945007251 - 0.1137443080529385002159133658573*I
```

A machine-real list runs through the packed kernel

```mathematica
In[10]:= Zeta[{0.5, 2.0, 3.0, -1.5}]
Out[10]= {-1.46035, 1.64493, 1.20206, -0.0254852}
```

The Laurent expansion at the pole defines the Stieltjes constants

```mathematica
In[11]:= Series[Zeta[s], {s, 1, 2}]
Out[11]= 1/(s - 1) + EulerGamma - StieltjesGamma[1] (s - 1) + 1/2 StieltjesGamma[2] (s - 1)^2 + O[s - 1]^3
```

Hurwitz form at an integer shift: Zeta[4] minus a finite power sum

```mathematica
In[12]:= Zeta[4, 5]
Out[12]= -22369/20736 + 1/90 Pi^4
```

The two-argument form uses ((k+a)^2)^(-s/2), so it is real for negative a

```mathematica
In[13]:= {Zeta[-2.5, -0.5], Zeta[3, -1.5]}
Out[13]= {0.169765, 16.7107}
```

## Algorithm

Mathilda -- the Riemann and Hurwitz zeta functions.

```text
  Zeta[s]      Riemann zeta      zeta(s)   = Sum_{k>=1} k^-s          (Re s > 1)
  Zeta[s, a]   Hurwitz zeta      zeta(s,a) = Sum_{k>=0} (k+a)^-s      (Re s > 1)
```

Both are defined elsewhere by analytic continuation; the evaluator routes each kind of argument to the cheapest exact or fastest numeric path:

```text
  exact integer s         ->  closed form:
                                s = 1        : ComplexInfinity (pole)
                                s = 0        : -1/2
                                s = 2n > 0   : rational * Pi^(2n)   (Bernoulli)
                                s = -m < 0   : rational            (Bernoulli)
                                s = 2n+1 > 0 : stays symbolic (no closed form)
  exact Hurwitz, integer a -> Zeta[s] - Sum_{k=1}^{a-1} k^-s
  machine / MPFR real s    -> MPFR mpfr_zeta (Riemann only)
  complex s, or any a != 1 -> Euler-Maclaurin complex-MPFR kernel
  everything else          -> stays symbolic (return NULL)
```

MPFR provides mpfr_zeta for real Riemann zeta only -- it has no Hurwitz and no complex zeta -- so the Hurwitz / complex kernel is implemented here from the Euler-Maclaurin summation formula (DLMF 25.11.5):

```text
  zeta(s,a) = Sum_{k=0}^{N-1} (a+k)^-s
            + (a+N)^(1-s)/(s-1)
            + 1/2 (a+N)^-s
            + Sum_{j>=1} B_{2j}/(2j)! (s)_{2j-1} (a+N)^(-s-2j+1)
```

with (s)_{2j-1} the rising factorial. N is chosen from the working precision and |s|; the correction series is truncated at its optimal (smallest) term.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_zeta` evaluates `Zeta[s]` (Riemann) and `Zeta[s, a]`
(Hurwitz) in a fixed dispatch order. (1) An `Interval` argument goes through
`interval_apply_function`. (2) **Machine-scalar gate**: when
`sf_machine_scalar_args` classifies the call as machine precision (some machine
`Real`, no MPFR leaf), the raw-`double` kernels of `sf_machine_zeta.c` answer —
`sf_zeta_r` / `sf_zeta_c` for one argument, `sf_hurwitz_r` / `sf_hurwitz_c` for
two (which only accept `Re a > 0`, where the symmetric and principal powers
coincide). `sf_zeta_r` reads tabulated ζ(2…64) and ζ(0…−60) at integers, sums the
Dirichlet series directly for `s ≥ 20`, uses Euler–Maclaurin for `s ≥ 1/2`, and
the functional equation (DLMF 25.4.1, with sin(πs/2) reduced exactly and the
rounding of `1 − s` corrected to first order in Γ and ζ) below 1/2; complex
`Re s < 1/2` declines. A decline (pole, domain, overflow) falls through. (3)
**Exact** one-argument closed forms: `1 → ComplexInfinity`, `0 → -1/2`, even
`2n > 0 →` rational·π²ⁿ from an exact Bernoulli cache, negative `−m →` rational
(zero at even `m`), odd positive integers stay symbolic, `Zeta[Infinity] = 1`.
(4) **Arbitrary precision** (an MPFR real, or a complex with an MPFR part) goes
to arb's `acb_dirichlet_zeta` via `flint_num_zeta`; then `mpfr_zeta` for any
remaining real `s`, and the complex-MPFR Euler–Maclaurin kernel
`zeta_hurwitz_cx` (with `a = 1`) for a complex `s` the gate declined. The
two-argument form reduces exactly at `a = 1`, `a = 1/2` (`(2^s − 1) Zeta[s]`)
and positive integer `a` (`Zeta[s] − Σ_{k<a} k⁻ˢ`); any other inexact input runs
`zeta_numeric` (complex-MPFR Euler–Maclaurin), which uses the **symmetric** power
`((a+k)²)^(−s/2)` so `Re a < 0` gives Mathematica's real value.

**Data structures.** `Expr`; an exact `mpq` Bernoulli cache for the closed forms;
the file-local `zcx` (`mpfr_t` re/im) toolkit for the Euler–Maclaurin kernel;
`double` tables (`ZM_ZPOS`, `ZM_ZNEG`, `ZM_B2K_FACT`) and TwoSum/double-double
helpers in the machine kernel; complex machine results share one `Complex` head
node (`sf_machine_complex_result`). ND: unary kernel `NDKU_Zeta` with both a real
arm (`sf_machine_zeta`) and a complex arm (`sf_zeta_c`), registered `REG_U`, so
packed and visible-`NDArray` arguments stay on the buffer. Attributes:
`Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Machine values cost well under a microsecond (one table
read or a short Euler–Maclaurin sum, `w ≥ max(10, |s|/2 + 4)`, correction cut at
2⁻⁵⁶). Arbitrary-precision cost is arb's. Exact-integer cap
`ZETA_EXACT_INT_CAP = 10000`; Hurwitz-integer-`a` cap `100000`. Complex
`Re s < 1/2` at machine precision is computed on the MPFR path (slower but
accurate). `Compile[]` lowers `Zeta[s]` at scalar and rank-1 array shapes, real
and complex (`Compiled -> True`); the two-argument `Zeta[s, a]` does **not**
compile — use `HurwitzZeta`, which does (identical for `Re a > 0`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [HurwitzZeta](../../special-functions/HurwitzZeta/)

- DLMF §25.2 — the Riemann zeta function; §25.4.1 — the functional equation; §25.11 — the Hurwitz zeta function.
- DLMF §25.11.5 / §25.2.9 — the Euler-Maclaurin summation formula used by the machine and complex-MPFR kernels.
- Source: [`src/special_functions/zeta.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/zeta.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dirichletl.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dirichletl.c)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)

## Notes & additional examples

### Notes

- `Zeta[s]` is the Riemann zeta function, `Zeta[s, a]` the generalized zeta
  Σ_{k≥0} ((k+a)²)^(−s/2). That is Mathematica's convention: for `Re a > 0` it
  equals `HurwitzZeta[s, a]`, and for negative real `a` it stays real where the
  principal-branch `HurwitzZeta` has poles or complex values.
- Exact closed forms: `Zeta[0] = -1/2`, `Zeta[1] = ComplexInfinity`, even positive
  integers give rational·π²ⁿ, negative integers give rationals, and
  `Zeta[Infinity] = 1`. Odd positive integers stay symbolic.
- Numeric surfaces: machine-precision real and complex `s` (with `Re s ≥ 1/2` for
  complex `s`) use one `double` kernel, shared bit-for-bit by the scalar,
  packed/`NDArray` and `Compile[]` paths. Arbitrary precision uses arb
  (`acb_dirichlet_zeta`) for `Zeta[s]` and a complex-MPFR Euler–Maclaurin kernel
  for `Zeta[s, a]`.
- `Compile[]` lowers `Zeta[s]` (scalar and array, real and complex) but not the
  two-argument form. Inside compiled code, use `HurwitzZeta[s, a]` instead.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `HurwitzZeta`, `LerchPhi`, `PolyLog`, `StieltjesGamma`, `BernoulliB`.
