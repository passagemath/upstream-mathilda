# HurwitzZeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HurwitzZeta[s, a]`**

is the Hurwitz zeta function zeta(s, a) = Sum\_{k\>=0} (k + a)^-s.

<details>
<summary>Notes</summary>

Identical to Zeta\[s, a\] for Re(a) \> 0, but built on the principal-branch power (k + a)^-s, so it differs from Zeta for non-positive real a and has poles at a = 0, -1, -2, ... . HurwitzZeta\[s, 1\] is Zeta\[s\], HurwitzZeta\[s, 1/2\] is (2^s - 1) Zeta\[s\], and a positive integer a reduces to Zeta\[s\] minus a finite power sum. A non-positive integer a gives ComplexInfinity for positive integer s and the Bernoulli-polynomial value for non-positive integer s. Real, complex, machine and arbitrary-precision (MPFR) arguments evaluate numerically via an Euler-Maclaurin kernel. Listable.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= HurwitzZeta[s, 1/2]
Out[1]= (-1 + 2^s) Zeta[s]

In[2]:= HurwitzZeta[3, -3.5]
Out[2]= 0.0307784

In[3]:= HurwitzZeta[-3, -2]
Out[3]= -1079/120
```

### Applications (9)

Reduces to Zeta[2] when the second argument is 1

```mathematica
In[4]:= HurwitzZeta[2, 1]
Out[4]= 1/6 Pi^2
```

The half-integer shift: (2^s - 1) Zeta[s]

```mathematica
In[5]:= HurwitzZeta[s, 1/2]
Out[5]= (-1 + 2^s) Zeta[s]
```

Positive integer second argument: Zeta[4] minus a finite power sum

```mathematica
In[6]:= HurwitzZeta[4, 5]
Out[6]= -22369/20736 + 1/90 Pi^4
```

At a non-positive integer a, -BernoulliB[4, -2]/4

```mathematica
In[7]:= HurwitzZeta[-3, -2]
Out[7]= -1079/120
```

Principal powers keep the pole at a = -1

```mathematica
In[8]:= HurwitzZeta[2, -1]
Out[8]= ComplexInfinity
```

Machine reals, threaded over a list

```mathematica
In[9]:= HurwitzZeta[3.0, {0.5, 1.0, 2.5}]
Out[9]= {8.4144, 1.20206, 0.118102}
```

Complex s and a on the machine kernel

```mathematica
In[10]:= HurwitzZeta[2.0 + 1.0 I, 0.5 + 0.5 I]
Out[10]= 1.65298 - 5.03825*I
```

Arbitrary precision via arb

```mathematica
In[11]:= N[HurwitzZeta[1/3, 2/7], 40]
Out[11]= 0.22385839328154108697904688142087606725148
```

Principal versus symmetric powers for negative a

```mathematica
In[12]:= {HurwitzZeta[3, -1.5], Zeta[3, -1.5]}
Out[12]= {0.118102, 16.7107}
```

## Algorithm

Mathilda -- the Hurwitz zeta function.

```text
  HurwitzZeta[s, a]   zeta(s,a) = Sum_{k>=0} (k+a)^-s        (Re s > 1)
```

defined elsewhere by analytic continuation. HurwitzZeta agrees with the two-argument Zeta for Re(a) > 0, but unlike Zeta it sums the *principal branch* powers (k+a)^-s rather than ((k+a)^2)^(-s/2). The consequences:

```text
  - the two functions disagree for non-positive real a, and
  - HurwitzZeta retains the singular summands that Zeta discards, so it has
    poles at a = 0, -1, -2, ... .
```

The evaluator routes each kind of argument to the cheapest exact or fastest numeric path:

```text
  s == 1 (exact)             ->  ComplexInfinity (pole, for any a)
  a == 1                     ->  Zeta[s]          (Riemann closed forms)
  a == 1/2                   ->  (2^s - 1) Zeta[s]
  a positive integer m >= 2  ->  Zeta[s] - Sum_{k=1}^{m-1} k^-s
  a non-positive integer:
      s positive integer     ->  ComplexInfinity (pole)
      s non-positive integer ->  -BernoulliB[1-s, a]/(1-s)   (polynomial)
  any inexact operand        ->  Euler-Maclaurin complex-MPFR kernel
  everything else            ->  stays symbolic (return NULL)
```

MPFR has no Hurwitz zeta, so the numeric kernel is implemented here from the Euler-Maclaurin summation formula (DLMF 25.11.5):

```text
  zeta(s,a) = Sum_{k=0}^{N-1} (a+k)^-s
            + (a+N)^(1-s)/(s-1)
            + 1/2 (a+N)^-s
            + Sum_{j>=1} B_{2j}/(2j)! (s)_{2j-1} (a+N)^(-s-2j+1)
```

with (s)_{2j-1} the rising factorial. N grows with the working precision and

```text
|s|; the correction series is truncated at its optimal (smallest) term. The
```

kernel uses the principal branch for every (a+k)^-s, which is exactly the HurwitzZeta convention. (The structure mirrors src/special_functions/zeta.c; the self-contained Bernoulli cache and complex-MPFR toolkit are replicated here so the two files stay independent.)

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_hurwitzzeta` handles `HurwitzZeta[s, a] = Sum_{k>=0} (k+a)^{-s}`, summing *principal-branch* powers (so, unlike two-argument `Zeta`, it keeps poles at `a = 0, -1, -2, …`). Dispatch order: (1) the **machine-scalar gate** `sf_machine_scalar_args` — for a machine-precision call, `sf_hurwitz_r` (real) or `sf_hurwitz_c` (complex `s` and/or `a`) from `sf_machine_zeta.c` answers when `s ≠ 1`, `Re s ≥ 1/2` and `Re a > 0`: Euler–Maclaurin with each `k + a` formed by an exact TwoSum and its power corrected to first order, declining a complex sum whose terms exceed 16× its value; `a = 1` is delegated to the Riemann kernel. (2) **Exact reductions** when neither operand is inexact: `s == 1 → ComplexInfinity`; `a == 1 → Zeta[s]`; `a == 1/2 → (2^s − 1) Zeta[s]`; positive integer `a = m → Zeta[s] − Σ_{k=1}^{m−1} k^{−s}` (cap `HZ_HURWITZ_A_CAP = 100000`); a non-positive integer `a` gives `ComplexInfinity` for positive integer `s` and the entire value `−BernoulliB[1−s, a]/(1−s)` for non-positive integer `s`. (3) **Numeric** (`hz_numeric`): above machine precision arb's `acb_dirichlet_hurwitz` via `flint_num_hurwitz_zeta`; otherwise (and on an arb decline) the complex-MPFR Euler–Maclaurin kernel `hz_hurwitz_cx` (head sum, `(a+N)^{1−s}/(s−1)` and `½(a+N)^{−s}` corrections, Bernoulli tail truncated at its optimal term). A head term `k + a = 0` with `Re s > 0` signals the pole and returns `ComplexInfinity`.

**Data structures.** `Expr`; a self-contained `mpq_t` Bernoulli cache and the alias-safe complex-MPFR toolkit `hcx`; the shared `double` kernels and `ZM_B2K_FACT` table in `sf_machine_zeta.c`. ND: binary kernel `NDKB_HurwitzZeta` (`REG_B`) whose real and complex arms both call `sf_hurwitz_c`, so packed and visible-`NDArray` arguments are evaluated on the buffer, bit-identical to the scalar path. `Compile[]` lowers `HurwitzZeta[s, a]` at scalar and rank-1 array shapes, real and complex (`Compiled -> True`).

**Complexity / limits.** Machine values ≈ 0.3 µs real, ≈ 1 µs complex. The Euler–Maclaurin head-term count and correction depth grow with precision and `|s|`; arbitrary-precision cost is arb's. `Re s < 1/2` at machine precision is computed on the MPFR path. Known gap: for real `s` and negative non-integer real `a` with non-integer `s` (where the principal value is complex) the MPFR path returns only the real part. Precision follows the minimum over inexact leaves, floored at machine. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Zeta](../../special-functions/Zeta/), [NDArray](../../linear-algebra/NDArray/)

- DLMF §25.11 — the Hurwitz zeta function ζ(s,a) = Sum_{k>=0} (k+a)^{-s} and the Euler–Maclaurin summation formula 25.11.5.
- Source: [`src/special_functions/hurwitzzeta.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hurwitzzeta.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_hurwitzzeta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hurwitzzeta.c)
- Tests: [`tests/test_integrate_intrep.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_intrep.c)

## Notes & additional examples

### Notes

- `HurwitzZeta[s, a] = Σ_{k≥0} (k+a)^(−s)`, continued to all `s ≠ 1`. Unlike
  two-argument `Zeta`, it uses *principal-branch* powers. It therefore keeps
  poles at `a = 0, -1, -2, …` and differs from `Zeta[s, a]` for negative real `a`.
- Exact reductions: `s == 1` gives `ComplexInfinity`, `a == 1` gives `Zeta[s]`,
  `a == 1/2` gives `(2^s - 1) Zeta[s]`, and a positive integer `a` gives
  `Zeta[s]` minus a finite power sum. At a non-positive integer `a`, the result
  is `ComplexInfinity` for positive integer `s` and the Bernoulli-polynomial
  value for non-positive integer `s`.
- Numeric surfaces: machine-precision real and complex arguments with
  `Re s ≥ 1/2` and `Re a > 0` use one Euler–Maclaurin `double` kernel, shared by
  the scalar, packed/`NDArray` and `Compile[]` paths (scalar and array, real and
  complex). Higher precision goes to arb's `acb_dirichlet_hurwitz`, and other
  inputs go to a complex-MPFR Euler–Maclaurin kernel.
- Known gap: for a negative non-integer real `a` with non-integer real `s`
  (e.g. `HurwitzZeta[-2.5, -0.5]`), the principal value is complex, but only its
  real part is returned.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `Zeta`, `LerchPhi`, `PolyGamma`, `BernoulliB`.
