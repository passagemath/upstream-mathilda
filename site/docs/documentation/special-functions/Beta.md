# Beta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Beta[a, b]`**

is the Euler beta function B(a, b) = Gamma(a) Gamma(b) / Gamma(a+b).

**`Beta[z, a, b]`**

is the incomplete beta function Integral\_0^z t^(a-1) (1-t)^(b-1) dt.

**`Beta[z0, z1, a, b]`**

is the generalized incomplete beta Beta\[z1, a, b\] - Beta\[z0, a, b\].

<details>
<summary>Notes</summary>

Exact for rational arguments (a positive integer gives a rational via Pochhammer); non-positive integer poles give ComplexInfinity. Machine and arbitrary-precision (MPFR) real and complex inputs evaluate numerically. The incomplete form reduces through Hypergeometric2F1. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Beta[3, 5]
Out[1]= 1/105

In[2]:= Beta[-2, 1]
Out[2]= -1/2

In[3]:= Beta[2.5 + I, 1 - I]
Out[3]= 0.0831078 + 0.142164*I
```

### Applications (10)

A positive-integer argument collapses to an exact rational

```mathematica
In[4]:= Beta[3, 5]
Out[4]= 1/105
```

Half-integers fold to Pi; other rationals keep a Gamma quotient

```mathematica
In[5]:= {Beta[1/2, 1/2], Beta[1/3, 1/3]}
Out[5]= {Pi, Gamma[1/3]^2/Gamma[2/3]}
```

A cancelling pair of gamma poles has a finite limit; a surviving pole is ComplexInfinity

```mathematica
In[6]:= {Beta[-2, 1], Beta[0, b]}
Out[6]= {-1/2, ComplexInfinity}
```

Machine real, computed directly in double

```mathematica
In[7]:= Beta[2.3, 3.2]
Out[7]= 0.0540298
```

Machine complex: the three Gammas under one exponential

```mathematica
In[8]:= Beta[2.5 + I, 1 - I]
Out[8]= 0.0831078 + 0.142164*I
```

Arbitrary precision through the Gamma ratio

```mathematica
In[9]:= N[Beta[22/10, 33/10], 50]
Out[9]= 0.0564856913732825668070517540044914293695377770152413
```

The incomplete beta terminates when b is a positive integer

```mathematica
In[10]:= Beta[1/2, 3, 4]
Out[10]= 7/640
```

A numeric incomplete beta through Hypergeometric2F1

```mathematica
In[11]:= Beta[0.3, 0.5, 2.5]
Out[11]= 0.938813
```

Rank-1 Compile[] runs the same kernel as the packed map

```mathematica
In[12]:= cb = Compile[{{v, _Real, 1}}, Beta[v, 2.5]]; cb[{1.5, 2.5}] == Beta[{1.5, 2.5}, 2.5]
Out[12]= True
```

Symbolic arguments stay unevaluated, as in Mathematica

```mathematica
In[13]:= Beta[a, b]
Out[13]= Beta[a, b]
```

## Algorithm

beta.c -- the Euler beta function and its incomplete / generalized forms.

```text
  Beta[a, b]         = Gamma(a) Gamma(b) / Gamma(a+b)
  Beta[z, a, b]      = Int_0^z t^(a-1) (1-t)^(b-1) dt          (incomplete)
  Beta[z0, z1, a, b] = Beta[z1, a, b] - Beta[z0, a, b]         (generalized)
```

A machine-precision Beta[a, b] (real or complex) is computed directly in double by sfg_beta / sfg_beta_c (sf_machine_gamma.c), the kernel the NDArray map and Compile[] share. Everything else is assembled from functions that already exist and are well tested:

```text
  - Beta[a, b] is reduced to a ratio of Gamma calls and handed back to the
    evaluator. The existing Gamma builtin already closes exact integers,
    half-integers, rationals (-> Sqrt[Pi] forms), machine / arbitrary-
    precision reals, and complex arguments, so Beta inherits all of that for
    free. Non-positive-integer poles (where a, b, or a+b hits a gamma pole)
    are detected up front: a surviving pole gives ComplexInfinity, while a
    cancelling pair of poles reduces by the finite limit of the gamma ratio.

  - Beta[z, a, b] is reduced through Hypergeometric2F1 via
      B_z(a, b) = z^a / a * 2F1(a, 1-b; a+1; z),
    which terminates to an exact closed form when b is a positive integer
    (1-b a non-positive integer) and evaluates numerically (real or complex)
    otherwise.
```

Memory: builtin_beta takes ownership of res. It returns a freshly built tree (the evaluator frees res) or NULL to leave the call unevaluated (the evaluator keeps res). It never frees res itself. Every intermediate tree is consumed by eval_and_free.

## Implementation notes

**Algorithm.** `builtin_beta` checks the machine case before anything else. When `sf_machine_scalar_args` accepts a two-argument call (machine `Real`, exact, or `Complex` leaves, with at least one inexact leaf and no MPFR), it computes the value directly in double:

- `sfg_beta` (real) takes `tgamma(a)/tgamma(a+b)·tgamma(b)`. Because `a + b` is rarely representable, it undoes that rounding to first order through `ψ(a+b)`, using `Γ(s+e) ≈ Γ(s)(1 + ψ(s) e)`.
- `sfg_beta_c` (complex) combines the three Gammas under one exponential, summing their double-double exponents from the `sf_machine_gamma.c` tiles or Stirling series, with the same ψ correction.

A decline falls through to the layered path below: a pole, `|a|`, `|b|` or `|a+b| > 170`, or underflow.

The layered path builds on the `Gamma` machinery. `beta_lattice_case` first sorts out the gamma poles on the integer lattice: a pole that survives gives `ComplexInfinity`, a cancelling pair gives its finite limit (`Beta[-2, 1] = -1/2`), and a pole only in the denominator gives `0`. Purely symbolic arguments stay unevaluated. If one argument is a positive integer `n`, the call collapses to `(n-1)!/Pochhammer[b, n]`, which is exact for rational `b` (`Beta[3, 1/3] = 27/14`). Otherwise the call is rebuilt as `Gamma[a] Gamma[b]/Gamma[a+b]` and evaluated. That keeps half-integers exact (`Beta[1/2, 1/2] = Pi`), keeps other rationals as a Gamma quotient, and sends MPFR and complex arguments to `Gamma`'s arbitrary-precision paths.

The incomplete form uses `Beta[z, a, b] = z^a/a · 2F1(a, 1-b; a+1; z)` through `Hypergeometric2F1`. The series terminates to a closed form when `b` is a positive integer, even for symbolic `z` and `a`. Otherwise it reduces only when every argument is numeric. Two special values are handled directly: `Beta[0, a, b] = 0` and `Beta[1, a, b] = Beta[a, b]`. The four-argument form `Beta[z0, z1, a, b]` evaluates as `Beta[z1, a, b] - Beta[z0, a, b]` when every argument is numeric.

**Data structures.** `Expr` trees driven through `eval_and_free`. The machine kernel works on `(re, im)` doubles in `sf_machine_gamma.c` and shares that file's Gamma and digamma tables. The ND kernel is the binary `NDKB_Beta`, registered with `REG_B`; its `ndk_Beta_c` calls `sfg_beta_c` for real and complex buffers, including the case where both arguments are arrays. So packed arrays and `NDArray`s are handled directly, and the scalar, packed, `NDArray` and `Compile[]` answers agree bit for bit. `Compile[]` lowers `Beta` through that kernel at real and complex scalar shapes, and at rank-1 shapes when one argument is an array and the other a scalar (`Beta[v, 2.5]`). An array in both arguments is not yet lowered.

**Complexity / limits.** The machine kernel is `O(1)`. Measured error: real median 0.8 ulp, max 4.1 ulp on `[0.3, 5]²`; complex median 1.3 ulp, max 6.3 ulp on `[0.5, 3] × ±2i`. The exact path costs what `Gamma`/`Pochhammer` cost, and the incomplete form costs one `2F1` evaluation. Symbolic two-argument calls, and the four-argument form with any symbolic argument, stay unevaluated; symbolic derivatives of the four-argument form are handled in `deriv.c`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Gamma](../../special-functions/Gamma/), [Pi](../../mathematical-constants/Pi/), [NDArray](../../linear-algebra/NDArray/), [Hypergeometric2F1](../../special-functions/Hypergeometric2F1/), [PolyGamma](../../special-functions/PolyGamma/)

- DLMF §5.12 — the Euler beta function B(a,b) = Gamma(a) Gamma(b)/Gamma(a+b).
- DLMF §5.15 — the digamma function psi, used for the first-order rounding correction of a+b.
- DLMF §8.17 — the incomplete beta function B_z(a,b) and its relation to the Gauss hypergeometric function.
- Source: [`src/special_functions/beta.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/beta.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_gamma_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gamma_machine.c)
- Tests: [`tests/test_integrate_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_beta.c)

## Notes & additional examples

### Notes

- `Beta[a, b] = Gamma[a] Gamma[b]/Gamma[a+b]`; `Beta[z, a, b]` is the incomplete beta `∫₀ᶻ t^(a-1) (1-t)^(b-1) dt` (the integration limit comes first), and `Beta[z0, z1, a, b]` the generalized form `Beta[z1, a, b] - Beta[z0, a, b]`.
- Machine-precision `Beta[a, b]` (real or complex) is computed in double by `sfg_beta`/`sfg_beta_c`, correcting the rounding of `a + b` through `PolyGamma`; poles and `|a|, |b|, |a+b| > 170` fall back to the exact/MPFR path.
- Packed and `NDArray` arguments run element-wise — including both arguments arrays — and `Compile[]` lowers scalar and array-with-scalar calls.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Gamma`, `Pochhammer`, `Hypergeometric2F1`.
