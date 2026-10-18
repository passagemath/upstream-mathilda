# LegendreP

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LegendreP[n, x]`**

gives the Legendre polynomial P\_n(x).

**`LegendreP[n, m, x]`**

gives the associated Legendre function P\_n^m(x).

**`LegendreP[n, m, a, x]`**

gives the Legendre function of type a (a in

<details>
<summary>Notes</summary>

{1, 2, 3}, default 1). Integer n yields the explicit polynomial; a non-integer order with an inexact argument evaluates numerically at machine or arbitrary (MPFR) precision, real or complex. Listable.

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= LegendreP[3, x]
Out[1]= -3/2 x + 5/2 x^3

In[2]:= LegendreP[10, 2, x]
Out[2]= (1 - x^2) (3465/128 - 45045/32 x^2 + 675675/64 x^4 - 765765/32 x^6 + 2078505/128 x^8)

In[3]:= LegendreP[5, 0.3 + 0.2 I]
Out[3]= 0.594236 - 0.0178925*I

In[4]:= LegendreP[2.5, 2]
Out[4]= 9.58312

In[5]:= D[LegendreP[n, x], x]
Out[5]= ((1 + n) LegendreP[1 + n, x] + (-1 - n) x LegendreP[n, x])/(-1 + x^2)
```

### Applications (11)

Exact polynomials with rational coefficients

```mathematica
In[6]:= Table[LegendreP[n, x], {n, 0, 4}]
Out[6]= {1, x, -1/2 + 3/2 x^2, -3/2 x + 5/2 x^3, 3/8 - 15/4 x^2 + 35/8 x^4}
```

Negative order by P_{-n-1} = P_n

```mathematica
In[7]:= LegendreP[-3, x]
Out[7]= -1/2 + 3/2 x^2
```

Every P_n equals 1 at the endpoint

```mathematica
In[8]:= LegendreP[n, 1]
Out[8]= 1
```

Orthogonality on the interval

```mathematica
In[9]:= Integrate[LegendreP[2, x] LegendreP[3, x], {x, -1, 1}]
Out[9]= 0
```

High degree by the machine recurrence

```mathematica
In[10]:= LegendreP[200, 0.7]
Out[10]= -0.00273489
```

Complex argument

```mathematica
In[11]:= LegendreP[5, 0.3 + 0.2 I]
Out[11]= 0.594236 - 0.0178925*I
```

The associated function P_2^1

```mathematica
In[12]:= LegendreP[2, 1, x]
Out[12]= -3 x Sqrt[1 - x^2]
```

Associated Legendre over a machine array, by Horner

```mathematica
In[13]:= LegendreP[30, 4, {0.1, 0.5}]
Out[13]= {123672.0, 132751.0}
```

Non-integer order through the Gauss 2F1 series

```mathematica
In[14]:= N[LegendreP[1/3, 1/5], 30]
Out[14]= 0.7848413653720337798085467296758
```

The derivative for symbolic order

```mathematica
In[15]:= D[LegendreP[n, x], x]
Out[15]= ((1 + n) LegendreP[1 + n, x] + (-1 - n) x LegendreP[n, x])/(-1 + x^2)
```

Non-integer order outside |(1-x)/2| < 1: stays unevaluated

```mathematica
In[16]:= LegendreP[0.5, 3.]
Out[16]= LegendreP[0.5, 3.0]
```

## Algorithm

Mathilda -- Legendre polynomials and associated Legendre functions.

```text
  LegendreP[n, x]        Legendre polynomial / function P_n(x).
  LegendreP[n, m, x]     associated Legendre function P_n^m(x) (type 1).
  LegendreP[n, m, a, x]  Legendre function of type a (a in {1, 2, 3}).
```

Evaluation is layered so each argument shape takes the cheapest exact or numeric route:

```text
  LegendreP[n, x]
    exact integer n            ->  the explicit degree-|n'| polynomial in x
                                   (n' = n, or -1-n for n < 0, since
                                   P_{-1-n} = P_n) with exact rational
                                   coefficients, built from the three-term
                                   recurrence; an inexact x then evaluates
                                   the monomials numerically.
    x == 1                     ->  1 (for any order n).
    non-integer n, some arg    ->  numeric Gauss series
       inexact                      P_n(x) = 2F1(-n, n+1; 1; (1-x)/2),
                                   real or complex, machine or MPFR
                                   precision (requires |(1-x)/2| < 1).
    everything else            ->  stays symbolic (return NULL).

  LegendreP[n, m, x] / [n, m, a, x]   (integer n, integer m >= 0)
    type 1 (default, a == 1)   ->  (-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x)
                                   (the Rodrigues derivative form; 0 when
                                   m > |n'|).
    types 2, 3                 ->  C(x) * R_a(x), where
                                   C(x) = 2F1Reg(-n, n+1, 1-m, (1-x)/2)
                                   is the (terminating, exact) regularized
                                   Gauss polynomial and the prefactor is
                                     R_2 = (1+x)^(m/2) (1-x)^(-m/2),
                                     R_3 = (1+x)^(m/2) (-1+x)^(-m/2).
    non-integer / negative m   ->  stays symbolic (return NULL).
```

Attributes: Listable, NumericFunction, Protected.

Deferred (left symbolic): symbolic Series / SeriesCoefficient, D[] rules, the non-integer associated and Legendre-function forms, and analytic continuation of the numeric series for |(1-x)/2| >= 1.

## Implementation notes

**Algorithm.** `builtin_legendre_p` handles `LegendreP[n, x]`,
`LegendreP[n, m, x]` and `LegendreP[n, m, a, x]`. For the two-argument form it
first tries the **machine gate** `legendre_machine`: when both arguments are
machine-scalar leaves with at least one machine `Real` (`sf_machine_scalar_args`),
the order is integer-valued and real, and `x` is not an exact `1`, it calls
`sf_machine_legendre_p_c` (`src/special_functions/sf_machine_plog.c`). That runs the
three-term recurrence `k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}` in `double`
(real `x`, via `sf_machine_legendre_p` in `sf_machine.c`) or in re/im doubles
(complex `x`), maps negative orders with `P_{-n-1} = P_n`, and declines above
degree 1000. Otherwise `legendre_two_arg` does the exact work: an exact integer
order builds the explicit degree-`|n'|` polynomial from the same recurrence with
GMP `mpq` coefficients (`legendre_p_coeffs`, cap `LEG_POLY_CAP = 2000`);
`x == 1 -> 1` for any order; a non-integer order with an inexact argument sums the
Gauss series `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)` in the `ncpx` MPFR-complex toolkit
(`leg_numeric_p`, real or complex, machine or arbitrary precision; requires
`|(1-x)/2| < 1`). The associated forms (integer `n`, integer `m >= 0`): type 1 is
the Rodrigues derivative `(-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x)` (0 when
`m > |n'|`); types 2 and 3 are the regularized Gauss polynomial
`2F1Reg(-n, n+1, 1-m, (1-x)/2)` (`legendre_2f1reg_core`) times
`(1+x)^(m/2) (1-x)^(-m/2)` or `(1+x)^(m/2) (-1+x)^(-m/2)`. For a machine buffer
argument, `leg_poly_from_coeffs` evaluates these polynomials by Horner in `double`
instead of as a sum of monomials, each of which would cost a `pow()` per element.
`D` for symbolic order lives in `calculus/deriv.c`.

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays for the polynomials;
`ncpx` (`mpfr_t` re/im) for the non-integer-order series; plain doubles (no
`<complex.h>`) in the machine kernels. ND: binary kernel
`NDKB_LegendreP = { ndk_LegendreP_c, ... }`, which calls the same
`sf_machine_legendre_p_c`, registered `REG_B` (so `packed_aware`); it threads `n`
and `x` element-wise. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** The machine recurrence is `O(n)` per point, up to degree
1000. The exact recurrence is `O(n^2)` in growing big rationals, capped at
`n = 2000`. The non-integer-order series costs `O(wp)` terms inside the
convergence disc. The Horner buffer path is `O(n)` fused multiply-adds per element
(degree 30 over 10^5 points about 10 ms). Left symbolic: non-integer or negative
`m`, symbolic `Series`, and non-integer order with `|(1-x)/2| >= 1`
(`LegendreP[0.5, 3.]`). `Compile[]` lowers `LegendreP[n, x]` at real and complex
scalar shapes and at rank-1 array shape (`Compiled -> True`), through the same
kernel. The associated three-argument form is not lowered.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/), [LegendreQ](../../special-functions/LegendreQ/), [Series](../../power-series/Series/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

- DLMF §14 — Legendre and related functions (§14.3 series, §14.7 integer degree, §14.10.3 the three-term recurrence).
- Source: [`src/special_functions/legendre.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/legendre.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve_m14_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m14_stress.c)
- Tests: [`tests/test_integer_assume_simplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_assume_simplify.c)
- Tests: [`tests/test_legendre.c`](https://github.com/stblake/mathilda/blob/main/tests/test_legendre.c)

## Notes & additional examples

### Notes

- For an exact integer order, `LegendreP[n, x]` is the explicit Legendre
  polynomial `P_n(x)` with exact rational coefficients (orders up to
  `n = 2000`). `LegendreP[n, 1]` is `1` for every order.
- A machine-precision call with an integer order, real or complex `x`, runs the
  three-term recurrence in doubles (degree up to 1000). That kernel is shared
  with packed arrays, `NDArray` arguments and `Compile[]` (scalar, complex and
  rank-1 array), so these all agree bit for bit.
- A non-integer order is evaluated numerically, when an argument is inexact,
  through `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)` at machine or arbitrary precision.
  This needs `|(1-x)/2| < 1`; outside that disc the call stays unevaluated.
- `LegendreP[n, m, x]` is the associated Legendre function, with the
  Condon–Shortley phase `(-1)^m`. `LegendreP[n, m, a, x]` selects type `a`
  (1, 2 or 3). These need integer `n` and integer `m >= 0`. Over a machine array
  the polynomial is evaluated by Horner, which keeps high degree fast.
- Attributes `Listable`, `NumericFunction`, `Protected`.
- See also: `LegendreQ`, `ChebyshevT`, `ChebyshevU`, `HermiteH`,
  `Hypergeometric2F1`.
