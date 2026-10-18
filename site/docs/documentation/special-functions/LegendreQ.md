# LegendreQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LegendreQ[n, x]`**

gives the Legendre function of the second kind Q\_n(x).

**`LegendreQ[n, m, x]`**

gives the associated Legendre function Q\_n^m(x).

**`LegendreQ[n, m, a, x]`**

gives the Legendre function of type a (a in

<details>
<summary>Notes</summary>

{1, 2, 3}, default 1). For integer n the explicit closed form P\_n(x) (Log\[1+x\] - Log\[1-x\])/2 + v\_n(x) is generated; a non-integer order with an inexact argument on the cut (|x| \< 1) evaluates numerically at machine or arbitrary (MPFR) precision, real or complex. D\[LegendreQ\[n,x\],x\] and the origin Series are supported. Listable.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= LegendreQ[2, x]
Out[1]= -3/2 x + (1/2 Log[1 + x] - 1/2 Log[1 - x]) (-1/2 + 3/2 x^2)

In[2]:= LegendreQ[2, 0.5]
Out[2]= -0.818663

In[3]:= LegendreQ[3, 2.]
Out[3]= 0.00487112 - 26.7035*I

In[4]:= LegendreQ[1/2, 0]
Out[4]= -(Gamma[3/4] Sqrt[Pi])/(2 Sqrt[2] Gamma[5/4])
```

### Applications (10)

The second-kind solution carries a logarithm

```mathematica
In[5]:= LegendreQ[0, x]
Out[5]= 1/2 Log[1 + x] - 1/2 Log[1 - x]
```

P_2(x) L(x) plus a polynomial

```mathematica
In[6]:= LegendreQ[2, x]
Out[6]= -3/2 x + (1/2 Log[1 + x] - 1/2 Log[1 - x]) (-1/2 + 3/2 x^2)
```

High precision from the closed form

```mathematica
In[7]:= N[LegendreQ[3, 1/3], 30]
Out[7]= 0.2476922409970481777113045678512
```

Listable, machine kernel

```mathematica
In[8]:= LegendreQ[3, {0.1, 0.5, 0.9}]
Out[8]= {0.626867, -0.198655, -0.66271}
```

Complex argument

```mathematica
In[9]:= LegendreQ[2, 0.3 + 0.4 I]
Out[9]= -0.755379 - 0.7532*I
```

Beyond x = 1 the type-1 value is complex

```mathematica
In[10]:= LegendreQ[3, 2.]
Out[10]= 0.00487112 - 26.7035*I
```

The associated function Q_2^1

```mathematica
In[11]:= LegendreQ[2, 1, 0.5]
Out[11]= 0.729806
```

An exact value at the origin for non-integer order

```mathematica
In[12]:= LegendreQ[1/2, 0]
Out[12]= -(Gamma[3/4] Sqrt[Pi])/(2 Sqrt[2] Gamma[5/4])
```

Non-integer order on the cut, by two Frobenius 2F1 series

```mathematica
In[13]:= N[LegendreQ[2.5, 0.4]]
Out[13]= -0.337214
```

Non-integer order off the cut: stays unevaluated

```mathematica
In[14]:= LegendreQ[1.5, 2.]
Out[14]= LegendreQ[1.5, 2.0]
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

**Algorithm.** `builtin_legendre_q` handles `LegendreQ[n, x]`,
`LegendreQ[n, m, x]` and `LegendreQ[n, m, a, x]`. For the two-argument form it
first tries the **machine gate** `legendre_machine` (machine-scalar arguments
with at least one machine `Real`, real integer-valued order, `x` not an exact
`1`), which calls `sf_machine_legendre_q_c` in
`src/special_functions/sf_machine_plog.c`. That kernel takes integer `n >= 0`
up to degree 1000 and uses Christoffel's form `Q_n = P_n Q_0 - W_{n-1}`, running
`P_n` and the polynomial `W_{n-1}` up the three-term recurrence. Seeding a
recurrence on `Q` itself with `Q_1 = x Q_0 - 1` would cancel near `|x| ~ 0.8`.
`Q_0` is `atanh(x)` on `(-1, 1)`. Off that interval it uses Kahan's
cancellation-free formulas, with the closed form's branch: a real `x > 1`
reads `Log[1 - x] = log(x - 1) + iπ`, so the value is complex. It declines at
`x = ±1`. Otherwise the exact path runs: an integer `n >= 0` (cap
`LEG_POLY_CAP = 2000`; `n < 0` is singular and stays symbolic) builds
`Q_n(x) = P_n(x) L(x) + v_n(x)` with `L(x) = (1/2)(Log[1+x] - Log[1-x])` and
`v_n` from the same recurrence seeded `v_0 = 0`, `v_1 = -1` (`legendre_q_closed`).
A non-integer order gives the exact special value
`Q_v(0) = -(Sqrt[Pi]/2) Sin[v Pi/2] Gamma[(v+1)/2] / Gamma[v/2+1]` at an exact
zero. For an inexact `x` on the cut `|x| < 1` it evaluates the two Frobenius
series `Q_v(0) 2F1(-v/2, (v+1)/2; 1/2; x^2) + Q_v'(0) x 2F1((1-v)/2, (v+2)/2; 3/2; x^2)`,
built as an `Expr` over `Hypergeometric2F1`/`Gamma` (`legendre_q_numeric`), so it
inherits their machine, MPFR and complex numerics. Associated `Q_n^m` (all three
types, integer `n, m >= 0`) differentiate `Q_n` in a fresh dummy variable, then
substitute. Over a machine buffer the polynomial parts are evaluated by Horner
(`leg_poly_from_coeffs`, shared with `LegendreP`).

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays; plain re/im doubles
in the machine kernel; the non-integer numeric path is an `Expr` over `2F1`/`Gamma`.
ND: binary kernel `NDKB_LegendreQ = { ndk_LegendreQ_c, ... }`, which calls the same
`sf_machine_legendre_q_c`, registered `REG_B` (so `packed_aware`). It is not
real-closed: a real element with `|x| > 1` returns complex. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** Machine kernel `O(n)` per point up to degree 1000;
exact recurrence `O(n^2)` big rationals, cap `n = 2000`. The non-integer series
converges only inside `|x| < 1` (the result is kept only if it collapsed to a
number), so `LegendreQ[1.5, 2.]` stays unevaluated. `LegendreQ[n, 1.]` is
`Indeterminate`. `D` for symbolic order is in `calculus/deriv.c`. `Compile[]`
lowers `LegendreQ[n, x]` at scalar and rank-1 array shapes (`Compiled -> True`)
through the same kernel.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/), [Log](../../elementary-functions/Log/), [LegendreP](../../special-functions/LegendreP/)

- DLMF §14 — Legendre functions of the second kind (§14.3, §14.5.22, §14.7, §14.10.3).
- W. Kahan, Branch cuts for complex elementary functions, in The State of the Art in Numerical Analysis (1987) — the cancellation-free Q_0.
- Source: [`src/special_functions/legendre.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/legendre.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_dsolve_m14_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m14_stress.c)
- Tests: [`tests/test_legendre_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_legendre_machine.c)
- Tests: [`tests/test_legendreq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_legendreq.c)

## Notes & additional examples

### Notes

- `LegendreQ[n, x]` is the Legendre function of the second kind `Q_n(x)`
  (type 1). Unlike `P_n` it carries a logarithm even at integer order: the
  result is `P_n(x) L(x) + v_n(x)` with `L(x) = (1/2)(Log[1+x] - Log[1-x])` and
  a polynomial `v_n` (orders up to `n = 2000`). Negative integer orders are
  singular and stay symbolic.
- A machine-precision call with an integer order, real or complex `x`, uses one
  double kernel (`Q_n = P_n Q_0 - W_{n-1}`, with a cancellation-free `Q_0`)
  shared with packed arrays, `NDArray` arguments and `Compile[]`. A real `x`
  with `|x| > 1` gives a complex value, following the branch of the closed form.
  `LegendreQ[n, 1.]` is `Indeterminate`.
- A non-integer order is evaluated through two Frobenius `2F1` series on the cut
  `|x| < 1`. The exact special value `Q_v(0)` is returned for an exact zero
  argument.
- `LegendreQ[n, m, x]` and `LegendreQ[n, m, a, x]` (type `a` in 1, 2, 3) are
  obtained by differentiating `Q_n`, for integer `n, m >= 0`.
- Attributes `Listable`, `NumericFunction`, `Protected`.
- See also: `LegendreP`, `ArcTanh`, `Hypergeometric2F1`.
