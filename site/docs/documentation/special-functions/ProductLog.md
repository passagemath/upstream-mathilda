# ProductLog

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ProductLog[z]`**

gives the principal solution w of z == w e^w (the Lambert W function).

**`ProductLog[k, z]`**

gives the k-th solution (k any integer, k == 0 the

**`ProductLog[E] = 1, ProductLog[-Pi/2] = I Pi/2 and ProductLog[k, 0] =`**

**`D[ProductLog[z], z] = ProductLog[z]/(z (1 + ProductLog[z])). Listable.`**

<details>
<summary>Notes</summary>

principal branch); branches are ordered by imaginary part. ProductLog\[z\] is real for z \>= -1/e and has a branch cut along (-Infinity, -1/e\]. Exact values include ProductLog\[0\] = 0, ProductLog\[-1/E\] = -1, -Infinity for k != 0. Inexact real or complex arguments evaluate numerically at machine or arbitrary (MPFR) precision. Satisfies

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ProductLog[1.0]
Out[1]= 0.567143

In[2]:= ProductLog[-1/E]
Out[2]= -1

In[3]:= ProductLog[-1, -0.2]
Out[3]= -2.54264

In[4]:= ProductLog[-0.5]
Out[4]= -0.794024 + 0.770112*I

In[5]:= Series[ProductLog[x], {x, 0, 4}]
Out[5]= x - x^2 + 3/2 x^3 - 8/3 x^4 + O[x]^5
```

### Applications (12)

Since e = 1 * e^1, the principal value is 1

```mathematica
In[6]:= ProductLog[E]
Out[6]= 1
```

The branch point where the k = 0 and k = -1 branches meet

```mathematica
In[7]:= ProductLog[-1/E]
Out[7]= -1
```

An exact complex value

```mathematica
In[8]:= ProductLog[-Pi/2]
Out[8]= (1/2*I) Pi
```

The omega constant

```mathematica
In[9]:= N[ProductLog[1], 30]
Out[9]= 0.56714329040978387299996866221
```

The other real branch, k = -1

```mathematica
In[10]:= ProductLog[-1, -0.2]
Out[10]= -2.54264
```

Below -1/e the principal value is complex

```mathematica
In[11]:= ProductLog[-0.5]
Out[11]= -0.794024 + 0.770112*I
```

Complex argument, machine complex kernel

```mathematica
In[12]:= ProductLog[1. + 2. I]
Out[12]= 0.823771 + 0.532929*I
```

A higher complex branch

```mathematica
In[13]:= ProductLog[3, 1.5]
Out[13]= -2.44586 + 17.137*I
```

Listable

```mathematica
In[14]:= ProductLog[{0.5, 1., 10.}]
Out[14]= {0.351734, 0.567143, 1.74553}
```

Arbitrary precision on the k = -1 branch

```mathematica
In[15]:= ProductLog[-1, -0.1`40]
Out[15]= -3.5771520639572972184093919635119948804016
```

The Lagrange-inversion series

```mathematica
In[16]:= Series[ProductLog[x], {x, 0, 4}]
Out[16]= x - x^2 + 3/2 x^3 - 8/3 x^4 + O[x]^5
```

```mathematica
In[17]:= D[ProductLog[x], x]
Out[17]= ProductLog[x]/(x (1 + ProductLog[x]))
```

## Algorithm

Mathilda -- ProductLog, the Lambert W function.

```text
  ProductLog[z]     principal branch W_0(z): the solution w of z = w e^w.
  ProductLog[k, z]  the k-th branch W_k(z), k any integer (k == 0 principal).
```

Evaluation is layered so each kind of argument takes the cheapest accurate route:

```text
  exact special values   ->  ProductLog[0] = 0, ProductLog[E] = 1,
                             ProductLog[-1/E] = -1, ProductLog[-Pi/2] = I Pi/2,
                             ProductLog[+-Infinity/ComplexInfinity] = Infinity,
                             ProductLog[k, 0] = -Infinity  (k != 0)
  numeric (real/complex)  ->  unified complex-MPFR Halley core; the result is
                             a Real / MPFR leaf when it is real-valued for the
                             chosen branch, otherwise Complex[..]
  everything else        ->  stays symbolic (return NULL)
```

The numeric core (pl_core) builds on the shared `ncpx` complex-MPFR toolkit (numeric_complex.h). It seeds an initial approximation by region --

```text
  - branch-point series in p = sqrt(2(e z + 1)) near z = -1/e (branches 0,-1);
  - the Maclaurin seed z(1 - z + 3/2 z^2) for the principal branch near 0;
  - otherwise the asymptotic L1 - L2 + L2/L1 with L1 = log z + 2 pi i k,
    L2 = log L1
```

-- and refines it with Halley's cubically-convergent iteration

```text
  w <- w - (w e^w - z) / (e^w (w+1) - (w+2)(w e^w - z)/(2w+2))   (Corless 1996).
```

A real seed keeps the whole iteration exactly real (every ncpx op preserves a zero imaginary part), so real-valued branches return a real leaf with no imaginary noise. Working precision carries guard bits above the requested output precision.

```text
D[ProductLog[z], z] = ProductLog[z] / (z (1 + ProductLog[z]))  (calculus/deriv.c).
```

Series at 0, at the branch point -1/E, and at Infinity live in calculus/series.c.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_productlog` evaluates the Lambert W function:
`ProductLog[z]` is the principal branch `W_0`, `ProductLog[k, z]` the `k`-th
branch (`k` an explicit machine `Integer`). It first tries the **machine gate**
`productlog_machine`: when every argument is a machine-scalar leaf with at least
one machine `Real` (`sf_machine_scalar_args`) and the branch is `0`, the raw-double
kernels in `src/special_functions/sf_machine_plog.c` answer. `sf_machine_productlog`
(real `x >= -1/e`) seeds with `x(1 - x + 3/2 x^2)`, `log(1 + x)` or the
asymptotic `L1 - L2 + L2/L1`, and runs Halley's iteration. Near the branch point it
forms `d = x + 1/e` with `1/e` split into a double head and tail (exact by
Sterbenz's lemma), so the real/complex boundary is decided correctly, and it
iterates in `u = w + 1` on `g(u) = (u-1)e^u + 1 - e d`, which cancels the two
`-1/e` terms analytically. `sf_machine_productlog_c` handles `x < -1/e` and complex
`z`: branch-point series, a Padé seed near the origin, or the asymptotic seed,
then complex Halley. Both kernels **verify** the result (`w + log w = log z`
modulo `2πi`, and that `w` lies in the range of `W_0`) and decline otherwise.
Anything the gate does not take goes to `productlog_eval`: exact special values
(`0`, `E -> 1`, `-1/E -> -1` on branches 0 and -1, `-Pi/2 -> I Pi/2`,
`Infinity`, and `ProductLog[k, 0] = -Infinity` for `k != 0`), then the
complex-MPFR Halley core `pl_core` (`ncpx`) for any inexact argument: other
branches at machine precision, and every branch at arbitrary precision, seeded by
the branch-point series in `p = sqrt(2(e z + 1))`, the Maclaurin seed, or
`L1 = log z + 2πik`. A real seed keeps the MPFR iteration exactly real, so a
real-valued branch returns a real leaf. `D` lives in `calculus/deriv.c`; Series at
0, `-1/E` and `Infinity` in `calculus/series.c`.

**Data structures.** `Expr`; raw `double` with complex values as re/im pairs (no
`<complex.h>`) on the machine path; the shared `ncpx` (`mpfr_t` re/im) toolkit with
guard bits on the MPFR path. ND: unary kernel
`NDKU_ProductLog = { sf_machine_productlog_c, ndk_ProductLog_r, ... }`, registered
`REG_U` (so `packed_aware`); a real buffer whose elements fall below `-1/e`
returns complex values. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Halley converges cubically (60 iterations at most on the
machine path, 100 plus a polishing step on the MPFR path). A real leaf is returned
only on the real-valued domains (`W_0` for `x >= -1/e`, `W_{-1}` for
`-1/e <= x < 0`). Exact non-special arguments stay symbolic. The branch index must
be an `Integer`: `N[ProductLog[-1, -1/10]]` numericalizes `k` to `-1.` and stays
unevaluated. `Compile[]` lowers the principal branch at real and complex scalar
shapes and at rank-1 array shape (`Compiled -> True`), calling the same
kernels.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/)

- DLMF §4.13 — the Lambert W-function.
- R. M. Corless, G. H. Gonnet, D. E. G. Hare, D. J. Jeffrey and D. E. Knuth, On the Lambert W function, Adv. Comput. Math. 5 (1996) 329-359.
- Source: [`src/special_functions/productlog.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/productlog.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

- `ProductLog[z]` is the Lambert W function, the principal branch `W_0(z)`: the
  solution `w` of `z = w e^w`. `ProductLog[k, z]` selects the `k`-th branch for
  an explicit integer `k`; `k = 0` is the principal branch and `k = -1` is the
  other real branch.
- A few exact values are returned directly (`0`, `E`, `-1/E`, `-Pi/2`). Every
  other exact argument stays symbolic.
- Machine-precision principal-branch values, real or complex, come from one
  raw-`double` kernel. It handles the branch point `-1/e` carefully (the shift
  `z + 1/e` is formed without cancellation) and checks every result. The same
  kernel serves packed arrays, `NDArray` arguments and `Compile[]` (scalar,
  complex and rank-1 array), so these agree bit for bit. Other branches, and
  arbitrary precision, use a complex-MPFR Halley iteration.
- The result is real on the real-valued domains (`W_0` for `x >= -1/e`, `W_{-1}`
  for `-1/e <= x < 0`) and `Complex` otherwise.
- The branch index must be an `Integer`. `N[ProductLog[-1, -1/10], 40]` turns
  the index into `-1.` and stays unevaluated. Write `ProductLog[-1, N[-1/10, 40]]`
  instead.
- Attributes `Listable`, `NumericFunction`, `Protected`.
- See also: `Log`, `Exp`, `Solve`.
