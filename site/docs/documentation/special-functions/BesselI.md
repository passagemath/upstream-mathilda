# BesselI

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselI[n, z]`**

gives the modified Bessel function of the first kind I\_n(z), the solution of z^2 y'' + z y' - (z^2 + n^2) y = 0 regular at the origin.

**`D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2. Listable.`**

<details>
<summary>Notes</summary>

I\_0(0) = 1, I\_n(0) = 0 for integer n != 0; I\_n grows like e^z as z -\> Inf and is even in n (I\_{-n} = I\_n). Has a branch cut along the negative real z axis for non-integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= BesselI[0, 2.0]
Out[1]= 2.27959

In[2]:= D[BesselI[n, x], x]
Out[2]= 1/2 (BesselI[-1 + n, x] + BesselI[1 + n, x])

In[3]:= BesselI[1, -z]
Out[3]= -BesselI[1, z]
```

### Applications (14)

Regular at the origin

```mathematica
In[4]:= BesselI[0, 0]
Out[4]= 1
```

A negative fractional order diverges there

```mathematica
In[5]:= BesselI[-1/2, 0]
Out[5]= ComplexInfinity
```

Listable over the order

```mathematica
In[6]:= BesselI[{0, 1, 2}, 1.]
Out[6]= {1.26607, 0.565159, 0.135748}
```

Half-integer orders give hyperbolic closed forms; I is not even at half-integer order

```mathematica
In[7]:= BesselI[-1/2, x]
Out[7]= Cosh[x] Sqrt[2/(Pi x)]
```

Even in an integer order

```mathematica
In[8]:= BesselI[-3, x]
Out[8]= BesselI[3, x]
```

Integer-order parity

```mathematica
In[9]:= BesselI[1, -z]
Out[9]= -BesselI[1, z]
```

All-positive Frobenius series

```mathematica
In[10]:= Series[BesselI[0, x], {x, 0, 8}]
Out[10]= 1 + 1/4 x^2 + 1/64 x^4 + 1/2304 x^6 + 1/147456 x^8 + O[x]^9
```

On the cut the principal value e^{i nu pi} I_nu(3) is purely imaginary

```mathematica
In[11]:= BesselI[2.5, -3.0]
Out[11]= 0.0 + 1.51534*I
```

Complex order and argument through the MPFR core

```mathematica
In[12]:= BesselI[3 + I, 1.5 - I]
Out[12]= -0.25665 + 0.0492771*I
```

Fifty digits

```mathematica
In[13]:= N[BesselI[0, 1], 50]
Out[13]= 1.26606587775200833559824462521471753760767031135496
```

Beyond the double range the value comes back at arbitrary precision

```mathematica
In[14]:= BesselI[0, 800.]
Out[14]= 3.8460816154192104336e+345
```

The Wronskian I_0 K_1 + I_1 K_0 = 1/z, here 1/2

```mathematica
In[15]:= N[BesselI[0, 2] BesselK[1, 2] + BesselI[1, 2] BesselK[0, 2], 30]
Out[15]= 0.5
```

From the integral tables

```mathematica
In[16]:= Integrate[BesselI[1, x], x]
Out[16]= BesselI[0, x]
```

Exponential growth

```mathematica
In[17]:= Limit[BesselI[0, x], x -> Infinity]
Out[17]= Infinity
```

## Implementation notes

**Algorithm.** `builtin_besseli` tries three layers in order. (1) **Machine scalar**: with a real order and machine arguments accepted by `sf_machine_scalar_args`, `sf_bessel_scalar(SF_BESSEL_I, …)` answers. On the real line `I_0` and `I_1` come from Cephes-style regions in `bes_i01`: the series, then piecewise Chebyshev fits of the exponentially scaled function, then Hankel. A negative real `x` uses parity at integer order and gives a `Complex` value otherwise, through the continuation `I_nu(z e^{±iπ}) = e^{±i nu π} I_nu(z)`. Other orders and complex `z` use `sf_bessel_i_c`. There the ascending series is used where `|z|^2/4 <= nu + 1`, so no cancellation is possible. Elsewhere Miller's backward recurrence runs, normalised by a Gegenbauer/Neumann sum (`i_miller`, Amos's CMLRI). That needs no `K`, so it costs about `|z|` complex multiplies. Hankel's expansion takes over for `|z| >= 21`. A negative order uses `I_{-nu} = I_nu + (2/π) sin(nu π) K_nu`. The kernel declines `|nu| > 1000` and overflow, and those calls fall through. (2) **Exact origin** (`bessel_origin_value`): `1` for `nu = 0`, `0` for integer `n != 0` or `Re nu > 0`, `ComplexInfinity` for non-integer `Re nu < 0`, and `Indeterminate` for a pure-imaginary order. (3) **MPFR** (`bi_eval`): MPFR has no modified-Bessel routine. `bi_core` uses the shared all-positive power series `bk_iseries` with about `2|z|/ln2` guard bits for complex cancellation. Past `|z| > 1.3·max(4, P ln2/2)` it switches to the two-exponential asymptotic form DLMF 10.40.5, valid on the whole principal sheet. Mathematically real results drop the sub-precision imaginary noise. Exact calls then reach `bessel_symbolic_rules` (formerly `bessel.m`). A half-integer order with a symbolic argument gives `Sinh`/`Cosh` closed forms; both `±1/2` base cases are needed, because `I` is not even at half-integer order. `I_{-n} = I_n` covers negative integer order. `bessel_z_parity_fold` applies `I_n(-z) = (-1)^n I_n(z)`.

**Data structures.** Scalars use `Expr` leaves. The machine kernel works on double `(re, im)` pairs with the coefficient tables in `sf_machine_bessel_tab.c`. The MPFR path uses the `ncpx` complex-MPFR toolkit. The ND kernel `NDKB_BesselI` (`REG_B` in `ndkernels.c`, through `sf_bessel_nd`) covers real and complex buffers with a real order. `Compile[]` lowers it at scalar real, scalar complex and rank-1 shapes through that kernel, bit-identical to the scalar. A machine result past the double range declines, so the value is returned at arbitrary precision (`BesselI[0, 800.] ≈ 3.846×10^345`).

**Complexity / limits.** At machine precision the cost is O(1) on the real axis for orders 0 and 1, and O(|z|) through Miller's recurrence otherwise. MPFR term counts and guard bits scale with `|z|` and `P`. Non-integer order has a branch cut along the negative real axis. Exact numeric arguments stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [BesselJ](../../special-functions/BesselJ/), [BesselK](../../special-functions/BesselK/), [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/), [Real](../../other-advanced/Real/), [Complex](../../arithmetic/Complex/), [Normal](../../data-structures/Normal/), [Plus](../../arithmetic/Plus/)

- DLMF §10.25.2 — power series I_nu(z) = Sum_k (z/2)^{nu+2k}/(k! Gamma(nu+k+1)).
- DLMF §10.40.5 — the large-argument asymptotic expansion of I_nu(z).
- DLMF §10.34.1 — I_nu(z e^{±i pi}) = e^{±i nu pi} I_nu(z).
- D. E. Amos, A portable package for Bessel functions of a complex argument and nonnegative order, ACM TOMS 12 (1986) 265-273 (Algorithm 644).
- Source: [`src/special_functions/bessel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/bessel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_besseli.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besseli.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)

## Notes & additional examples

### Notes

- `BesselI[n, z]` is the modified Bessel function of the first kind. It grows like `e^z`, is even in an integer order, and has a cut along the negative real axis for non-integer order.
- Machine-number calls with a real order use the raw-double kernel in `sf_machine_bessel.c`: Cephes-style `I_0`/`I_1` on the real line, and otherwise the ascending series or Miller's backward recurrence. Packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1) call the same helper. Overflow, complex order and arbitrary precision fall through to the MPFR power or asymptotic series.
- `D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselK`, `BesselJ`, `BesselY`.
