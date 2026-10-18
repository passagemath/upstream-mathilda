# BesselY

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselY[n, z]`**

gives the Bessel function of the second kind Y\_n(z), the solution of z^2 y'' + z y' + (z^2 - n^2) y = 0 singular at the origin.

<details>
<summary>Notes</summary>

Y\_0(0) = -Infinity, Y\_n(0) = ComplexInfinity for integer n != 0; Y\_n has a logarithmic branch point at 0 and a branch cut along the negative real z axis, with Y\_{-n} = (-1)^n Y\_n for integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision; D\[BesselY\[n, z\], z\] = (BesselY\[n-1, z\] - BesselY\[n+1, z\])/2. Listable.

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= BesselY[0, 2.5]
Out[1]= 0.49807

In[2]:= D[BesselY[n, x], x]
Out[2]= 1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])

In[3]:= BesselY[0, -1.0]
Out[3]= 0.088257 + 1.5304*I
```

### Applications (13)

Logarithmic singularity at the origin

```mathematica
In[4]:= BesselY[0, 0]
Out[4]= -Infinity
```

Negative half-odd orders vanish there: Y_{-1/2} = J_{1/2}

```mathematica
In[5]:= BesselY[-1/2, 0]
Out[5]= 0
```

Listable, machine kernel

```mathematica
In[6]:= BesselY[0, {1.0, 2.0, 3.0}]
Out[6]= {0.088257, 0.510376, 0.37685}
```

Half-integer orders close in elementary functions

```mathematica
In[7]:= BesselY[3/2, x]
Out[7]= -(Cos[x] Sqrt[2/(Pi x)])/x - Sin[x] Sqrt[2/(Pi x)]
```

Negative integer order reflects: Y_{-n} = (-1)^n Y_n

```mathematica
In[8]:= BesselY[-3, x]
Out[8]= -BesselY[3, x]
```

The logarithmic series DLMF 10.8.1

```mathematica
In[9]:= Series[BesselY[0, x], {x, 0, 3}]
Out[9]= (2 EulerGamma - 2 Log[2] + 2 Log[x])/Pi - ((-2 + 2 EulerGamma - 2 Log[2] + 2 Log[x]) x^2)/(4 Pi) + O[x]^4
```

Across the branch cut the value is complex

```mathematica
In[10]:= BesselY[0, -1.0]
Out[10]= 0.088257 + 1.5304*I
```

Fractional order, complex argument

```mathematica
In[11]:= BesselY[2.5, 4.0 + 1.0 I]
Out[11]= 0.0616085 + 0.389628*I
```

Complex order goes to the MPFR core

```mathematica
In[12]:= BesselY[0.5 I, 3 - I]
Out[12]= 1.04686 + 0.884784*I
```

Fifty digits

```mathematica
In[13]:= N[BesselY[0, 1], 50]
Out[13]= 0.0882569642156769579829267660235151628278175230906755
```

The Wronskian J_1 Y_0 - J_0 Y_1 = 2/(Pi z) holds to working precision

```mathematica
In[14]:= N[BesselJ[1, 5] BesselY[0, 5] - BesselJ[0, 5] BesselY[1, 5] - 2/(5 Pi), 30]
Out[14]= 2.138211768073756516912429173721e-50
```

No parity fold: the negative axis is a branch cut

```mathematica
In[15]:= BesselY[1, -z]
Out[15]= BesselY[1, -z]
```

Exact numeric arguments stay unevaluated

```mathematica
In[16]:= BesselY[11/2, 1]
Out[16]= BesselY[11/2, 1]
```

## Implementation notes

**Algorithm.** `builtin_bessely` tries three layers in order. (1) **Machine scalar**: when the shared gate `sf_machine_scalar_args` accepts the call (machine and exact numbers, at least one machine number, real order), `sf_bessel_scalar(SF_BESSEL_Y, …)` in `sf_machine_bessel.c` answers. For real `x > 0` and integer order `|n| <= 200` that means libm `yn`, reflected by hand. Everything else (fractional order, `x < 0`, complex `z`) goes to Amos's construction in `jy_core`: `Y_nu(z) = -(2/π) e^{-i nu π/2} K_nu(-iz) + i e^{i nu π/2} I_nu(-iz)` for `Im z >= 0`, and conjugate symmetry below. `K` comes from Temme's series, tabulated Taylor nodes or Steed's CF2, plus Hankel's expansion for `|u| >= 21`, then upward recurrence. `I` comes from the ascending series or Miller's recurrence. On the cut `x < 0` the result is a `Complex`. A decline (`|nu| > 1000`, overflow) passes the call on. (2) **Exact origin**: `bessel_origin_value` returns `-Infinity` for `nu = 0` and `0` for a negative half-odd-integer order (where `cos(nu π) = 0` cancels the divergence, `Y_{-1/2} = J_{1/2}`). It returns `Indeterminate` for a pure-imaginary order and `ComplexInfinity` otherwise. (3) **MPFR** (`by_eval`, at the least input precision): an integer order with real `z > 0` uses the correctly rounded `mpfr_yn`. Otherwise `by_core` picks one of three methods. For `|z| > 1.3·max(4, P ln2/2)` and `|arg z| < 0.9π` it uses the asymptotic series DLMF 10.17.4. For smaller `|z|` it uses the logarithmic series DLMF 10.8.1 at integer order, or the connection `(J_nu cos nuπ - J_{-nu})/sin nuπ` at non-integer order, with `-log2|sin nuπ|` extra guard bits. An exact call then reaches `bessel_symbolic_rules`, the former `bessel.m` DownValues, now in C. A half-integer order with a symbolic argument gives the spherical-Bessel closed form by recurrence. A negative integer order is rewritten by `Y_{-n} = (-1)^n Y_n`. There is deliberately no argument-parity fold, because of the branch cut.

**Data structures.** Scalars use `Expr` leaves. The machine kernel works on double `(re, im)` pairs with generated tables (`sf_machine_bessel_tab.c`, `sf_machine_bessel_ktab.c`). The MPFR path uses the `ncpx` complex-MPFR toolkit and reuses `BesselJ`'s power series. The ND kernel `NDKB_BesselY` (`REG_B` in `ndkernels.c`, via `sf_bessel_nd`) runs over real and complex buffers with a real order. `Compile[]` lowers it at scalar real, scalar complex and rank-1 shapes through that kernel, and every surface agrees bit for bit with the scalar.

**Complexity / limits.** At machine precision the cost is O(1) on the real axis at integer order and O(|z|) complex operations otherwise. The documented accuracy is ≤ 3.6 ulp of the local magnitude on the real axis and median ≈ 1 for complex arguments. MPFR term counts and guard bits scale with `|z|` and `P`. `Y_nu` has a logarithmic branch point at 0 and a cut along the negative real axis. Exact numeric arguments stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [BesselJ](../../special-functions/BesselJ/), [BesselK](../../special-functions/BesselK/), [BesselI](../../special-functions/BesselI/), [NDArray](../../linear-algebra/NDArray/), [Real](../../other-advanced/Real/), [Complex](../../arithmetic/Complex/)

- DLMF §10.8.1 — the logarithmic series for integer-order Y_n(z) near the origin.
- DLMF §10.2.3 — Y_nu = (J_nu cos(nu pi) - J_{-nu})/sin(nu pi).
- DLMF §10.17.4 — the large-argument asymptotic expansion of Y_nu(z).
- D. E. Amos, A portable package for Bessel functions of a complex argument and nonnegative order, ACM TOMS 12 (1986) 265-273 (Algorithm 644).
- N. M. Temme, On the numerical evaluation of the modified Bessel function of the third kind, J. Comput. Phys. 19 (1975) 324-337.
- Source: [`src/special_functions/bessel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/bessel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_bessely.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bessely.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m61_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m61_stress.c)

## Notes & additional examples

### Notes

- `BesselY[n, z]` is the Bessel function of the second kind. It is singular at the origin (`Y_0(0) = -Infinity`, `ComplexInfinity` for most other orders) and has a cut along the negative real axis.
- Machine-number calls with a real order use the shared raw-double kernel (`sf_machine_bessel.c`): libm `yn` for integer order with `x > 0`, and Amos's construction from `K_nu(-iz)` and `I_nu(-iz)` otherwise. The same helper serves packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1), bit for bit. Complex order, arbitrary precision and `|nu| > 1000` use MPFR: `mpfr_yn`, the logarithmic series, the connection formula or Hankel's expansion.
- Half-integer orders with a symbolic argument give spherical-Bessel closed forms.
- `D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselJ`, `BesselK`, `BesselI`.
