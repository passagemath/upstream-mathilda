# BesselJ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselJ[n, z]`**

gives the Bessel function of the first kind J\_n(z), a solution of z^2 y'' + z y' + (z^2 - n^2) y = 0 regular at the origin.

**`D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2. Listable.`**

<details>
<summary>Notes</summary>

J\_0(0) = 1, J\_n(0) = 0 for integer n != 0. Has a branch cut along the negative real z axis for non-integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= BesselJ[0, 5.2]
Out[1]= -0.11029

In[2]:= D[BesselJ[n, x], x]
Out[2]= 1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])

In[3]:= BesselJ[1, -2 x]
Out[3]= -BesselJ[1, 2 x]
```

### Applications (13)

Regular at the origin: J_0(0) = 1

```mathematica
In[4]:= BesselJ[0, 0]
Out[4]= 1
```

Listable over the order, machine kernel

```mathematica
In[5]:= BesselJ[{0, 1, 2}, 2.5]
Out[5]= {-0.0483838, 0.497094, 0.446059}
```

Half-integer orders close in elementary functions

```mathematica
In[6]:= BesselJ[5/2, x]
Out[6]= -Sin[x] Sqrt[2/(Pi x)] + (3 (-Cos[x] Sqrt[2/(Pi x)] + (Sin[x] Sqrt[2/(Pi x)])/x))/x
```

Integer-order parity J_n(-z) = (-1)^n J_n(z)

```mathematica
In[7]:= BesselJ[1, -2 x]
Out[7]= -BesselJ[1, 2 x]
```

The alternating Frobenius series

```mathematica
In[8]:= Series[BesselJ[0, x], {x, 0, 6}]
Out[8]= 1 - 1/4 x^2 + 1/64 x^4 - 1/2304 x^6 + O[x]^7
```

On the branch cut a real argument gives a complex value

```mathematica
In[9]:= BesselJ[1/2, -2.]
Out[9]= 0.0 + 0.513016*I
```

Complex argument at machine precision

```mathematica
In[10]:= BesselJ[1.5, 2.0 + 3.0 I]
Out[10]= 2.97761 + 1.4353*I
```

Complex order takes the MPFR core

```mathematica
In[11]:= BesselJ[7/3 + I, 4.5 - I]
Out[11]= 1.18908 + 0.715653*I
```

Fifty digits from the arbitrary-precision path

```mathematica
In[12]:= N[BesselJ[0, 4], 50]
Out[12]= -0.39714980986384737228659076845169804197561868528939
```

An NDArray runs through the element-wise kernel

```mathematica
In[13]:= BesselJ[1, NDArray[{0.5, 1., 1.5}]]
Out[13]= NDArray[{0.242268, 0.440051, 0.557937}]
```

Compiled, bit-identical to the scalar values

```mathematica
In[14]:= cf = Compile[{{v, _Real, 1}}, BesselJ[1, v]]; cf[{0.5, 1., 1.5}]
Out[14]= {0.242268, 0.440051, 0.557937}
```

An exact numeric argument stays unevaluated, as in Mathematica

```mathematica
In[15]:= BesselJ[11/2, 1]
Out[15]= BesselJ[11/2, 1]
```

Numericalise explicitly instead

```mathematica
In[16]:= N[BesselJ[11/2, 1]]
Out[16]= 7.38531e-05
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| BesselJ[0, .] over 10^6 | 3.46e+03 s | 1.74e+03 s | 53 s |
| Zeta over 10^6 | 8.17 s | 4.31e+03 s | 5.7 s |
| AiryAi over 10^6 | 3.14 s | 107 s | 58.9 s |
| PolyGamma[0, .] over 10^6 | 1.55 s | 151 s | 8.19 s |
| Gamma over 10^6 | 1.16 s | 1.35 s | 7.34 s |
| Erf over 10^6 | 0.968 s | 1.23 s | 7.43 s |

## Implementation notes

**Algorithm.** `builtin_besselj` tries three layers in order. (1) **Machine scalar**: `bessel_machine` asks the shared gate `sf_machine_scalar_args` (`sf_scalar.c`) whether every argument is a machine or exact number with at least one machine `Real`/`Complex`. If the order is real, the call goes to `sf_bessel_scalar(SF_BESSEL_J, …)` in `sf_machine_bessel.c`. On the real line, an integer order `|n| <= 200` uses libm `jn`, reflected by hand for negative `n`. Every other case goes through Amos's scheme (TOMS 644): `jy_core` evaluates `J_nu(z) = e^{i nu pi/2} I_nu(-iz)` in the right half plane. There `K_mu, K_{mu+1}` (`|mu| <= 1/2`) come from Temme's series (`|u| <= 1.25`), from Taylor continuation off tabulated nodes (`ktab_eval`, orders 0 and ±1/3), or from Steed's algorithm for Temme's CF2, and Hankel's expansion takes over for `|u| >= 21`. `I_nu` comes from the ascending series or from Miller's backward recurrence (`i_miller`), and conjugate symmetry covers `Im z < 0`. A real argument whose value is complex (`BesselJ[1/2, -2.]`) returns a `Complex`. The kernel declines `|nu| > 1000` (it has no Debye expansion), overflow and underflow, and passes them on. (2) **Exact origin**: `bessel_origin_value` classifies the order. It returns `1` (`nu = 0`), `0` (integer `n != 0` or `Re nu > 0`), `ComplexInfinity` (non-integer `Re nu < 0`) or `Indeterminate` (pure-imaginary `nu`). (3) **MPFR**: when an argument is inexact (arbitrary precision, a complex order, or a machine decline), `bj_eval` runs at the least input precision. An integer order with a real argument uses MPFR's correctly rounded `mpfr_jn`. Otherwise the complex-MPFR core `bj_core` runs. It uses the power series DLMF 10.2.2 (one `Gamma(nu+1)`, about `2|z|/ln2` guard bits for the alternating cancellation), or the asymptotic series DLMF 10.17.3 summed to optimal truncation when `|z| > 1.3·max(4, P ln2/2)` and `|arg z| < 0.9π`. Whatever stays exact goes to `bessel_symbolic_rules`, the former `src/internal/bessel.m` DownValues now in C so numeric calls no longer pay for the pattern matching. A half-integer order with a non-numeric argument becomes an elementary closed form through the up/down recurrence. A negative integer order is rewritten by `J_{-n} = (-1)^n J_n`. `bessel_z_parity_fold` applies `J_n(-z) = (-1)^n J_n(z)` for integer order. `D` lives in `calculus/deriv.c` and `Series` at 0/∞ in `calculus/series.c`.

**Data structures.** Scalars use `Expr` leaves. The machine kernel works on `(re, im)` pairs of doubles with Smith division, double-double error-free transforms and generated coefficient tables (`sf_machine_bessel_tab.c`, `sf_machine_bessel_ktab.c`). The MPFR path uses the `ncpx` complex-MPFR toolkit (`numeric_complex.h`). The ND kernel `NDKB_BesselJ` (`REG_B` in `ndkernels.c`) calls the same `sf_bessel_nd` dispatch over real **and** complex buffers with a real order, so packed lists and visible `NDArray`s run element-wise. A real element with a complex value declines to the `List` path. `Compile[]` lowers `BesselJ[n, z]` at scalar real, scalar complex and rank-1 array shapes through that kernel. The scalar, List, packed, `NDArray` and compiled answers agree bit for bit.

**Complexity / limits.** At machine precision the cost is O(1) for real integer order and O(|z|) complex multiplies for Miller recurrence otherwise. The documented accuracy (`tools/sf_accuracy/bessel.py`) is ≤ 3.6 ulp on the real axis, median ≈ 1 and max ≤ 8 for complex arguments, and max ≈ 20 for fractional orders. The MPFR path's term counts and guard bits grow with `|z|` and `P`. Non-integer order has a branch cut along the negative real axis. An exact numeric argument (`BesselJ[11/2, 1]`) stays unevaluated, as in Mathematica. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [BesselY](../../special-functions/BesselY/), [BesselI](../../special-functions/BesselI/), [BesselK](../../special-functions/BesselK/), [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/), [Real](../../other-advanced/Real/), [Complex](../../arithmetic/Complex/), [Gamma](../../special-functions/Gamma/)

- DLMF §10.2.2 — power series J_nu(z) = Sum_k (-1)^k (z/2)^{nu+2k}/(k! Gamma(nu+k+1)).
- DLMF §10.17.3 — the large-argument asymptotic expansion of J_nu(z).
- DLMF §10.27.6, §10.34 — J_nu(z) = e^{i nu pi/2} I_nu(-i z) and the analytic-continuation formulas used by the machine kernel.
- D. E. Amos, A portable package for Bessel functions of a complex argument and nonnegative order, ACM TOMS 12 (1986) 265-273 (Algorithm 644).
- N. M. Temme, On the numerical evaluation of the modified Bessel function of the third kind, J. Comput. Phys. 19 (1975) 324-337.
- I. J. Thompson and A. R. Barnett, Modified Bessel functions I_nu(z) and K_nu(z) of real order and complex argument, to selected accuracy, Comput. Phys. Commun. 47 (1987) 245-257.
- Source: [`src/special_functions/bessel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/bessel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m61_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m61_stress.c)

## Notes & additional examples

### Notes

- `BesselJ[n, z]` is the Bessel function of the first kind, regular at the origin. `J_0(0) = 1`, `J_n(0) = 0` for integer `n != 0`, and a non-integer order with `Re nu < 0` gives `ComplexInfinity` at 0.
- Machine-number calls with a **real** order go to the raw-double kernel in `sf_machine_bessel.c`: libm `jn` for real integer orders, and Amos's I/K scheme otherwise. A packed list, a visible `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1) call the same helper, so all five surfaces agree bit for bit. A complex order, an argument with more than machine precision, or `|nu| > 1000` uses the MPFR core (`mpfr_jn`, power series or Hankel asymptotic).
- Half-integer orders with a symbolic argument rewrite to `Sin`/`Cos` forms, and negative integer orders reflect by `J_{-n} = (-1)^n J_n`. These are C rules now, no longer `bessel.m` DownValues.
- Non-integer order has a branch cut along the negative real axis.
- `D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselY`, `BesselI`, `BesselK`, `BesselJZero`, `AiryAi`.
