# BesselK

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselK[n, z]`**

gives the modified Bessel function of the second kind K\_n(z), a solution of z^2 y'' + z y' - (z^2 + n^2) y = 0.

**`D[BesselK[n, z], z] = -(BesselK[n-1, z] + BesselK[n+1, z])/2. Listable.`**

<details>
<summary>Notes</summary>

K\_n(z) is even in n (K\_{-n} = K\_n) and decays like e^{-z} as z -\> Inf. K\_0(0) = Infinity, K\_n(0) = ComplexInfinity. Has a branch cut along the negative real z axis. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= BesselK[0, 4.0]
Out[1]= 0.0111597

In[2]:= D[BesselK[n, x], x]
Out[2]= -1/2 (BesselK[-1 + n, x] + BesselK[1 + n, x])

In[3]:= BesselK[0, -1.0]
Out[3]= 0.421024 - 3.97746*I
```

### Applications (14)

A directed divergence at the origin

```mathematica
In[4]:= BesselK[0, 0]
Out[4]= Infinity
```

Pure-imaginary order oscillates boundedly, so there is no limit

```mathematica
In[5]:= BesselK[I, 0]
Out[5]= Indeterminate
```

Listable over the order

```mathematica
In[6]:= BesselK[{1, 2, 3}, 1.0]
Out[6]= {0.601907, 1.62484, 7.10126}
```

Half-integer orders give decaying elementary forms

```mathematica
In[7]:= BesselK[3/2, x]
Out[7]= E^(-x) Sqrt[Pi/(2 x)] + (E^(-x) Sqrt[Pi/(2 x)])/x
```

Even in the order

```mathematica
In[8]:= BesselK[-2, x]
Out[8]= BesselK[2, x]
```

Logarithmic expansion at the origin

```mathematica
In[9]:= Series[BesselK[0, x], {x, 0, 3}]
Out[9]= 1/2 (-2 EulerGamma + 2 Log[2] - 2 Log[x]) + 1/8 (2 - 2 EulerGamma + 2 Log[2] - 2 Log[x]) x^2 + O[x]^4
```

Exponential decay

```mathematica
In[10]:= Limit[BesselK[1, x], x -> Infinity]
Out[10]= 0
```

Fractional order, complex argument, machine kernel

```mathematica
In[11]:= BesselK[1/3, 2.0 + 1.0 I]
Out[11]= 0.0378793 - 0.104026*I
```

Across the cut: K_0(1) - i Pi I_0(1)

```mathematica
In[12]:= BesselK[0, -1.0]
Out[12]= 0.421024 - 3.97746*I
```

Complex order takes the MPFR core

```mathematica
In[13]:= BesselK[1 + I, 3.0 - 2 I]
Out[13]= -0.0225108 + 0.0169607*I
```

Fifty digits

```mathematica
In[14]:= N[BesselK[0, 4], 50]
Out[14]= 0.0111596760858530242697451959798334892250090238884744
```

An underflowing double is promoted rather than flushed to zero

```mathematica
In[15]:= BesselK[0, 800.]
Out[15]= 1.625030835448122e-349
```

An NDArray runs through the element-wise kernel

```mathematica
In[16]:= BesselK[0, NDArray[{0.5, 1., 2.}]]
Out[16]= NDArray[{0.924419, 0.421024, 0.113894}]
```

No parity fold across the branch cut

```mathematica
In[17]:= BesselK[1, -z]
Out[17]= BesselK[1, -z]
```

## Implementation notes

**Algorithm.** `builtin_besselk` tries three layers in order. (1) **Machine scalar**: with a real order and arguments accepted by `sf_machine_scalar_args`, `sf_bessel_scalar(SF_BESSEL_K, …)` answers. For real `x > 0` and integer order, `K_0` and `K_1` come from `bes_k01` (series and Chebyshev fits of the scaled function), followed by **upward** recurrence `K_{k+1} = K_{k-1} + (2k/x) K_k`, the stable direction for `K`. Fractional orders and complex `z` use `sf_bessel_k_c`. With `nu = n + mu`, `|mu| <= 1/2`, it computes `K_mu, K_{mu+1}` by one of three methods, then recurses upward. The methods are Temme's series for `|z| <= 1.25`, Taylor continuation of the modified Bessel equation from tabulated nodes (`ktab_eval`, orders 0 and ∓1/3) or Steed's algorithm for Temme's CF2 for `1.25 < |z| < 21`, and Hankel's expansion beyond. The left half plane is reached by `K_nu(z e^{±iπ}) = e^{∓i nu π} K_nu(z) ∓ iπ I_nu(z)`, so `BesselK[0, -1.]` is complex. The kernel declines underflow to `0`, overflow and `|nu| > 1000`. (2) **Exact origin**: `K_0(0) = Infinity`, `Indeterminate` for a pure-imaginary order (bounded oscillation), and `ComplexInfinity` for every other order. (3) **MPFR** (`bk_eval`): `bk_core` uses the asymptotic series DLMF 10.40.2 for `|z| > 1.3·max(4, P ln2/2)`. Below that it uses the logarithmic series DLMF 10.31.1 at integer order, or the connection `(π/2)(I_{-nu} - I_nu)/sin(nu π)` at non-integer order, with `-log2|sin nu π|` guard bits. `builtin_besselk` wraps the result in `numeric_promote_result_if_degenerate`, so a flush to zero is re-evaluated with an extended exponent (`BesselK[0, 800.] ≈ 1.625×10^-349`). Exact calls go to `bessel_symbolic_rules` (formerly `bessel.m`). A half-integer order with a symbolic argument gives `Sqrt[Pi/(2 z)] E^(-z)` times a polynomial in `1/z`, and `K_{-nu} = K_nu`. There is no argument-parity fold, because of the branch cut.

**Data structures.** Scalars use `Expr` leaves. The machine kernel works on double `(re, im)` pairs. Its K-node tables (`sf_machine_bessel_ktab.c`, about 7k lines of generated Taylor nodes) and its scaled-function fits (`sf_machine_bessel_tab.c`) are generated. The MPFR path uses the `ncpx` toolkit. The ND kernel `NDKB_BesselK` (`REG_B` in `ndkernels.c`, through `sf_bessel_nd`) covers real and complex buffers with a real order. `Compile[]` lowers it at scalar real, scalar complex and rank-1 shapes through that kernel, bit-identical to the scalar.

**Complexity / limits.** At machine precision the cost is O(1) plus `n` recurrence steps at integer order, and a bounded table or CF2 evaluation at fractional order. MPFR term counts and guard bits scale with `|z|` and `P`. `K` diverges at the origin for every order and has a cut along the negative real axis. Exact numeric arguments stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [BesselJ](../../special-functions/BesselJ/), [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/), [Real](../../other-advanced/Real/), [Complex](../../arithmetic/Complex/)

- DLMF §10.27.4 — K_nu(z) = (π/2)(I_{-nu}(z) - I_nu(z))/sin(nu π).
- DLMF §10.31.1 — the logarithmic series for integer-order K_n(z).
- DLMF §10.40.2 — the large-argument asymptotic expansion K_nu(z) ~ sqrt(π/(2z)) e^{-z} Sum a_k/z^k.
- DLMF §10.34.2 — K_nu(z e^{±i pi}) = e^{∓i nu pi} K_nu(z) ∓ i pi I_nu(z).
- D. E. Amos, A portable package for Bessel functions of a complex argument and nonnegative order, ACM TOMS 12 (1986) 265-273 (Algorithm 644).
- N. M. Temme, On the numerical evaluation of the modified Bessel function of the third kind, J. Comput. Phys. 19 (1975) 324-337.
- I. J. Thompson and A. R. Barnett, Modified Bessel functions I_nu(z) and K_nu(z) of real order and complex argument, to selected accuracy, Comput. Phys. Commun. 47 (1987) 245-257.
- Source: [`src/special_functions/bessel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/bessel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_besselk.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselk.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)

## Notes & additional examples

### Notes

- `BesselK[n, z]` is the modified Bessel function of the second kind. It decays like `e^{-z}`, is even in the order, diverges at the origin, and has a cut along the negative real axis.
- Machine-number calls with a real order use the raw-double kernel in `sf_machine_bessel.c`. Integer orders take table-fitted `K_0`/`K_1` and upward recurrence. Fractional orders take Temme's series, tabulated Taylor nodes or Steed's CF2, plus Hankel. Packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1) call the same helper, bit for bit. Complex order, arbitrary precision and underflow go to the MPFR core: the logarithmic series, the connection formula or the asymptotic series.
- `D[BesselK[n, z], z] = -(BesselK[n-1, z] + BesselK[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselI`, `BesselJ`, `BesselY`, `AiryAi`.
