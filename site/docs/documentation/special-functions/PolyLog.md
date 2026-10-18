# PolyLog

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolyLog[n, z]`**

gives the polylogarithm Li\_n(z) = Sum\_{k\>=1} z^k/k^n.

**`PolyLog[n, p, z]`**

gives the Nielsen generalized polylogarithm S\_{n,p}(z) (accepted but

**`PolyLog[1, z] = -Log[1-z], PolyLog[0, z] = z/(1-z), negative integer`**

**`PolyLog[n, -1] = (2^(1-n)-1) Zeta[n] for integer n >= 2, with exact forms`**

<details>
<summary>Notes</summary>

left unevaluated). Special arguments reduce in closed form: PolyLog\[n, 0\] = 0, orders give rational functions, PolyLog\[n, 1\] = Zeta\[n\] and for PolyLog\[2, 1/2\] and PolyLog\[3, 1/2\]. Inexact real or complex arguments evaluate numerically at machine or arbitrary (MPFR) precision via a power series or the Jonquiere/zeta expansion. There is a branch cut from 1 to Infinity. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= PolyLog[3, 1/2]
Out[1]= 1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]

In[2]:= PolyLog[2, 0.9]
Out[2]= 1.29971

In[3]:= PolyLog[2, 3.0]
Out[3]= 2.32018 - 3.45139*I
```

### Applications (10)

At z = 1 the polylogarithm is a zeta value

```mathematica
In[4]:= PolyLog[2, 1]
Out[4]= 1/6 Pi^2
```

Euler's closed form mixing Log[2], Pi^2 and Zeta[3]

```mathematica
In[5]:= PolyLog[3, 1/2]
Out[5]= 1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]
```

The dilogarithm at I: Catalan's constant appears

```mathematica
In[6]:= PolyLog[2, I]
Out[6]= I Catalan - 1/48 Pi^2
```

Negative integer orders are rational functions

```mathematica
In[7]:= PolyLog[-2, z]
Out[7]= (z + z^2)/(1 - z)^3
```

A complex argument at 30 digits

```mathematica
In[8]:= N[PolyLog[3, 1/2 + I/2], 30]
Out[8]= 0.48615953708556007896672148708 + 0.5700774070887689781956097575898*I
```

A fractional order at high precision via arb

```mathematica
In[9]:= N[PolyLog[1/3, 9/10], 40]
Out[9]= 5.1131317264553797957509609457956961662068
```

On the cut the value is taken from below: Im = -Pi Log[3]

```mathematica
In[10]:= PolyLog[2, 3.0]
Out[10]= 2.32018 - 3.45139*I
```

Real for z < -1, through the 1/z inversion

```mathematica
In[11]:= PolyLog[2, -1000.0]
Out[11]= -25.5025
```

A machine list runs through the packed dilogarithm kernel

```mathematica
In[12]:= PolyLog[2, {0.25, 0.5, 0.75}]
Out[12]= {0.267653, 0.582241, 0.978469}
```

A non-integer order far from the origin is outside every regime, so it stays unevaluated

```mathematica
In[13]:= PolyLog[2.5, 1000.0]
Out[13]= PolyLog[2.5, 1000.0]
```

## Algorithm

Mathilda -- the polylogarithm PolyLog.

```text
  PolyLog[n, z]     Li_n(z) = Sum_{k>=1} z^k / k^n      (|z| < 1; analytic
                    continuation elsewhere, branch cut [1, Infinity))
  PolyLog[n, p, z]  Nielsen generalized polylogarithm S_{n,p}(z).  Accepted
                    for surface compatibility but left symbolic -- there is
                    no closed-form / numeric engine for it here.
```

Evaluation is layered so each kind of argument takes the cheapest exact or fastest numeric route, mirroring src/gamma.c and src/zeta.c:

```text
  exact special values (any z) ->  closed forms
      PolyLog[n, 0]  = 0
      PolyLog[1, z]  = -Log[1 - z]
      PolyLog[0, z]  = z/(1 - z)
      PolyLog[-m, z] = Eulerian-number rational function   (m >= 1)
      PolyLog[n, 1]  = Zeta[n]                              (integer n >= 2)
      PolyLog[n, -1] = (2^(1-n) - 1) Zeta[n]                (integer n >= 2)
      PolyLog[2, 1/2] = Pi^2/12 - Log[2]^2/2
      PolyLog[3, 1/2] = Log[2]^3/6 - Pi^2 Log[2]/12 + 7 Zeta[3]/8
  numeric (>= 1 inexact operand, all numeric):
      real s, real -1 < z < 1   -> direct real MPFR power series (fast path)
      |z| <= 1/2                -> direct power series (complex MPFR)
      1/2 < |z|, |ln z| < 2 Pi  -> Jonquiere / zeta expansion (DLMF 25.12.11/12)
      otherwise                 -> stays symbolic
  everything else -> stays symbolic (return NULL)

The zeta expansion needs zeta(s - k) and Gamma(1 - s).  For real order these
```

are evaluated directly with MPFR (mpfr_zeta / mpfr_gamma); for the rare complex-order case they are obtained by evaluating the Zeta / Gamma builtins (which carry their own arbitrary-precision complex kernels).

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_polylog` evaluates `PolyLog[n, z] = Li_n(z)`. (1)
**Machine-scalar gate** (any `z` except the exact `0`): `sf_machine_scalar_args`
routes a machine-precision call to the `double` kernels of `sf_machine_zeta.c` —
`sf_machine_polylog` (= `sf_polylog_r`) for a real order and real `z ∈ [−1, 1]`
(every real `z` for the dilogarithm), `sf_polylog_c` otherwise. `PolyLog[2, z]`
uses `sf_dilog_r` / `sf_dilog_c`: the Bernoulli series in `u = −log(1−z)` after
the maps `z → 1/z` (`|z| > 1`) and `z → 1−z` (`Re z > 1/2`), so `|u| ≤ 1.25`.
Other real orders: closed forms at `s = 0, 1`, the power series for `|z| ≤ 1/2`,
the duplication formula for real `x ∈ [−1, −1/2)` with `s ≥ 1`, the Jonquière
expansion in `μ = log z` for `|μ| < 3.6`, and for integer `s ≥ 2` and `|z| > 1`
the Bernoulli-polynomial inversion formula. A real `z > 1` answers from below the
cut (`Im z = −0`); a real `z < −1` keeps only the (real) real part. A complex
order declines. (2) **Exact closed forms**: `Li_n(0) = 0`; `Li_1(z) = −Log[1−z]`;
`Li_0(z) = z/(1−z)`; `Li_{−m}(z)` the Eulerian-number rational function (cap
`POLYLOG_NEGINT_CAP = 400`); for integer `n ≥ 2`, `Li_n(1) = Zeta[n]`,
`Li_n(−1) = (2^(1−n)−1) Zeta[n]`, `Li_2(1/2)`, `Li_3(1/2)`, and `Li_2(±I) =
−π²/48 ± I Catalan`. (3) **Numeric** (`polylog_numeric`, at least one inexact
operand): real order with real `−1 < z < 1` uses arb's `acb_polylog`
(`flint_num_polylog`) above machine precision, else the real-MPFR series
`pl_real_series`; everything else runs the complex-MPFR kernel — power series
for `|z| ≤ 1/2`, the Jonquière/zeta expansion (`ζ(s−k)` and `Γ(1−s)`, reflected
through the functional equation in the left half-plane) for `|ln z| < 2π − 0.05`,
and the `Li_2` inversion `pl_li2_inversion` for `|z| > 1`; otherwise symbolic.
`PolyLog[n, p, z]` (Nielsen) is accepted and left symbolic.

**Data structures.** `Expr`; the local `pcx` (`mpfr_t` re/im) toolkit; real paths
use `mpfr_zeta`/`mpfr_gamma`, complex-order paths reuse the `Zeta`/`Gamma`
builtins. The machine kernel carries `static const` tables of ζ at the integers
(`ZM_ZPOS`, `ZM_ZNEG`), Bernoulli coefficients (`ZM_DILOG_C`, `ZM_BN`) and π²/6
as a two-double constant. ND: a hand-written binary kernel `NDKB_PolyLog` with a
complex arm (`ndk_PolyLog_c` → `sf_polylog_c`), registered `REG_B`, so packed
and visible-`NDArray` arguments, real or complex, stay on the buffer.
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Machine values ≈ 0.15 µs real, ≈ 0.4 µs complex; a
10⁶-element `PolyLog[2, v]` takes a few ms. Each MPFR regime is gated by
`|z|`/`|ln z|`; outside them (e.g. non-integer order with large `|z|`, or an
order the machine kernel declines) the call stays symbolic. `Compile[]` lowers
`PolyLog[n, z]` at scalar and rank-1 array shapes, real and complex
(`Compiled -> True`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Gamma](../../special-functions/Gamma/), [Zeta](../../special-functions/Zeta/), [NDArray](../../linear-algebra/NDArray/)

- DLMF §25.12 — the polylogarithm and Jonquiere's function (§25.12.11/12 for the zeta expansion, §25.12.4 and §25.12.13 for the inversion formulas).
- G. 't Hooft and M. Veltman, Scalar one-loop integrals, Nucl. Phys. B153 (1979) 365-401 — the Bernoulli series for Li_2 in u = -log(1-z).
- R. E. Crandall, Note on fast polylogarithm computation (2006) — the dilogarithm series and the integer-order inversion formula.
- Source: [`src/special_functions/polylog.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/polylog.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_dilog_exp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog_exp.c)
- Tests: [`tests/test_cherry_polylog_exp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_polylog_exp.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

- `PolyLog[n, z] = Li_n(z) = Σ_{k≥1} zᵏ/kⁿ`, with a branch cut along
  `[1, ∞)`. Values on the cut are taken from below, as in Mathematica.
- Exact closed forms: `PolyLog[n, 0] = 0`, `PolyLog[1, z] = -Log[1 - z]`, and
  `PolyLog[0, z] = z/(1 - z)`. Negative integer orders give rational functions.
  `PolyLog[n, ±1]` reduces to `Zeta[n]`. `Li_2(1/2)`, `Li_3(1/2)` and `Li_2(±I)`
  have closed forms.
- Numeric surfaces: a real order at machine precision runs on one `double`
  kernel (the Bernoulli-series dilogarithm, plus power series, Jonquière and
  inversion for other orders). The scalar, packed/`NDArray` and `Compile[]`
  paths (scalar and array, real and complex `z`) all share it. A complex order
  and arbitrary precision use complex MPFR; real `-1 < z < 1` at high precision
  uses arb's `acb_polylog`.
- `PolyLog[n, p, z]` (Nielsen) is accepted but left unevaluated.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `LerchPhi`, `Zeta`, `Log`.
