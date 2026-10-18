# MeijerG

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MeijerG[{{a1, ..., an}, {an+1, ..., ap}}, {{b1, ..., bm}, {bm+1, ..., bq}}, z]`**

is the Meijer G-function G^{m,n}\_{p,q}(z | a; b), the Mellin-Barnes integral

**`MeijerG[alist, blist, z, r]`**

is the generalized form, s -\> r s with a factor r

<details>
<summary>Notes</summary>

(1/(2 Pi I)) Int prod\_{j\<=m} Gamma\[bj+s\] prod\_{j\<=n} Gamma\[1-aj-s\] / (prod\_{j\>n} Gamma\[aj+s\] prod\_{j\>m} Gamma\[1-bj-s\]) z^-s ds. (equal to MeijerG\[alist, blist, z^(1/r)\] on the principal branch). Parameters are sorted within each group. Simple cases reduce automatically to Exp, powers, Bessel functions, HypergeometricU, EllipticK/E and HypergeometricPFQ (Slater's theorem); z = 0 gives 0, a finite value or ComplexInfinity. Evaluates to machine, arbitrary-precision and complex numbers (output precision = lowest input precision), including the log case of integer-spaced parameters; threads over a list or packed/NDArray z and compiles. Inconsistent parameters (m = n = 0, or ak - bj a positive integer for k \<= n, j \<= m) give MeijerG::hdiv.

</details>

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= N[MeijerG[{{},{2}},{{1/2,3/2},{}},3],50]
Out[1]= 0.0415911598753969836703303907497787586265297151590396

In[2]:= MeijerG[{{},{}},{{v},{-v}},z]
Out[2]= BesselJ[2 v, 2 Sqrt[z]]

In[3]:= MeijerG[{{1/2},{}},{{0},{0}},x]
Out[3]= Sqrt[Pi] BesselI[0, 1/2 x] E^(-1/2 x)

In[4]:= MeijerG[{{3+I},{2.}},{{2/3,1/5},{4}},3]
Out[4]= -0.664578 - 3.02966*I

In[5]:= MeijerG[{{1/3},{}},{{1/2},{1}},{0.1, 0.3, 0.5},1/2]
Out[5]= {0.0511285, 0.126112, 0.133903}

In[6]:= D[MeijerG[{{a1,a2},{a3,a4}},{{v},{-v}},z],z]
Out[6]= MeijerG[{{-1, -1 + a1, -1 + a2}, {-1 + a3, -1 + a4}}, {{-1 + v}, {0, -1 - v}}, z]
```

### Worked examples (2)

```mathematica
In[7]:= MeijerG[{{1},{2}},{{2/3,1/5},{4}},z]
Out[7]= MeijerG[{{1}, {2}}, {{1/5, 2/3}, {4}}, z]

In[8]:= MeijerG[{{31,4},{2,3}},{{5,3/2},{1/10}},3]
Out[8]= -MeijerG[{{2, 4, 31}, {3}}, {{3/2}, {1/10, 5}}, 3]
```

### Applications (10)

The simplest G-function, G^{1,0}_{0,1}, is the exponential

```mathematica
In[9]:= MeijerG[{{},{}},{{0},{}},z]
Out[9]= E^(-z)
```

A one-term Slater sum that names its 2F1: the logarithm

```mathematica
In[10]:= MeijerG[{{1,1},{}},{{1},{0}},z]
Out[10]= Log[1 + z]
```

A 2F1 with parameters 1/2, 1/2; 1 is named as EllipticK

```mathematica
In[11]:= MeijerG[{{1/2,1/2},{}},{{0},{0}},z]
Out[11]= 2 EllipticK[-z]
```

Parameters sorted within each group; two 2F2 terms do not reduce

```mathematica
In[12]:= MeijerG[{{1},{2}},{{2/3,1/5},{4}},z]
Out[12]= MeijerG[{{1}, {2}}, {{1/5, 2/3}, {4}}, z]
```

The log case, b's an integer apart: machine kernel

```mathematica
In[13]:= MeijerG[{{},{2}},{{1/2,3/2},{}},3.]
Out[13]= 0.0415912
```

Checked against G^{2,0}_{1,2}(z|a;b,c) = z^b e^(-z) U(a-c, 1+b-c, z)

```mathematica
In[14]:= N[Sqrt[3] Exp[-3] HypergeometricU[1/2, 0, 3], 30] - N[MeijerG[{{},{2}},{{1/2,3/2},{}},3], 30]
Out[14]= 0.0
```

Complex parameters

```mathematica
In[15]:= MeijerG[{{3+I},{2.}},{{2/3,1/5},{4}},3]
Out[15]= -0.664578 - 3.02966*I
```

The generalized form: MeijerG[a, b, z, r] = MeijerG[a, b, z^(1/r)], so e^-2

```mathematica
In[16]:= MeijerG[{{},{}},{{0},{}},4., 2]
Out[16]= 0.135335
```

The derivative is another MeijerG

```mathematica
In[17]:= D[MeijerG[{{1},{2}},{{1/5,2/3},{4}},z],z]
Out[17]= MeijerG[{{-1, 0}, {1}}, {{-4/5, -1/3}, {0, 3}}, z]
```

Pole collision, a_1 - b_1 a positive integer: MeijerG::hdiv, unevaluated

```mathematica
In[18]:= MeijerG[{{2},{}},{{1},{}},z]
Out[18]= MeijerG[{{2}, {}}, {{1}, {}}, z]
```

## Algorithm

meijer_g.c -- MeijerG, the Meijer G-function, for Mathilda.

```text
  MeijerG[{{a1..an},{an+1..ap}}, {{b1..bm},{bm+1..bq}}, z]
  MeijerG[{{a1..an},{an+1..ap}}, {{b1..bm},{bm+1..bq}}, z, r]
```

defined by the Mellin-Barnes integral

```text
           r    /  prod_{j<=m} G(b_j+rs) prod_{j<=n} G(1-a_j-rs)
  G  =  ------  |  ------------------------------------------------ z^-s ds
        2 pi i  /  prod_{j>n} G(a_j+rs) prod_{j>m} G(1-b_j-rs)
```

(G = Gamma; r = 1 is the ordinary G-function). Substituting t = r s shows MeijerG[a, b, z, r] = MeijerG[a, b, w] with w = z^(1/r) on the principal branch, which is how every path below treats r.

Evaluation order in builtin_meijer_g:

```text
  1. validate: 3 or 4 arguments (else hdiv + argt), each parameter argument a
     pair of Lists (else hdiv); r == 1 drops to the 3-argument form;
  2. buffer z (NDArray / packed List): the machine kernel over the buffer;
     a plain List z threads element-wise (MeijerG is not Listable: its first
     two arguments ARE lists);
  3. canonical order: each of the four parameter groups sorted (the function
     is symmetric within a group), as Mathematica does;
  4. consistency: m = n = 0, or some a_k - b_j (k <= n, j <= m) a positive
     integer (the contour cannot separate the two pole families): hdiv;
  5. machine-precision call: sf_mg_eval (sf_machine_meijerg.c), the helper
     shared with the buffer and Compile[] surfaces, before any exact rewrite;
  6. arbitrary precision (any MPFR leaf, all leaves numeric): the arb Slater
     engine (flint_num_meijerg), at the minimum input precision;
  7. z = 0: 0, a finite value, or ComplexInfinity where the expansion decides;
  8. automatic reduction to simpler functions (see mg_try_reduce);
  9. otherwise unevaluated.
```

Memory contract: the builtin takes ownership of res but must not free it.

## Implementation notes

**Algorithm.** `builtin_meijer_g` runs the following steps in order.

1. **Validate.** The call must have 3 or 4 arguments, and each parameter argument must be a pair of lists. `MeijerG::hdiv` and `MeijerG::argt` are emitted through `mth_message`. When r == 1 the call drops to the 3-argument form.
2. **Buffer z.** An `NDArray` or packed z goes to `mg_ndarray`, which builds the plan once and runs `nd_parallel_for`. A plain `List` z is threaded element by element.
3. **Canonical order.** `mg_canonical` sorts each of the four parameter groups.
4. **Consistency.** `mg_collision` rejects m = n = 0, or `a_k − b_j ∈ ℤ⁺` with k ≤ n and j ≤ m. Before that rejection, `mg_move_pair` moves a cancelling Gamma pair from one group to the other.
5. **Machine precision.** `mg_machine` calls `sf_mg_eval`. The evaluation is a Slater sum in w = z^(1/r) (when p < q, or p = q with |w| < 1) or in 1/w. Each term is a regularized pFq from the shared `sf_machine_pfq` / `sf_machine_pfq_c` kernels, multiplied by the z-independent coefficients `∏ π/sin(π(b_j−b_h)) ∏Γ(1+b_h−a_j)/∏Γ(a_j−b_h)`.
   - **Log case.** When two of b_1…b_m are an integer apart, `build_log_table` tabulates the residue series `Σ w^(β+k) Σ_j c[k][j] (log w)^j` up to `SF_MG_LOG_KMAX = 256` terms. It takes the Laurent seeds from arb and steps them exactly in double-double.
   - **Declines.** Too much cancellation between terms, or a kernel decline, falls back to arb at the same point.
6. **Arbitrary precision.** `mg_arbitrary` calls `flint_num_meijerg`, which uses regularized `acb_hypgeom_pfq`. The log case is handled by perturbing `b_j → b_j + (j+1)ε`.
7. **z = 0.** `mg_at_zero` handles this point.
8. **Automatic reduction.** `mg_try_reduce` applies the rules in this order:
   - `G^{2,0}_{0,2}` → `BesselK`.
   - `G^{2,1}_{1,2}` → `HypergeometricU`.
   - Any one-term Slater expansion reduces.
   - A two-term expansion reduces when its terms are 0F1, 1F1 or 2F1 and there is no log case.
   - Each pFq is then named by `mg_term_expr`: `Exp`, a power, `BesselJ`/`BesselI`, 1F1(u;2u) → `BesselI`, or `EllipticK`/`EllipticE`.

`D` returns another MeijerG. `meijerg_series_terms` provides the series at 0 and at Infinity, including the log case through exact Laurent expansion of the Gamma factors.

**Data structures.**
- `MgView` is a flat view of the a and b parameters with the group sizes n, p, m and q.
- `MgExp` holds one expansion: in z, or the dual with A = 1 − b and B = 1 − a.
- `MgTerm` holds one term's coefficient, its exponent, its pFq parameters and the flags for regularized lower parameters.
- On the machine side, an `SfMgPlan` stores per-term coefficients and parameters in fixed arrays of at most `SF_MG_MAXP = 12`, plus the optional log-case table.
- Plans are memoised in a 16-slot table, guarded by a mutex in threaded builds.
- `Compile[]` lowers MeijerG through a flattened n-ary layout (`ndk_MeijerG_n`, `OP_KERNN`, and the strip-mined `OP_VKERNN_*`) at both scalar and rank-1 shapes.

**Complexity / limits.** Building a plan costs O(m·(p+q)) Gamma evaluations, and each point after that is m pFq evaluations. A log-case table is 256·(class size) coefficients and is built once per parameter set. The following are not implemented:
- the exponential (Stokes) asymptotics at infinity for p < q;
- p = q ≥ 3 exactly on |w| = 1, which stays unevaluated;
- the general-term `SeriesCoefficient` at z = 1.

Reduced forms carry `HypergeometricPFQ[...]/Gamma[...]` because Mathilda has no regularized-pFq heads. On a negative real z where G is real, the machine path can leave an imaginary residue of about 1e-15.

**Attributes:** `NumericFunction`, `Protected`.

## References

**See also:** [Quiet](../../control-flow/Quiet/), [Check](../../control-flow/Check/), [Exp](../../elementary-functions/Exp/), [BesselJ](../../special-functions/BesselJ/), [BesselI](../../special-functions/BesselI/), [EllipticK](../../special-functions/EllipticK/), [EllipticE](../../special-functions/EllipticE/), [HypergeometricPFQ](../../special-functions/HypergeometricPFQ/)

- L. J. Slater, *Generalized Hypergeometric Functions* (Cambridge University Press, 1966) — Slater's theorem, the residue-sum expansion of the G-function.
- DLMF §16.17 — the Meijer G-function, its Mellin–Barnes definition and expansions in pFq.
- Source: [`src/special_functions/meijer_g.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/meijer_g.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_meijerg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_meijerg.c)

## Notes & additional examples

### Notes

- Attributes: `NumericFunction`, `Protected`. It is not `Listable`, because its first two arguments are lists, but it threads over a `List` given as z.
- Canonical form: the parameters are sorted within each of the four groups. A pole collision (`a_k − b_j` a positive integer, k ≤ n, j ≤ m), or m = n = 0, emits `MeijerG::hdiv` and the expression stays unevaluated. A cancelling Gamma pair is moved between groups, as Mathematica does.
- Automatic reduction comes from Slater's theorem. `G^{2,0}_{0,2}` becomes `BesselK`, and `G^{2,1}_{1,2}` becomes `HypergeometricU`. A one-term expansion always reduces. A two-term expansion reduces when its terms are 0F1, 1F1 or 2F1 and no log case arises. Each pFq is then named where possible (`Exp`, a power, `BesselJ`/`BesselI`, `EllipticK`/`EllipticE`).
- Machine numerics use `sf_machine_meijerg.c`: a Slater sum of regularized pFq terms in w = z^(1/r) or in 1/w, with the coefficients precomputed once per parameter set. When two b's are an integer apart (the log case), a tabulated residue series is used, with a double-double re-sum. Cancellation past the budget sends the point to arb. Packed lists, `NDArray` and `Compile[]` (scalar and rank-1) share the plan.
- With precision above machine, `flint_num_meijerg` (arb) is used. The log case is handled by perturbing the parameters.
- A real z < 0 on a function that is entire there can come back with a spurious imaginary part of about 1e-15, e.g. `MeijerG[{{1/2},{}},{{0},{0}},-2.]`.
- See also: `FunctionExpand`, `HypergeometricPFQ`, `HypergeometricU`.
