---
source:
  - src/special_functions/meijer_g.c
  - src/special_functions/sf_machine_meijerg.c
  - src/special_functions/sf_machine_meijerg.h
  - src/flint_num_bridge.c
  - src/compile/compile_infer.c
  - src/ndkernels.c
references:
  - "L. J. Slater, *Generalized Hypergeometric Functions* (Cambridge University Press, 1966) — Slater's theorem, the residue-sum expansion of the G-function."
  - "DLMF §16.17 — the Meijer G-function, its Mellin–Barnes definition and expansions in pFq."
---
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
