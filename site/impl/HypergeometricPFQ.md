---
source: src/special_functions/hypergeopfq.c
references:
  - "DLMF §16 — generalized hypergeometric functions pFq."
  - "T. J. Dekker, A floating-point technique for extending the available precision, Numer. Math. 18 (1971) 224-242 — the double-double re-sum."
  - "Y. Hida, X. S. Li and D. H. Bailey, Algorithms for quad-double precision floating point arithmetic, ARITH-15 (2001)."
---
**Algorithm.** `builtin_hypergeometric_pfq` handles `pFq(a;b;z)`, the series
`Sum_k (prod_i (a_i)_k / prod_j (b_j)_k) z^k / k!`. Dispatch order: validate
(arity 3, first two arguments `List`); a packed / `NDArray` `z` goes to
`hg_pfq_ndarray` (else `ndarray_delist_and_reeval`); a plain `List` `z` is
threaded by hand; then the **machine gate** `hg_machine_pfq` — every parameter and
`z` a machine-scalar leaf (Integer, Rational, Real, or `Complex` of those) with at
least one machine `Real` — answers from the raw-`double` kernel *before* any exact
rewrite, so exact parameters with a machine `z` take the same path and bits as the
buffer and `Compile[]` surfaces. Only after that come the exact steps: **generic
cancellation** of a parameter common to both lists (`try_cancel`), `z == 0 -> 1`,
**termination** to the explicit polynomial for a non-positive integer upper
parameter (`try_terminate`), **reduction** to closed forms (`try_reduce`: `0F0 =
Exp`, `1F0 = (1-z)^-a`, `0F1(1/2)`/`0F1(3/2)` as `Cosh`/`Sinh`, `1F1(1;2;z)`, the
`2F1` logarithm and `ArcSin` families, and a table-backed very-well-poised `4F3`
Ramanujan `1/Pi` value), and finally `try_numeric` for non-machine inexact input.
There the governing precision is the *minimum* inexact precision of `z` and every
parameter; above 53 bits a `2F1` with `|z| < 1` goes to arb
(`flint_num_hyp2f1`, `acb_hypgeom_2f1`), anything else to the MPFR `(re, im)`
series `hgpfq_mpfr_sum`. A machine partial sum that overflows promotes the whole
call to MPFR (`numeric_promote_machine_call`).

**Data structures.** The machine kernel lives in
`src/special_functions/sf_machine_hyper.c` (`sf_machine_pfq` / `sf_machine_pfq_c`,
plain C99 `double` with complex values as explicit re/im pairs, no
`<complex.h>`). `sf_machine_pfq` sends `{}/{b}`, `{a}/{b}` and `{a,b}/{c}` to the
`0F1` / `1F1` / `2F1` kernels, so every spelling of one function reaches the same
code; other shapes use a **loss-measuring** series: summed once in double while
tracking its peak term, and re-summed in double-double (`hdd`, Dekker /
Hida–Li–Bailey) when `peak/|sum|` exceeds 2 or the peak lies past the 6th term,
declining beyond `2^45` of cancellation (`HY_LOSS_MAX`). Parameters are collected
in `HgParams` (at most 16 upper and 16 lower). The buffer path runs the vector
driver `sf_machine_pfq_vec` for real parameters over a float64 buffer (the `2F1`
coefficients are tabulated once per buffer and four elements are summed in
lock-step, bit-identical to the scalar loop), and otherwise a general loop that
widens the output dtype to complex when a parameter is complex or a `2F1` element
falls on the cut. ND: n-ary kernel `NDKN_HypergeometricPFQ` (`REG_N`), which
receives `p` as a leading operand; all four hypergeometric heads are on `pack.c`'s
`AWARE` list. Attributes: `NumericFunction`, `Protected` — not `Listable`, since
the first two arguments are themselves lists.

**Complexity / limits.** Series term caps: 20000 (double), 6000 (double-double),
`HGPFQ_MACHINE_MAX_TERMS = 200000` / `HGPFQ_MPFR_MAX_TERMS = 1000000` on the
generic path. The general series is summed only where it converges — `p <= q`
for any `z`, `p == q+1` for `|z| < 1`; other shapes stay unevaluated, except that
the `0F1`/`1F1`/`2F1` shapes have their own continuations (`2F1` is continued to
the whole plane at machine precision; see `Hypergeometric2F1`). `Compile[]`
lowers the head at real and complex scalar shapes (the two literal lists are
destructured in `compile_emit_mathfn.c` into one `OP_KERNN` / `OP_KERNN_C`
operand block, capped at `KERNN_MAX_OPERANDS = 24`); the rank-1 array shape
reports `Compiled -> False`, while the packed / `NDArray` buffer path still runs
the kernel at the REPL.
