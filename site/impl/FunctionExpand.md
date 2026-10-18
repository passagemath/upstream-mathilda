---
source:
  - src/special_functions/functionexpand.c
  - src/special_functions/meijer_g.c
references:
  - "L. J. Slater, *Generalized Hypergeometric Functions* (Cambridge University Press, 1966) — the Slater expansion of the G-function in pFq."
---
**Algorithm.** `builtin_function_expand` accepts one or two arguments and calls `fe_walk` on the first. `fe_walk` copies the tree bottom-up. Wherever a child changed, it re-evaluates the rebuilt node with `eval_and_free`. Then, if the node has head `MeijerG`, it replaces the node with `meijerg_slater_expand(node, named = true)`. That function works as follows:

- `mg_pick` chooses the expansion. It uses z when p < q and 1/z when p > q. When p = q it uses the shorter expansion, or the one that converges at a numeric z, and it declines when m + n ≤ p.
- `mg_log_case` declines when two of the b's are an integer apart.
- Otherwise `mg_sum_expr` builds the sum of terms `C_h · z^(exponent) · pFq`. Each pFq is named where a classical form exists: `Exp`, a power, Bessel functions, `EllipticK`/`EllipticE`. Anything else stays `HypergeometricPFQ[...]` divided by the Gamma of its regularized lower parameters.

A node that does not expand comes back unchanged. The assumptions argument is accepted and ignored.

**Data structures.** The walk is a plain recursive copy of the `Expr` tree. It allocates a temporary child array per node and threads a `changed` flag upward. The expansion uses the `MgView`, `MgExp` and `MgTerm` structures from `meijer_g.c`, the same ones the automatic reduction in `MeijerG` uses.

**Complexity / limits.** The walk is linear in tree size, plus one evaluation per changed node. MeijerG is the only head with an expansion rule. `Gamma`, `PolyGamma`, Bessel functions at half-integer order, `StruveH`, `WhittakerM`/`WhittakerW` and `HypergeometricU` all pass through unchanged. So does the log case of `MeijerG`, and so does p = q with m + n ≤ p. A packed or `NDArray` buffer passes straight through.
