# FunctionExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FunctionExpand[expr]`**

expands special functions in expr into simpler ones. Currently each

**`FunctionExpand[expr, assum]`**

accepts assumptions (not yet used).

<details>
<summary>Notes</summary>

MeijerG is rewritten to its Slater expansion: a sum of HypergeometricPFQ terms in z (p \<= q) or 1/z (p \> q), with classical functions (Exp, powers, Bessel, EllipticK/E) where a term has one. The log case (parameters an integer apart) and p == q with m + n == p are left unchanged.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FunctionExpand[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},z]]
Out[1]= 1/12 Gamma[1/5] Sqrt[Pi] HypergeometricPFQ[{1/5}, {1/2, 3, 4}, 1/z] - (256 Gamma[7/10] HypergeometricPFQ[{7/10}, {3/2, 7/2, 9/2}, 1/z])/(1575 Sqrt[Pi] Sqrt[z])

In[2]:= FunctionExpand[Sin[x] + MeijerG[{{},{}},{{0},{}},z]]
Out[2]= Sin[x] + E^(-z)

In[3]:= FunctionExpand[MeijerG[{{},{2}},{{1/2,3/2},{}},z]]
Out[3]= MeijerG[{{}, {2}}, {{1/2, 3/2}, {}}, z]
```

### Applications (7)

A MeijerG that does not reduce on its own becomes its Slater sum of 2F2s

```mathematica
In[4]:= FunctionExpand[MeijerG[{{1},{2}},{{1/5,2/3},{4}},z]]
Out[4]= (Gamma[-7/15] Gamma[2/3] z^(2/3) HypergeometricPFQ[{2/3, -1/3}, {22/15, -7/3}, -z])/(Gamma[-7/3] Gamma[4/3]) + (Gamma[1/5] Gamma[7/15] z^(1/5) HypergeometricPFQ[{1/5, -4/5}, {8/15, -14/5}, -z])/(Gamma[-14/5] Gamma[9/5])
```

The expansion agrees with the arb value of the MeijerG

```mathematica
In[5]:= {N[FunctionExpand[MeijerG[{{1},{2}},{{1/5,2/3},{4}},z]] /. z -> 2, 20], N[MeijerG[{{1},{2}},{{1/5,2/3},{4}},2], 20]}
Out[5]= {-2.0758904683947045901, -2.0758904683947045901}
```

With p > q the expansion is in 1/z

```mathematica
In[6]:= FunctionExpand[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},z]]
Out[6]= 1/12 Gamma[1/5] Sqrt[Pi] HypergeometricPFQ[{1/5}, {1/2, 3, 4}, 1/z] - (256 Gamma[7/10] HypergeometricPFQ[{7/10}, {3/2, 7/2, 9/2}, 1/z])/(1575 Sqrt[Pi] Sqrt[z])
```

Acts on subexpressions; other heads are left alone

```mathematica
In[7]:= FunctionExpand[Sin[x] + MeijerG[{{},{}},{{0},{}},z]]
Out[7]= Sin[x] + E^(-z)
```

The log case has no pFq Slater sum and stays

```mathematica
In[8]:= FunctionExpand[MeijerG[{{},{2}},{{1/2,3/2},{}},z]]
Out[8]= MeijerG[{{}, {2}}, {{1/2, 3/2}, {}}, z]
```

Only MeijerG is expanded so far

```mathematica
In[9]:= FunctionExpand[Gamma[1/2 + n]]
Out[9]= Gamma[1/2 + n]
```

Likewise: no Struve elementary forms yet

```mathematica
In[10]:= FunctionExpand[StruveH[1/2, z]]
Out[10]= StruveH[1/2, z]
```

## Algorithm

functionexpand.c -- FunctionExpand for Mathilda (minimal).

```text
  FunctionExpand[expr]
  FunctionExpand[expr, assum]
```

expands special functions in expr into simpler ones. In this first version the only rule is MeijerG -> its Slater expansion (src/special_functions/ meijer_g.c, meijerg_slater_expand): a sum of HypergeometricPFQ terms in z (p < q, or p == q with the shorter expansion) or in 1/z (p > q), each pFq rewritten to a classical function where one exists. A MeijerG with no such expansion -- the log case of integer-spaced parameters, or p == q with m + n == p, where the two expansions are not continuations of each other -- is left as it is, as is everything else. The assumptions argument is accepted for compatibility and not yet used.

Memory contract: the builtin takes ownership of res but must not free it.

## Implementation notes

**Algorithm.** `builtin_function_expand` accepts one or two arguments and calls `fe_walk` on the first. `fe_walk` copies the tree bottom-up. Wherever a child changed, it re-evaluates the rebuilt node with `eval_and_free`. Then, if the node has head `MeijerG`, it replaces the node with `meijerg_slater_expand(node, named = true)`. That function works as follows:

- `mg_pick` chooses the expansion. It uses z when p < q and 1/z when p > q. When p = q it uses the shorter expansion, or the one that converges at a numeric z, and it declines when m + n ≤ p.
- `mg_log_case` declines when two of the b's are an integer apart.
- Otherwise `mg_sum_expr` builds the sum of terms `C_h · z^(exponent) · pFq`. Each pFq is named where a classical form exists: `Exp`, a power, Bessel functions, `EllipticK`/`EllipticE`. Anything else stays `HypergeometricPFQ[...]` divided by the Gamma of its regularized lower parameters.

A node that does not expand comes back unchanged. The assumptions argument is accepted and ignored.

**Data structures.** The walk is a plain recursive copy of the `Expr` tree. It allocates a temporary child array per node and threads a `changed` flag upward. The expansion uses the `MgView`, `MgExp` and `MgTerm` structures from `meijer_g.c`, the same ones the automatic reduction in `MeijerG` uses.

**Complexity / limits.** The walk is linear in tree size, plus one evaluation per changed node. MeijerG is the only head with an expansion rule. `Gamma`, `PolyGamma`, Bessel functions at half-integer order, `StruveH`, `WhittakerM`/`WhittakerW` and `HypergeometricU` all pass through unchanged. So does the log case of `MeijerG`, and so does p = q with m + n ≤ p. A packed or `NDArray` buffer passes straight through.

**Attributes:** `Protected`.

## References

**See also:** [MeijerG](../../special-functions/MeijerG/), [HypergeometricPFQ](../../special-functions/HypergeometricPFQ/)

- L. J. Slater, *Generalized Hypergeometric Functions* (Cambridge University Press, 1966) — the Slater expansion of the G-function in pFq.
- Source: [`src/special_functions/functionexpand.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/functionexpand.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_meijerg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_meijerg.c)

## Notes & additional examples

### Notes

- Attribute: `Protected` only.
- So far there is one rule: each `MeijerG` subexpression is replaced by its Slater expansion (`meijerg_slater_expand`, with classical functions named where possible). The expansion is in z when p < q and in 1/z when p > q. When p = q it is in z unless m + n ≤ p; in that case the two expansions are not continuations of each other and the node is left alone. The tree is walked bottom-up and re-evaluated where it changed.
- A log-case `MeijerG` (two b's an integer apart) is returned unchanged, and so is every other head: `Gamma`, Bessel half-integer orders, `StruveH`, `WhittakerM` and so on.
- `FunctionExpand[expr, assum]` accepts assumptions but does not use them yet.
- Unlike Mathematica's `FunctionExpand`, it does not rewrite `Gamma`, `PolyGamma`, `Zeta` or the elementary special values.
- See also: `MeijerG`, `HypergeometricPFQ`, `Simplify`.
