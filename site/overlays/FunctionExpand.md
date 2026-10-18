### Worked examples

```mathematica
In[1]:= FunctionExpand[MeijerG[{{1},{2}},{{1/5,2/3},{4}},z]]  (* a MeijerG that does not reduce on its own becomes its Slater sum of 2F2s *)
```

```mathematica
In[1]:= {N[FunctionExpand[MeijerG[{{1},{2}},{{1/5,2/3},{4}},z]] /. z -> 2, 20], N[MeijerG[{{1},{2}},{{1/5,2/3},{4}},2], 20]}  (* the expansion agrees with the arb value of the MeijerG *)
```

```mathematica
In[1]:= FunctionExpand[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},z]]  (* with p > q the expansion is in 1/z *)
```

```mathematica
In[1]:= FunctionExpand[Sin[x] + MeijerG[{{},{}},{{0},{}},z]]  (* acts on subexpressions; other heads are left alone *)
```

```mathematica
In[1]:= FunctionExpand[MeijerG[{{},{2}},{{1/2,3/2},{}},z]]  (* the log case has no pFq Slater sum and stays *)
```

```mathematica
In[1]:= FunctionExpand[Gamma[1/2 + n]]  (* only MeijerG is expanded so far *)
```

```mathematica
In[1]:= FunctionExpand[StruveH[1/2, z]]  (* likewise: no Struve elementary forms yet *)
```

### Notes

- Attribute: `Protected` only.
- So far there is one rule: each `MeijerG` subexpression is replaced by its Slater expansion (`meijerg_slater_expand`, with classical functions named where possible). The expansion is in z when p < q and in 1/z when p > q. When p = q it is in z unless m + n ≤ p; in that case the two expansions are not continuations of each other and the node is left alone. The tree is walked bottom-up and re-evaluated where it changed.
- A log-case `MeijerG` (two b's an integer apart) is returned unchanged, and so is every other head: `Gamma`, Bessel half-integer orders, `StruveH`, `WhittakerM` and so on.
- `FunctionExpand[expr, assum]` accepts assumptions but does not use them yet.
- Unlike Mathematica's `FunctionExpand`, it does not rewrite `Gamma`, `PolyGamma`, `Zeta` or the elementary special values.
- See also: `MeijerG`, `HypergeometricPFQ`, `Simplify`.
