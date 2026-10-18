### Worked examples

```mathematica
In[1]:= MeijerG[{{},{}},{{0},{}},z]  (* the simplest G-function, G^{1,0}_{0,1}, is the exponential *)
```

```mathematica
In[1]:= MeijerG[{{1,1},{}},{{1},{0}},z]  (* a one-term Slater sum that names its 2F1: the logarithm *)
```

```mathematica
In[1]:= MeijerG[{{1/2,1/2},{}},{{0},{0}},z]  (* a 2F1 with parameters 1/2, 1/2; 1 is named as EllipticK *)
```

```mathematica
In[1]:= MeijerG[{{1},{2}},{{2/3,1/5},{4}},z]  (* parameters sorted within each group; two 2F2 terms do not reduce *)
```

```mathematica
In[1]:= MeijerG[{{},{2}},{{1/2,3/2},{}},3.]  (* the log case, b's an integer apart: machine kernel *)
```

```mathematica
In[1]:= N[Sqrt[3] Exp[-3] HypergeometricU[1/2, 0, 3], 30] - N[MeijerG[{{},{2}},{{1/2,3/2},{}},3], 30]  (* checked against G^{2,0}_{1,2}(z|a;b,c) = z^b e^(-z) U(a-c, 1+b-c, z) *)
```

```mathematica
In[1]:= MeijerG[{{3+I},{2.}},{{2/3,1/5},{4}},3]  (* complex parameters *)
```

```mathematica
In[1]:= MeijerG[{{},{}},{{0},{}},4., 2]  (* the generalized form: MeijerG[a, b, z, r] = MeijerG[a, b, z^(1/r)], so e^-2 *)
```

```mathematica
In[1]:= D[MeijerG[{{1},{2}},{{1/5,2/3},{4}},z],z]  (* the derivative is another MeijerG *)
```

```mathematica
In[1]:= MeijerG[{{2},{}},{{1},{}},z]  (* pole collision, a_1 - b_1 a positive integer: MeijerG::hdiv, unevaluated *)
```

### Notes

- Attributes: `NumericFunction`, `Protected`. It is not `Listable`, because its first two arguments are lists, but it threads over a `List` given as z.
- Canonical form: the parameters are sorted within each of the four groups. A pole collision (`a_k − b_j` a positive integer, k ≤ n, j ≤ m), or m = n = 0, emits `MeijerG::hdiv` and the expression stays unevaluated. A cancelling Gamma pair is moved between groups, as Mathematica does.
- Automatic reduction comes from Slater's theorem. `G^{2,0}_{0,2}` becomes `BesselK`, and `G^{2,1}_{1,2}` becomes `HypergeometricU`. A one-term expansion always reduces. A two-term expansion reduces when its terms are 0F1, 1F1 or 2F1 and no log case arises. Each pFq is then named where possible (`Exp`, a power, `BesselJ`/`BesselI`, `EllipticK`/`EllipticE`).
- Machine numerics use `sf_machine_meijerg.c`: a Slater sum of regularized pFq terms in w = z^(1/r) or in 1/w, with the coefficients precomputed once per parameter set. When two b's are an integer apart (the log case), a tabulated residue series is used, with a double-double re-sum. Cancellation past the budget sends the point to arb. Packed lists, `NDArray` and `Compile[]` (scalar and rank-1) share the plan.
- With precision above machine, `flint_num_meijerg` (arb) is used. The log case is handled by perturbing the parameters.
- A real z < 0 on a function that is entire there can come back with a spurious imaginary part of about 1e-15, e.g. `MeijerG[{{1/2},{}},{{0},{0}},-2.]`.
- See also: `FunctionExpand`, `HypergeometricPFQ`, `HypergeometricU`.
