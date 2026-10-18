### Worked examples

```mathematica
In[1]:= Residue[1/z, {z, 0}]  (* a simple pole at the origin *)
In[2]:= Residue[1/(z^2 + 1), {z, I}]  (* simple pole at a complex location, via P(z0)/Q'(z0) *)
In[3]:= Residue[1/(z^4 + 1), {z, Exp[I Pi/4]}]  (* an algebraic pole location, kept in one radical spelling *)
In[4]:= Residue[1/Sin[z]^5, {z, 0}]  (* fifth-order pole of a transcendental integrand *)
In[5]:= Residue[Pi Cot[Pi z]/(z^2 + 1), {z, I}]  (* the summation kernel Pi Cot[Pi z] at a pole of the rational factor *)
In[6]:= Residue[Zeta[s], {s, 1}]  (* residue 1 at the only pole of Zeta *)
In[7]:= Residue[Exp[I z]/(z^2 + a^2)^2, {z, I a}]  (* symbolic double pole, self-validated across expansion orders *)
In[8]:= Residue[z^(1/3)/(1 + z^2)^2, {z, I}]  (* fractional power times a double pole, via the dropped-pole guard *)
In[9]:= Residue[1/Sqrt[z], {z, 0}]  (* a branch point has no residue, so the call stays unevaluated *)
In[10]:= Residue[Gamma[z], {z, -3}]  (* limitation: Series has no expansion of Gamma at its poles, so this stays unevaluated; NResidue gives -0.166667 *)
```

### Notes

- `Residue[f, {z, z0}]` is the coefficient of `(z - z0)^-1` in the Laurent
  expansion of `f` at `z0`. An analytic point gives `0`. Poles of any order, at
  numeric or symbolic centres, are handled.
- For higher-order poles the coefficient is accepted only after two consecutive
  expansion orders agree. This closes the silent wrong answers at symbolic
  double and triple poles that existed before v0.268. Those wrong residues also
  fed `Integrate[..., Method -> "Residue"]`.
- A branch point (a fractional-power expansion) is left unevaluated, matching
  Mathematica. An essential singularity (`E^(1/z)` at `0`) and the poles of
  `Gamma` are also left unevaluated, because `Series` has no expansion there.
  Mathematica returns `1` and `(-1)^n/n!` for these. Use `NResidue` for a
  numerical value.
- Only the two-argument form is supported. Attributes: `Protected`.
- See also: `NResidue`, `Series`, `SeriesCoefficient`, `Integrate` (contour /
  residue-theorem definite integration).
