---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 3."
---
### Worked examples

```mathematica
In[1]:= Series[Sin[x], {x, 0, 7}]  (* a Taylor series, returned as SeriesData with an O-term *)
In[2]:= Series[1/(Exp[x] - 1), {x, 0, 4}]  (* a Laurent series; the coefficients are Bernoulli numbers over factorials *)
In[3]:= Series[Sqrt[Sin[x]], {x, 0, 5}]  (* a Puiseux series in half-integer powers *)
In[4]:= Series[Sin[x] - x, x -> 0]  (* leading-term form *)
In[5]:= Series[Zeta[s], {s, 1, 1}]  (* Laurent expansion at the pole, with Stieltjes constants *)
In[6]:= Series[LogGamma[x], {x, Infinity, 3}]  (* Stirling's expansion, keeping the additive growth head symbolic *)
In[7]:= Series[ArcSin[x], {x, 1, 1}]  (* a branch point, with Mathematica's (-1)^Floor[...] sheet discriminator *)
In[8]:= Series[Sqrt[a^2 + x], {x, 0, 2}, Assumptions -> a > 0]  (* assumptions simplify Sqrt[a^2] to a in every coefficient *)
In[9]:= Series[ArcTan[a x], {x, Infinity, 1}, Assumptions -> a < 0]  (* the sign of the scale selects the -Pi/2 branch *)
In[10]:= Series[Gamma[x], {x, 0, 1}]  (* limitation: the pole of Gamma at 0 is not expanded, so the call stays unevaluated *)
```

### Notes

- `Series[f, {x, x0, n}]` returns a `SeriesData` object covering Taylor, Laurent
  and Puiseux expansions. At infinity it also gives asymptotic forms, in which a
  symbolic prefactor (`E^(-x^2)`) or an additive head (`Log[x]`) multiplies or
  adds to a series in `1/x`. Apply `Normal` to drop the O-term.
- Attributes: `Protected` only. `Series` does **not** hold its arguments, as in
  Mathematica, so a named spec `s = {x, 0, 3}` works. The expansion variable
  must be a symbol; anything else emits `Series::ivar`.
- `Assumptions -> assm` (or an ambient `Assuming` / `$Assumptions`) selects
  branches. It decides the `Log[x]` vs `Log[-x]` form at `0`, resolves `Abs`,
  `Sign`, `UnitStep` and `Sqrt[g^2]` by the sign of their argument at `x0`, picks
  the sign of even-order roots under `x < 0`, and picks the branch for scaled
  arguments `c x` at infinity. Uninformative assumptions reproduce the
  assumption-free result exactly.
- `{x, -Infinity, n}` is computed by reflecting to `+Infinity`.
- Dedicated asymptotic expansions at infinity exist for `ExpIntegralEi`,
  `ExpIntegralE[1, x]`, `Erf`/`Erfc`/`Erfi`, `LogGamma`, `PolyGamma`, `Zeta`, the
  sine/cosine/hyperbolic integrals, `FresnelS`/`FresnelC`, the Airy and Bessel
  families, `HermiteH` and `ProductLog`.
- Limitation: `Gamma` is not expanded at its poles or at infinity. Expand
  `Gamma[x + 1]/x` instead of `Gamma[x]` at `0`.
- Limitation: `StieltjesGamma[n]` (from the `Zeta` expansion) has no numerical
  value in this build.
- See also: `Normal`, `SeriesCoefficient`, `SeriesData`, `Limit`, `Residue`,
  `NSeries`.
