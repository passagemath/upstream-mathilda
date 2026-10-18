---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 3."
---
### Worked examples

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]  (* the classic removable singularity *)
In[2]:= Limit[(1 + a/x)^x, x -> Infinity]  (* a symbolic 1^Infinity form via the f^g = Exp[g Log f] reduction *)
In[3]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]  (* one-sided; the two-sided limit is ComplexInfinity *)
In[4]:= Limit[Floor[x], x -> 2, Direction -> "FromBelow"]  (* one-sided limit of a step function at an integer *)
In[5]:= Limit[x^a Log[x], x -> 0, Direction -> "FromAbove", Assumptions -> a > 0]  (* assumptions decide that the power dominates the logarithm *)
In[6]:= Limit[x (Gamma[a, x] Exp[x] x^(1 - a) - 1), x -> Infinity]  (* sub-leading incomplete-Gamma asymptotics *)
In[7]:= Limit[(Zeta[x] - 1) 2^x, x -> Infinity]  (* the Dirichlet tail of Zeta at infinity *)
In[8]:= Limit[Exp[x] (Exp[1/x - Exp[-x]] - Exp[1/x]), x -> Infinity, Method -> "Gruntz"]  (* nested-exponential cancellation by Gruntz's mrv algorithm *)
In[9]:= Limit[Cos[x] + Sin[Sqrt[2] x], x -> Infinity]  (* two incommensurate oscillations, proved to have no limit *)
In[10]:= Limit[f[x], x -> 0]  (* limitation: no continuity is assumed for an undefined head *)
```

### Notes

- `Limit[f, x -> a]` can return a finite value, `±Infinity`, `ComplexInfinity`,
  `Indeterminate` (for example, a bounded oscillation with no limit), or the
  unevaluated call. The default is two-sided.
  `Direction -> "FromAbove"`/`"FromBelow"` (or `-1`/`+1`, following Mathematica's
  sign convention) selects one side. `Direction -> Reals`, `Complexes` or a
  complex number such as `I` give the other approaches.
- Calling forms: nested `Limit[f, {x -> a, y -> b}]` (the last rule is applied
  innermost) and joint `Limit[f, {x, y} -> {a, b}]`. The joint form returns
  `Indeterminate` when directional limits disagree, as with `x y/(x^2 + y^2)` at
  the origin.
- `Assumptions -> cond`, or an ambient `Assuming` / `$Assumptions`, supplies the
  signs of parameters. It is used for exponents, rates, dominant products and
  `Re`/`Im` decay. Without it the cascade behaves exactly as it does with no
  assumptions.
- `Method` values: `Automatic` (the full cascade), `"Substitution"`,
  `"RationalFunction"`, `"Series"`, `"LHospital"`, `"Asymptotic"` (the
  Log/dominant-term, `f^g` and special-function asymptotics layers), `"Bounded"`
  (squeeze and `Interval`), `"Oscillatory"` (the phase normal form at infinity)
  and `"Gruntz"`. A named method applies only to the outermost call and never
  falls back to the full cascade. Each method is also callable as a head, for
  example ``Limit`Gruntz[f, x -> a]``.
- Inexact input is rationalised, solved exactly, and numericalised back.
- Discontinuous heads at a jump are left unevaluated when the limit is
  two-sided, for example `Limit[Floor[x], x -> 2]`. Ask for a one-sided limit
  instead.
- See also: `Series`, `D`, `Residue`, `NLimit`.
