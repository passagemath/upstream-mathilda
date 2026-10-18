# PossibleZeroQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PossibleZeroQ[expr]`**

gives True if symbolic and numerical methods suggest that expr has value zero, and False otherwise.

**`PossibleZeroQ[expr, Assumptions -> assum]`**

tests under the assumptions assum; it also respects an ambient Assuming\[\] / $Assumptions scope. Assumptions restrict the numeric sampler to the assumed region (integer, real or complex domain, sign, range, and Re/Im-part constraints), so identities that hold only there are recognised.

<details>
<summary>Notes</summary>

The general problem of deciding whether an expression is zero is undecidable; PossibleZeroQ is a quick but not always accurate test.

</details>

## Examples (23)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= PossibleZeroQ[E^(I Pi/4) - (-1)^(1/4)]
Out[1]= True

In[2]:= PossibleZeroQ[(x + 1)(x - 1) - x^2 + 1]
Out[2]= True

In[3]:= PossibleZeroQ[(E + Pi)^2 - E^2 - Pi^2 - 2 E Pi]
Out[3]= True

In[4]:= PossibleZeroQ[E^Pi - Pi^E]
Out[4]= False

In[5]:= PossibleZeroQ[2^(2 I) - 2^(-2 I) - 2 I Sin[Log[4]]]
Out[5]= True

In[6]:= PossibleZeroQ[Sqrt[x^2] - x]
Out[6]= False

In[7]:= PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]
Out[7]= True

In[8]:= Assuming[Re[x] > 0, PossibleZeroQ[-x + Sqrt[x^2]]]
Out[8]= True
```

### Options (4)

```mathematica
In[9]:= PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x >= 0]
Out[9]= True

In[10]:= PossibleZeroQ[Sin[n Pi], Assumptions -> Element[n, Integers]]
Out[10]= True

In[11]:= PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> -Pi < Im[z] <= Pi]
Out[11]= True

In[12]:= PossibleZeroQ[Cos[n Pi], Assumptions -> Element[n, Integers]]
Out[12]= False
```

### Worked examples (1)

```mathematica
In[13]:= Gamma[x] - Gamma[x]
Out[13]= 0
```

### Applications (10)

A rational identity, decided exactly by Together and Cancel

```mathematica
In[14]:= PossibleZeroQ[(x - 1) (x + 1) - (x^2 - 1)]
Out[14]= True
```

A nested-radical denesting, confirmed on the MPFR ladder

```mathematica
In[15]:= PossibleZeroQ[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]
Out[15]= True
```

Constant-phase angle expansion in the exact trig/exp test

```mathematica
In[16]:= PossibleZeroQ[Tan[Pi/2 - x] - Cot[x]]
Out[16]= True
```

Huge operands: the cancellation is judged by its trend up to 1000 bits

```mathematica
In[17]:= PossibleZeroQ[Cosh[14]^30 - (1 - Tanh[14]^2)^(-15)]
Out[17]= True
```

The same scale with a genuine +1 is rejected

```mathematica
In[18]:= PossibleZeroQ[Cosh[14]^30 - (1 - Tanh[14]^2)^(-15) + 1]
Out[18]= False
```

Numeric-function declines are re-drawn, not taken as undecidable

```mathematica
In[19]:= PossibleZeroQ[PolyLog[2, x^2] + 1]
Out[19]= False
```

Structural certificate: an exponential never vanishes

```mathematica
In[20]:= PossibleZeroQ[-E^(-a x)]
Out[20]= False
```

Sampling is confined to the assumed region

```mathematica
In[21]:= PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x > 0]
Out[21]= True
```

A coupling inequality is enforced by rejection sampling

```mathematica
In[22]:= PossibleZeroQ[{Max[x, y] - x, Max[x, y] - y}, Assumptions -> x > y]
Out[22]= {True, False}
```

Zero on the integer lattice only

```mathematica
In[23]:= PossibleZeroQ[Exp[2 Pi I k] - 1, Assumptions -> Element[k, Integers]]
Out[23]= True
```

## Algorithm

zero_test.c — PossibleZeroQ: hybrid symbolic-numeric zero recognition.

Pipeline (early exit at any stage that yields a definite verdict):

```text
  Stage 0.5 (normalisation) — when the input carries a constant-base
            exponential with a NON-LINEAR symbol-dependent exponent
            (Power[E, f] / Exp[f], f super-linear in a free symbol), ExpandAll
            distributes the sums so same-base exponentials
            become adjacent Times factors and collapse (E^a * E^(-a) -> 1).
            Value-preserving, so the verdict is unchanged; it spares the numeric
            ladder the tiny*huge catastrophic cancellation of a same-base
            exponential product split across summands (a Gaussian * Erf residual
            of an exact-ODE integrating-factor solution otherwise hangs the
            ladder).

  Stage 0 — O(1) structural shortcuts: literal Integer/Real/BigInt/MPFR
            zero, Complex[0, 0], List of zeros, unbound symbol, …

  Stage 1 — Rational normalisation via Together ∘ Cancel + Expand,
            then is_zero_poly. Decides every identity in Q(x_1,…,x_n).

  Stage 2 — Closed-form numeric path: numericalize at machine precision,
            compare |z| against the IEEE catastrophic-cancellation
            ambiguity threshold scale * 2^(-p/2 + 4). If ambiguous, bump
            precision (machine → 200 → 500 → 1000 bits) and retry.
            A non-zero result stays roughly constant across precisions;
            a true zero shrinks geometrically. Surviving the full ladder
            implies "True".

  Stage 3 — Schwartz–Zippel. For inputs with free symbols, substitute
            each free symbol with a random REAL rational of moderate
            magnitude, recurse into Stage 2, and require independent
            confirmations. Sampling is real-line only: an analytic identity
            holding on a real interval holds on a complex neighbourhood
            (identity theorem), so real points confirm it, while complex
            samples needlessly cross branch cuts (Log/ArcTan/Sqrt) and blow
            up special functions (Gamma), manufacturing false negatives.
            The draw stream is seeded deterministically from the input's
            structural hash, so the verdict is a pure function of the input
            (no run-to-run flakiness) and the user's RNG stream is untouched.
```

## Implementation notes

**Algorithm.** `builtin_possible_zero_q` (`src/zero_test.c`) threads a `List`
first argument by hand and broadcasts the remaining arguments unchanged, so a
list of assumptions is never threaded against the expressions. The head is
not `Listable`. Next it forms the effective assumption: an
`Assumptions -> X` option overrides `$Assumptions`, and a positional
assumption is `And`-ed with it. With no usable facts it calls
`zero_test_decide`; otherwise it builds an `AssumeCtx` and calls
`zero_test_decide_assuming`. Both entries first apply
`zt_normalize_exp_kernels`, which combines exponentials so that
Gaussian × `Erf` style products cannot hang the later stages. They then run
a staged pipeline that exits on the first definite verdict:

- *Stage 0, structural* (`decide_structural`): literal zeros, `Complex[0,0]`
  and nonzero constants are settled at once. `provably_nonzero` then
  certifies shapes such as `c E^f`, which can never vanish, so `-E^(-a x)`
  is `False` even though sampling sees only underflow.
- *Guards.* A free symbol mixed with an algebraic constant, or a symbolic
  exponent, skips straight to sampling, because `Together ∘ Cancel` over
  `Q(α)` can fail to terminate.
- *Stage 1, rational* (`decide_rational`): `Together ∘ Cancel` and `Expand`
  decide identities in `Q(x_1, …, x_n)`. `False` is trusted only on the pure
  rational-function path.
- *Stage 1b, exact trig/exp* (`trigexp_rational_is_zero`): a rational
  function of one kernel `t = E^(i x)` is decided by exact evaluation on a
  Nullstellensatz grid. A constant phase `k x + c Pi` is angle-expanded
  first (`expand_pi_phase`), so `Tan[Pi/2 - x] - Cot[x]` comes out `True`
  rather than a spurious `False`.
- *Stage 2, numeric* (`decide_numeric`): a symbol-free input is evaluated at
  machine precision and checked against a cancellation threshold, then up an
  MPFR ladder (53 → 200 → 500 → 1000 bits) while the verdict is ambiguous.
  Once the operand scale passes `ZT_HUGE_SCALE_BITS` (2^40), the machine
  rejection is not trusted and `decide_numeric_huge_scale` reads the trend
  of the residual across the ladder: a plateau means nonzero, shrinking to
  the top rung means zero.
- *Stage 3, Schwartz–Zippel* (`decide_schwartz_zippel_core`): each free
  symbol is replaced by a random real point with |v| ≥ 1, or by a point
  drawn from a per-symbol `SampleSpec` (real, integer or complex channels,
  with bounds) when assumptions apply. `ZT_SCREEN_SAMPLES` = 24 screening
  points are followed by `ZT_CONFIRM_SAMPLES` = 4 confirmations. IEEE
  overflow and numeric-function declines
  (`residue_is_numeric_decline`: `PolyLog[2, z]` returning itself for a
  large argument) are re-drawn on smaller magnitude shells, not
  treated as undecidable. An inequality that couples two or more symbols
  (`x > y`) is enforced by rejection sampling (`assignment_conforms`, up to
  `ZT_CONFORM_MAX_TRIES` = 128 tries), so every test point lies inside the
  region. A measure-zero fact (`a == b`, `Or`) cannot be sampled into, so
  under such a fact a `False` is downgraded to unknown.

Under assumptions, only `True` from the exact stages is trusted: an
unconditional `False` can be wrong on an integer lattice, as with
`E^(2 Pi I k) - 1`. Constrained sampling alone settles `False`.

**Result mapping.** `False` is returned only for a definite `False`. Both
`True` and unknown give `True`, the documented "assume zero when
uncertain" behaviour, and no `PossibleZeroQ::ztest1` message is emitted.
The attribute set is just `Protected`.

**Complexity / limits.** The cost is a few symbolic normalisations plus at
most about 28 sampled evaluations, with MPFR only on ambiguous rungs.
Without assumptions the free symbols are sampled on the real line. A
function that is not analytic in its argument (`Conjugate`, `Re`, `Im`)
therefore looks like the identity there: `PossibleZeroQ[Conjugate[x] - x]`
returns `True`. Supply an assumption such as `Element[x, Complexes]`, or
expand with `ComplexExpand`, before relying on that verdict.

**Attributes:** `Protected`.

## References

**See also:** [ExpandAll](../../structural-manipulation/ExpandAll/), [Times](../../arithmetic/Times/), [Erf](../../special-functions/Erf/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [Together](../../algebra/Together/), [Cancel](../../algebra/Cancel/), [Expand](../../algebra/Expand/)

- J. T. Schwartz, "Fast probabilistic algorithms for verification of polynomial identities", JACM 27 (1980).
- R. Zippel, "Probabilistic algorithms for sparse polynomials", EUROSAM 1979.
- Source: [`src/zero_test.c`](https://github.com/stblake/mathilda/blob/main/src/zero_test.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_divisorsigma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisorsigma.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m5_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m5_stress.c)
- Tests: [`tests/test_dsolve_m62_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m62_stress.c)

## Notes & additional examples

### Notes

- **Heuristic, as in Mathematica.** `PossibleZeroQ` runs exact symbolic
  stages first (rational normalisation, the exact trig/exp kernel test,
  structural non-zero certificates). Then come numeric evaluation on an MPFR
  precision ladder and Schwartz–Zippel sampling at random points. `False`
  means a definite non-zero was found. `True` means the expression was
  proved zero, or behaved like zero at every point tested; when unsure it
  answers `True` and prints no message.
- **Assumptions.** `PossibleZeroQ[e, Assumptions -> a]` and an enclosing
  `Assuming` restrict the sample points to the assumed region. That covers
  sign and interval bounds, `Element[x, Integers]`,
  `Element[x, Complexes]`, and inequalities coupling several variables. An
  option replaces `$Assumptions`, and a positional assumption is conjoined
  with it.
- **Real sampling by default.** Without assumptions, free symbols are
  sampled as reals. Functions that are not analytic, such as `Conjugate`,
  `Re` and `Im`, then look like the identity:
  `PossibleZeroQ[Conjugate[x] - x]` is `True`. Pass
  `Assumptions -> Element[x, Complexes]` to sample off the real line, which
  gives `False`.
- **Threading.** A list in the first argument threads with the options
  broadcast. The head is `Protected` and not `Listable`, so a list
  assumption is read as one conjunction.
- See also: `Simplify`, `Refine`, `Equal`, `ComplexExpand`.
