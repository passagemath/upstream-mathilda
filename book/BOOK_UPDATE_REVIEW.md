# The Mathilda Book — Update Review (git history since the last book pass)

> **STATUS: IMPLEMENTED 2026-10-20** (book now current through v0.481; 560 pages). Window `70e10fc0..ed2d3bd9`, Mathilda
> **v0.290 → v0.480**. The previous review (v0.113–v0.189, implemented 2026-09-25)
> is in git history (`git log -- book/BOOK_UPDATE_REVIEW.md`).

This is a **triage backlog**, not prose. It maps every *book-relevant* capability
that landed after the book's last substantive pass onto the chapter that must cover
it. Read `CONTEXT.md` first; this file is scoped to *the delta*.

## Anchor and window

- **Book anchor:** `70e10fc0` (2026-10-05, v0.290), "§5.3.2 Limits — vertical
  asymptote and finite-limit-at-infinity now prove". Later commits touching `book/`
  (v0.301 Locked removal, v0.367–v0.372 Notebook) were narrow edits, not a pass.
- **Window:** 267 commits, v0.290 → v0.480, 2026-10-05 → 2026-10-20. Sources: the git
  log and `docs/spec/changelog/2026-10-05.md`, `2026-10-12.md`, `2026-10-19.md`.
- **Method:** three triage passes, one per version slice (0.290–0.350, 0.351–0.420,
  0.421–0.480). Each grepped `book/` for every head it named, and all 23 example
  sessions were re-run on v0.480 and diffed against the v0.367 transcripts.

| Tag | Meaning |
|-----|---------|
| **REVIEW** | A behaviour change that falsifies existing prose or an example. Fix first. |
| **INCLUDE** | New capability with no book coverage. |
| **UPDATE** | Extends a section that already covers the area. |

## REVIEW — things the book said that were no longer true

| # | Where | What changed | Fix |
|---|-------|--------------|-----|
| R1 | §5.2 PossibleZeroQ callout, `possiblezeroq-limits/3` | `Simplify` now proves `ArcTan[1/3]+ArcTan[1/2]-Pi/4 == 0` (v0.327, exactly via e^{ie}=1 over the algebraic numbers) | Callout and pair rewritten; timings re-measured |
| R2 | §5.3 `dsolve-special/5` + Kovacic Theory callout | Series answers are opt-in since v0.361; `y''+Sin[x]y==0` returns unevaluated | New pair 6 with `Method -> "FrobeniusSeries"`; prose rewritten |
| R3 | §5.3 `feynman/1` | `[0,1]` Feynman integral now needs `a > -1` proved (FeynmanParameter stage, v0.405) | `Assumptions -> a > -1` added; prose states the condition |
| R4 | Ch. 13 `cliff/4` | `MoebiusMu` compiles since v0.471 | Now `EulerPhi` (compiles) vs `PrimeOmega + EulerPhi` (does not) |
| R5 | §5.7 Gamma / Erf / Ei Under-the-hood notes | Machine complex Gamma, machine complex Erf, machine Ei family now use Mathilda's own double kernels shared with packed + Compile (v0.437–0.441) | Three notes rewritten (real Gamma still `tgamma`) |
| R6 | §5.6 FactorInteger note | Cascade rebuilt (v0.447–0.473) | Folded into INCLUDE A |
| R7 | Ch. 15 evaluation-order figure | UpValues (v0.295) fire before DownValues; mixed exact/inexact retry (v0.477) | Figure now 15 steps; prose updated |
| R8 | Ch. 15 matcher | Flat-head partial rule matching (v0.299, v0.462); OneIdentity only via defaults (v0.475) | Paragraphs rewritten, verified on binary |
| R9 | Ch. 7 Optimization opening + LP paragraph | Exact `Minimize`/`Maximize` (v0.304+) | With INCLUDE B |
| R10 | all transcripts | `Series` at ∞ prints `1/(3 x^3)` (v0.464); version strings, timings | Regenerated; 27 files changed, none other than R3/R4 changing meaning |

Kernel bug found during the pass: a successful definite `Integrate` (e.g.
`∫₀^∞ log(1+a²x²)/(1+x²) dx`, `a>0`) emits a spurious `Integrate::nonelem`, so
`Check[]` returns its failure branch. Fixed separately (see the changelog).

## INCLUDE / UPDATE — by chapter

| Item | Target | Capabilities (versions) |
|------|--------|-------------------------|
| **A** | §5.6 Number theory | FactorInteger rebuild: sub-2^64 engine, Lehman, CFRAC, SIQS + block Lanczos, named method heads, deterministic parallel work, `FactorInteger[n,k]`, Gaussian integers (v0.447–0.473); **Blake Rational Base Descent** and its six variants (v0.452/0.457); `RandomPrime`, `PrimePowerQ`; PrimePi Deleglise–Rivat, `Prime[n]`; partitions engines; NT review batches (Gaussian GCD/PowerMod, `Mod[7,0]`, `CoprimeQ[3]`, Roman `FromDigits`) |
| **B** | Ch. 7 Optimization (+ Ch. 3 pointer) | Exact `Minimize`/`Maximize`: KKT + `Reduce` certificate, `::natt`/`::infeas`, Integers, LP/IP, Abs objectives, separable objectives, `TimeConstraint` (v0.304–0.318) |
| **C** | §5.7 Special functions | Orthogonal polynomials (ChebyshevT/U, HermiteH, LaguerreL); Whittaker M/W, HypergeometricU; MeijerG + FunctionExpand; DirichletCharacter / DirichletL; Nielsen `PolyLog[n,p,z]`; ExpIntegralE, StruveH; Zeta/LerchPhi options; mixed exact+inexact numericalization; machine-precision performance standing (v0.296–0.479) |
| **D** | §5.3 Calculus | Classical definite integrals (King, Frullani/FeynmanParameter, log-sine, series term-by-term, polylog ladders, Holonomic, Meijer-G convolution, integral representations, Cosh^n residues, Bose/Fermi → HurwitzZeta); piecewise integrands; symbolic endpoints under Assumptions; PV extensions; cascade reorder + soundness corpora; Limit and Series under Assumptions; DSolve forcing/Whittaker/InverseOrientation/Riccati/corpus; linear-argument substitution; LogByParts |
| **E** | Ch. 12 Programming | UpValues family (`^=`, `^:=`, `/:`, `=.`, `UpValues`, `Definition`); Flat-head partial matching; `NonCommutativeMultiply`; `Unprotect`/`ClearAttributes` on builtins (v0.480); matcher corpus; TimeConstrained wall-clock + Windows |
| **F** | Ch. 11 Data structures | Strings and string patterns (named captures, Unicode, `\[Name]`, v0.455–0.456); complex64 packing + RandomComplex (v0.402); n-ary special functions over packed arrays |
| **G** | §5.2, Ch. 2, Ch. 13, Ch. 15 | Simplify under assumptions campaign (v0.329–0.350), exact trig at rational π multiples, FullSimplify log-trig, Refine radicals, `Reduce[expr]`, certified transcendental `Reduce` with `Root[{f,x0}]`, SOS certificates; embedding API (`make install`, pkg-config) and native Windows; compile frontier (NT heads, n-ary SF, Chop/Clip); `Hash`, shortest-round-trip `InputForm`, empty slot = `Null` |

## Skipped (no book impact)

Roughly 45 pure bug-fix commits (portability, hangs, edge cases) and roughly 60
internal, tooling, docs, test or benchmark commits. The evaluator speedups (about 10%
on numeric calls) are not worth a sentence.

## Kernel findings from the writing pass (follow-ups, not fixed here)

These were found while building examples and written around. Each was checked on
v0.480/0.481.

**Wrong or misleading output**
- `Limit[x^s Zeta[s,x], x->Infinity, Assumptions->s>1]` returns
  `DirectedInfinity[1/(-1+s)]`, not `Infinity`.
- `PolyLog[0, 2, 2.]` has a spurious imaginary part of about 6e-16.
- `DSolve[{y'[t]+y[t]==UnitStep[t-a], y[0]==0}, y[t], t]` leaks the internal
  ``DSolve`vpS`` into an inert integral.
- `StringCases[..., RegularExpression[...] :> {"$1","$2"}]` does not expand `$n`
  inside a list right-hand side.
- `Definition` shows an `^=` upvalue as `:=`.
- `ClearAttributes[Plus, Orderless]` has no visible effect, because Plus/Times sort
  themselves.

**Budget overruns and hangs**
- `Minimize` of a halved 3-variable Styblinski–Tang runs about 147 s, past its 30 s
  `TimeConstraint`.
- `Integrate[Exp[-x t]/t^2, {t,1,Infinity}, Assumptions->x>0]` hangs (over 120 s)
  with repeated `::idiv`.

**Correct but unsimplified**
- `Minimize` optima such as `14-2Sqrt[13]` come back unsimplified.
- Several definite integrals return MeijerG/PolyLog forms that are correct but
  unreduced: `ArcTan[x]/x` on [0,1] (Catalan), `x Log[2Cos x]`,
  `E^-x Sin[x]/x` (π/4, plus a spurious `Infinity::indet`).

**Declines and gaps**
- 3-variable integer programs decline (`MZ_MAX_INEQ = 6`).
- `D[StruveH[nu,x],x]` is not implemented.
- `DirichletL` with an exact-rational `Interval` stays unevaluated.
- `Boole` does not compile.
- `Chop[x, 10^-6]` with an exact tolerance does not compile.
- `IntegerString[n, "Roman"]` stays unevaluated.
- `ToCharacterCode[s, "UTF8"]` stays unevaluated.
- The `?Minimize` docstring is stale.
