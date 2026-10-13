(* ExtremeDefiniteIntegralsCorpus.m
 *
 * Second, harder 50-case definite-integration stress corpus
 * (battery `extremeAssumedIntegralStressTest`): Frullani integrals & advanced
 * rational poles (1-10), fractional branch cuts / Mellin kernels (11-20),
 * essential singularities / oscillatory / hyperbolic kernels (21-30),
 * complex-contour deflections & Cauchy principal values (31-40), and
 * polylogarithms / special removals / Gamma-Zeta bounds (41-50).
 *
 * Each record is  {Hold[rawIntegral], Hold[expected], {parameter -> value, ...}}.
 * The integral is Held so Get[] does not evaluate all 50 at load; the third
 * element instantiates the free parameters at a generic point satisfying the
 * case's Assumptions, so a runner can verify  N[(result - expected) /. inst] ~ 0
 * without relying on Simplify collapsing the symbolic form.
 *
 * FOUR BATTERY ORACLES WERE CORRECTED HERE (all verified by NIntegrate to >=12
 * digits; the user-supplied battery had them wrong).  DO NOT "re-fix" Mathilda
 * to the old values:
 *   - Case  6  (e^{-a x^2}-e^{-b x^2})/x^2   = Sqrt[Pi](Sqrt[b]-Sqrt[a])
 *              (battery had an extra factor of 2)
 *   - Case 24  x Sin[a x]/Cosh[b x]^2        = (Pi/(2 b^2)) Csch[a Pi/(2b)]
 *                                              (a Pi/(2b) Coth[a Pi/(2b)] - 1)
 *              (battery gave the formula for Integrate[Cos[a x]/Cosh[b x]^2])
 *   - Case 36  Sqrt[x]/(x^2-a^2)  (PV)       = +Pi/(2 Sqrt[a])
 *              (battery had the wrong sign, -Pi/(2 Sqrt[a]))
 *   - Case 37  x Log[x]/(x^2-a^2)^2 (PV)     = double-pole Hadamard finite part,
 *              a regularization-convention-dependent value; the battery oracle
 *              -1/(4 a^2) matches no standard regularization.  This case is
 *              EXCLUDED from scoring in the runner (never WRONG, never CORRECT).
 *   - Case 19  x^a/(1-x)^2 (PV)              = ALSO a double-pole "PV": the Cauchy
 *              PV (symmetric deletion) diverges as 2A/eps, and the Hadamard finite
 *              part is -Pi a Cot[Pi a] (verified by NIntegrate with the 2A/eps
 *              subtraction), the OPPOSITE sign of the battery's +Pi a Cot[Pi a].
 *              Convention-dependent -> EXCLUDED from scoring, like case 37.
 *
 * TWO INSTANTIATIONS were repointed off degenerate values where the expected
 * form collapses to 0 (which would let an eval-to-0 score a false CORRECT):
 *   - Case 12  a->1, b->3   (a/b=1/2 gives Tan[Pi/2]=Infinity -> expected 0)
 *   - Case 19  a->1/3       (a=1/2   gives Cot[Pi/2]=0        -> expected 0)
 *)

{
  (* ---- 1-10: Frullani integrals and advanced rational poles ---- *)
  (*  1 *) {Hold[Integrate[(Exp[-a x] - Exp[-b x])/x, {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], Hold[Log[b/a]], {a -> 2, b -> 3}},
  (*  2 *) {Hold[Integrate[(Cos[a x] - Cos[b x])/x, {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], Hold[Log[b/a]], {a -> 2, b -> 3}},
  (*  3 *) {Hold[Integrate[(ArcTan[a x] - ArcTan[b x])/x, {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], Hold[(Pi Log[a/b])/2], {a -> 2, b -> 3}},
  (*  4 *) {Hold[Integrate[Sin[a x]^3/x^3, {x, 0, Infinity}, Assumptions -> a > 0]], Hold[(3 a^2 Pi)/8], {a -> 2}},
  (*  5 *) {Hold[Integrate[Sin[a x]^4/x^4, {x, 0, Infinity}, Assumptions -> a > 0]], Hold[(a^3 Pi)/3], {a -> 2}},
  (*  6 *) {Hold[Integrate[(Exp[-a x^2] - Exp[-b x^2])/x^2, {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], Hold[Sqrt[Pi] (Sqrt[b] - Sqrt[a])], {a -> 2, b -> 3}},
  (*  7 *) {Hold[Integrate[1/((x^2 + a^2) (x^2 + b^2)), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], Hold[Pi/(2 a b (a + b))], {a -> 2, b -> 3}},
  (*  8 *) {Hold[Integrate[x^2/((x^2 + a^2) (x^2 + b^2)), {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]], Hold[Pi/(2 (a + b))], {a -> 2, b -> 3}},
  (*  9 *) {Hold[Integrate[1/((x^2 + a^2)^3), {x, 0, Infinity}, Assumptions -> a > 0]], Hold[(3 Pi)/(16 a^5)], {a -> 2}},
  (* 10 *) {Hold[Integrate[x^4/((x^2 + a^2)^3), {x, 0, Infinity}, Assumptions -> a > 0]], Hold[(3 Pi)/(16 a)], {a -> 2}},

  (* ---- 11-20: fractional branch cuts and Mellin-type kernels ---- *)
  (* 11 *) {Hold[Integrate[x^(a - 1)/(1 + x^b), {x, 0, Infinity}, Assumptions -> b > a > 0]], Hold[Pi/(b Sin[(a Pi)/b])], {a -> 3/2, b -> 3}},
  (* 12 *) {Hold[Integrate[x^(a - 1)/(1 - x^b), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> b > a > 0]], Hold[Pi/(b Tan[(a Pi)/b])], {a -> 1, b -> 3}},
  (* 13 *) {Hold[Integrate[1/(1 + x^a), {x, 0, Infinity}, Assumptions -> a > 1]], Hold[Pi/(a Sin[Pi/a])], {a -> 3}},
  (* 14 *) {Hold[Integrate[x^2/(1 + x^a), {x, 0, Infinity}, Assumptions -> a > 3]], Hold[Pi/(a Sin[(3 Pi)/a])], {a -> 4}},
  (* 15 *) {Hold[Integrate[x^(a - 1) Exp[-b x^c], {x, 0, Infinity}, Assumptions -> a > 0 && b > 0 && c > 0]], Hold[Gamma[a/c]/(c b^(a/c))], {a -> 2, b -> 1, c -> 2}},
  (* 16 *) {Hold[Integrate[x^(a - 1)/(Exp[x] - 1), {x, 0, Infinity}, Assumptions -> a > 1]], Hold[Gamma[a] Zeta[a]], {a -> 2}},
  (* 17 *) {Hold[Integrate[x^(a - 1)/(Exp[x] + 1), {x, 0, Infinity}, Assumptions -> a > 0]], Hold[Gamma[a] Zeta[a] (1 - 2^(1 - a))], {a -> 2}},
  (* 18 *) {Hold[Integrate[x^a/(1 + x)^2, {x, 0, Infinity}, Assumptions -> -1 < a < 1]], Hold[Pi a Csc[Pi a]], {a -> 1/2}},
  (* 19 *) {Hold[Integrate[x^a/(1 - x)^2, {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> -1 < a < 1]], Hold[Pi a Cot[Pi a]], {a -> 1/3}},
  (* 20 *) {Hold[Integrate[1/(x^a (1 + x)), {x, 0, Infinity}, Assumptions -> 0 < a < 1]], Hold[Pi Csc[a Pi]], {a -> 1/2}},

  (* ---- 21-30: essential singularities, phase jumps, highly oscillatory ---- *)
  (* 21 *) {Hold[Integrate[Exp[-I a x^2], {x, -Infinity, Infinity}, Assumptions -> a > 0]], Hold[Sqrt[Pi/a] Exp[-I Pi/4]], {a -> 2}},
  (* 22 *) {Hold[Integrate[Exp[-a x^2] Cos[b x], {x, 0, Infinity}, Assumptions -> a > 0 && Element[b, Reals]]], Hold[(Sqrt[Pi] Exp[-b^2/(4 a)])/(2 Sqrt[a])], {a -> 2, b -> 3}},
  (* 23 *) {Hold[Integrate[Cos[a x]/Cosh[b x], {x, 0, Infinity}, Assumptions -> Element[a, Reals] && b > 0]], Hold[Pi/(2 b Cosh[(a Pi)/(2 b)])], {a -> 2, b -> 3}},
  (* 24 *) {Hold[Integrate[x Sin[a x]/Cosh[b x]^2, {x, 0, Infinity}, Assumptions -> Element[a, Reals] && b > 0]], Hold[(Pi/(2 b^2)) Csch[(a Pi)/(2 b)] ((a Pi)/(2 b) Coth[(a Pi)/(2 b)] - 1)], {a -> 2, b -> 3}},
  (* 25 *) {Hold[Integrate[Sin[a x]/Sinh[b x], {x, 0, Infinity}, Assumptions -> Element[a, Reals] && b > 0]], Hold[(Pi Tanh[(a Pi)/(2 b)])/(2 b)], {a -> 2, b -> 3}},
  (* 26 *) {Hold[Integrate[BesselJ[0, a x] Cos[b x], {x, 0, Infinity}, Assumptions -> a > b > 0]], Hold[1/Sqrt[a^2 - b^2]], {a -> 3, b -> 2}},
  (* 27 *) {Hold[Integrate[BesselJ[0, a x] Cos[b x], {x, 0, Infinity}, Assumptions -> b > a > 0]], Hold[0], {a -> 2, b -> 3}},
  (* 28 *) {Hold[Integrate[BesselJ[0, a x] Sin[b x], {x, 0, Infinity}, Assumptions -> b > a > 0]], Hold[1/Sqrt[b^2 - a^2]], {a -> 2, b -> 3}},
  (* 29 *) {Hold[Integrate[x/Sinh[a x], {x, 0, Infinity}, Assumptions -> a > 0]], Hold[Pi^2/(4 a^2)], {a -> 2}},
  (* 30 *) {Hold[Integrate[x^2/Sinh[a x]^2, {x, 0, Infinity}, Assumptions -> a > 0]], Hold[Pi^2/(6 a^3)], {a -> 2}},

  (* ---- 31-40: complex contour deflections & Cauchy principal values ---- *)
  (* 31 *) {Hold[Integrate[1/(x^3 - a^3), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[(-Pi)/(3 Sqrt[3] a^2)], {a -> 2}},
  (* 32 *) {Hold[Integrate[x/(x^3 - a^3), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[Pi/(3 Sqrt[3] a)], {a -> 2}},
  (* 33 *) {Hold[Integrate[1/(x^4 - a^4), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[(-Pi)/(4 a^3)], {a -> 2}},
  (* 34 *) {Hold[Integrate[x^2/(x^4 - a^4), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[Pi/(4 a)], {a -> 2}},
  (* 35 *) {Hold[Integrate[Log[x]/(x^2 - a^2), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[Pi^2/(4 a)], {a -> 2}},
  (* 36 *) {Hold[Integrate[Sqrt[x]/(x^2 - a^2), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[Pi/(2 Sqrt[a])], {a -> 2}},
  (* 37 *) {Hold[Integrate[x Log[x]/(x^2 - a^2)^2, {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[(-1)/(4 a^2)], {a -> 2}},
  (* 38 *) {Hold[Integrate[Log[Abs[x - a]]/(x^2 + b^2), {x, -Infinity, Infinity}, Assumptions -> Element[a, Reals] && b > 0]], Hold[(Pi Log[a^2 + b^2])/(2 b)], {a -> 2, b -> 3}},
  (* 39 *) {Hold[Integrate[Log[x]/(x^2 - 1), {x, 0, Infinity}, PrincipalValue -> True]], Hold[Pi^2/4], {}},
  (* 40 *) {Hold[Integrate[Log[x]^2/(x^2 - 1), {x, 0, Infinity}, PrincipalValue -> True]], Hold[0], {}},

  (* ---- 41-50: polylogarithms, special removals, Gamma/Zeta bounds ---- *)
  (* 41 *) {Hold[Integrate[Log[x]^3/(1 - x), {x, 0, 1}]], Hold[(-Pi^4)/15], {}},
  (* 42 *) {Hold[Integrate[Log[1 + x]/x, {x, 0, 1}]], Hold[Pi^2/12], {}},
  (* 43 *) {Hold[Integrate[Log[1 - x]/x, {x, 0, 1}]], Hold[(-Pi^2)/6], {}},
  (* 44 *) {Hold[Integrate[Log[x^2 + a^2]/(x^2 + b^2), {x, -Infinity, Infinity}, Assumptions -> a > 0 && b > 0]], Hold[(2 Pi Log[a + b])/b], {a -> 2, b -> 3}},
  (* 45 *) {Hold[Integrate[Log[x]^2/(1 + x^2), {x, 0, Infinity}]], Hold[Pi^3/8], {}},
  (* 46 *) {Hold[Integrate[Log[1 + x^2]/(1 + x^2), {x, 0, Infinity}]], Hold[Pi Log[2]], {}},
  (* 47 *) {Hold[Integrate[Log[x]/Sqrt[1 - x^2], {x, 0, 1}]], Hold[(-Pi Log[2])/2], {}},
  (* 48 *) {Hold[Integrate[x/(1 + Exp[x]), {x, 0, Infinity}]], Hold[Pi^2/12], {}},
  (* 49 *) {Hold[Integrate[x^3/(Exp[x] - 1), {x, 0, Infinity}]], Hold[Pi^4/15], {}},
  (* 50 *) {Hold[Integrate[Log[Gamma[x]], {x, 0, 1}]], Hold[Log[2 Pi]/2], {}}
}
