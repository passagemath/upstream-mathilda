(* SingularDefiniteIntegralsCorpus.m
 *
 * 50-case definite-integration stress corpus targeting the Newton-Leibniz /
 * Fundamental-Theorem-of-Calculus boundaries: integrable endpoint singularities,
 * singularities/kinks inside the domain (Abs/Sign/CubeRoot/Min/Max), Cauchy
 * principal values (PrincipalValue -> True), branch cuts in the antiderivative
 * (trig-rational over a full period, Poisson kernels), and essential/oscillatory
 * endpoints (reciprocal-argument, Fresnel, Dirichlet, log-rational half-line).
 *
 * Each record is  {Hold[rawIntegral], Hold[expected], {parameter -> value, ...}}.
 * The integral is Held so Get[] does not evaluate all 50 at load (one process);
 * the third element instantiates the free parameters at a generic point (chosen
 * to satisfy the case's Assumptions and avoid special values) so a runner can
 * verify  N[(result - expected) /. instantiation] ~ 0  without relying on
 * Simplify collapsing the symbolic form.
 *
 * TWO BATTERY ORACLES WERE CORRECTED HERE (the user-supplied battery had them
 * wrong; both verified by NIntegrate). DO NOT "re-fix" Mathilda to the old
 * values:
 *   - Case 15  Integrate[Abs[x^2-a^2],{x,-2a,2a}]  = 4 a^3   (battery said 6 a^3)
 *   - Case 41  Integrate[Sin[a/x]/x^2,{x,a,Inf}]   = (1-Cos[1])/a  (battery had
 *              the sign flipped: (Cos[1]-1)/a)
 *)

{
  (* ---- 1-10: integrable singularities at the limits of integration ---- *)
  (*  1 *) {Hold[Integrate[x^(a - 1), {x, 0, 1}, Assumptions -> a > 0]], Hold[1/a], {a -> 17/10}},
  (*  2 *) {Hold[Integrate[(1 - x)^(a - 1), {x, 0, 1}, Assumptions -> a > 0]], Hold[1/a], {a -> 17/10}},
  (*  3 *) {Hold[Integrate[x^(a - 1) Log[x], {x, 0, 1}, Assumptions -> a > 0]], Hold[-1/a^2], {a -> 17/10}},
  (*  4 *) {Hold[Integrate[1/Sqrt[a^2 - x^2], {x, 0, a}, Assumptions -> a > 0]], Hold[Pi/2], {a -> 17/10}},
  (*  5 *) {Hold[Integrate[1/Sqrt[x (a - x)], {x, 0, a}, Assumptions -> a > 0]], Hold[Pi], {a -> 17/10}},
  (*  6 *) {Hold[Integrate[Log[x]/Sqrt[x], {x, 0, 1}]], Hold[-4], {}},
  (*  7 *) {Hold[Integrate[Log[x]/(1 - x), {x, 0, 1}]], Hold[-Pi^2/6], {}},
  (*  8 *) {Hold[Integrate[Log[x] Log[1 - x], {x, 0, 1}]], Hold[2 - Pi^2/6], {}},
  (*  9 *) {Hold[Integrate[Log[Sin[x]], {x, 0, Pi}]], Hold[-Pi Log[2]], {}},
  (* 10 *) {Hold[Integrate[Log[Cos[x]], {x, 0, Pi/2}]], Hold[-Pi/2 Log[2]], {}},

  (* ---- 11-20: singularities and kinks WITHIN the domain ---- *)
  (* 11 *) {Hold[Integrate[1/Sqrt[Abs[x]], {x, -a, a}, Assumptions -> a > 0]], Hold[4 Sqrt[a]], {a -> 17/10}},
  (* 12 *) {Hold[Integrate[1/Abs[x]^b, {x, -a, a}, Assumptions -> a > 0 && 0 < b < 1]], Hold[2 a^(1 - b)/(1 - b)], {a -> 17/10, b -> 2/5}},
  (* 13 *) {Hold[Integrate[Log[Abs[x - a]], {x, 0, 2 a}, Assumptions -> a > 0]], Hold[2 a (Log[a] - 1)], {a -> 17/10}},
  (* 14 *) {Hold[Integrate[Sign[x] Abs[x]^a, {x, -1, 1}, Assumptions -> a > -1]], Hold[0], {a -> 3/5}},
  (* 15 *) {Hold[Integrate[Abs[x^2 - a^2], {x, -2 a, 2 a}, Assumptions -> a > 0]], Hold[4 a^3], {a -> 17/10}},
  (* 16 *) {Hold[Integrate[1/CubeRoot[x], {x, -a, a}, Assumptions -> a > 0]], Hold[0], {a -> 17/10}},
  (* 17 *) {Hold[Integrate[Abs[Sin[x]], {x, -Pi, Pi}]], Hold[4], {}},
  (* 18 *) {Hold[Integrate[Sqrt[1 - Cos[2 x]], {x, 0, 2 Pi}]], Hold[4 Sqrt[2]], {}},
  (* 19 *) {Hold[Integrate[Max[x^2, a^2], {x, -2 a, 2 a}, Assumptions -> a > 0]], Hold[20 a^3/3], {a -> 17/10}},
  (* 20 *) {Hold[Integrate[Min[x^2, a^2], {x, -2 a, 2 a}, Assumptions -> a > 0]], Hold[8 a^3/3], {a -> 17/10}},

  (* ---- 21-30: Cauchy principal values ---- *)
  (* 21 *) {Hold[Integrate[1/x, {x, -a, a}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[0], {a -> 17/10}},
  (* 22 *) {Hold[Integrate[1/x, {x, -a, b}, PrincipalValue -> True, Assumptions -> a > 0 && b > 0]], Hold[Log[b/a]], {a -> 17/10, b -> 23/10}},
  (* 23 *) {Hold[Integrate[1/(x - a), {x, 0, 2 a}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[0], {a -> 17/10}},
  (* 24 *) {Hold[Integrate[Tan[x], {x, 0, Pi}, PrincipalValue -> True]], Hold[0], {}},
  (* 25 *) {Hold[Integrate[Cot[x], {x, -Pi/2, Pi/2}, PrincipalValue -> True]], Hold[0], {}},
  (* 26 *) {Hold[Integrate[x/(x^2 - a^2), {x, 0, 2 a}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[Log[3]/2], {a -> 17/10}},
  (* 27 *) {Hold[Integrate[Sin[x]/x^2, {x, -a, a}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[0], {a -> 17/10}},
  (* 28 *) {Hold[Integrate[Cos[x]/x, {x, -a, a}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[0], {a -> 17/10}},
  (* 29 *) {Hold[Integrate[Exp[x]/x, {x, -a, a}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[2 SinhIntegral[a]], {a -> 17/10}},
  (* 30 *) {Hold[Integrate[Cos[x]/(x^2 - a^2), {x, -Infinity, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]], Hold[-Pi Sin[a]/a], {a -> 17/10}},

  (* ---- 31-40: branch cuts in the antiderivative & trig phase jumps ---- *)
  (* 31 *) {Hold[Integrate[1/(a + b Cos[x]), {x, 0, 2 Pi}, Assumptions -> a > b > 0]], Hold[2 Pi/Sqrt[a^2 - b^2]], {a -> 3, b -> 2}},
  (* 32 *) {Hold[Integrate[1/(a - b Cos[x]), {x, 0, 2 Pi}, Assumptions -> a > b > 0]], Hold[2 Pi/Sqrt[a^2 - b^2]], {a -> 3, b -> 2}},
  (* 33 *) {Hold[Integrate[1/(1 - 2 a Cos[x] + a^2), {x, 0, 2 Pi}, Assumptions -> 0 < a < 1]], Hold[2 Pi/(1 - a^2)], {a -> 2/5}},
  (* 34 *) {Hold[Integrate[1/(1 - 2 a Cos[x] + a^2), {x, 0, 2 Pi}, Assumptions -> a > 1]], Hold[2 Pi/(a^2 - 1)], {a -> 5/2}},
  (* 35 *) {Hold[Integrate[x Sin[x]/(1 + Cos[x]^2), {x, 0, Pi}]], Hold[Pi^2/4], {}},
  (* 36 *) {Hold[Integrate[Log[1 - 2 a Cos[x] + a^2], {x, 0, Pi}, Assumptions -> 0 < a < 1]], Hold[0], {a -> 2/5}},
  (* 37 *) {Hold[Integrate[Log[1 - 2 a Cos[x] + a^2], {x, 0, Pi}, Assumptions -> a > 1]], Hold[2 Pi Log[a]], {a -> 5/2}},
  (* 38 *) {Hold[Integrate[1/(1 + a^2 Tan[x]^2), {x, 0, Pi/2}, Assumptions -> a > 0]], Hold[Pi/(2 (1 + a))], {a -> 17/10}},
  (* 39 *) {Hold[Integrate[Cos[x]^2/(1 + a Sin[x]^2), {x, 0, Pi}, Assumptions -> a > 0]], Hold[Pi (Sqrt[1 + a] - 1)/(a Sqrt[1 + a])], {a -> 17/10}},
  (* 40 *) {Hold[Integrate[x/Sqrt[a^2 - x^2], {x, 0, a}, Assumptions -> a > 0]], Hold[a], {a -> 17/10}},

  (* ---- 41-50: essential singularities, oscillatory endpoints, conditionally convergent ---- *)
  (* 41 *) {Hold[Integrate[Sin[a/x]/x^2, {x, a, Infinity}, Assumptions -> a > 0]], Hold[(1 - Cos[1])/a], {a -> 17/10}},
  (* 42 *) {Hold[Integrate[Cos[a/x]/x^2, {x, a, Infinity}, Assumptions -> a > 0]], Hold[Sin[1]/a], {a -> 17/10}},
  (* 43 *) {Hold[Integrate[Exp[-I a/x]/x^2, {x, a, Infinity}, Assumptions -> a > 0]], Hold[I (Exp[-I] - 1)/a], {a -> 17/10}},
  (* 44 *) {Hold[Integrate[Sin[a x^2], {x, 0, Infinity}, Assumptions -> a > 0]], Hold[1/2 Sqrt[Pi/(2 a)]], {a -> 17/10}},
  (* 45 *) {Hold[Integrate[Cos[a x^2], {x, 0, Infinity}, Assumptions -> a > 0]], Hold[1/2 Sqrt[Pi/(2 a)]], {a -> 17/10}},
  (* 46 *) {Hold[Integrate[Sin[x]/Sqrt[x], {x, 0, Infinity}]], Hold[Sqrt[Pi/2]], {}},
  (* 47 *) {Hold[Integrate[Cos[x]/Sqrt[x], {x, 0, Infinity}]], Hold[Sqrt[Pi/2]], {}},
  (* 48 *) {Hold[Integrate[Log[x]/(a^2 + x^2), {x, 0, Infinity}, Assumptions -> a > 0]], Hold[Pi Log[a]/(2 a)], {a -> 17/10}},
  (* 49 *) {Hold[Integrate[Log[x]^2/(a^2 + x^2), {x, 0, Infinity}, Assumptions -> a > 0]], Hold[Pi (Pi^2 + 4 Log[a]^2)/(8 a)], {a -> 17/10}},
  (* 50 *) {Hold[Integrate[Sin[a x]/x, {x, -Infinity, Infinity}, Assumptions -> a > 0]], Hold[Pi], {a -> 17/10}}
}
