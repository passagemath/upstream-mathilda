(* MseTrickIntegralsCorpus.m
 *
 * Definite integrals from the MSE thread "Tough integrals that can be easily
 * beaten by using simple techniques" (math.stackexchange.com/q/1821075), framed
 * as a campaign for fully-general, algorithmic definite-integration methods.
 * The thread's "simple techniques" reduce to three general method families:
 *   A. General King / reflection averaging   (x -> a+b-x ; f(x)+f(a+b-x) simpler)
 *   B. Series expand -> term-by-term -> resum (R(x) Log[x]^k on [0,1])
 *   C. Feynman on a parameter already present (differentiate under the integral)
 *
 * Each record is  {Hold[rawIntegral], Hold[expected], {parameter -> value, ...}}
 * (an empty instantiation list for a parameter-free case).  The integral is
 * Held so Get[] does not evaluate all cases at load; the runner verifies
 * N[(result - expected) /. instantiation] ~ 0, so an unreduced-but-correct
 * symbolic answer still scores CORRECT.  All closed forms are NIntegrate-
 * verified to >= 6 digits.
 *
 * INVARIANT: zero WRONG answers.  Open cases (a method not yet landed) are
 * UNEVALUATED, never WRONG.  The runner gates WRONG == 0 (hard) and
 * CORRECT >= a baseline that ratchets UP as each method lands.
 *)

{
  (* ---- B: series term-by-term on [0,1], alternating (1 + x^d) family ---- *)
  (*  1 *) {Hold[Integrate[Log[x]/(1 + x^2), {x, 0, 1}]], Hold[-Catalan], {}},
  (*  2 *) {Hold[Integrate[Log[x]^2/(1 + x^2), {x, 0, 1}]], Hold[Pi^3/16], {}},
  (*  3 *) {Hold[Integrate[Log[x]/(1 + x), {x, 0, 1}]], Hold[-Pi^2/12], {}},
  (*  4 *) {Hold[Integrate[Log[x]^2/(1 + x), {x, 0, 1}]], Hold[(3/2) Zeta[3]], {}},
  (*  5 *) {Hold[Integrate[Log[x]^3/(1 + x), {x, 0, 1}]], Hold[-7 Pi^4/120], {}},
  (*  6 *) {Hold[Integrate[x Log[x]/(1 + x^2), {x, 0, 1}]], Hold[-Pi^2/48], {}},
  (*  7 *) {Hold[Integrate[x^2 Log[x]/(1 + x^2), {x, 0, 1}]], Hold[Catalan - 1], {}},

  (* ---- B: regression -- integrate_ratlog01_try owns these (1 - x^m);
   *         the series stage must NOT perturb them ---- *)
  (*  8 *) {Hold[Integrate[Log[x]/(1 - x), {x, 0, 1}]], Hold[-Pi^2/6], {}},
  (*  9 *) {Hold[Integrate[Log[x]^2/(1 - x), {x, 0, 1}]], Hold[2 Zeta[3]], {}},
  (* 10 *) {Hold[Integrate[Log[x]^3/(1 - x), {x, 0, 1}]], Hold[-Pi^4/15], {}},

  (* ---- A: general King / reflection averaging ---- *)
  (* 11 *) {Hold[Integrate[1/(1 + Tan[x]^Sqrt[2]), {x, 0, Pi/2}]], Hold[Pi/4], {}},
  (* 12 *) {Hold[Integrate[x Sin[x]/(1 + Cos[x]^2), {x, 0, Pi}]], Hold[Pi^2/4], {}},

  (* ---- C: Feynman on a parameter already present ---- *)
  (* 13 *) {Hold[Integrate[(x^a - 1)/Log[x], {x, 0, 1}, Assumptions -> a > -1]], Hold[Log[a + 1]], {a -> 2}},

  (* ---- Negative controls: must give the right value via OTHER stages and NOT
   *      be misfired by the series stage ---- *)
  (* 14 *) {Hold[Integrate[Log[1 - x]/x, {x, 0, 1}]], Hold[-Pi^2/6], {}},
  (* 15 *) {Hold[Integrate[1/(1 + x^2), {x, 0, 1}]], Hold[Pi/4], {}}
}
