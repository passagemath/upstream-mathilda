(* ====================================================================== *)
(* gamma.m -- FullSimplify identities for the gamma-function family.       *)
(*                                                                          *)
(* Loaded lazily by FullSimplify when one of the heads Gamma, LogGamma,     *)
(* PolyGamma, Pochhammer, Beta, or Factorial appears in the input. Each      *)
(* identity is a unary function applied (by Simplify's TransformationFunc-   *)
(* tions machinery) to every subexpression; it uses top-level Replace since  *)
(* Simplify already walks every node. Identities that do not apply return    *)
(* their input unchanged.                                                    *)
(* ====================================================================== *)

(* Gamma[z+1] == z Gamma[z]  (recurrence; the engine reuses it both to       *)
(* collapse Gamma[x+1]/Gamma[x] -> x and to lower integer-shifted gammas).   *)
(* Recurrence Gamma[z+1] == z Gamma[z]: collapses Gamma[x+1]/Gamma[x] -> x and *)
(* shifts integer-offset gammas.                                              *)
(*                                                                            *)
(* The reflection / conjugate-pair identities below fire ONLY as PAIR rules   *)
(* (both gammas must be present in the same Times), so they never raise the   *)
(* gamma count, and each is guarded to avoid the spurious ComplexInfinity at  *)
(* integer z that made a naive lone-gamma reflection unsafe during the search.*)
RegisterTransforms[Gamma, {
    Function[e, Replace[e, Gamma[z_ + 1] :> z Gamma[z]]],
    (* Raising / consolidation: z Gamma[z] == Gamma[z+1]. Absorbs a factor   *)
    (* that matches the gamma argument (e.g. x Gamma[x] -> Gamma[1+x]). The   *)
    (* Orderless Times match lets the factor and the Gamma appear in any      *)
    (* position, and the leading r___ carries any unrelated cofactors along.  *)
    Function[e, Replace[e,
        HoldPattern[Times[r___, z_, Gamma[z_]]] :> Gamma[z + 1] Times[r]]],
    (* Reflection PAIR: Gamma[z] Gamma[1-z] == Pi / Sin[Pi z].  Fires only    *)
    (* when the two arguments sum to 1 and z is not an integer (guard against *)
    (* Sin[Pi n] == 0 -> ComplexInfinity).  e.g. Gamma[1/4] Gamma[3/4] ->     *)
    (* Pi Sqrt[2].                                                            *)
    Function[e, Replace[e,
        HoldPattern[Times[r___, Gamma[z_], Gamma[w_]]] /;
            (! IntegerQ[z]) && PossibleZeroQ[z + w - 1] :>
            Pi / Sin[Pi z] Times[r]]],
    (* Shifted reflection PAIR (u + v == 2): with d = (u - v)/2,               *)
    (* Gamma[1+d] Gamma[1-d] == d Gamma[d] Gamma[1-d] == Pi d / Sin[Pi d].     *)
    (* Holds for EVERY complex d (an analytic identity), so it needs no        *)
    (* conjugate-pair or reality guard: Gamma[1+I b] Gamma[1-I b] collapses to *)
    (* Pi b / Sinh[Pi b] for any b.  Integer d is excluded (Sin[Pi d] == 0).   *)
    Function[e, Replace[e,
        HoldPattern[Times[r___, Gamma[u_], Gamma[v_]]] /;
            (! IntegerQ[Together[(u - v)/2]]) && PossibleZeroQ[u + v - 2] :>
            (Pi ((u - v)/2) / Sin[Pi (u - v)/2]) Times[r]]],
    (* Symmetric reflection PAIR (u + v == 1): with d = (u - v)/2,             *)
    (* Gamma[1/2+d] Gamma[1/2-d] == Pi / Cos[Pi d] (the reflection formula     *)
    (* rewritten about 1/2), again valid for every complex d.  It gives the    *)
    (* compact Pi Sech[Pi b] for Gamma[1/2+I b] Gamma[1/2-I b], where the      *)
    (* Pi/Sin[Pi z] form above leaves Csc[Pi (1/2 + I b)].  (These replace the *)
    (* former conjugate-pair rules, whose PossibleZeroQ[u - Conjugate[v]]       *)
    (* guard was decided at real points only and so fired for complex b,       *)
    (* returning Pi Sech[Pi Im[u]] -- wrong off the real line.)                 *)
    Function[e, Replace[e,
        HoldPattern[Times[r___, Gamma[u_], Gamma[v_]]] /;
            (! IntegerQ[u]) && PossibleZeroQ[u + v - 1] :>
            (Pi / Cos[Pi (u - v)/2]) Times[r]]]
}];

(* LogGamma[z+1] == Log[z] + LogGamma[z]. *)
RegisterTransforms[LogGamma, {
    Function[e, Replace[e, LogGamma[z_ + 1] :> Log[z] + LogGamma[z]]]
}];

(* PolyGamma[n, z+1] == PolyGamma[n, z] + (-1)^n n! / z^(n+1). *)
RegisterTransforms[PolyGamma, {
    Function[e, Replace[e,
        PolyGamma[n_, z_ + 1] :> PolyGamma[n, z] + (-1)^n n! / z^(n + 1)]],
    (* PolyGamma[0, z+1] == PolyGamma[0, z] + 1/z. *)
    Function[e, Replace[e,
        PolyGamma[0, z_ + 1] :> PolyGamma[0, z] + 1/z]]
}];

(* Pochhammer[a, n] == Gamma[a+n]/Gamma[a]. *)
RegisterTransforms[Pochhammer, {
    Function[e, Replace[e, Pochhammer[a_, n_] :> Gamma[a + n]/Gamma[a]]]
}];

(* Beta[a, b] == Gamma[a] Gamma[b] / Gamma[a+b]. *)
RegisterTransforms[Beta, {
    Function[e, Replace[e, Beta[a_, b_] :> Gamma[a] Gamma[b] / Gamma[a + b]]]
}];

(* Factorial[n] == Gamma[n+1]; lets gamma identities reach factorials. *)
RegisterTransforms[Factorial, {
    Function[e, Replace[e, Factorial[n_] :> Gamma[n + 1]]]
}];
