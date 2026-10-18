(* ====================================================================== *)
(* logexp.m -- FullSimplify identities for logarithms of trig/hyperbolic   *)
(*             combinations that collapse to a single inverse-hyperbolic,   *)
(*             inverse-trigonometric, or linear form.                        *)
(*                                                                          *)
(* Loaded lazily by FullSimplify when the head Log appears in the input     *)
(* (every identity below has Log on the left). Each identity is a unary      *)
(* function applied (by Simplify's TransformationFunctions machinery) to     *)
(* every subexpression; it uses a top-level Replace since Simplify already   *)
(* walks every node. Identities that do not apply return their input         *)
(* unchanged.                                                                *)
(*                                                                          *)
(* THE CLASS. These are the real-log antiderivative (Gudermannian) family:   *)
(* the logarithm of a two-term +/- combination of cofunctions, or of a       *)
(* half-angle tangent/cotangent, collapses to a single inverse function.     *)
(* They have no general algorithmic route (a generic log/exp canonicaliser   *)
(* does not find them), which is exactly why they live in FullSimplify's     *)
(* rule tier rather than in Simplify. All signs were verified numerically.   *)
(*                                                                          *)
(*   Log[Sec[x] +/- Tan[x]]   == +/- ArcTanh[Sin[x]]                         *)
(*   Log[Csc[x] +/- Cot[x]]   == +/- ArcTanh[Cos[x]]                         *)
(*   Log[Cosh[x] +/- Sinh[x]] == +/- x                                       *)
(*   Log[Coth[x] +/- Csch[x]] == +/- ArcCoth[Cosh[x]]                        *)
(*   Log[Tan[x/2]]  == -ArcTanh[Cos[x]]   Log[Cot[x/2]]  == ArcTanh[Cos[x]]  *)
(*   Log[Tanh[x/2]] == -ArcCoth[Cosh[x]]  Log[Coth[x/2]] == ArcCoth[Cosh[x]] *)
(*                                                                          *)
(* TWO SURFACE FORMS PER TERM. A coefficient may reach us in either shape:   *)
(* the Simplify pipeline FACTORS a positive sum (2 Sec + 2 Tan becomes       *)
(* Times[2, Plus[Sec, Tan]]) but leaves a difference DISTRIBUTED             *)
(* (2 Sec - 2 Tan stays Plus[Times[2, Sec], Times[-2, Tan]]). So every       *)
(* two-term identity ships FOUR rules: a factored '+' and '-' with one       *)
(* optional coefficient a_. on the whole group, and a distributed '+' and    *)
(* '-' with two INDEPENDENT optional coefficients a_. (first cofunction) and *)
(* b_. (second) guarded by b === a / b === -a. The pair covers whichever     *)
(* form the pipeline hands us, and the bare case (a == 1) through either. A  *)
(* single optional coefficient on both distributed terms cannot express the  *)
(* '-' case at a non-unit coefficient, because -2 is one integer atom and    *)
(* will not match a literal -1 times c.                                      *)
(*                                                                          *)
(* POSITIVITY GUARD. Log[a z] == Log[a] + Log[z] holds on the same branch    *)
(* only for a > 0, so every rule requires TrueQ[Positive[a]]. For the bare   *)
(* forms the optional default a == 1 is positive, so they fire; a symbolic   *)
(* or negative coefficient is left untouched (sound: no branch is crossed).  *)
(* Log[a] vanishes to 0 on evaluation when a == 1. FreeQ[{a, b}, y] keeps    *)
(* the coefficients independent of the angle.                                *)
(*                                                                          *)
(* DOMAIN GUARD. Each identity is a REAL-LOG identity: it holds only where  *)
(* the log argument is positive (otherwise the left side carries +/- I Pi   *)
(* and the right side does not -- Log[Sec[3] + Tan[3]] is 0.142 + 3.14 I     *)
(* but ArcTanh[Sin[3]] is 0.142). So every rule additionally requires that   *)
(* the assumptions in force (FullSimplify's positional assumption, which     *)
(* the driver binds to $Assumptions) PROVE the condition under which the     *)
(* argument is positive (by interval arithmetic on the angle, below):         *)
(*   Sec +/- Tan          (1 +/- Sin)/Cos    > 0  <=  Cos[y] > 0             *)
(*   Csc +/- Cot, Tan[y/2], Cot[y/2]  (1 -/+ Cos)/Sin > 0  <=  Sin[y] > 0     *)
(*   Cosh +/- Sinh        E^(+/- y)          Log[E^y] == y  <=  y real        *)
(*   Coth +/- Csch, Tanh[y/2], Coth[y/2]  (Cosh +/- 1)/Sinh > 0  <=  y > 0     *)
(* With no assumption the input is left alone (sound); under e.g.            *)
(* -Pi/2 < x < Pi/2 the Sec/Tan rule fires.                                  *)
(*                                                                          *)
(* DIRECTION / TERMINATION. The rules only collapse (Log[...] -> inverse);   *)
(* no right-hand side contains a head that re-matches any left-hand side, so *)
(* there is no rewrite loop. Each collapse strictly lowers SimplifyCount, so *)
(* Simplify's search keeps it. Every statement is ';'-terminated as the file *)
(* loader requires (newlines are not statement separators).                  *)
(* ====================================================================== *)

Begin["FullSimplify`Private`"];

(* ---------------------------------------------------------------------- *)
(* Domain-guard helpers (see DOMAIN GUARD above). These PROVE a sign        *)
(* condition from the inequality facts in $Assumptions by plain interval     *)
(* arithmetic; they never call Reduce/Refine on a trigonometric inequality   *)
(* (those decisions are not reliable enough to gate a branch identity on).   *)
(* Every helper answers True only on proof; anything it cannot read --      *)
(* a non-affine angle, several variables, an unbounded range -- is False,   *)
(* which merely leaves the input unsimplified.                               *)
(* ---------------------------------------------------------------------- *)

(* Flatten $Assumptions (And / List nesting) to a list of facts. *)
LEFacts[a_And]  := Flatten[Map[LEFacts, List @@ a]];
LEFacts[a_List] := Flatten[Map[LEFacts, a]];
LEFacts[a_]     := {a};

(* A fact as a list of {lhs, op, rhs} comparisons (chains are split). *)
LEPairs[Inequality[s__]] := Module[{l = {s}},
    Table[{l[[i]], l[[i + 1]], l[[i + 2]]}, {i, 1, Length[l] - 2, 2}]];
LEPairs[(h : Less | LessEqual | Greater | GreaterEqual)[s__]] := Module[{l = {s}},
    Table[{l[[i]], h, l[[i + 1]]}, {i, 1, Length[l] - 1}]];
LEPairs[_] := {};

LEFlip[Less] = Greater; LEFlip[Greater] = Less;
LEFlip[LessEqual] = GreaterEqual; LEFlip[GreaterEqual] = LessEqual;

LERealNum[c_] := NumericQ[c] && TrueQ[Im[N[c]] == 0];

(* {lo, loOpen, hi, hiOpen} for the symbol v from $Assumptions; an absent   *)
(* bound is None. *)
LEVarBounds[v_Symbol] := Module[{pairs, lo = None, loO = True, hi = None,
                                 hiO = True, l, op, r, c, open},
    pairs = Flatten[Map[LEPairs, LEFacts[$Assumptions]], 1];
    Do[
        l = pairs[[i, 1]]; op = pairs[[i, 2]]; r = pairs[[i, 3]];
        If[l =!= v && r === v && LERealNum[l],
            {l, op, r} = {r, LEFlip[op], l}];
        If[l === v && LERealNum[r],
            c = r; open = (op === Less || op === Greater);
            If[op === Less || op === LessEqual,
                If[hi === None || TrueQ[c < hi] || (TrueQ[c == hi] && open),
                    hi = c; hiO = open],
                If[lo === None || TrueQ[c > lo] || (TrueQ[c == lo] && open),
                    lo = c; loO = open]]],
        {i, Length[pairs]}];
    {lo, loO, hi, hiO}];
LEVarBounds[_] := {None, True, None, True};

(* Bounds of an affine real angle y = c0 + c1 v (one free symbol v, real     *)
(* numeric c0, c1 != 0); $Failed when y is not of that shape.                 *)
LEAngleBounds[y_] := Module[{vs, v, c0, c1, b},
    vs = Union[Cases[{y}, _Symbol?(! NumericQ[#] &), Infinity]];
    If[Length[vs] =!= 1, Return[$Failed]];
    v = First[vs];
    c1 = D[y, v];
    c0 = y /. v -> 0;
    If[! (LERealNum[c1] && LERealNum[c0] && TrueQ[c1 != 0] &&
          PossibleZeroQ[y - c0 - c1 v]), Return[$Failed]];
    b = LEVarBounds[v];
    If[TrueQ[c1 > 0],
        {If[b[[1]] === None, None, c0 + c1 b[[1]]], b[[2]],
         If[b[[3]] === None, None, c0 + c1 b[[3]]], b[[4]]},
        {If[b[[3]] === None, None, c0 + c1 b[[3]]], b[[4]],
         If[b[[1]] === None, None, c0 + c1 b[[1]]], b[[2]]}]];

(* lo/hi strictly inside (a, b), an endpoint touching only when it is open. *)
LEInside[{lo_, loO_, hi_, hiO_}, a_, b_] :=
    lo =!= None && hi =!= None &&
    (TrueQ[lo > a] || (TrueQ[lo == a] && loO)) &&
    (TrueQ[hi < b] || (TrueQ[hi == b] && hiO));
LEInside[_, _, _] := False;

(* Cos[y] > 0 for every y the assumptions allow: the y-range sits inside     *)
(* one window (-Pi/2 + 2 Pi k, Pi/2 + 2 Pi k). *)
CosPositiveQ[y_] := Module[{b = LEAngleBounds[y], k},
    If[b === $Failed || b[[1]] === None, Return[False]];
    k = Floor[N[(b[[1]] + Pi/2)/(2 Pi)]];
    If[! IntegerQ[k], Return[False]];
    LEInside[b, -Pi/2 + 2 Pi k, Pi/2 + 2 Pi k]];

(* Sin[y] > 0  <=>  Cos[y - Pi/2] > 0. *)
SinPositiveQ[y_] := CosPositiveQ[y - Pi/2];

(* y > 0: an affine y with a positive lower bound. *)
AnglePositiveQ[y_] := Module[{b = LEAngleBounds[y]},
    b =!= $Failed && b[[1]] =!= None &&
    (TrueQ[b[[1]] > 0] || (TrueQ[b[[1]] == 0] && b[[2]]))];

(* y real under the assumptions (element decision, not an inequality). *)
AngleRealQ[y_] := TrueQ[Refine[Element[y, Reals]]];

RegisterTransforms[Log, {

    (* ------- Sec/Tan: Log[a (Sec[y] +/- Tan[y])] == Log[a] +/- ArcTanh[Sin[y]] *)
    Function[e, Replace[e,
        Log[a_. (Sec[y_] + Tan[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            CosPositiveQ[y]) :> Log[a] + ArcTanh[Sin[y]]]],
    Function[e, Replace[e,
        Log[a_. (Sec[y_] - Tan[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            CosPositiveQ[y]) :> Log[a] - ArcTanh[Sin[y]]]],
    Function[e, Replace[e,
        Log[a_. Sec[y_] + b_. Tan[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a &&
            CosPositiveQ[y]) :>
            Log[a] + ArcTanh[Sin[y]]]],
    Function[e, Replace[e,
        Log[a_. Sec[y_] + b_. Tan[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a &&
            CosPositiveQ[y]) :>
            Log[a] - ArcTanh[Sin[y]]]],

    (* ------- Csc/Cot: Log[a (Csc[y] +/- Cot[y])] == Log[a] +/- ArcTanh[Cos[y]] *)
    Function[e, Replace[e,
        Log[a_. (Csc[y_] + Cot[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            SinPositiveQ[y]) :> Log[a] + ArcTanh[Cos[y]]]],
    Function[e, Replace[e,
        Log[a_. (Csc[y_] - Cot[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            SinPositiveQ[y]) :> Log[a] - ArcTanh[Cos[y]]]],
    Function[e, Replace[e,
        Log[a_. Csc[y_] + b_. Cot[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a &&
            SinPositiveQ[y]) :>
            Log[a] + ArcTanh[Cos[y]]]],
    Function[e, Replace[e,
        Log[a_. Csc[y_] + b_. Cot[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a &&
            SinPositiveQ[y]) :>
            Log[a] - ArcTanh[Cos[y]]]],

    (* ------- Cosh/Sinh: Log[a (Cosh[y] +/- Sinh[y])] == Log[a] +/- y        *)
    (* (Cosh +/- Sinh == E^(+/- y); base Simplify already collapses the bare   *)
    (* case, these add the coefficient-carrying forms and insure the class.)   *)
    Function[e, Replace[e,
        Log[a_. (Cosh[y_] + Sinh[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            AngleRealQ[y]) :> Log[a] + y]],
    Function[e, Replace[e,
        Log[a_. (Cosh[y_] - Sinh[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            AngleRealQ[y]) :> Log[a] - y]],
    Function[e, Replace[e,
        Log[a_. Cosh[y_] + b_. Sinh[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a &&
            AngleRealQ[y]) :> Log[a] + y]],
    Function[e, Replace[e,
        Log[a_. Cosh[y_] + b_. Sinh[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a &&
            AngleRealQ[y]) :> Log[a] - y]],

    (* ------- Coth/Csch: Log[a (Coth[y] +/- Csch[y])] == Log[a] +/- ArcCoth[Cosh[y]] *)
    Function[e, Replace[e,
        Log[a_. (Coth[y_] + Csch[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            AnglePositiveQ[y]) :> Log[a] + ArcCoth[Cosh[y]]]],
    Function[e, Replace[e,
        Log[a_. (Coth[y_] - Csch[y_])] /;
            (FreeQ[a, y] && TrueQ[Positive[a]] &&
            AnglePositiveQ[y]) :> Log[a] - ArcCoth[Cosh[y]]]],
    Function[e, Replace[e,
        Log[a_. Coth[y_] + b_. Csch[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === a &&
            AnglePositiveQ[y]) :>
            Log[a] + ArcCoth[Cosh[y]]]],
    Function[e, Replace[e,
        Log[a_. Coth[y_] + b_. Csch[y_]] /;
            (FreeQ[{a, b}, y] && TrueQ[Positive[a]] && b === -a &&
            AnglePositiveQ[y]) :>
            Log[a] - ArcCoth[Cosh[y]]]],

    (* ------- Half-angle tangents/cotangents (single term; Log[k] -> 0 at k == 1). *)
    (* Log[k Tan[y/2]]  == Log[k] - ArcTanh[Cos[y]]   (real-log form of Int Csc). *)
    Function[e, Replace[e,
        Log[k_. Tan[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]] &&
            SinPositiveQ[y]) :>
            Log[k] - ArcTanh[Cos[y]]]],
    (* Log[k Cot[y/2]]  == Log[k] + ArcTanh[Cos[y]].  *)
    Function[e, Replace[e,
        Log[k_. Cot[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]] &&
            SinPositiveQ[y]) :>
            Log[k] + ArcTanh[Cos[y]]]],
    (* Log[k Tanh[y/2]] == Log[k] - ArcCoth[Cosh[y]]  (real-log form of Int Csch). *)
    Function[e, Replace[e,
        Log[k_. Tanh[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]] &&
            AnglePositiveQ[y]) :>
            Log[k] - ArcCoth[Cosh[y]]]],
    (* Log[k Coth[y/2]] == Log[k] + ArcCoth[Cosh[y]]. *)
    Function[e, Replace[e,
        Log[k_. Coth[y_/2]] /; (FreeQ[k, y] && TrueQ[Positive[k]] &&
            AnglePositiveQ[y]) :>
            Log[k] + ArcCoth[Cosh[y]]]]
}];

End[];
