#!/usr/bin/env wolframscript
(* ::Package:: *)
(*
  verify.wls -- numeric oracle for the Mathar integral entries.

  Usage:
    wolframscript -file tools/mathar/verify.wls <entries.json> <verified.json> [nkernels] [overrides|id,id,..]

  For every entry with status "ok" (translated by extract.py) this decides,
  independently of Mathilda, whether the translated integral and the paper's
  value agree.  It is the ONLY source of the reference values stored in the
  corpus, so a translation slip or a typo in the paper can never become an
  expected answer.

  Per entry:
   1. Expand notational ambiguities into variants and keep the first one that
      verifies: elliptic K(k)/E(k) as modulus (EllipticK[k^2]) or parameter;
      incomplete F(.,.)/E(.,.) amplitude-first or modulus-first; K'(k).
   2. Isolate the integral: the side holding it may be c*Integrate[..] + d.
   3. Instantiate the free parameters at admissible values (integer-named
      parameters at small positive integers first, then reals), filtered by the
      translated side conditions when there are any.
   4. Definite: compare the value with NIntegrate (no convergence messages
      allowed) at two distinct instantiations; relative error < 1e-7.
      Indefinite: D[F, x] - f must vanish at three points per instantiation.
   5. Record the reference value N[value /. inst, 30] for each instantiation.

  Output statuses: VERIFIED, MISMATCH (integral converged but disagrees with
  every variant), NOCONV (no instantiation gave a trustworthy numeric
  integral), BADEXPR (does not parse/evaluate), NOVALUE (value not numeric).
*)

args = $ScriptCommandLine;
If[Length[args] < 3, Print["usage: verify.wls in.json out.json [nk]"]; Exit[1]];
inFile = args[[2]]; outFile = args[[3]];
nk = If[Length[args] >= 4, ToExpression[args[[4]]], 8];

entries = Select[Import[inFile, "RawJSON"], #["status"] === "ok" &];
(* optional 5th argument: "overrides" (only hand-repaired entries) or a
   comma-separated id list *)
If[Length[args] >= 5,
  entries = If[args[[5]] === "overrides",
    Select[entries, TrueQ[Lookup[#, "override", False]] &],
    Select[entries, MemberQ[StringSplit[args[[5]], ","], #["id"]] &]]];
Print["entries to verify: ", Length[entries]];

intNames = {"n", "m", "k", "l", "j", "N", "nu0"};
realPool = {3/8, 5/8, 5/6, 5/4, 7/4, 7/3, 1/5, 3};
intPool = {1, 2, 3, 4, 0};

(* notational variants *)
variantRules = {
  (* modulus convention, amplitude first *)
  {MatharK[k_] :> EllipticK[k^2], MatharE[k_] :> EllipticE[k^2],
   MatharK1[k_] :> 2/(1 + k) EllipticK[((1 - k)/(1 + k))^2],
   MatharE1[k_] :> (1 + k) EllipticE[((1 - k)/(1 + k))^2] - 2 k/(1 + k) EllipticK[((1 - k)/(1 + k))^2],
   MatharF2[a_, b_] :> EllipticF[a, b^2], MatharE2[a_, b_] :> EllipticE[a, b^2],
   MatharPi[n_, k_] :> EllipticPi[n, k^2]},
  (* modulus convention, modulus first *)
  {MatharK[k_] :> EllipticK[k^2], MatharE[k_] :> EllipticE[k^2],
   MatharK1[k_] :> 2/(1 + k) EllipticK[((1 - k)/(1 + k))^2],
   MatharE1[k_] :> (1 + k) EllipticE[((1 - k)/(1 + k))^2] - 2 k/(1 + k) EllipticK[((1 - k)/(1 + k))^2],
   MatharF2[a_, b_] :> EllipticF[b, a^2], MatharE2[a_, b_] :> EllipticE[b, a^2],
   MatharPi[n_, k_] :> EllipticPi[n^2, k^2]},
  (* parameter convention *)
  {MatharK[k_] :> EllipticK[k], MatharE[k_] :> EllipticE[k],
   MatharK1[k_] :> EllipticK[1 - k], MatharE1[k_] :> EllipticE[1 - k],
   MatharF2[a_, b_] :> EllipticF[a, b], MatharE2[a_, b_] :> EllipticE[a, b],
   MatharPi[n_, k_] :> EllipticPi[n, k]},
  {MatharK[k_] :> EllipticK[k], MatharE[k_] :> EllipticE[k],
   MatharK1[k_] :> EllipticK[1 - k], MatharE1[k_] :> EllipticE[1 - k],
   MatharF2[a_, b_] :> EllipticF[b, a], MatharE2[a_, b_] :> EllipticE[b, a],
   MatharPi[n_, k_] :> EllipticPi[n, k]}
};
hasMathar[e_] := ! FreeQ[e, MatharK | MatharE | MatharK1 | MatharE1 | MatharF2 | MatharE2 | MatharPi];

(* Free Global` symbols, excluding the dummy indices of Sum/Product/Table
   (an index name that is also used freely elsewhere is rare enough to accept). *)
symNames[e_] := Complement[
  Union[Cases[e, s_Symbol /; Context[s] === "Global`" :> SymbolName[s], {0, Infinity}, Heads -> True]],
  Union[Cases[e, (Sum | Product | Table)[_, its__List] :>
      Sequence @@ Cases[{its}, {v_Symbol, ___} :> SymbolName[v]], {0, Infinity}, Heads -> True]]];

relErr[a_, b_] := Abs[a - b]/Max[1, Abs[b]];

(* Trustworthy numeric integral: no messages at all. *)
nint[f_, {x_, lo_, hi_}] := Module[{r},
  r = Check[TimeConstrained[
      NIntegrate[f, {x, lo, hi}, PrecisionGoal -> 10, AccuracyGoal -> 12,
        MaxRecursion -> 40, WorkingPrecision -> MachinePrecision], 25, $Failed],
    $Failed];
  If[NumericQ[r], r, $Failed]];

cmplxOK[v_] := NumericQ[v] && Abs[v] < 10^30;

(* Parameters the side conditions declare integer (Element[p, Integers]) are
   always drawn from the integer pool; conventional integer names (n, m, k, ...)
   are tried as integers first and as reals after. *)
intDeclared[conds_] := Union[Cases[conds, Element[s_Symbol, Integers] :> SymbolName[s], {0, Infinity}]];

instances[params_, conds_, seed_] := Module[{tries = {}, a, i, decl = intDeclared[conds]},
  SeedRandom[seed];
  Do[
    a = Table[
      p -> Which[
        MemberQ[decl, p], RandomChoice[Join[intPool[[1 ;; 4]], {5, 6}]],
        MemberQ[intNames, p] && i <= 24, RandomChoice[intPool[[1 ;; 4]]],
        True, RandomChoice[realPool]], {p, params}];
    AppendTo[tries, a], {i, 64}];
  DeleteDuplicates[tries]];

condsHold[conds_, inst_] := AllTrue[conds, TrueQ[(# /. inst)] &];

verifyEntry[e_] := Module[
  {lhs, rhsList, id = e["id"], res, conds, condsOK},
  res = Catch[
    lhs = Quiet@Check[ToExpression[e["lhs"], InputForm, Hold], Throw[<|"status" -> "BADEXPR", "why" -> "lhs"|>]];
    rhsList = Quiet@DeleteCases[Check[ToExpression[#, InputForm, Hold], $Failed] & /@ e["rhs"], $Failed];
    If[rhsList === {}, Throw[<|"status" -> "BADEXPR", "why" -> "rhs"|>]];
    conds = Quiet@DeleteCases[Check[ToExpression[#, InputForm], $Failed] & /@ e["conds"], $Failed];
    tryVariants[e, lhs, rhsList, conds]
  ];
  Join[<|"id" -> id|>, res]];

(* isolate Integrate in held lhs: lhs == rhs  -> integral value *)
(* Returns {{integrand, limits}, value}; the Integrate itself is never evaluated. *)
isolate[lhsH_, rhsH_] := Module[{ints, ii, eqn, sol, l, r},
  ints = Cases[lhsH, HoldPattern[Integrate[a_, b_]] :> {a, b}, {0, Infinity}];
  If[Length[ints] =!= 1, Return[$Failed]];
  l = ReleaseHold[lhsH /. HoldPattern[Integrate[_, _]] -> ii];
  r = ReleaseHold[rhsH];
  If[FreeQ[l, ii], Return[$Failed]];
  If[l === ii, Return[{ints[[1]], r}]];
  sol = Quiet@Solve[l == r, ii];
  If[! MatchQ[sol, {{ii -> _}}], Return[$Failed]];
  {ints[[1]], ii /. sol[[1]]}];

tryVariants[e_, lhs_, rhsList_, conds_] := Module[
  {vrs, best = <|"status" -> "NOCONV"|>, out},
  vrs = If[hasMathar[{lhs, rhsList}], variantRules, {{}}];
  Do[
    Do[
      out = tryOne[e, lhs /. vr, rhs /. vr, conds, vr];
      If[out["status"] === "VERIFIED", Return[out, Module]];
      If[out["status"] === "MISMATCH" && best["status"] =!= "MISMATCH", best = out];
      If[best["status"] === "NOCONV" && out["status"] =!= "NOCONV", best = out],
      {rhs, rhsList}],
    {vr, vrs}];
  best];

tryOne[e_, lhsH_, rhsH_, conds_, vr_] := Module[
  {iso, intg, val, f, var, lim, params, insts, good = {}, bad = {}, numv, refv, inst,
   fi, vi, kind, xs, okInst, extra},
  iso = isolate[lhsH, rhsH];
  If[iso === $Failed, Return[<|"status" -> "BADEXPR", "why" -> "isolate"|>]];
  {intg, val} = iso;
  f = intg[[1]]; lim = intg[[2]];
  kind = If[ListQ[lim], "definite", "indefinite"];
  var = If[ListQ[lim], lim[[1]], lim];
  params = DeleteCases[symNames[{f, If[ListQ[lim], Rest[lim], {}]}], SymbolName[var]];
  extra = Complement[symNames[val], params, {SymbolName[var]}];
  If[kind === "indefinite", extra = Complement[extra, {SymbolName[var]}]];
  If[extra =!= {}, Return[<|"status" -> "BADEXPR", "why" -> "value has extra symbols " <> ToString[extra]|>]];
  insts = instances[params, conds, Hash[e["id"]]];
  insts = (Thread[Symbol /@ params -> (# /. Rule[_, v_] :> v)] &) /@ insts;
  If[params === {}, insts = {{}}];
  Do[
    If[Length[good] >= 2, Break[]];
    If[conds =!= {} && ! condsHold[conds, inst] && Length[insts] > 1, Continue[]];
    If[kind === "definite",
      fi = f /. inst; vi = Quiet[N[val /. inst, 20]];
      If[! cmplxOK[vi], AppendTo[bad, {inst, "novalue"}]; Continue[]];
      numv = Quiet@nint[fi, {var, lim[[2]] /. inst, lim[[3]] /. inst}];
      If[numv === $Failed, Continue[]];
      (* A value near zero makes the relative test vacuous (relErr divides by
         Max[1, |b|]); such an instantiation proves nothing -- skip it. *)
      If[Abs[vi] < 10^-6 || Abs[numv] < 10^-6, Continue[]];
      If[relErr[numv, vi] < 10^-7,
        refv = Quiet[N[val /. inst, 30]];
        AppendTo[good, <|"inst" -> ToString[inst, InputForm], "ref" -> ToString[refv, InputForm],
          "num" -> ToString[numv, InputForm]|>],
        AppendTo[bad, {inst, numv, vi}]],
      (* indefinite: derivative test at three points *)
      xs = {0.31, 0.77, 1.43};
      okInst = Quiet@Check[
          Max[Table[
              Abs[N[(D[val /. Abs[u_] :> Sqrt[u^2], var] - f) /. inst /. var -> xx, 20]]/Max[1, Abs[N[f /. inst /. var -> xx, 20]]],
              {xx, xs}]], $Failed];
      If[! NumericQ[okInst], Continue[]];
      If[okInst < 10^-8,
        AppendTo[good, <|"inst" -> ToString[inst, InputForm], "ref" -> "", "num" -> ""|>],
        AppendTo[bad, {inst, okInst}]]
    ],
    {inst, If[params === {}, {{}}, insts]}];
  Which[
    Length[good] >= Min[2, Max[1, Length[insts]]],
      <|"status" -> "VERIFIED", "kind" -> kind, "var" -> ToString[var, InputForm],
        "integrand" -> ToString[f, InputForm],
        "limits" -> If[kind === "definite", ToString[Rest[lim], InputForm], ""],
        "value" -> ToString[val, InputForm], "params" -> params,
        "variant" -> ToString[vr, InputForm], "good" -> good|>,
    bad =!= {} && good === {},
      <|"status" -> "MISMATCH", "kind" -> kind, "integrand" -> ToString[f, InputForm],
        "value" -> ToString[val, InputForm], "bad" -> ToString[Take[bad, UpTo[2]], InputForm]|>,
    Length[good] == 1,
      <|"status" -> "VERIFIED1", "kind" -> kind, "var" -> ToString[var, InputForm],
        "integrand" -> ToString[f, InputForm],
        "limits" -> If[kind === "definite", ToString[Rest[lim], InputForm], ""],
        "value" -> ToString[val, InputForm], "params" -> params,
        "variant" -> ToString[vr, InputForm], "good" -> good|>,
    True, <|"status" -> "NOCONV", "kind" -> kind|>]];

LaunchKernels[nk];
DistributeDefinitions[verifyEntry, tryVariants, tryOne, isolate, instances, intDeclared, condsHold,
  nint, cmplxOK, relErr, symNames, hasMathar, variantRules, intNames, realPool, intPool];
t0 = AbsoluteTime[];
results = ParallelMap[
  Function[en, Module[{r},
    r = TimeConstrained[verifyEntry[en], 300, <|"id" -> en["id"], "status" -> "TIMEOUT"|>];
    r]], entries, Method -> "FinestGrained"];
Print["elapsed s: ", Round[AbsoluteTime[] - t0]];
Print[Counts[#["status"] & /@ results]];
Export[outFile, results, "RawJSON"];
