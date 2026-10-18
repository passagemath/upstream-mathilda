#!/usr/bin/env wolframscript
(* ::Package:: *)
(*
  mma_integrate.wls -- the "algorithmically solvable" signal for triage.

  Usage:
    wolframscript -file tools/mathar/mma_integrate.wls <verified.json> <out.json> [nk] [timeout]

  For every VERIFIED entry, run Mathematica's own Integrate on the symbolic
  integral (with the parameters assumed positive / integer as in the corpus)
  under TimeConstrained.  Status:
    CLOSED   result has no Integrate head (closed form, possibly special functions)
    OPEN     Integrate returned unevaluated
    TIMEOUT  exceeded the budget
    ABORTED  the kernel aborted the evaluation
  This never feeds the corpus' reference values (those come from verify.wls);
  it only ranks Mathilda's failures: an integral Mathematica closes is
  evidence a general algorithm reaches it.
*)
args = $ScriptCommandLine;
rs = Select[Import[args[[2]], "RawJSON"], MemberQ[{"VERIFIED", "VERIFIED1"}, #["status"]] &];
nk = If[Length[args] >= 4, ToExpression[args[[4]]], 6];
tmo = If[Length[args] >= 5, ToExpression[args[[5]]], 60];

one[r_] := Module[{f, var, lim, ps, assum, res, t},
  f = ToExpression[r["integrand"]];
  var = ToExpression[r["var"]];
  ps = Symbol /@ r["params"];
  assum = And @@ (# > 0 & /@ ps);
  lim = If[r["kind"] === "definite", Prepend[ToExpression[r["limits"]], var], var];
  {t, res} = AbsoluteTiming[TimeConstrained[
      Quiet@Integrate[f, lim, Assumptions -> assum], tmo, $TimedOut]];
  <|"id" -> r["id"],
    "status" -> Which[res === $TimedOut, "TIMEOUT", ! FreeQ[res, Integrate], "OPEN", True, "CLOSED"],
    "secs" -> Round[t, 0.1], "result" -> StringTake[ToString[res, InputForm], UpTo[300]]|>];

(* Crash-tolerant: each worker appends one JSON line per entry to its own
   file, so a kernel crash or a global abort loses only the entry in flight. *)
safe[v_] := If[StringQ[v] || (NumberQ[v] && Head[v] =!= Complex), v, ToString[v, InputForm]];
emit[a_Association] := Module[{fn = args[[3]] <> ".part" <> ToString[$KernelID], s},
  s = OpenAppend[fn];
  WriteString[s, ExportString[Map[safe, a], "RawJSON", "Compact" -> True] <> "\n"];
  Close[s]];
LaunchKernels[nk];
DistributeDefinitions[one, tmo, safe, emit, args];
ParallelDo[
  emit[Quiet@CheckAbort[MemoryConstrained[one[r], 4*^9,
      <|"id" -> r["id"], "status" -> "MEMORY", "secs" -> 0, "result" -> ""|>],
    <|"id" -> r["id"], "status" -> "ABORTED", "secs" -> 0, "result" -> ""|>]],
  {r, rs}, Method -> "FinestGrained"];
out = ImportString[#, "RawJSON"] & /@ Select[Flatten[StringSplit[Import[#, "Text"], "\n"] & /@ FileNames[args[[3]] <> ".part*"]], StringLength[#] > 2 &];
Print[Counts[#["status"] & /@ out]];
Export[args[[3]], out, "RawJSON"];
