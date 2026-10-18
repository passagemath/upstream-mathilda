#!/usr/bin/env wolframscript
(* instance.wls -- add the numeric-instance integrand/limits (first verified
   instantiation substituted) to verify.wls output, so build_corpus.py never
   has to substitute textually.
   Usage: wolframscript -file tools/mathar/instance.wls <verified.json> <out.json> *)
args = $ScriptCommandLine;
rs = Import[args[[2]], "RawJSON"];
add[r_] := Module[{f, lim, inst},
  If[! MemberQ[{"VERIFIED", "VERIFIED1"}, r["status"]] || r["good"] === {}, Return[r]];
  inst = ToExpression[r["good"][[1]]["inst"]];
  f = ToExpression[r["integrand"], InputForm, Hold];
  lim = If[r["limits"] === "", Hold[{}], ToExpression[r["limits"], InputForm, Hold]];
  Join[r, <|"inst_integrand" -> ToString[ReleaseHold[f] /. inst, InputForm],
    "inst_limits" -> ToString[ReleaseHold[lim] /. inst, InputForm]|>]];
Export[args[[3]], add /@ rs, "RawJSON"];
