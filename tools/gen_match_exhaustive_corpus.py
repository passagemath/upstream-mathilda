#!/usr/bin/env python3
"""gen_match_exhaustive_corpus.py -- build tests/match_exhaustive_corpus.m.

The hand-written tests/match_stress_corpus.m checks a few hundred chosen
cases. This generator covers the pattern language systematically instead:
every construct family is crossed with a shared set of subject expressions,
and the expected value of every case is computed by the reference kernel,
never written from memory.

  * a pattern with no names is asked as  MatchQ[subject, pattern]
  * a pattern with names is asked as     ReplaceList[subject, pattern :> {names}]
    which pins the bindings AND the backtracking order AND completeness,
    not merely yes/no.
  * further families exercise the query heads (Cases / Count / Position /
    DeleteCases / Replace* at level specs), attributes on user heads
    (Orderless / Flat / OneIdentity), definitions and their specificity
    ordering, and the parse shape of every pattern operator
    (ToString[FullForm[Hold[...]]]).

Each case runs in its own context in the reference kernel (and in its own forked
process in the Mathilda runner), so definitions and attributes never leak
between cases. Cases the reference answers with a message, a timeout, or a value that
does not survive a round trip as text (temporary $-symbols, contexts) are
dropped.

The output is checked in, so neither CI nor a contributor needs the reference;
re-run this only to add families:

    python3 tools/gen_match_exhaustive_corpus.py            # regenerate
    python3 tools/gen_match_exhaustive_corpus.py --kernel PATH

The kernel is the reference's command-line script runner (invoked as
`PATH -file script`); pass it with --kernel or set REFERENCE_KERNEL.

The runner is tests/test_match_stress_corpus.c (ctest:
match_exhaustive_corpus_tests).
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "tests", "match_exhaustive_corpus.m")
OBS_OUT = os.path.join(ROOT, "tests", "match_exhaustive_observe.m")
DEFAULT_KERNEL = os.environ.get("REFERENCE_KERNEL")

# --------------------------------------------------------------------------
# Subjects. A shared core crossed with every pattern family; families that
# need a particular shape (g[]-carrying products for #102, sums for the
# Optional defaults, patterns-as-data for Verbatim) add their own.
# --------------------------------------------------------------------------

CORE = [
    "0", "5", "-3", "2/3", "1.5", "3 + 4 I", "x", "\"s\"",
    "f[]", "f[a]", "f[a, b]", "f[a, b, c]", "f[a, a]", "f[a, b, a]",
    "f[1, 2, 3]", "f[1, a, 2]", "f[a, f[b]]", "f[g[a]]", "g[a, b]",
    "f[a, b][c]",
    "{}", "{a}", "{a, b, c}", "{1, 2, 3}", "{1, {2}}", "{a, a, b}",
    "a + b", "a + b + c", "a b", "2 a", "a^2", "a^b", "Sin[x]",
]

G_SUBJ = ["g[]", "g[a]", "2 g[]", "a g[]", "3 g[a]", "g[] h[]",
          "f[1, g[]]", "f[1, 2, g[]]", "f[g[]]", "f[1, 2 g[]]",
          "f[1, 2, 3 g[a]]", "x + g[]", "x + 2 g[a]", "g[]^2", "x^g[]"]

POLY_SUBJ = ["x", "2 x", "a + x", "a + b x", "x^2", "a x^3", "1 + x",
             "(a + b x)^3", "(1 + x)^2", "a", "f[]", "f[1]", "f[1, 2]",
             "x + y", "3", "b x + a + c"]

LITERAL_SUBJ = ["1/2", "-1/3", "3 + 4 I", "I", "a + b", "a - b", "-a",
                "1/a", "a/b", "Sqrt[2]", "Sqrt[a]", "E^a", "Exp[2]",
                "a b c", "a^-2", "2^(1/3)", "1/2 + I/3", "x", "5", "2.5"]

PATDATA_SUBJ = ["_", "x_", "f[_]", "f[x_]", "f[a]", "HoldPattern[a]",
                "a", "f[_, _]", "_Integer", "x__"]

ASSOC_SUBJ = ["<|a -> 1|>", "<|a -> 1, b -> 2|>", "<|b -> 2|>", "<||>",
              "{a -> 1, b -> 2}", "{b -> 2}", "a -> 1", "{a -> 1}",
              "f[a -> 1]", "{}"]

ATTR_SUBJ = ["h[a]", "h[a, b]", "h[b, a]", "h[a, b, c]", "h[1, a]",
             "h[a, a]", "h[a, b, a]", "h[h[a, b], c]", "h[]", "a"]

# --------------------------------------------------------------------------
# Pattern families: (title, subjects, patterns). H(...) builds both the f[..]
# and the {..} spelling of a sequence pattern.
# --------------------------------------------------------------------------

def H(*seqs):
    out = []
    for s in seqs:
        out.append("f[%s]" % s)
        out.append("{%s}" % s)
    return out


BLANK = [
    "_", "_Integer", "_Real", "_Rational", "_Complex", "_String", "_Symbol",
    "_f", "_g", "_List", "_Plus", "_Times", "_Power", "_Sin",
    "f[_]", "f[_, _]", "f[_f]", "f[_, _, _]", "f[]", "_[_]", "_[]",
    "_[_, _]", "_f[_]", "f[_][_]", "_[_][_]", "_?AtomQ", "_[___]",
    "f[_Integer, _]", "{_, _}", "{_Integer..}",
]

NAMED_BLANK = [
    "x_", "x_Integer", "f[x_]", "f[x_, y_]", "h_[x_]", "h_[x_, y_]",
    "h_[]", "f[x_][y_]", "h_[x_][y_]", "x_f", "x_List",
]

SEQ = H("__", "___", "x__", "x___", "x__, y__", "x___, y___", "x__, y___",
        "x___, y__", "x__, a", "a, x__", "x___, a, y___", "x__Integer",
        "x___Symbol", "x__, x__", "x___, b, x___", "x_, y__, z_",
        "__, _", "___, _, ___", "x___, 2, y___", "x__f", "x___, y_, z___",
        "x__, y_Integer", "x_, y___", "x___, y_")

NAMED = H("x_, x_", "x_, y_", "x_, f[x_]", "x_, g[x_]", "x_, _, x_",
          "x_Symbol, x_", "x_, {x_}", "x_, y_, x_") + [
    "x_ + x_", "x_ x_", "x_^x_", "x_ + y_", "x_ y_", "f[x_, f[y_]]",
    "f[g[x_]]", "f[x : g[_]]", "x : f[_, _]", "f[x : _, x_]",
]

# The issue #102 family: `sym:patt` composed with every operator.
COLON = [
    "n_*a:_g", "n_.*a:_g", "f[__, n_*a:_g]", "f[_, n_.*a:_g]",
    "f[__, n_.*a:_g]", "f[_, _, n_.*a:_g]", "x_ + a:_g", "x_. + a:_g",
    "a:_g^n_", "a:_g^n_.", "x_^a:_g", "f[a:_g]", "f[a:_g..]",
    "f[a:(_Integer | _g)]", "f[a:_, n_.*b:_g]", "a:_g", "f[a:_:0]",
    "f[x_, a:_g:g[0]]", "-a:_g", "f[_, n_Integer*a:_g]",
    "n_*a:_g /; n > 2", "f[x_, n_.*a:_g]", "a:_g + b:_", "f[a:_g, ___]",
    "f[___, n_.*a:_g, ___]", "n_.*a:_g^m_.",
]

OPTIONAL = H("x_, y_:0", "x_:1, y_", "x_:1, y_:2", "x_:0, y___",
             "x___, y_:0", "x_Integer:5", "x_:0, y_:0, z_:0", "x_:1") + [
    "x_. + y_", "x_.*y_", "x_^n_.", "a_. + b_.*x", "(a_. + b_.*x)^n_.",
    "Optional[x_, 0] + y_", "x_.*y_^n_.", "a_. + b_.*x^2", "c_.*x^n_.",
    "p:x_.*y_", "x_ + y_.", "x_ y_.",
]

ALTERNATIVES = H("a | b", "x : (a | b)", "x_Integer | x_Symbol",
                 "(a | b)..", "x__ | y_", "a | f[_]", "x : (_Integer | _Real)",
                 "___, a | b, ___") + [
    "a | b", "_Integer | _Real", "_List | _f", "x_Integer | x_String",
    "f[a] | f[b] | g[_, _]", "(_ | _) + _",
]

TESTCOND = H("x_?EvenQ", "x__?IntegerQ", "x_ /; x > 1", "x_, y_ /; y > x",
             "x__ /; Length[{x}] == 2", "x_?(# > 2 &)", "x_Integer?Positive",
             "x_?NumericQ, y_", "x___ /; Length[{x}] > 1, y_") + [
    "_?IntegerQ", "_?(# === a &)", "x_ /; IntegerQ[x]",
    "f[x_, y_] /; x < y", "x_ /; MemberQ[{a, b}, x]", "_?AtomQ",
]

REPEATED = [
    "{a..}", "{a...}", "{_Integer..}", "{x:_Integer..}", "{(a | b)..}",
    "f[a..]", "f[x:a..]", "{Repeated[_, {2}]}", "{Repeated[_, {1, 2}]}",
    "{Repeated[_, {0, 1}]}", "{Repeated[a, 2]}", "{_.., _}", "{{_, _}..}",
    "{x_..}", "{RepeatedNull[a], b}", "f[PatternSequence[a, b]..]",
    "{x:{_, _}..}", "f[x__..]", "{_Symbol..}", "{a.., b...}", "{x : a..., y___}",
    "f[Repeated[_Integer, {2, 3}]]", "f[g[_]..]", "{Repeated[x_, {2}]}",
]

EXCEPT = [
    "Except[a]", "Except[_Integer]", "f[Except[a]]", "f[Except[a]..]",
    "Except[_Integer, _?NumericQ]", "f[x:Except[0]]", "f[Except[_f], ___]",
    "{Except[a]...}", "Except[_List, _[___]]", "f[Except[a, _Symbol]..]",
    "HoldPattern[_ + _]", "HoldPattern[f[x_]]", "HoldPattern[x_]",
    "Verbatim[_]", "Verbatim[x_]", "f[Verbatim[_]]", "Verbatim[f][_]",
    "f[Verbatim[x_]]", "HoldPattern[Verbatim[_]]", "Verbatim[f[_]]",
]

PATSEQ = [
    "f[PatternSequence[a, b]]", "f[PatternSequence[a, b], c]",
    "f[x:PatternSequence[_, _], ___]", "{PatternSequence[a, b]..}",
    "f[OrderlessPatternSequence[a, b]]", "f[OrderlessPatternSequence[a, b], ___]",
    "f[OrderlessPatternSequence[x_Integer, y_Symbol]]",
    "f[PatternSequence[]]", "f[___, PatternSequence[a, b], ___]",
    "{x:PatternSequence[_, _]..}", "f[PatternSequence[x_, x_]]",
    "f[OrderlessPatternSequence[x_, a]]",
]

LONGSHORT = H("Longest[x__], y___", "Shortest[x__], y__",
              "x___, Longest[y__]", "Longest[x___], y___",
              "Shortest[x___], Shortest[y___]", "Shortest[x__], y___",
              "x___, Shortest[y__]", "Longest[x:a...], y___")

LITERAL = [
    "Rational[a_, b_]", "Complex[a_, b_]", "x_ + y_", "x_*y_", "x_^y_",
    "Plus[x__]", "Times[-1, x_]", "Sqrt[x_]", "1/x_", "-x_", "x_/y_",
    "E^x_", "Exp[x_]", "x_^Rational[1, n_]", "Complex[_, b_]",
    "Times[x__]", "x_Rational", "x_Complex", "Power[x_, -1]",
    "x_Integer/y_Integer", "a_ + b_ I",
]

KVP = [
    "KeyValuePattern[{a -> _}]", "KeyValuePattern[{a -> x_, b -> y_}]",
    "KeyValuePattern[a -> x_]", "KeyValuePattern[{}]", "<|a -> _|>",
    "_Association", "{___Rule}", "{(_ -> _)..}", "x_Rule", "<|x_ -> y_|>",
    "KeyValuePattern[{b -> 2}]", "_ -> _Integer",
]

PATTERN_FAMILIES = [
    ("Blank and head-typed Blank (_, _h, head patterns)", CORE, BLANK),
    ("Named blanks and head-binding patterns", CORE, NAMED_BLANK),
    ("BlankSequence / BlankNullSequence", CORE, SEQ),
    ("Named patterns and nonlinear reuse", CORE, NAMED),
    ("sym:patt composed with operators (issue #102)", G_SUBJ + CORE[:10], COLON),
    ("Optional and default values", POLY_SUBJ + CORE[8:20], OPTIONAL),
    ("Alternatives", CORE, ALTERNATIVES),
    ("PatternTest and Condition", CORE, TESTCOND),
    ("Repeated / RepeatedNull", CORE + ["{a, a, a}", "{1, 2}", "{a, b, b}",
                                       "{{1, 2}, {3, 4}}", "f[a, b, a, b]",
                                       "f[g[1], g[2]]", "{a, a, b, b}"],
     REPEATED),
    ("Except, HoldPattern, Verbatim", CORE + PATDATA_SUBJ, EXCEPT),
    ("PatternSequence / OrderlessPatternSequence",
     CORE + ["f[b, a]", "f[a, b, a, b]", "f[x, 1]", "f[1, x]", "f[a, b, b]"],
     PATSEQ),
    ("Longest / Shortest (ReplaceList order)",
     ["f[a, b, c]", "f[a, a, b]", "f[a]", "f[]", "{a, b, c, d}", "{a, a, a}"],
     LONGSHORT),
    ("Literal structure of numbers and arithmetic", LITERAL_SUBJ, LITERAL),
    # Replace, not ReplaceList: the reference's ReplaceList over a
    # KeyValuePattern drops bindings ({{1}} for {a -> x_, b -> y_} :> {x, y}).
    ("KeyValuePattern and Association/Rule patterns", ASSOC_SUBJ, KVP, "Replace"),
]

NAME_RE = re.compile(r"(?<![A-Za-z0-9`$])([a-z][A-Za-z0-9]*)(?=_|\s*:(?![:=>]))")


def pattern_names(p):
    """Names a pattern binds, in order of first appearance."""
    seen = []
    for m in NAME_RE.finditer(p):
        n = m.group(1)
        if n not in seen:
            seen.append(n)
    return seen


def pattern_cases(subjects, patterns, query="ReplaceList"):
    out = []
    for p in patterns:
        names = pattern_names(p)
        for s in subjects:
            if names:
                out.append("%s[%s, %s :> {%s}]" % (query, s, p, ", ".join(names)))
            else:
                out.append("MatchQ[%s, %s]" % (s, p))
    return out

# --------------------------------------------------------------------------
# Attributes on a user head h, declared inside each case.
# --------------------------------------------------------------------------

ATTR_SETS = ["{Orderless}", "{Flat}", "{Flat, Orderless}", "{Flat, OneIdentity}",
             "{OneIdentity}", "{Flat, Orderless, OneIdentity}"]
ATTR_PATS = ["h[a, x_]", "h[x_, a]", "h[x_, y_]", "h[x__, y_]", "h[a, b]",
             "h[x_Integer, y_]", "h[a, x__]", "h[x_, x_]", "h[x_]",
             "h[x_, h[y_]]", "h[x___, b]", "h[x_Symbol, y___]"]


def attr_cases():
    out = []
    for a in ATTR_SETS:
        for p in ATTR_PATS:
            names = pattern_names(p)
            for s in ATTR_SUBJ:
                q = ("ReplaceList[%s, %s :> {%s}]" % (s, p, ", ".join(names))
                     if names else "MatchQ[%s, %s]" % (s, p))
                out.append("SetAttributes[h, %s]; %s" % (a, q))
    return out

# --------------------------------------------------------------------------
# Query heads over fixed data at level specs.
# --------------------------------------------------------------------------

QDATA = ["{1, a, f[2], {3, b}, f[g[4]]}", "f[a, g[b, h[c]], 1, 2.5]",
         "{{1, 2}, {3, 4}, {5}}", "{a, {a, {a}}, b}"]
QPATS = ["_Integer", "_Symbol", "f[_]", "_List", "_f", "_?NumericQ",
         "_[_]", "_", "a", "{_, _}"]
QLEVELS = ["{1}", "{2}", "Infinity", "{0, Infinity}", "{-1}", "{1, 2}", "{2, -1}"]


def query_cases():
    out = []
    for d in QDATA:
        for p in QPATS:
            out.append("Cases[%s, %s]" % (d, p))
            out.append("Count[%s, %s]" % (d, p))
            out.append("Position[%s, %s]" % (d, p))
            out.append("DeleteCases[%s, %s]" % (d, p))
            out.append("FreeQ[%s, %s]" % (d, p))
            out.append("MemberQ[%s, %s]" % (d, p))
            out.append("FirstCase[%s, %s]" % (d, p))
            for lv in QLEVELS:
                out.append("Cases[%s, %s, %s]" % (d, p, lv))
                out.append("Count[%s, %s, %s]" % (d, p, lv))
                out.append("Position[%s, %s, %s]" % (d, p, lv))
            out.append("Cases[%s, %s, Infinity, Heads -> True]" % (d, p))
            out.append("Position[%s, %s, Infinity, Heads -> False]" % (d, p))
    rules = ["x_Integer :> x + 1", "f[x_] :> g[x]", "{x_, y_} :> x + y",
             "a -> z", "x_Symbol :> s[x]", "_List :> L"]
    for d in QDATA:
        for r in rules:
            out.append("Replace[%s, %s]" % (d, r))
            out.append("Replace[%s, %s, {1}]" % (d, r))
            out.append("Replace[%s, %s, Infinity]" % (d, r))
            out.append("ReplaceAll[%s, %s]" % (d, r))
            out.append("Cases[%s, %s, Infinity]" % (d, r))
    out += [
        "ReplaceRepeated[f[f[f[a]]], f[x_] :> x]",
        "ReplaceRepeated[{3, 1, 2}, {x___, a_, b_, y___} /; a > b :> {x, b, a, y}]",
        "ReplaceRepeated[x^2 + x^3, x^n_ :> n x]",
        "ReplaceAll[f[a, b], {f[x_, y_] :> {x, y}, a -> z}]",
        "ReplaceAll[{a, b, c}, {a -> b, b -> c}]",
        "Replace[f[a], {{f[x_] :> 1}, {f[x_] :> 2}}]",
        "ReplaceList[f[a, b, c], f[x__, y__] :> {{x}, {y}}]",
        "ReplaceList[{a, b, c}, {x___, y_, z___} :> y]",
        "ReplaceList[a + b + c, x_ + y_ :> {x, y}]",
        "ReplaceList[a b c, x_ y_ :> {x, y}]",
        "ReplaceList[f[a, b, c, d], f[x__, y__, z__] :> {{x}, {y}, {z}}, 3]",
        "ReplaceAll[f[a, g[a]], a -> b]",
        "ReplaceAll[Hold[a + a], a -> 1]",
        "ReplaceAll[a + b, Plus -> List]",
        "ReplaceAll[{f[1], f[2]}, f[x_] /; x > 1 :> g[x]]",
        "Cases[{1, 2, 3, 4}, x_ /; OddQ[x] :> x^2]",
        "Cases[{f[1, 2], f[3], g[4]}, f[x__] :> {x}]",
        "DeleteCases[{a, b, a, c}, a, {1}, 1]",
        "Cases[{1, 2, 3, 4, 5}, _Integer, {1}, 2]",
        "Position[{a, b, a, a}, a, {1}, 2]",
        "Count[{{a}, a, {{a}}}, a, Infinity]",
        "ReplaceList[{1, 2, 3}, {___, x_, ___, y_, ___} /; x < y :> {x, y}]",
    ]
    return out

# --------------------------------------------------------------------------
# Definitions: DownValue specificity, Optional/sequence arguments, Condition
# and recursion. Hand-written programs, still graded by the reference kernel.
# --------------------------------------------------------------------------

DEFS = [
    "q[x_] := 1; q[x_Integer] := 2; q[0] := 3; {q[0], q[5], q[a]}",
    "q[x_Integer] := 2; q[x_] := 1; q[0] := 3; {q[0], q[5], q[a]}",
    "q[x_, y_:3] := x + y; {q[1], q[1, 2], q[]}",
    "q[x__] := {x}; q[x_] := single; {q[1], q[1, 2], q[]}",
    "q[x___] := {x}; {q[], q[1], q[1, 2]}",
    "q[x_, x_] := same; q[x_, y_] := diff; {q[1, 1], q[1, 2]}",
    "q[x_] /; x > 0 := pos; q[x_] := other; {q[1], q[-1], q[a]}",
    "q[x_ /; x > 0] := pos; {q[1], q[-1]}",
    "q[x_?EvenQ] := even; q[x_?OddQ] := odd; {q[2], q[3], q[a]}",
    "q[f[x_]] := x; q[g[x_, y_]] := x y; {q[f[2]], q[g[2, 3]], q[h[1]]}",
    "q[0] = 1; q[n_Integer] := n q[n - 1]; q[5]",
    "q[n_] := q[n] = If[n < 2, n, q[n - 1] + q[n - 2]]; q[20]",
    "q[{x___, y_}] := y; {q[{1, 2, 3}], q[{}]}",
    "q[x : {___Integer}] := Total[x]; {q[{1, 2}], q[{1, a}]}",
    "q[x_. + y_] := {x, y}; {q[a], q[a + b]}",
    "q[n_.*x_Symbol] := {n, x}; {q[x], q[3 x]}",
    "q[n_.*a : _g] := {n, a}; {q[g[]], q[2 g[]], q[x]}",
    "q[a : _g, n_.*b : _g] := {a, n, b}; {q[g[1], g[2]], q[g[1], 3 g[2]]}",
    "q[a_, b_:Automatic] := {a, b}; {q[1], q[1, 2]}",
    "Default[q] = 7; q[x_.] := x; {q[], q[1]}",
    "Default[q, 2] = 9; q[x_, y_.] := {x, y}; {q[1], q[1, 2]}",
    "SetAttributes[q, Orderless]; q[a, x_] := x; {q[b, a], q[a, b]}",
    "SetAttributes[q, Flat]; q[x_, y_] := {x, y}; q[1, 2, 3]",
    "SetAttributes[q, HoldAll]; q[x_ + y_] := {Hold[x], Hold[y]}; q[1 + 2]",
    "q[x_] := x^2; q[y_] := y^3; {q[2], DownValues[q] // Length}",
    "q[x_Integer | x_Real] := num; q[_] := other; {q[1], q[1.5], q[a]}",
    "q[x : (_Integer | _Real)] := x; {q[1], q[a]}",
    "q[Except[0, x_Integer]] := 1/x; {q[2], q[0], q[a]}",
    "q[x_List?(Length[#] > 1 &)] := many; q[_List] := few; {q[{1}], q[{1, 2}]}",
    "q[HoldPattern[Plus[x__]]] := {x}; q[a + b + c]",
    "q[Verbatim[_]] := blank; q[_] := other; {q[_], q[x_], q[a]}",
    "q[x_, OptionsPattern[{o -> 1}]] := {x, OptionValue[o]}; {q[1], q[1, o -> 2]}",
    "Options[q] = {o -> 1, p -> 2}; q[x_, OptionsPattern[]] := {x, OptionValue[o], OptionValue[p]}; {q[1], q[1, p -> 5], q[1, {o -> 3}]}",
    "q[x_][y_] := x + y; q[1][2]",
    "q[n_Integer?Positive] := n; q[_] := bad; {q[3], q[-3], q[3.]}",
    "q[{a_, b_}, {c_, d_}] := a d - b c; q[{1, 2}, {3, 4}]",
    "q[x_] := Module[{y = x}, y + 1] /; x > 0; {q[1], q[-1]}",
    "q[x : f[___]] := Length[x]; {q[f[]], q[f[a, b]]}",
    "q[PatternSequence[a_, b_], c_] := {a, b, c}; q[1, 2, 3]",
    "q[x__Integer] := Plus[x]; {q[1, 2, 3], q[1, a]}",
    "q[s_String] := StringLength[s]; {q[\"abc\"], q[abc]}",
]

# --------------------------------------------------------------------------
# Parse shapes of the pattern operators (the #102 table and more).
# --------------------------------------------------------------------------

PARSE = [
    "n_.*a:_g", "n_*a:_g", "x+a:_g", "x + a:_g == b", "x + a:_g + c",
    "n a:_g", "n a:_g b", "a:_g b", "a.b:_g", "x_:0", "a:b:c", "a:b:c:d",
    "x_:a:b", "a:x_:0", "-a:_g", "!a:_", "2 a:_g", "b/a:_g", "x^a:_g",
    "f@a:_g", "a:_g->1", "a:_/;c", "a:_g..", "a:b|c", "x|a:b", "x|a:b|c",
    "a ~~ x:b", "x:a ~~ b", "a:_+b", "x->a:_g", "{a:_g, b:_h}",
    "f[a:_g, n_.*b:_h]", "a:_g:1", "n_.*a:_g:1", "x__*a:_g", "a:_g?f",
    "x?f:_", "a:_g&&b", "x&&a:_g", "x==a:_g", "x>a:_g", "a:_g;b",
    "x: _", "x : _g * y", "x_^n_.*a:_", "Times[a:_g]", "p:{___}",
    "f[x:_:0]", "n x_:0", "n_.*a:_g:1 + c", "a:_g:1+c", "x + a:b:c",
    "a:_:1 b", "x + a:_g:1 == b", "x + a:_g:1 -> b", "x /; a:_g",
    "x:a:b:c:d:e", "a:_g:b:_h", "2 a:_g:1", "f[x] := n_*a:_g",
    "a:_g:1?f", "x^a:_:1", "{n_*a:_g, x}", "a:_g + x_:0", "a:f[b:c:d]",
    "a:{b:c, d:e:f}", "a:(b:c)", "a:b|c:d", "a:_?(x:_)", "x_:a:b:c",
    "x_:a:b:c:d", "a:b..:c", "x:_ ~~ y:_", "x_:0 ~~ y", "a:x_+y_:1",
    "f[x_, y_:0 z]", "f[x_, y_:0, z_:1]", "a:b->c:d", "a:b:>c:d",
    "x_:a+b:c", "x_h:0", "x__:0", "x___h:0", "_:0", "__:0", "x_:0?f",
    "x:_:0 /; y", "h[a:_:0, b:_.]", "x_:a:b|c", "a:-b:c", "a:!b:c",
    "f[_, n_.*a:_g:2]", "p:x_*y_", "p:x_*y_:1", "x:a.b:c",
    "_", "__", "___", "_h", "x_", "x__h", "x___h", "_.", "x_.", "x_h.",
    "_..", "_...", "x_..", "a..", "a...", "_Integer..", "x__..",
    "a|b|c", "a|b..", "(a|b)..", "x_?f", "x_?f..", "x_?(#>1&)",
    "x_ /; x>1", "f[x_] /; c", "f[x_ /; c]", "x_ -> y", "x_ :> y /; c",
    "f[x_] := y /; c", "_?f /; c", "x:_?f", "x_h?f", "a|x_?f",
    "Except[a]|b", "x__ ~~ y__", "x_ ~~ \"a\" ~~ y_", "a ~~ b | c",
    "f[x_, y_] :> g[x] /; y > 0", "x_ -> y_ -> z", "a:b->c",
    "Verbatim[_]", "HoldPattern[x_+y_]", "f[__, _.]",
]

# --------------------------------------------------------------------------

def build_cases():
    sections = []
    for fam in PATTERN_FAMILIES:
        title, subjects, patterns = fam[:3]
        query = fam[3] if len(fam) > 3 else "ReplaceList"
        sections.append((title, pattern_cases(subjects, patterns, query)))
    sections.append(("Attributes on a user head (Orderless / Flat / OneIdentity)",
                     attr_cases()))
    sections.append(("Query heads and level specifications", query_cases()))
    sections.append(("Definitions: specificity, defaults, conditions", DEFS))
    sections.append(("Parse shapes of the pattern operators",
                     ["ToString[FullForm[Hold[%s]]]" % s for s in PARSE]))
    return sections


KERNEL_DRIVER = r'''
cases = Import["%(inp)s", "RawJSON"];
res = Table[
  Block[{$Context = "MthCase" <> ToString[i] <> "`", $ContextPath = {"System`"}},
    ToString[Quiet[Check[TimeConstrained[ToExpression[cases[[i]]], 2,
                                         "MTH$TIMEOUT"], "MTH$MESSAGE"]],
             InputForm]],
  {i, Length[cases]}];
Export["%(out)s", res, "RawJSON"];
'''

# A result that cannot be checked by re-reading its text: a message or a
# timeout, a temporary or context-qualified symbol, a held/opaque object.
REJECT = re.compile(r"MTH\$|\$|`|\\\[|Removed\[|Unique|Graphics|"
                    r"InterpolatingFunction|CompiledFunction|SparseArray|"
                    r"Dispatch|\bTemporary\b")


def kernel_eval(cases, kernel):
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "cases.json")
        out = os.path.join(d, "out.json")
        drv = os.path.join(d, "drive.wl")
        with open(inp, "w") as f:
            json.dump(cases, f)
        with open(drv, "w") as f:
            f.write(KERNEL_DRIVER % {"inp": inp, "out": out})
        subprocess.run([kernel, "-file", drv], check=True, timeout=3600)
        with open(out) as f:
            return json.load(f)


def m_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


# --------------------------------------------------------------------------
# Known divergences. A case matching one of these goes to the report-only
# observe corpus (match_exhaustive_observe.m, run with --observe) instead of
# the gate, with the reason recorded beside it. `sorted_gate` marks an
# ordering-only difference: the same matches, listed in another order. Such a
# case also gets a GATING twin, Sort[ReplaceList[...]], so the set of matches
# stays pinned even though its order is not. The twin compares against the reference's
# list sorted by the RUNNER's own Sort (Sort[...] === Sort[{wl...}], expected
# True), so it is an exact multiset check that does not depend on the two
# systems agreeing on canonical order (they do not, e.g. {2, x} vs {1, 2 x}).
# --------------------------------------------------------------------------

OBSERVE = [
    (re.compile(r"^SetAttributes\[h, \{Flat(, Orderless)?\}\]; ReplaceList"),
     "Flat without OneIdentity: the reference binds a one-operand slice both as h[b] and "
     "as b, and repeats matches; Mathilda binds b once", False),
    (re.compile(r"^SetAttributes\[q, Flat\]"),
     "Flat without OneIdentity in a definition: the reference wraps single operands "
     "(x -> q[1])", False),
    (re.compile(r"^SetAttributes\[h, \{Flat, Orderless, OneIdentity\}\]; ReplaceList"),
     "Flat+Orderless enumeration order", True),
    (re.compile(r"^ReplaceList\[.*, (x_ \+ y_\.|x_ y_\.|x_\.\*y_\^n_\.) :>"),
     "Orderless enumeration order: the reference's canonical order puts y_. before x_ "
     "(patterns sort after other expressions), Mathilda's does not", True),
    (re.compile(r"^q\[x_\]\[y_\] :="),
     "SubValues (f[x_][y_] := ...) are not implemented", False),
]


def classify(case):
    for rx, why, sorted_gate in OBSERVE:
        if rx.search(case):
            return why, sorted_gate
    return None, False


def write_corpus(path, header, sections):
    """sections: [(title, [(input, expected, note-or-None)])]"""
    lines = header + ["{"]
    blocks = []
    for idx, (title, entries) in enumerate(sections, 1):
        if not entries:
            continue
        rows = []
        for inp, exp, note in entries:
            row = "    {%s, %s}" % (m_string(inp), exp)
            if note:
                row = "    (* %s *)\n%s" % (note, row)
            rows.append(row)
        blocks.append("    (* ---- %d. %s (%d cases) ---- *)\n%s"
                      % (idx, title, len(entries), ",\n".join(rows)))
    lines.append(",\n".join(blocks))
    lines.append("}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kernel", default=DEFAULT_KERNEL,
                    help="reference kernel script runner (or set REFERENCE_KERNEL)")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--observe-out", default=OBS_OUT)
    args = ap.parse_args()

    # Duplicates across families add runtime, not coverage.
    seen = set()
    sections = []
    for title, cs in build_cases():
        keep = [c for c in cs if not (c in seen or seen.add(c))]
        sections.append((title, keep))
    queries = [c for _, cs in sections for c in cs]
    if not args.kernel:
        ap.error("no reference kernel: pass --kernel PATH or set REFERENCE_KERNEL")
    print("evaluating %d cases in the reference kernel..." % len(queries), file=sys.stderr)
    expect = dict(zip(queries, kernel_eval(queries, args.kernel)))

    gate, observe = [], []
    n_gate = n_obs = dropped = 0
    for title, cs in sections:
        g, o = [], []
        for c in cs:
            r = expect[c]
            if REJECT.search(r):
                dropped += 1
                continue
            why, _ = classify(c)
            if why:
                o.append((c, r, why))
                if classify(c)[1]:      # ordering-only: gate the multiset
                    g.append(("Sort[%s] === Sort[%s]" % (c, r), "True", None))
            else:
                g.append((c, r, None))
        gate.append((title, g)); observe.append((title, o))
        n_gate += len(g); n_obs += len(o)

    common = [
        " * tools/gen_match_exhaustive_corpus.py; do not edit by hand, re-run the",
        " * generator. Every expected value was computed by the reference",
        " * kernel, each case in a fresh context. Same {inputString, expected}",
        " * format as match_stress_corpus.m; run by test_match_stress_corpus.c.",
    ]
    write_corpus(args.out, ["(* match_exhaustive_corpus.m -- GENERATED by"] + common + [
        " * Gating: every case must match.",
        " *)"], gate)
    write_corpus(args.observe_out, ["(* match_exhaustive_observe.m -- GENERATED by"] + common + [
        " * REPORT-ONLY (run with --observe): known divergences from the reference,",
        " * each annotated with its reason. Ordering-only ones have a gating",
        " * Sort[...] twin in match_exhaustive_corpus.m.",
        " *)"], observe)
    print("wrote %d gating + %d observe cases (%d dropped)" % (n_gate, n_obs, dropped),
          file=sys.stderr)


if __name__ == "__main__":
    main()
