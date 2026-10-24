#!/usr/bin/env python3
"""
sparse_head_sweep.py -- which builtins treat a SparseArray differently from the
reference system?

tools/sparse_stress.py checks the invariant Normal[f[s]] == f[Normal[s]] over a
curated list of operations, so it cannot see a head nobody listed. This sweep
names nobody: it takes every System symbol both systems define, calls it on
SparseArray arguments in a handful of call shapes, and runs the same calls
through the reference (wolframscript). A row is reported when the two
disagree, after Normal, on the answer or on whether the call evaluated -- an
unevaluated call, a SparseArray treated as an atom, or a silent no-op.

The report needs triage, not a blind fix: the reference iterates a SparseArray
in stored-values-then-background order (Scan / Select / Do differ by design;
see docs/spec/builtins/data-structures.md), refuses some heads (Sort, Fold,
Union) that Mathilda answers, and prints numbers differently (1.0 vs 1.).

Usage:
    python3 tools/sparse_head_sweep.py [--mathilda ./Mathilda]
        [--wolframscript PATH] [--only Head1,Head2] [--out report.tsv]

The reference binary is --wolframscript, else $WOLFRAMSCRIPT, else
`wolframscript` on PATH.

Exits 0; the report is written to --out (default: sparse_head_sweep.tsv in the
current directory) and the mismatch count is printed.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

# Heads with side effects, interactive behaviour, or no meaning on an array.
BLACKLIST_RE = re.compile(
    r"(File|Director|Dialog|Export|Import|Put|Get|Run|Install|Link|Quit|Exit|"
    r"Abort|Pause|Print|Message|Interrupt|Uninstall|Remove|Clear|Unset|Protect|"
    r"Needs|Begin|End|Dynamic|Notebook|Front|Environment|Process|Socket|URL|"
    r"Cloud|Database|Speak|Beep|Echo|Monitor|Reset|Kernel|Parallel|Launch|"
    r"Close|Open|Read|Write|Stream|Create|Copy|Rename|Splice|Save|Dump|"
    r"Timing|Time|Date|Random|Seed|Share|Memory|Graphics|Plot|Show|Image|"
    r"Sound|Play|Audio|Video|Package|Context|Attributes|Options|Values|"
    r"Definition|Information|Names|Trace|Stack|Hold|Unevaluated|Evaluate|"
    r"Compile|Function|Set|Throw|Catch|Return|Break|Continue|Goto|Label|"
    r"Sow|Reap|Check|Quiet|Off|On$|Block|Module|With|Do$|For$|While|Table|"
    r"Sum$|Product$|Array$|Assert|Abs.*Time|Constrained|Interval|Locked|"
    r"Pipe|Encode|Decode|Hash|System|Language|Version|Machine|Precision$)"
)

SHAPES = [
    "h[v]", "h[m]", "h[v, 2]", "h[m, 2]", "h[ff, v]", "h[ff, m]", "h[v, ff]",
    "h[v, v]", "h[m, m]", "h[m, v]",
]

# The fixtures are rebuilt before every call: a mutating head (AppendTo,
# AddTo, ...) would otherwise change them for every later call.
PRELUDE = r"""
fix[] := (ClearAll[ff];
  v = SparseArray[{1 -> 2, 3 -> -1}, 3];
  m = SparseArray[{{1, 1} -> 1, {2, 3} -> 2, {3, 2} -> -1}, {3, 3}];
  If[TrueQ[$denseFixtures], v = Normal[v]; m = Normal[m]]);
nrm[r_] := StringTake[ToString[r /. a_SparseArray :> Normal[a], InputForm], UpTo[400]];
R[id_, e_] := (fix[]; Module[{r = TimeConstrained[Quiet[e], 5, $TimedOut]},
  WriteString["stdout", id, "\t", ToString[Head[r]], "\t", nrm[r], "\n"]]);
SetAttributes[R, HoldRest];
"""


def run(cmd, text, timeout):
    try:
        p = subprocess.run(cmd, input=text, capture_output=True, text=True, timeout=timeout)
        return p.stdout
    except subprocess.TimeoutExpired as e:
        return e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")


def names(mathilda, ws):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        f.write('Scan[WriteString["stdout", #, "\\n"] &, Names["*"]]\n')
        path = f.name
    mn = set(run([mathilda, "-file", path], "", 120).split())
    os.unlink(path)
    wn = set(run([ws, "-code", 'Scan[WriteString["stdout", #, "\\n"] &, Names["System`*"]]'],
                 "", 300).split())
    return sorted(n for n in mn & wn
                  if n[:1].isupper() and not n.startswith("$") and not BLACKLIST_RE.search(n))


def calls(heads):
    out = []
    for h in heads:
        for k, s in enumerate(SHAPES):
            out.append((f"{h}#{k}", s.replace("h[", h + "[", 1)))
    return out


def program(cs, dense=False):
    pre = "$denseFixtures = True;\n" if dense else ""
    return pre + PRELUDE + "".join(f'R["{i}", {c}];\n' for i, c in cs)


def parse(text):
    res = {}
    for line in text.splitlines():
        parts = line.split("\t", 2)
        if len(parts) == 3 and "#" in parts[0]:
            res[parts[0]] = (parts[1], parts[2])
    return res


def norm(s):
    s = s.replace(" ", "").replace("*", "")
    s = re.sub(r"(\d)\.0(?![0-9])", r"\1.", s)
    s = re.sub(r"(\d)\.(?![0-9])", r"\1", s)        # 1. vs 1 is reported separately
    return s


def run_mathilda(mathilda, cs, chunk=150, dense=False):
    res = {}
    for i in range(0, len(cs), chunk):
        part = cs[i:i + chunk]
        with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
            f.write(program(part, dense))
            path = f.name
        got = parse(run([mathilda, "-file", path], "", 900))
        os.unlink(path)
        res.update(got)
        missing = [c for c in part if c[0] not in got]
        # a crash ends the chunk: re-run the rest one by one
        if missing and len(part) > 1:
            for c in missing:
                res.update(run_mathilda(mathilda, [c], 1, dense) or {c[0]: ("CRASH", "")})
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mathilda", default="./Mathilda")
    ap.add_argument("--wolframscript",
                    default=os.environ.get("WOLFRAMSCRIPT") or shutil.which("wolframscript"))
    ap.add_argument("--only", default="")
    ap.add_argument("--out", default="sparse_head_sweep.tsv")
    a = ap.parse_args()

    heads = a.only.split(",") if a.only else names(a.mathilda, a.wolframscript)
    cs = calls(heads)
    print(f"{len(heads)} heads, {len(cs)} calls", file=sys.stderr)
    with tempfile.NamedTemporaryFile("w", suffix=".wl", delete=False) as f:
        f.write(program(cs))
        wpath = f.name
    ref = parse(run([a.wolframscript, "-file", wpath], "", 7200))
    os.unlink(wpath)
    got = run_mathilda(a.mathilda, cs)
    dense = run_mathilda(a.mathilda, cs, dense=True)

    # Kinds, most important first:
    #   CRASH    the sparse call crashed Mathilda
    #   SPARSE   Mathilda's sparse answer differs from its own dense answer,
    #            which the reference agrees with (or the reference refuses):
    #            a sparse bug
    #   HEAD     same Normal everywhere, but SparseArray vs List differs from
    #            the reference
    #   DENSE    Mathilda's dense and sparse answers agree and the reference
    #            differs: a dense-side difference, not a sparse one
    rows = []
    for cid, call in cs:
        w = ref.get(cid)
        g = got.get(cid, ("MISSING", ""))
        d = dense.get(cid, ("MISSING", ""))
        if not w:
            continue
        h = cid.split("#")[0]
        if g[0] in ("CRASH", "MISSING"):
            rows.append(("CRASH", call, g[0], g[1], w[0], w[1]))
            continue
        same_gd = norm(g[1]) == norm(d[1]) or (g[0] == h and d[0] == h)
        same_gw = norm(g[1]) == norm(w[1]) or (g[0] == h and w[0] == h)
        if not same_gd and not same_gw:
            # SPARSE: the dense answer is the reference's (or the reference
            # refused); SPARSE?: all three differ -- triage by hand
            kind = "SPARSE" if norm(d[1]) == norm(w[1]) or w[0] == h else "SPARSE?"
            rows.append((kind, call, g[0], g[1], d[0], d[1]))
            continue
        if not same_gd:
            continue      # the sparse answer is the reference's (Head[s], ...)
        if not same_gw:
            rows.append(("DENSE", call, g[0], g[1], w[0], w[1]))
            continue
        sa = {g[0], w[0]} & {"SparseArray"}
        if sa and g[0] != w[0] and h not in (g[0], w[0]):
            rows.append(("HEAD", call, g[0], g[1], w[0], w[1]))

    with open(a.out, "w") as f:
        f.write("kind\tcall\tmathilda_head\tmathilda\tother_head\tother\n")
        for r in sorted(rows):
            f.write("\t".join(r) + "\n")
    counts = {}
    for r in rows:
        counts[r[0]] = counts.get(r[0], 0) + 1
    print(f"{len(rows)} mismatches {counts} -> {a.out}")


if __name__ == "__main__":
    main()
