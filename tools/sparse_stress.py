#!/usr/bin/env python3
"""
sparse_stress.py -- differential stress test for SparseArray (src/sparse/).

The one invariant every SparseArray operation promises is

    Normal[op[s1, s2, ...]]  ==  op[Normal[s1], Normal[s2], ...]

(element-wise functions, reductions, structural operations, linear algebra,
Part and Part assignment). This tool generates random arrays -- rank 1..3,
dimensions down to 1, every value kind (machine Integer / Real / Complex,
Rational, big integers, symbolic, mixed), zero / nonzero / symbolic
backgrounds, explicit entries equal to the background, repeated and unsorted
rules -- applies ~120 operations to them, and checks the invariant inside one
Mathilda process per chunk. A chunk that crashes is bisected down to the case
that crashed it.

The pattern heads (Count, Position, Cases, MemberQ, ReplaceAll, Level) are NOT
in the table: the language treats a SparseArray as an atom there, so the
invariant deliberately does not hold (see src/sparse/sparse_gate.c).

Usage:
    python3 tools/sparse_stress.py                 # default: 3000 cases, seed 1
    python3 tools/sparse_stress.py -n 20000 --seed 7
    python3 tools/sparse_stress.py --valgrind -n 200
    python3 tools/sparse_stress.py --keep FILE      # write the generated script

Exit status 0 when every case holds (or is EXEMPT), 1 otherwise.
"""

import argparse
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BINARY = os.path.join(ROOT, "Mathilda")

# Operations by arity class. `S`, `T` are sparse operands of the case's shape;
# `v` a dense vector matching the last dimension, `k` a small positive int,
# `x` a symbol. Each entry is a Mathilda template over those names.
ELEMENTWISE = [
    "S + T", "S - T", "S T", "S / 2", "2 S", "S + 1", "S - x", "S^2", "S^3",
    "-S", "S + x", "x S", "S + 2.5", "1.5 S", "S^0", "2^S", "S T + S",
    "Abs[S]", "Sign[S]", "Floor[S]", "Ceiling[S]", "Round[S]", "Re[S]", "Im[S]",
    "Conjugate[S]", "Sin[S]", "Cos[S]", "Exp[S]", "N[S]", "Chop[N[S]]",
    "Mod[S, 3]", "Max[S, 0]", "Min[S, 1]", "Clip[S]", "Unitize[S]", "Boole[Positive[S]]",
    "Positive[S]", "UnitStep[S]", "Sqrt[S^2]", "Log[1 + S^2]",
]
REDUCTION = [
    "Total[S]", "Total[S, 2]", "Total[S, -1]", "Total[S, {2}]", "Max[S]", "Min[S]",
    "MinMax[S]", "Mean[S]", "Tr[S]", "Accumulate[S]", "Differences[S]",
    "Length[S]", "Dimensions[S]", "ArrayDepth[S]", "Depth[S]", "Norm[Flatten[S]]",
    "Total[Flatten[S]]", "Variance[Flatten[N[S]]]", "Median[Flatten[N[S]]]",
]
STRUCTURAL = [
    "Transpose[S]", "Flatten[S]", "Reverse[S]", "Take[S, 1]", "Take[S, -1]",
    "Drop[S, 1]", "Rest[S]", "Most[S]", "First[S]", "Last[S]", "Join[S, T]",
    "Join[S, S]", "Append[S, First[S]]", "Prepend[S, Last[S]]", "RotateLeft[S]",
    "RotateRight[S, 2]", "ArrayReshape[S, {Times @@ Dimensions[S]}]",
    "Partition[Flatten[S], 1]", "S[[1]]", "S[[-1]]", "S[[All, 1]]",
    "S[[1 ;; -1]]", "S[[{1, 1}]]", "S[[All, -1]]", "S[[2 ;; All]]",
    "Map[f, S]", "Apply[List, S]", "Apply[Plus, S]", "Map[Total, S]",
    "PadRight[S, Dimensions[S] + 1]", "ArrayPad[S, 1]", "Thread[S]",
    "Outer[Times, Flatten[S], {1, 2}]", "Equal[S, Normal[S]]",
    "S == T", "S === S",
]
MATRIX = [  # rank-2 square
    "Transpose[S]", "S . S", "S . T", "S . v", "v . S", "S . Normal[T]",
    "Normal[S] . T", "Det[S]", "Tr[S]", "Diagonal[S]", "Diagonal[S, 1]",
    "MatrixRank[S]", "Inverse[N[S] + 10 IdentityMatrix[Length[S]]]",
    "LinearSolve[N[S] + 10 IdentityMatrix[Length[S]], v]",
    "MatrixPower[S, 2]", "KroneckerProduct[S, T]", "ConjugateTranspose[S]",
    "LowerTriangularize[S]", "UpperTriangularize[S]", "S . S . v",
    "ArrayFlatten[{{S, T}, {T, S}}]", "Inner[Times, S, v, Plus]",
    "Eigenvalues[N[S + Transpose[S]]]",
]
VECTOR = [  # rank-1
    "S . T", "S . v", "Norm[S]", "Total[S]", "Accumulate[S]", "Reverse[S]",
    "Map[f, S]", "Select[S, # != 0 &]", "Outer[Times, S, T]",
    "KroneckerProduct[S, T]", "Riffle[S, T]", "DiagonalMatrix[S]",
    "Tally[S]", "Sort[Normal[S]]", "Ordering[Normal[S]]",
]

# Known, accepted differences (template -> reason). Keep this short and
# commented: every line is a semantics decision, not a bug to hide.
EXEMPT = {
}

# Heads whose cases only use machine operands (see gen_cases).
NUMERIC_ONLY = ("LinearSolve", "Inverse", "Eigenvalues")

KINDS = ["int", "real", "complex", "rational", "big", "symbolic", "mixed"]


def rand_value(rng, kind):
    if kind == "mixed":
        kind = rng.choice(KINDS[:-1])
    if kind == "int":
        return str(rng.randint(-9, 9))
    if kind == "real":
        return repr(round(rng.uniform(-5, 5), 3)).replace("e", "*^")
    if kind == "complex":
        return f"Complex[{round(rng.uniform(-3, 3), 2)}, {round(rng.uniform(-3, 3), 2)}]"
    if kind == "rational":
        return f"{rng.randint(-9, 9)}/{rng.randint(1, 7)}"
    if kind == "big":
        return str(rng.choice([-1, 1]) * rng.randint(10**20, 10**25))
    return rng.choice(["a", "b", "y", "a + 1", "Sqrt[2]", "Pi", "z^2"])


def rand_dims(rng, rank, square):
    if square and rank == 2:
        n = rng.randint(1, 5)
        return [n, n]
    return [rng.randint(1, 5) for _ in range(rank)]


def rand_array(rng, dims, kind, bg, density):
    total = 1
    for d in dims:
        total *= d
    nrules = max(0, int(total * density))
    rules = []
    for _ in range(nrules):
        pos = [rng.randint(1, d) for d in dims]
        if rng.random() < 0.1:   # negative indices count from the end
            j = rng.randrange(len(pos))
            pos[j] = pos[j] - dims[j] - 1
        val = bg if (rng.random() < 0.08) else rand_value(rng, kind)
        p = "{" + ", ".join(map(str, pos)) + "}" if len(dims) > 1 else str(pos[0])
        rules.append(f"{p} -> {val}")
    dims_s = "{" + ", ".join(map(str, dims)) + "}"
    return f"SparseArray[{{{', '.join(rules)}}}, {dims_s}, {bg}]"


def gen_cases(n, seed):
    rng = random.Random(seed)
    cases = []
    for i in range(n):
        rank = rng.choice([1, 1, 2, 2, 2, 3])
        family = rng.choice(["elem", "red", "struct", "lin"])
        square = family == "lin" and rank == 2
        dims = rand_dims(rng, rank, square)
        kind = rng.choice(KINDS)
        bg = rng.choice(["0", "0", "0", "0.", "1", "x"]) if rng.random() < 0.3 else "0"
        dens = rng.choice([0.0, 0.2, 0.4, 0.7, 1.0])
        s = rand_array(rng, dims, kind, bg, dens)
        t = rand_array(rng, dims, rng.choice(KINDS), "0" if rng.random() < .8 else bg, rng.choice([0.2, 0.5]))
        v = "{" + ", ".join(rand_value(rng, "int") for _ in range(dims[-1])) + "}"
        if family == "elem":
            op = rng.choice(ELEMENTWISE)
        elif family == "red":
            op = rng.choice(REDUCTION)
        elif family == "struct":
            op = rng.choice(STRUCTURAL)
        else:
            op = rng.choice(MATRIX if rank == 2 else VECTOR if rank == 1 else STRUCTURAL)
            if any(h in op for h in NUMERIC_ONLY) and kind not in ("int", "real", "rational"):
                # a numeric solver over symbolic / 25-digit entries is slow on
                # both sides alike; regenerate the operands numerically
                s = rand_array(rng, dims, "real", "0", dens)
                t = rand_array(rng, dims, "int", "0", 0.3)
        cases.append((i, op, s, t, v))
    return cases


PRELUDE = r"""
(* Normal of every array; a machine zero is the exact 0 (the two are conflated
   on purpose: a sparse result's background positions are exact even where the
   dense computation produced 0.) *)
deN[e_] := (e /. a_SparseArray :> Normal[a]) /. r_Real /; r == 0 :> 0;
numClose[a_, b_] := Module[{d = Quiet[Check[N[deN[a] - deN[b]], $Failed]]},
   d =!= $Failed && VectorQ[Flatten[{d}], NumericQ] &&
   Quiet[Check[Max[Abs[Flatten[{d}]]] < 1.*^-8 (1 + Max[Abs[N[Flatten[{deN[b]}]]]]), False]]];
same[a_, b_] := deN[a] === deN[b] || numClose[a, b] ||
   Quiet[Check[Simplify[deN[a] - deN[b]] === 0 || Union[Flatten[{Simplify[deN[a] - deN[b]]}]] === {0}, False]];
C[id_, op_, s_, t_, v_] := Module[{S = s, T = t, r1, r2},
   r1 = Quiet[TimeConstrained[op[S, T, v], 10, $TIMEOUT]];
   r2 = Quiet[TimeConstrained[op[Normal[S], Normal[T], v], 10, $TIMEOUT]];
   If[r1 === $TIMEOUT || r2 === $TIMEOUT, Print["TIMEOUT\t", id],
     If[! same[r1, r2],
       Print["FAIL\t", id, "\t", ToString[deN[r1], InputForm], "\t", ToString[deN[r2], InputForm]]]]];
"""


def case_line(c):
    i, op, s, t, v = c
    # the op as a pure function of (S, T, v); x and f stay symbolic
    return f'C[{i}, Function[{{S, T, v}}, {op}], {s}, {t}, {v}];'


def run_chunk(cases, valgrind=False, timeout=900):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write(PRELUDE)
        for c in cases:
            fh.write(case_line(c) + "\n")
        path = fh.name
    cmd = [BINARY, "-file", path]
    if valgrind:
        cmd = ["valgrind", "--error-exitcode=99", "--leak-check=full",
               "--errors-for-leak-kinds=definite", "-q"] + cmd
    try:
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        rc, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        rc, out, err = -1, "", "chunk timeout"
    os.unlink(path)
    return rc, out, err


def run_all(cases, chunk, valgrind):
    results = []   # (id, status, detail)
    stack = [cases[i:i + chunk] for i in range(0, len(cases), chunk)]
    while stack:
        part = stack.pop()
        rc, out, err = run_chunk(part, valgrind)
        for line in out.splitlines():
            f = line.split("\t")
            if f[0] in ("FAIL", "TIMEOUT"):
                results.append((int(f[1]), f[0], f[2:] if len(f) > 2 else []))
        if rc != 0:
            if len(part) == 1:
                kind = "VALGRIND" if rc == 99 else "CRASH"
                results.append((part[0][0], kind, [f"rc={rc}", err.strip()[-400:]]))
            else:
                done = {r[0] for r in results}
                rest = [c for c in part if c[0] not in done]
                mid = len(rest) // 2
                stack += [rest[:mid], rest[mid:]] if mid else [rest]
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=250)
    ap.add_argument("--valgrind", action="store_true")
    ap.add_argument("--keep", metavar="FILE")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args()

    cases = gen_cases(a.n, a.seed)
    if a.keep:
        with open(a.keep, "w") as fh:
            fh.write(PRELUDE + "\n".join(case_line(c) for c in cases) + "\n")
    byid = {c[0]: c for c in cases}
    results = run_all(cases, a.chunk if not a.valgrind else max(1, a.chunk // 10), a.valgrind)
    bad = []
    for cid, status, detail in sorted(results):
        op = byid[cid][1]
        if op in EXEMPT:
            continue
        bad.append((cid, status, detail))
    from collections import Counter
    ops = Counter(byid[cid][1] for cid, _, _ in bad)
    for cid, status, detail in bad:
        _, op, s, t, v = byid[cid]
        print(f"{status} #{cid}: {op}")
        if a.verbose or status != "FAIL":
            print(f"    S = {s}\n    T = {t}\n    v = {v}")
        for d in detail:
            print(f"    {d[:300]}")
    print(f"\n{len(cases)} cases, {len(bad)} failing; by op:")
    for op, k in ops.most_common():
        print(f"  {k:5d}  {op}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
