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
    "Total[S, {1}]", "Total[S, Infinity]", "Total[S, {-1}]", "Total[S, {1, 2}]",
    "Total[S, 0]", "Total[N[S]]", "Total[N[S], {2}]", "Mean[N[S]]", "Variance[S]",
    "StandardDeviation[N[S]]", "Max[S, x]", "Min[S, -1, T]", "MinMax[S]",
    "Accumulate[N[S]]", "S == T", "S != T", "Equal[S, T, S]", "S == Normal[T]",
    "Unequal[S, Normal[S]]", "S == N[S]",
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
    # native structural kernels (src/sparse/sparse_struct.c)
    "Flatten[S, 1]", "Take[S, {1, -1, 2}]",
    "Take[S, {-1, 1, -2}]", "Take[S, 0]", "Drop[S, {1, -1, 2}]", "Drop[S, {2}]",
    "Drop[S, -1]", "Join[S, T, 1]", "Partition[S, 2]", "Partition[S, 2, 1]",
    "ArrayReshape[S, {2, 3}]", "ArrayReshape[S, {3, 3}, x]", "ArrayPad[S, {1, 2}]",
    "ArrayPad[S, 1, x]", "PadLeft[S, Dimensions[S] + 2, 7]", "PadRight[S, Dimensions[S] + 1, x]",
    "Transpose[S, Range[ArrayDepth[S]]]", "Flatten[S, {Reverse[Range[ArrayDepth[S]]]}]",
    "Append[S, S[[-1]]]", "Prepend[S, x]",
]
MATRIX = [  # rank-2 square
    "Transpose[S]", "S . S", "S . T", "S . v", "v . S", "S . Normal[T]",
    "Normal[S] . T", "Det[S]", "Tr[S]", "Diagonal[S]", "Diagonal[S, 1]",
    "MatrixRank[S]", "Inverse[N[S] + 10 IdentityMatrix[Length[S]]]",
    "LinearSolve[N[S] + 10 IdentityMatrix[Length[S]], v]",
    "MatrixPower[S, 2]", "KroneckerProduct[S, T]", "ConjugateTranspose[S]",
    "LowerTriangularize[S]", "UpperTriangularize[S]", "S . S . v",
    "ArrayFlatten[{{S, T}, {T, S}}]", "Inner[Times, S, v, Plus]",
    "MatrixPower[S, 3]", "MatrixPower[S, 2, v]", "MatrixPower[S, 1]",
    "Dot[S, T, S]", "v . S . T", "S . T . v", "N[S] . T", "S . N[T]",
    "N[S] . N[T]", "S . N[v]", "N[v] . S", "S . Transpose[Normal[T]]",
    "Normal[T] . S . v", "Inner[Times, S, T, Plus]", "(S - x) . v",
    "Eigenvalues[N[S + Transpose[S]]]",
    "Norm[S]", "Norm[N[S]]", "Norm[S, 1]", "Norm[S, Infinity]", "Norm[S, \"Frobenius\"]",
    "Norm[N[S], \"Frobenius\"]", "Tr[S]", "Tr[S, List]", "Accumulate[S]", "Mean[S]",
    "Variance[N[S]]",
    # native structural kernels (src/sparse/sparse_struct.c)
    "Transpose[S, {2, 1}]", "Reverse[S, 2]", "Reverse[S, {1, 2}]", "Diagonal[S, -1]", "UpperTriangularize[S, 1]",
    "LowerTriangularize[S, -1]", "Join[S, T, 2]", "KroneckerProduct[S, Normal[T]]",
    "KroneckerProduct[S, v]", "ArrayFlatten[{{S, 0}, {0, T}}]", "ArrayFlatten[{{S, Normal[T]}}]",
    "Partition[S, 2]", "DiagonalMatrix[Diagonal[S]]", "Take[S, -1, {1, -1, 2}]",
    "Drop[S, {1}, {-1}]", "Append[S, v]", "Map[f, S]",
]
VECTOR = [  # rank-1
    "S . T", "S . v", "Norm[S]", "Total[S]", "Accumulate[S]", "Reverse[S]",
    "v . S", "N[S] . T", "S . N[v]", "S . Normal[T]", "Normal[S] . T",
    "S . Outer[Times, v, v]", "Outer[Times, v, v] . S",
    "Map[f, S]", "Select[S, # != 0 &]", "Outer[Times, S, T]",
    "KroneckerProduct[S, T]", "Riffle[S, T]", "DiagonalMatrix[S]",
    "Tally[S]", "Sort[Normal[S]]", "Ordering[Normal[S]]",
    "Norm[S, 1]", "Norm[S, Infinity]", "Norm[S, 3]", "Norm[N[S], 3]", "Norm[S, 2.5]",
    "Norm[N[S]]", "Max[S]", "Mean[S]", "Variance[S]", "StandardDeviation[N[S]]",
    "Total[S]", "Tr[S]", "Accumulate[S]", "S == Normal[T]", "Norm[S, \"Frobenius\"]",
    # native structural kernels (src/sparse/sparse_struct.c)
    "DiagonalMatrix[S, 1]", "DiagonalMatrix[S, -1, Length[S] + 2]", "Partition[S, 2, 1]",
    "Map[Sin, S]", "Map[{#} &, S]", "Append[S, x]", "Prepend[S, 0]",
    "PadLeft[S, Length[S] + 2, x]", "KroneckerProduct[S, v]", "Take[S, {-1, 1, -2}]",
    "Join[S, T]", "Rest[S]", "Most[S]", "Last[S]",
]

# Functional heads: level specs below level 1, folds, ordering and set heads
# answered as on Normal, threading and nesting (src/sparse/sparse_func.c).
FUNCTIONAL = [
    "Map[f, S, {ArrayDepth[S]}]", "Map[f, S, {-1}]", "Map[f, S, {1, 2}]", "Map[f, S, All]",
    "Map[f, S, {0, 1}]", "Map[f, S, 2]", "Map[f, S, Infinity]", "Map[f, S, {-2}]",
    "MapIndexed[f, S]", "MapIndexed[f, S, {-1}]", "MapIndexed[f, S, {ArrayDepth[S]}]",
    "Apply[f, S, {1}]", "Apply[f, S, {0, 1}]", "f @@ S", "Apply[Plus, S, {1}]",
    "Reap[Scan[Sow, S, {-1}]]", "Reap[Scan[Sow, S]]", "Reap[Scan[Sow, S, {ArrayDepth[S]}]]",
    "MapAt[f, S, 1]", "MapAt[f, S, -1]", "Fold[Plus, 0, S]", "FoldList[Plus, S]",
    "Fold[f, x, S]", "Fold[Plus, S, {S, T}]", "Nest[Reverse, S, 3]", "NestList[f /@ # &, S, 2]",
    "FixedPoint[Clip[# - 1, {0, 3}] &, S, 4]", "Sort[S]", "Union[S]", "DeleteDuplicates[S]",
    "Split[S]", "Gather[S]", "Ordering[S]", "Counts[S]", "Tally[S]", "Reverse[Sort[S]]",
    "AnyTrue[S, # === 0 &]", "AllTrue[S, # =!= x &]", "NoneTrue[S, ListQ]",
    "Select[S, # =!= 0 &]", "Thread[f[S, T]]", "Thread[f[S]]", "Thread[{S, T}]",
    "Flatten[{S, T}]", "Flatten[{S, {x}}, 1]", "Outer[f, S, {1, 2}]",
    "Outer[Times, S, {1, 2}]", "Outer[f, {x}, S]", "S Normal[T]", "Normal[T] S",
    "S + Normal[T]", "x S Normal[T]", "Clip[S, {-1, 1}]", "Chop[N[S]]",
    "Composition[f, Reverse][S]", "Through[{Length, First}[S]]", "MapThread[f, {S, T}]",
    "S + First[S]", "S First[T]",
]
# Procedural programming over a SparseArray: iterators, loops, scoping and
# in-place mutation.
PROCEDURAL = [
    "Table[i, {i, S}]", "Table[f[i], {i, S}]", "Sum[i, {i, S}]", "Product[f[i], {i, S}]",
    "Module[{l = {}}, Do[AppendTo[l, i], {i, S}]; l]",
    "Module[{c = 0}, Do[If[i === 0, Continue[]]; c++, {i, S}]; c]",
    "Catch[Do[If[i =!= 0, Throw[i]], {i, S}]; none]",
    "Module[{c = 0}, Do[c++; If[c == 2, Break[]], {i, S}]; c]",
    "Module[{w = S}, AppendTo[w, First[S]]; w]", "Module[{w = S}, PrependTo[w, Last[S]]; w]",
    "Module[{w = S}, w[[1]] = x; w]", "Module[{w = S}, w[[-1]] = First[T]; w]",
    "Module[{w = S}, w += T; w]", "Module[{w = S}, w *= 2; w]",
    "Module[{w = S, i}, For[i = 1, i <= Length[w], i++, w[[i]] = f[w[[i]]]]; w]",
    "Module[{i = 1, acc = 0}, While[i <= Length[S], acc += Total[Flatten[{S[[i]]}]]; i++]; acc]",
    "With[{w = S}, Reverse[w]]", "Block[{w = S}, Length[w] + Total[Flatten[{w}]]]",
    "Reap[Do[Sow[i], {i, S}]]", "Module[{w = S}, w[[All]] = 0; w]",
]

# Heads (or call forms) with no dense builtin to compare against: the dense
# side of these templates is the stated formula instead (a function of S, T,
# v as Normal).
DENSE_REF = {
    "KroneckerProduct[S, T]":
        "If[ArrayDepth[S] == 1, Outer[Times, S, T], ArrayFlatten[Outer[Times, S, T]]]",
    "KroneckerProduct[S, Normal[T]]": "ArrayFlatten[Outer[Times, S, T]]",
    "KroneckerProduct[S, v]":
        "If[ArrayDepth[S] == 1, Outer[Times, S, v], ArrayFlatten[Outer[Times, S, {v}]]]",
    # the dense Transpose leaves a vector unevaluated; a vector's transpose
    # is itself
    "Transpose[S]": "If[ArrayDepth[S] == 1, S, Transpose[S]]",
    "Transpose[S, Range[ArrayDepth[S]]]": "S",
    "Flatten[S, {Reverse[Range[ArrayDepth[S]]]}]":
        "If[ArrayDepth[S] == 1, S, Flatten[Transpose[S, Reverse[Range[ArrayDepth[S]]]]]]",
    "UpperTriangularize[S]":
        "Table[If[j >= i, S[[i, j]], 0], {i, Length[S]}, {j, Length[S[[1]]]}]",
    "UpperTriangularize[S, 1]":
        "Table[If[j >= i + 1, S[[i, j]], 0], {i, Length[S]}, {j, Length[S[[1]]]}]",
    "LowerTriangularize[S]":
        "Table[If[j <= i, S[[i, j]], 0], {i, Length[S]}, {j, Length[S[[1]]]}]",
    "LowerTriangularize[S, -1]":
        "Table[If[j <= i - 1, S[[i, j]], 0], {i, Length[S]}, {j, Length[S[[1]]]}]",
}

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


FAMILIES = ["elem", "red", "struct", "lin", "func", "proc"]


def gen_cases(n, seed, families=FAMILIES):
    rng = random.Random(seed)
    cases = []
    for i in range(n):
        rank = rng.choice([1, 1, 2, 2, 2, 3])
        family = rng.choice(families)
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
        elif family == "func":
            op = rng.choice(FUNCTIONAL)
        elif family == "proc":
            op = rng.choice(PROCEDURAL)
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
(* Numerically equal, with values put in for the free symbols: a product over
   symbolic entries can carry an exact coefficient on one side and its machine
   value on the other (the dense sum also adds 0 * 1.5 -> 0. terms, the sparse
   one never forms them), and the two are the same number. *)
symSub = {a -> 0.3711, b -> 1.1353, y -> -0.7137, z -> 0.5311, x -> 0.2917};
numClose[p_, q_] := Module[{d = Quiet[Check[N[(deN[p] - deN[q]) /. symSub], $Failed]]},
   d =!= $Failed && VectorQ[Flatten[{d}], NumericQ] &&
   Quiet[Check[Max[Abs[Flatten[{d}]]] < 1.*^-8 (1 + Max[Abs[N[Flatten[{deN[q] /. symSub}]]]]), False]]];
same[a_, b_] := deN[a] === deN[b] || numClose[a, b] ||
   Quiet[Check[Simplify[deN[a] - deN[b]] === 0 || Union[Flatten[{Simplify[deN[a] - deN[b]]}]] === {0}, False]];
C[id_, op_, s_, t_, v_] := C[id, op, op, s, t, v];
C[id_, op_, ref_, s_, t_, v_] := Module[{S = s, T = t, r1, r2},
   r1 = Quiet[TimeConstrained[op[S, T, v], 10, $TIMEOUT]];
   r2 = Quiet[TimeConstrained[ref[Normal[S], Normal[T], v], 10, $TIMEOUT]];
   If[r1 === $TIMEOUT || r2 === $TIMEOUT, Print["TIMEOUT\t", id],
     If[! same[r1, r2],
       Print["FAIL\t", id, "\t", ToString[deN[r1], InputForm], "\t", ToString[deN[r2], InputForm]]]]];
"""


def case_line(c):
    i, op, s, t, v = c
    # the op as a pure function of (S, T, v); x and f stay symbolic
    if op in DENSE_REF:
        return (f'C[{i}, Function[{{S, T, v}}, {op}], Function[{{S, T, v}}, {DENSE_REF[op]}], '
                f'{s}, {t}, {v}];')
    return f'C[{i}, Function[{{S, T, v}}, {op}], {s}, {t}, {v}];'


def run_chunk(cases, valgrind=False, timeout=900):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write(PRELUDE)
        for c in cases:
            fh.write(case_line(c) + "\n")
        path = fh.name
    cmd = [BINARY, "-file", path]
    if valgrind:
        # Memory errors (invalid reads/writes, uninitialised use) fail the run
        # through the exit code. Leaks are compared against an empty script's
        # (the startup baseline), so they are measured, not counted as errors.
        cmd = ["valgrind", "--error-exitcode=99", "--leak-check=full",
               "--errors-for-leak-kinds=none"] + cmd
    try:
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        rc, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        rc, out, err = -1, "", "chunk timeout"
    os.unlink(path)
    if valgrind:
        lost = definitely_lost(err)
        base = leak_baseline()
        if lost is not None and base is not None and lost > base:
            out += f"LEAK\t{cases[0][0]}\tdefinitely lost {lost} bytes (startup baseline {base})\n"
    return rc, out, err


def definitely_lost(valgrind_err):
    import re
    m = re.search(r"definitely lost: ([0-9,]+) bytes", valgrind_err)
    return int(m.group(1).replace(",", "")) if m else 0


_LEAK_BASE = []


def leak_baseline():
    # what an empty script leaks (one-time initialisation), measured once
    if not _LEAK_BASE:
        rc, out, err = run_chunk_raw(PRELUDE, valgrind=True)
        _LEAK_BASE.append(definitely_lost(err))
    return _LEAK_BASE[0]


def run_chunk_raw(text, valgrind=False, timeout=900):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write(text)
        path = fh.name
    cmd = [BINARY, "-file", path]
    if valgrind:
        cmd = ["valgrind", "--error-exitcode=99", "--leak-check=full",
               "--errors-for-leak-kinds=none"] + cmd
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    os.unlink(path)
    return p.returncode, p.stdout, p.stderr


def run_all(cases, chunk, valgrind):
    results = []   # (id, status, detail)
    stack = [cases[i:i + chunk] for i in range(0, len(cases), chunk)]
    while stack:
        part = stack.pop()
        rc, out, err = run_chunk(part, valgrind)
        for line in out.splitlines():
            f = line.split("\t")
            if f[0] in ("FAIL", "TIMEOUT", "LEAK"):
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
    ap.add_argument("--families", default=",".join(FAMILIES),
                    help="comma-separated subset of " + ",".join(FAMILIES))
    a = ap.parse_args()

    cases = gen_cases(a.n, a.seed, a.families.split(","))
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
