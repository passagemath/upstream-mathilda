#!/usr/bin/env python3
"""Experiment 120 -- SparseArray against scipy.sparse.  See sparse_arrays.m.

Baseline: scipy.sparse CSR arrays (csr_array), the standard compressed sparse
row container; element-wise operations are scipy's own (`+`, `.multiply`,
`.power`, `.sin`, `.sqrt`), all of which preserve sparsity as SparseArray does.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, check, require

import numpy as np
import scipy.sparse as sp

require(["scipy.sparse:csr_array", "numpy"])

n, nnz = 100000, 1000000
kk = np.arange(1, nnz + 1, dtype=np.int64)
ii = (48271 * kk) % 2147483647 % n
jj = (69621 * kk) % 2147483647 % n
jj2 = (69621 * kk + 7 * (kk % 2)) % 2147483647 % n
vv = (kk % 97) / 10.0 + 0.1


def digest(s):
    return int(round(1000 * float(s.data.sum())))


res = {}

# ---- build ------------------------------------------------------------------

bench("build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)",
      lambda: res.__setitem__("A", sp.csr_array((vv, (ii, jj)), shape=(n, n))))
A = res["A"]
check("build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)", digest(A))
B = sp.csr_array((vv, (ii, jj2)), shape=(n, n))

# ---- arithmetic -------------------------------------------------------------

for label, fn in [("arith: A + B", lambda: A + B),
                  ("arith: A B (element-wise)", lambda: A.multiply(B)),
                  ("arith: 2.5 A", lambda: 2.5 * A),
                  ("arith: A^2 (element-wise)", lambda: A.power(2)),
                  ("arith: A + A", lambda: A + A),
                  ("fn: Sin[A]", lambda: A.sin()),
                  ("fn: Sqrt[A]", lambda: A.sqrt())]:
    bench(label, lambda fn=fn: res.__setitem__("r", fn()))
    check(label, digest(sp.csr_array(res["r"])))

# ---- access -----------------------------------------------------------------

bench("access: row 500 x100", lambda: [A[[499], :] for _ in range(100)])
check("access: row 500 x100", int(round(1000 * float(A[[499], :].sum()))))
bench("access: element (500, 600) x1000", lambda: [A[499, 599] for _ in range(1000)])


def explicit():
    c = A.tocoo()
    return np.stack([c.row + 1, c.col + 1], axis=1), c.data


bench("access: explicit positions and values", lambda: res.__setitem__("e", explicit()))
check("access: explicit positions and values", int(res["e"][0].sum()))

# ---- dense <-> sparse -------------------------------------------------------

i = np.arange(1, 2001).reshape(-1, 1)
j = np.arange(1, 2001).reshape(1, -1)
dense = np.where((i + 3 * j) % 97 == 0, 1.5 + (i * j) % 7, 0.0)
bench("dense: 2000 x 2000 dense -> sparse", lambda: res.__setitem__("D", sp.csr_array(dense)))
D = res["D"]
check("dense: 2000 x 2000 dense -> sparse", digest(D))
bench("dense: 2000 x 2000 sparse -> dense", lambda: res.__setitem__("dd", D.toarray()))
check("dense: 2000 x 2000 sparse -> dense", int(round(1000 * float(res["dd"].sum()))))

# ---- products ---------------------------------------------------------------

x = (np.arange(1, n + 1) % 13) / 7.0
bench("dot: A . x (SpMV) x100", lambda: [A @ x for _ in range(100)])
check("dot: A . x (SpMV) x100", int(round(1000 * float((A @ x).sum()))))
bench("dot: x . A x10", lambda: [x @ A for _ in range(10)])
check("dot: x . A x10", int(round(1000 * float((x @ A).sum()))))
bench("dot: A . A (SpGEMM)", lambda: res.__setitem__("r", A @ A))
check("dot: A . A (SpGEMM)", digest(res["r"]))
bench("struct: Transpose[A]", lambda: res.__setitem__("At", A.T.tocsr()))
At = res["At"]
check("struct: Transpose[A]", digest(At))
bench("dot: A . Transpose[A] (SpGEMM)", lambda: res.__setitem__("r", A @ At))
check("dot: A . Transpose[A] (SpGEMM)", digest(res["r"]))

# ---- reductions -------------------------------------------------------------

import scipy.sparse.linalg as spla

for label, fn, dig in [
        ("reduce: Total[A] (column sums)", lambda: A.sum(axis=0), lambda r: int(round(1000 * float(r.sum())))),
        ("reduce: Total[A, {2}] (row sums)", lambda: A.sum(axis=1), lambda r: int(round(1000 * float(r.sum())))),
        ("reduce: Total[A, 2]", lambda: A.sum(), lambda r: int(round(1000 * float(r)))),
        ("reduce: Max[A]", lambda: A.max(), lambda r: int(round(1000 * float(r)))),
        ("reduce: Norm[A, Frobenius]", lambda: spla.norm(A, "fro"), lambda r: int(round(1000 * float(r))))]:
    bench(label, lambda fn=fn: res.__setitem__("r", fn()))
    check(label, dig(res["r"]))

# ---- assembly ---------------------------------------------------------------

bench("struct: Join[A, B] (vstack)", lambda: res.__setitem__("r", sp.vstack([A, B], format="csr")))
check("struct: Join[A, B] (vstack)", digest(res["r"]))
bench("struct: ArrayFlatten 2x2 blocks", lambda: res.__setitem__("r", sp.block_array([[A, B], [B, A]], format="csr")))
check("struct: ArrayFlatten 2x2 blocks", digest(res["r"]))
bench("struct: IdentityMatrix[10^6, SparseArray]", lambda: res.__setitem__("r", sp.identity(10**6, format="csr")))
check("struct: IdentityMatrix[10^6, SparseArray]", int(round(float(res["r"].data.sum()))))

# ---- 2-D Poisson: 5-point Laplacian on a 300 x 300 grid ---------------------

m = 300


def build_l():
    t1 = sp.diags_array([-1.0, 2.0, -1.0], offsets=[-1, 0, 1], shape=(m, m), format="csr")
    i1 = sp.identity(m, format="csr")
    return (sp.kron(t1, i1) + sp.kron(i1, t1)).tocsr()


bench("poisson: build L (Band, KroneckerProduct)", lambda: res.__setitem__("L", build_l()))
L = res["L"]
check("poisson: build L (Band, KroneckerProduct)", digest(L))
bb = (np.arange(1, m * m + 1) % 7) - 3.0
bench("poisson: LinearSolve[L, b] direct", lambda: res.__setitem__("sol", spla.spsolve(L.tocsc(), bb)))
check("poisson: LinearSolve[L, b] direct", int(round(1000 * float(res["sol"].sum()))))
bench("poisson: LinearSolve[L, b] CG 1e-10", lambda: res.__setitem__("cg", spla.cg(L, bb, rtol=1e-10)[0]))
check("poisson: LinearSolve[L, b] CG 1e-10", int(round(float(res["cg"].sum()))))
bench("poisson: 6 largest eigenvalues", lambda: res.__setitem__("ev", spla.eigsh(L, k=6, which="LM", return_eigenvectors=False)))
check("poisson: 6 largest eigenvalues", int(round(1e6 * float(res["ev"].sum()))))
bench("poisson: 6 smallest eigenvalues", lambda: res.__setitem__("ev", spla.eigsh(L, k=6, sigma=0, which="LM", return_eigenvectors=False)))
check("poisson: 6 smallest eigenvalues", int(round(1e6 * float(res["ev"].sum()))))
