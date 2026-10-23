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
