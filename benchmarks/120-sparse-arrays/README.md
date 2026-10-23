# 120 — SparseArray against scipy.sparse

The sparse-array substrate under sparse linear algebra and graph workloads:
building a 10^5 × 10^5 matrix from 10^6 coordinates, element-wise arithmetic
and functions that preserve sparsity, element / row / entry access, and a
2000 × 2000 dense ⇄ sparse round trip. Positions come from a Lehmer LCG spelled
identically in both languages (`B` shares half of `A`'s positions); no position
repeats, so first-rule-wins and scipy's sum-duplicates build the same matrix.
Checks are `Round[1000 Σ explicit values]` of each result (Integers, compared
exactly).

| case | Mathilda | scipy |
|---|---|---|
| build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows) | `SparseArray[pos -> vals, {n, n}]` | `csr_array((v, (i, j)))` |
| arith: A + B / A B / 2.5 A / A^2 / A + A | `Plus` / `Times` / `Power` | `+`, `.multiply`, `*`, `.power` |
| fn: Sin[A] / Sqrt[A] | Listable heads over the values | `.sin()`, `.sqrt()` |
| access: row 500 x100 / element (500, 600) x1000 | `A[[500]]`, `A[[500, 600]]` | `A[[499], :]`, `A[499, 599]` |
| access: explicit positions and values | `A["NonzeroPositions"]`, `A["NonzeroValues"]` | `tocoo()` |
| dense: 2000 x 2000 dense -> sparse / sparse -> dense | `SparseArray[dense]`, `Normal` | `csr_array(dense)`, `toarray()` |

## Findings

Measured (ms, minimum of 3, i9-9880H), Mathilda / scipy 1.18:

| case | Mathilda | scipy | ratio |
|---|---|---|---|
| build: 10^6 entries, 10^5 x 10^5 | 29.6 | 46.0 | 0.64× |
| arith: A + B | 7.0 | 11.5 | 0.61× |
| arith: A B (element-wise) | 5.2 | 10.3 | 0.50× |
| arith: 2.5 A | 1.7 | 4.0 | 0.41× |
| arith: A^2 (element-wise) | 1.7 | 4.0 | 0.43× |
| arith: A + A | 1.7 | 6.8 | 0.25× |
| fn: Sin[A] | 2.8 | 15.2 | 0.18× |
| fn: Sqrt[A] | 2.0 | 4.1 | 0.49× |
| access: row 500 x100 | 5.5 | 7.3 | 0.76× |
| access: element (500, 600) x1000 | 0.33 | 16.7 | 0.02× |
| access: explicit positions and values | 7.1 | 13.7 | 0.52× |
| dense: 2000 x 2000 dense -> sparse | 7.3 | 23.6 | 0.31× |
| dense: 2000 x 2000 sparse -> dense | 5.4 | 5.2 | 1.03× |

Every case is at parity or ahead. What makes the difference:

1. **Shared index structure.** `2.5 A`, `A^2`, `Sin[A]`, `Sqrt[A]` and `A + A`
   reuse `A`'s row pointers and column indices (reference-counted) and allocate
   only the value vector; scipy copies its index arrays for each result.
2. **The same kernel as the packed surface.** `Sin`, `Sqrt`, `Exp`, ... run the
   head's registered ND kernel over the value vector, threaded where the
   kernel is. The cheap real powers (`x^2`, `x^0.5`, `1/x`) are hoisted out of the
   per-element complex-power path, which made `Sqrt` 3× faster on packed arrays
   too.
3. **Construction** is one stable LSD radix sort on a composite
   `(row << kbits) | column` key with the input index packed into the low bits,
   so each pass streams a single 8-byte array. Three 12-bit passes measured
   about 30% faster than four 9-bit ones. The packed `positions -> values` rule
   reaches the constructor still packed: the evaluator evaluates the sides of
   SparseArray's Rule argument directly, because Rule itself is not
   packed-aware.
4. **sparse -> dense** is at parity. It is bounded by writing the 32 MB dense
   buffer; a zero background starts from `calloc` pages, as `np.zeros` does.

Written as `r = {positions, values}`, the explicit-entries case would cost
~700 ms: a List materialises its packed elements (the packing layer's
no-nesting invariant), so the case assigns the two buffers separately, as the
Python returns a tuple of arrays.
