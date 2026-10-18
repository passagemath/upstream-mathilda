---
source: src/strings/stringpartition.c
---
**Algorithm.** `builtin_stringpartition` takes 2–3 arguments and gives `::argtu`/`::argt` otherwise. The block length `n` is an integer, in which case every block is exactly `n` characters, or `UpTo[n]`, which clips each block to what remains. The offset `d` defaults to `n`. Block starts run over `0, d, 2d, …` in characters. A fixed `n` stops when a full block no longer fits. `UpTo[n]` takes every start before the end, so with `d < n` the tail repeats shorter and shorter blocks, as in Mathematica. With an explicit `d`, `n = 0` is allowed, and every block is then `""`.

**Data structures.** A `U8View` over the subject, so lengths and offsets count Unicode code points. Each block is created by `expr_new_string_n` straight from the byte span. The block array is bounded by `len/d + 1`.

**Complexity / limits.** `O(output)`. As in Mathematica, the first argument must be a single string. A list of strings is not threaded (a v0.456 change) and gives `StringPartition::string`. A bad block length gives `::intpm`, or `::intnm` once an offset is present. A non-positive or non-integer offset gives `::intp`. Each of these leaves the call unevaluated, and `UpTo[negative]` stays unevaluated without a message.
