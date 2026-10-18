---
source: src/strings/stringtakedrop.c, src/strings/stringtake.c, src/strings/stringdrop.c
---
**Algorithm.** `builtin_stringtakedrop` requires exactly two arguments (else `::argr` via `strings_argr_error`). The subject must be a `String` or a `List` of strings (else `::strse`), and the spec is parsed by `strings_seq_parse` into one `StrSeq` — `n`, `-n`, `{n}`, `{m, n}`, `{m, n, s}`, `UpTo[n]`, `All`, `None` (else `::seqs`). Before producing anything, the spec is resolved against *every* subject with `strings_chars_init` and `strings_seq_resolve`, so an out-of-range position reports `StringTakeDrop::take` rather than the underlying kernels' own messages. Each pair is then `{strings_take_seq(s, q), strings_drop_seq(s, q)}` — exactly the scalar kernels of `StringTake` and `StringDrop` — so the two halves always partition the string's characters.

**Data structures.** `StrSeq` (a kind tag plus `m`, `n`, step `s`) and a `StrChars` view that indexes character positions over the UTF-8 bytes, built only for the prefix the spec needs (`strings_seq_prefix`).

**Complexity / limits.** `O(B)` per subject for the take and the drop. Only a single sequence spec is accepted (no list of specs), matching Mathematica.
