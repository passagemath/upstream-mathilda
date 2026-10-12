# Task: Implement `Hash` builtin (Expression hash)

Scope: `"Expression"` hash type only; string formats that need no new type.

## Implementation
- [x] `src/core.c`: add `hash_make_integer`, `hash_format`, `builtin_hash` (near `builtin_bytecount`)
- [x] `src/core.c`: ensure `<inttypes.h>` is included
- [x] `src/core.h`: declare `builtin_hash`
- [x] `src/core.c` `core_init()`: register `Hash` + `ATTR_PROTECTED`
- [x] `src/info.c` `info_init()`: `symtab_set_docstring("Hash", ...)`
- [x] `src/sym_names.h` + `src/sym_names.c`: add `SYM_Hash` (3 spots)
- [x] `src/version.h`: bump 0.362 -> 0.363 (NUMBER + STRING)

## Tests
- [x] `tests/test_hash.c`: invariant + direct-C cross-check cases (11 groups)
- [x] `tests/CMakeLists.txt`: register `hash_tests`

## Docs
- [x] `docs/spec/builtins/expression-information.md`: `## Hash` section
- [x] `docs/spec/changelog/2026-10-12.md`: `## Hash` entry (`; v0.363`)
- [x] `Mathilda_spec.md`: changelog row already present (no change)

## Verification
- [x] `make -j` clean build (no errors/warnings)
- [x] REPL smoke via -file script (all behaviours correct; `$VersionNumber = 0.363`)
- [x] `?Hash` docstring confirmed via pty
- [x] `hash_tests` pass (11/11) + `expr_tests` no regression
- [x] valgrind clean (leak profile identical to sibling; only macOS platform baseline)
- [x] `make check-messages`, `make check-c99` clean (exit 0)
- [x] refresh code-review graph

## Review

**Done.** `Hash` implemented for the `"Expression"` type across the no-new-type string
formats. Key design choices:
- Surfaces the existing `expr_hash` (64-bit FNV-1a) rather than adding new hashing code —
  the structural hash was already the system's Association/memo keying primitive.
- `uint64` result returned as `EXPR_INTEGER` when it fits int64, else a bigint via
  `mpz_import` (portable, no `unsigned long` width assumption) + `expr_bigint_normalize`.
- Unsupported types/formats + wrong arg counts stay unevaluated with no message — keeps
  `Quiet[]`/`Check[]` correct and `make check-messages` clean.
- `ATTR_PROTECTED` only (deliberately not `Listable`).
- Structural/non-numeric head, so naturally outside the packed/NDArray/Compile audit surfaces.

**Not done (deferred to a commit step, per CLAUDE.md release policy):** the commit itself and
the `v0.363` git tag — pending the user's go-ahead to commit.
