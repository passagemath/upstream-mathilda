---
source: src/loadmodule.c, src/loadmodule.h
---
**Algorithm.** `builtin_loadmodule` takes one `EXPR_STRING` relpath, calls
`mathilda_load_module`, and returns `True` or `False`. Any other argument shape
leaves the call unevaluated. `mathilda_load_module` returns `True` at once if
the relpath is already in the load-once table (`lm_already_loaded`). Otherwise it
resolves the path, runs the file with `mathilda_run_file`, and on a successful
open records the relpath (`lm_mark_loaded`), so a repeat call never re-reads it.
A relpath that resolves nowhere emits `LoadModule::nofile` through `mth_message`,
so `Quiet` and `Check` both see it, and returns `False`.

`mathilda_resolve_internal` is deliberately **independent of the working
directory**. It first tries the cached base directory (`lm_base`) from an
earlier hit, then searches in this order:

1. `$MATHILDA_HOME/<relpath>`, then `$MATHILDA_HOME/internal/<relpath>`. The
   parent form was added for issue #96, where a kernelspec pointing at the
   directory containing `internal/` silently missed late-loaded tables.
2. Relative to the running executable (`lm_exe_dir`, using
   `readlink("/proc/self/exe")` on Linux, `_NSGetExecutablePath` on macOS and
   `GetModuleFileNameA` on Windows): `<exe_dir>/src/internal/` and
   `<exe_dir>/../share/mathilda/internal/`.
3. Relative to the image that holds this code (`lm_self_dir`, using `dladdr` on
   POSIX and `GetModuleHandleExA` on Windows): `<lib_dir>/../share/mathilda/internal/`.
   An installed shared `libmathilda` loaded into Python or Jupyter therefore
   finds its modules without `$MATHILDA_HOME` (issue #88). This step is skipped
   when that directory is the executable's own.
4. The compile-time `MATHILDA_PREFIX/share/mathilda/internal/`.
5. A CWD ladder: `src/internal/` up to three levels up, for the test binaries.

`mathilda_run_file` is the file-reading core that `Get` also wraps. It slurps the
file and loops `parse_next_expression` / `evaluate`, stopping at end of file or at
the first unparseable expression, and keeps the last value. `LoadModule`
discards that value and reports only whether the file was opened.

**Data structures.** The bookkeeping uses fixed-size tables: `lm_loaded[256][256]`
holds the relpaths (past 256 modules new entries are silently not recorded), and
`lm_base` caches the base directory. The file is read into one transient buffer.
On Linux `_GNU_SOURCE` is defined before any include, so `readlink` and `dladdr`
are visible under `-std=c99` (SPEC.md §10).

**Complexity / limits.** A first load is `O(file size)` to read and evaluate. A
repeat is a linear scan of the table, which is effectively `O(1)` for the small
fixed module set. Each module is loaded **at most once**, keyed by the literal
relpath string. That is what lets callers load lazily without re-registering
rules: `core.c` loads `simp/FullSimplify.m` at startup, `FullSimplify` loads its
`simp/transforms/*.m` families, `init.m` loads `distributions.m`, and
`Integrate` loads `PolyLogLadder.m`, `mixed/ParallelMixed.m` and
`ParallelMixedSpecial.m`. A file read with `Get` is not recorded. A syntax error
part-way through a module still counts as a successful load. Attributes:
`ATTR_PROTECTED`.
