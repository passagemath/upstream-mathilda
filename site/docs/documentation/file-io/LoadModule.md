# LoadModule

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LoadModule["relpath"]`**

loads the internal Mathilda source module at relpath (relative to src/internal), resolving the location independently of the current working directory. Each module is loaded at most once. Returns True if the module was located and loaded, False otherwise.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (3)

Load an internal module

```mathematica
In[1]:= LoadModule["simp/FullSimplify.m"]
Out[1]= True
```

Already loaded: still True, not re-read

```mathematica
In[2]:= LoadModule["simp/FullSimplify.m"]
Out[2]= True
```

Not found -> False

```mathematica
In[3]:= Quiet[LoadModule["no/such/module.m"]]
Out[3]= False
```

### Applications (7)

Already loaded by the kernel at startup, so True without re-reading the file

```mathematica
In[4]:= LoadModule["simp/FullSimplify.m"]
Out[4]= True
```

One of the identity libraries FullSimplify loads on demand

```mathematica
In[5]:= LoadModule["simp/transforms/gamma.m"]
Out[5]= True
```

Load-once: a repeat is a table hit, not a re-read

```mathematica
In[6]:= LoadModule["simp/transforms/gamma.m"]; LoadModule["simp/transforms/gamma.m"]
Out[6]= True
```

Not found: False, with the LoadModule::nofile message silenced

```mathematica
In[7]:= Quiet[LoadModule["no/such/module.m"]]
Out[7]= False
```

The message routes through the funnel, so Check sees it

```mathematica
In[8]:= Check[LoadModule["no/such/module.m"], "failed"]
Out[8]= "failed"
```

The guard idiom internal modules use for an optional dependency

```mathematica
In[9]:= If[TrueQ[LoadModule["no/such/module.m"]], "loaded", "fallback"]
Out[9]= "fallback"
```

A non-string argument is left unevaluated

```mathematica
In[10]:= LoadModule[42]
Out[10]= LoadModule[42]
```

## Algorithm

loadmodule.c -- runtime loading of Mathilda (.m) source modules.

See loadmodule.h for the contract. The file-reading loop here is the same one Get[] uses (builtin_get in readwrite.c is now a thin wrapper over mathilda_run_file); mathilda_load_module adds working-directory-independent path resolution and load-once bookkeeping on top.

Path resolution (mathilda_resolve_internal) is deliberately CWD-independent: a relocated or installed binary must still find its bundled src/internal tree — and so must an installed shared libmathilda loaded into some other program (a Python extension, a Jupyter kernel), which is located through dladdr. On glibc, readlink("/proc/self/exe") needs the POSIX feature-test macro and dladdr/Dl_info need _GNU_SOURCE to be visible under -std=c99; define it before any system header.

## Implementation notes

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

- `Protected`. Returns `True` if the module was located and loaded (or had already
  been loaded), `False` otherwise.
- Resolution is independent of the current working directory and tries, in
  order: `$MATHILDA_HOME/<relpath>` (then `$MATHILDA_HOME/internal/<relpath>`, so
  pointing it at the directory *containing* `internal/` also works);
  `<exe_dir>/src/internal/<relpath>` and
  `<exe_dir>/../share/mathilda/internal/<relpath>` (relative to the running
  binary, so a relocated or installed executable still finds its modules);
  `$(PREFIX)/share/mathilda/internal/<relpath>` when built with a compile-time
  `MATHILDA_PREFIX`; and finally a CWD ladder (`src/internal/`,
  `../src/internal/`, `../../src/internal/`, `../../../src/internal/`). This
  works from the REPL (run at the repo root or from anywhere with the binary in
  place), from the test binaries (run from `tests/build/`), and from a binary
  copied to a `bin` directory with `MATHILDA_HOME` pointing at `src/internal`.
  The winning base directory is cached after the first successful lookup.
- Each module is loaded **at most once**, so repeated calls — and the lazy
  per-family loading used by [`FullSimplify`](../simplification/index.md) —
  never re-register rules.
- Generalises the bespoke fallback previously hard-coded for the CRC integral
  tables; `Get` (above) shares its file-reading core.

**Attributes:** `Protected`.

## References

**See also:** [Get](../../file-io/Get/)

- Source: [`src/loadmodule.c`](https://github.com/stblake/mathilda/blob/main/src/loadmodule.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)

## Notes & additional examples

### Notes

- `LoadModule["relpath"]` loads a Mathilda `.m` source module, where `relpath` is relative to `src/internal` (or the installed `share/mathilda/internal`). It returns `True` if the module was found and loaded, or had already been loaded, and `False` otherwise.
- Each module is loaded at most once per session, keyed by the exact `relpath` string. This makes it safe to call from inside other modules and from lazy loaders: `FullSimplify` loads its `simp/transforms/*.m` families on demand, and `Integrate` loads `PolyLogLadder.m` and `mixed/ParallelMixed.m` this way. A module first loaded with `Get` is not in the table, so `LoadModule` would read it again.
- Resolution does not depend on the working directory. The search order is: `$MATHILDA_HOME/` (or its `internal/` subdirectory, since v0.449); the running executable's `src/internal/` and `../share/mathilda/internal/`; the directory of the image that holds the code, so an installed shared `libmathilda` inside Python or Jupyter finds its modules; a compile-time `PREFIX`; and finally a CWD ladder. The winning base directory is cached.
- A missing file emits `LoadModule::nofile` through `mth_message`, so both `Quiet` and `Check` work. A file that opens but stops at a syntax error still counts as loaded (`True`), because reading simply ends at the first unparseable expression, as with `Get`.
- Attributes: `Protected`. The value of the module's last expression is discarded. Use `Get` to get it back.
- See also: `Get`, `FullSimplify`.
