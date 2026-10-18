### Worked examples

```mathematica
In[1]:= LoadModule["simp/FullSimplify.m"]  (* already loaded by the kernel at startup, so True without re-reading the file *)
```

```mathematica
In[1]:= LoadModule["simp/transforms/gamma.m"]  (* one of the identity libraries FullSimplify loads on demand *)
```

```mathematica
In[1]:= LoadModule["simp/transforms/gamma.m"]; LoadModule["simp/transforms/gamma.m"]  (* load-once: a repeat is a table hit, not a re-read *)
```

```mathematica
In[1]:= Quiet[LoadModule["no/such/module.m"]]  (* not found: False, with the LoadModule::nofile message silenced *)
```

```mathematica
In[1]:= Check[LoadModule["no/such/module.m"], "failed"]  (* the message routes through the funnel, so Check sees it *)
```

```mathematica
In[1]:= If[TrueQ[LoadModule["no/such/module.m"]], "loaded", "fallback"]  (* the guard idiom internal modules use for an optional dependency *)
```

```mathematica
In[1]:= LoadModule[42]  (* a non-string argument is left unevaluated *)
```

### Notes

- `LoadModule["relpath"]` loads a Mathilda `.m` source module, where `relpath` is relative to `src/internal` (or the installed `share/mathilda/internal`). It returns `True` if the module was found and loaded, or had already been loaded, and `False` otherwise.
- Each module is loaded at most once per session, keyed by the exact `relpath` string. This makes it safe to call from inside other modules and from lazy loaders: `FullSimplify` loads its `simp/transforms/*.m` families on demand, and `Integrate` loads `PolyLogLadder.m` and `mixed/ParallelMixed.m` this way. A module first loaded with `Get` is not in the table, so `LoadModule` would read it again.
- Resolution does not depend on the working directory. The search order is: `$MATHILDA_HOME/` (or its `internal/` subdirectory, since v0.449); the running executable's `src/internal/` and `../share/mathilda/internal/`; the directory of the image that holds the code, so an installed shared `libmathilda` inside Python or Jupyter finds its modules; a compile-time `PREFIX`; and finally a CWD ladder. The winning base directory is cached.
- A missing file emits `LoadModule::nofile` through `mth_message`, so both `Quiet` and `Check` work. A file that opens but stops at a syntax error still counts as loaded (`True`), because reading simply ends at the first unparseable expression, as with `Get`.
- Attributes: `Protected`. The value of the module's last expression is discarded. Use `Get` to get it back.
- See also: `Get`, `FullSimplify`.
