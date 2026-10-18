# Compile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Compile[{x, ...}, expr]`**

or Compile\[{{x, \_Real}, ...}, expr\] builds a CompiledFunction that evaluates expr over machine numbers (types \_Real, \_Integer, \_Complex; default \_Real), falling back to the interpreter for symbolic arguments or non-compilable bodies. With RuntimeAttributes -\> Listable the object threads over List arguments; the default is RuntimeAttributes -\> {}. RuntimeOptions -\> {"CatchMachineIntegerOverflow" -\> False} (or the shorthand RuntimeOptions -\> "Speed") lets machine-integer arithmetic wrap instead of falling back to the interpreter, which is faster and gives a different answer from the interpreter once a result leaves the machine-integer range; the default True never does. WorkingPrecision -\> n compiles real/complex arithmetic in MPFR at n decimal digits (one fixed precision for the whole function), for the straight-line arithmetic + elementary-function subset; MachinePrecision (the default) keeps the machine path unchanged. "BigIntegers" -\> True makes integer arithmetic exact (GMP) instead of int64.

## Examples (22)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= f = Compile[{{x, _Real}}, x^2 + 1]
Out[1]= CompiledFunction[{x}, x^2 + 1]

In[2]:= f[3.0]
Out[2]= 10.0
```

Symbolic argument -> interpreter fallback

```mathematica
In[3]:= f[a]
Out[3]= 1 + a^2
```

```mathematica
In[4]:= g = Compile[{{n, _Integer}}, Module[{s = 0.}, Do[s = s + 1/i^2, {i, 1, n}]; s]]; g[100]
Out[4]= 1.63498

In[5]:= Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I]
Out[5]= -3.0 + 4.0*I

In[6]:= Compile[{{m, _Real, 2}}, m[[All, 1]]][{{1., 2.}, {3., 4.}}]
Out[6]= {1.0, 3.0}
```

A 5-point stencil: read an argument grid, write a local copy

```mathematica
In[7]:= Compile[{{a, _Real, 2}}, Module[{n = Length[a], b = a}, Do[b[[i, j]] = (a[[i - 1, j]] + a[[i + 1, j]] + a[[i, j - 1]] + a[[i, j + 1]])/4, {i, 2, n - 1}, {j, 2, n - 1}]; b]][Table[1.0 (10 i + j), {i, 1, 3}, {j, 1, 3}]]
Out[7]= {{11.0, 12.0, 13.0}, {21.0, 22.0, 23.0}, {31.0, 32.0, 33.0}}
```

### Scope (1)

```mathematica
In[8]:= Compile[{{n, _Integer}}, n^20, "BigIntegers" -> True][99]
Out[8]= 8179069375972308708891986605443361898001
```

### Options (4)

RuntimeAttributes -> Listable: the object threads over lists

```mathematica
In[9]:= h = Compile[{{x, _Real}}, If[x > 0, 1., -1.], RuntimeAttributes -> Listable]; h[{1., -2., 3.}]
Out[9]= {1.0, -1.0, 1.0}
```

A rank-1 parameter consumes one level, so this maps over the rows

```mathematica
In[10]:= Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable][ {{1., 2.}, {3., 4.}}]
Out[10]= {3.0, 7.0}
```

```mathematica
In[11]:= f = Compile[{{x, _Real}}, Sin[x] Cos[x] + x^3, WorkingPrecision -> 40]; f[N[7/5, 40]]
Out[11]= 2.9114940750779524597719268763562110530148

In[12]:= Compile[{{z, _Complex}}, Exp[z] + z^2, WorkingPrecision -> 45][N[1/2 + I/3, 45]]
Out[12]= 1.696859506173835946311071541529964410492750645 + 0.8727861896014286091867555415942467494395774456*I
```

### Applications (10)

A machine-number function applied to a Real

```mathematica
In[13]:= f = Compile[{{x, _Real}}, x^2 + 1]; f[3.0]
Out[13]= 10.0
```

A symbolic argument transparently falls back to the interpreter

```mathematica
In[14]:= f[a]
Out[14]= 1 + a^2
```

Integer arithmetic stays exact past the int64 range: overflow defers to the interpreter

```mathematica
In[15]:= Compile[{{n, _Integer}}, n^3][3000000]
Out[15]= 27000000000000000000
```

Opting out of overflow checking keeps the wrapped int64, 27*10^18 mod 2^64

```mathematica
In[16]:= Compile[{{n, _Integer}}, n^3, RuntimeOptions -> "Speed"][3000000]
Out[16]= 8553255926290448384
```

Special functions with a complex machine kernel lower directly

```mathematica
In[17]:= Compile[{{z, _Complex}}, Gamma[z] + Erf[z]][0.5 + 1.0 I]
Out[17]= 1.50554 + 0.599433*I
```

Whole-array Chop and Clip over a rank-1 parameter

```mathematica
In[18]:= Compile[{{v, _Real, 1}}, Clip[Chop[v - 1., 0.000001], {-1., 1.}]][{0.5, 1.0000000001, 3., -4.}]
Out[18]= {-0.5, 0.0, 1.0, -1.0}
```

A rank-1 parameter consumes one level, so a Listable object maps over the rows

```mathematica
In[19]:= Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable][{{1., 2.}, {3., 4.}}]
Out[19]= {3.0, 7.0}
```

One unlowered head sends the whole body to the interpreter; the diagnostics name it

```mathematica
In[20]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + myf[x]]
Out[20]= {"Compiled" -> False, "Reason" -> "no machine lowering for this head at these argument types", "Subexpression" -> "myf[x]"}
```

Opt-in MPFR arithmetic at a fixed 40 digits

```mathematica
In[21]:= Compile[{{x, _Real}}, Sin[x] Cos[x] + x^3, WorkingPrecision -> 40][N[7/5, 40]]
Out[21]= 2.9114940750779524597719268763562110530148
```

Opt-in exact GMP integers inside the compiled program

```mathematica
In[22]:= Compile[{{n, _Integer}}, n^20, "BigIntegers" -> True][99]
Out[22]= 8179069375972308708891986605443361898001
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NDSolve Van der Pol mu=10 | 0.441 s | 0.566 s | 3.32 s |
| NDSolve harmonic oscillator | 0.216 s | 0.214 s | 7.49 s |
| NDSolve long horizon, t to 200 | 0.111 s | 0.152 s | 2.89 s |
| NDSolve y'=-y on [0,10] | 0.065 s | 0.156 s | 2.41 s |
| NDSolve y'=y^2 t nonlinear | 0.034 s | 0.135 s | 0.3 s |

## Implementation notes

**Algorithm.** `Compile` is `HoldAll, Protected`, so the body is compiled in its
raw unevaluated form and the argument symbols stay local (a global value for one
does not leak in). `builtin_compile` first seeds `RuntimeAttributes` and
`RuntimeOptions` from the registered defaults, so `SetOptions[Compile, …]` takes
effect. It then strips trailing options right to left (`RuntimeAttributes`,
`RuntimeOptions`, `WorkingPrecision`, `"BigIntegers"`), and the rightmost setting
wins. An option it has no meaning for leaves `Compile[…]` unevaluated rather than
being silently ignored.

`compiled_function_new` parses the argument spec. `compile_infer.c` then
type-infers the body over the lattice Boolean / Integer / Real / Complex / array
(plus BigInteger / MPFRReal / MPFRComplex in managed mode). The per-category
emitters (`compile_emit_arith/int/logic/mathfn/ctrl/iter/array.c`,
`compile_assoc.c`) lower it to typed bytecode. Heads with an existing `NDArray`
entry point (structural, reductions, linalg, Fourier) are delegated through the
dispatch tables in `compile_ndtables.c`, so for those heads the compiled subset is
the interpreted one by construction. Element-wise special functions lower to
`KERN_*` opcodes over the machine kernel registry in `ndkernels.c`, and
`REG_N` N-ary kernels cover the hypergeometric/Whittaker family. `optimize.c`
then runs per-block value numbering (constant folding, CSE, copy propagation),
liveness-driven dead-code elimination and loop-invariant code motion. Its passes
are restricted to transforms that keep results bitwise identical (no FP
reassociation, no FMA contraction). A program with `WorkingPrecision -> n`
(converted to MPFR bits) or `"BigIntegers" -> True` instead goes through
`compile_mgd.c`, a separate managed-scalar emitter (MPFR real/complex, GMP
integers), with the optimiser and CSE forced off, so the machine path is
provably unchanged. A body that fails anywhere records a bail reason and the
offending subexpression, which `CompileDiagnostics` reports.

**Data structures.** The result is an `EXPR_COMPILED` atom holding a
reference-counted, immutable `CompiledFunction`. It stores the argument names and
types, the held body (for printing and the interpreter fallback), the program
(`NULL` when the body is outside the subset), and the runtime attributes. The VM
(`compile_vm.c`) is a register machine with three banks: scalar `R`, array handle
`V`, and strip-mining tile `T`. It uses a reusable per-program frame (no per-call
`malloc`). Element-wise array expressions are fused into one strip-mined pass,
and large ones are split across threads with `nd_parallel_for` when threads are
enabled. Managed registers hold pointers to per-call `mpfr_t` / `mpz_t` / complex
containers that are bound at frame entry. A `List` argument is packed into a
machine buffer at the boundary (an `int64` buffer reaching a `_Real` slot costs
one cast, and the call declines past 2^53). A packed or `NDArray` argument is
borrowed in place, and the result inherits its presentation. `autocompile.c`
reuses the same compiler for `Plot`, `Table`, `NIntegrate`, `NDSolve` and loop
bodies, gated by `$AutoCompilation`. `MATHILDA_COMPILE_DIAG=1` names every
automatic fallback.

**Complexity / limits.** Applying the object to numeric arguments runs the
bytecode with no `Expr` allocation. A symbolic argument, or a body outside the
compilable subset, transparently falls back to the interpreter, so the value is
always what the uncompiled body would give. Integer arithmetic is exact: an
`int64` overflow abandons the compiled call and the interpreter promotes to a
bignum, unless `RuntimeOptions -> "Speed"` opts into wrapping. The check is
compiled into each instruction, so leaving it on costs nothing measurable.
`Divide`, `Sqrt` and the transcendentals widen integers to `Real`, as in the
Wolfram Language. The managed (arbitrary-precision) subset is straight-line
arithmetic plus the unary elementary functions. Control flow, arrays and other
special functions fall back to the interpreter at that precision. Precision is
held fixed for the whole call, so the last digit can differ from `N[expr, n]`.
Arbitrary-precision arrays are not implemented yet.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Mod](../../arithmetic/Mod/), [Quotient](../../arithmetic/Quotient/), [Power](../../arithmetic/Power/), [Gamma](../../special-functions/Gamma/), [Erf](../../special-functions/Erf/), [BesselJ](../../special-functions/BesselJ/), [Zeta](../../special-functions/Zeta/)

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bitwise.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bitwise.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_arbprec.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_arbprec.c)

## Notes & additional examples

### Notes

- `Compile[{arg, …}, expr]` builds a `CompiledFunction` over raw machine numbers. Argument types are `_Real` (the default for a bare symbol), `_Integer`, `_Complex`, rank-`r` arrays (`{v, _Real, 1}`, `{m, _Integer, 2}`) and read-only `_Association` parameter bags.
- Attributes `HoldAll`, `Protected`: the body compiles in its raw form and argument symbols stay local, so a global value for one does not leak in.
- Fallback is part of the contract. A symbolic argument, or a body outside the compilable subset, runs through the interpreter, so a `CompiledFunction` returns what the uncompiled body would. The subset is a cliff: one unlowered head (here `myf`) costs the whole body. `CompileDiagnostics` names the subexpression that bailed, and `CompilePrint` shows the bytecode.
- Integers: an `int64` overflow abandons the compiled call and the interpreter promotes to a bignum. `RuntimeOptions -> "Speed"` (that is, `"CatchMachineIntegerOverflow" -> False`) keeps the wrapped value instead. That changes the answer, not just the speed, so it is opt-in. Integer-closed heads (`Factorial`, `Binomial`, `Power` with a non-negative exponent, …) return `Integer`s. `Divide`, `Sqrt` and the transcendentals widen to `Real`, so `n/2` gives `1.5`, as in the Wolfram Language.
- Every head with a machine kernel in the shared registry lowers through the same helper as the scalar and packed surfaces, so results agree bit for bit. Since v0.266 the registry gained the hypergeometric family, `Erf`/`Erfc`/`Erfi`/Fresnel, `Zeta`/`HurwitzZeta`/`LerchPhi`/`PolyLog`, `StruveH` and complex arms, plus whole-array `Chop`/`Clip`.
- `Chop` compiles only with a machine-real tolerance (`0.000001`). An exact one (`10^-6`) does not lower, and the interpreted `Chop` then answers with an exact `0` where the compiled one gives `0.0`.
- `WorkingPrecision -> n` and `"BigIntegers" -> True` are off by default and use a separate managed emitter, so the machine path is untouched. Precision is held fixed for the whole call, so the last digit or two can differ from `N[expr, n]`. Only straight-line arithmetic and the elementary functions lower in this mode; anything else falls back at that precision.
- `$AutoCompilation = False` switches off the hidden compilation inside `Plot`, `Table`, `NIntegrate`, `NDSolve`, loops and so on. It does not affect an explicit `Compile`.
- See also: `CompileDiagnostics`, `CompilePrint`, `$AutoCompilation`, `NDArray`, `Function`.
