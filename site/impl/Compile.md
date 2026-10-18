---
source: src/compile/compiled_function.c, src/compile/compile.c, src/compile/compile_infer.c, src/compile/compile_emit_mathfn.c, src/compile/compile_mgd.c, src/compile/optimize.c, src/compile/compile_vm.c, src/compile/autocompile.c
---
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
