### Worked examples

```mathematica
In[1]:= f = Compile[{{x, _Real}}, x^2 + 1]; f[3.0]  (* a machine-number function applied to a Real *)
```

```mathematica
In[1]:= f[a]  (* a symbolic argument transparently falls back to the interpreter *)
```

```mathematica
In[1]:= Compile[{{n, _Integer}}, n^3][3000000]  (* integer arithmetic stays exact past the int64 range: overflow defers to the interpreter *)
```

```mathematica
In[1]:= Compile[{{n, _Integer}}, n^3, RuntimeOptions -> "Speed"][3000000]  (* opting out of overflow checking keeps the wrapped int64, 27*10^18 mod 2^64 *)
```

```mathematica
In[1]:= Compile[{{z, _Complex}}, Gamma[z] + Erf[z]][0.5 + 1.0 I]  (* special functions with a complex machine kernel lower directly *)
```

```mathematica
In[1]:= Compile[{{v, _Real, 1}}, Clip[Chop[v - 1., 0.000001], {-1., 1.}]][{0.5, 1.0000000001, 3., -4.}]  (* whole-array Chop and Clip over a rank-1 parameter *)
```

```mathematica
In[1]:= Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable][{{1., 2.}, {3., 4.}}]  (* a rank-1 parameter consumes one level, so a Listable object maps over the rows *)
```

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + myf[x]]  (* one unlowered head sends the whole body to the interpreter; the diagnostics name it *)
```

```mathematica
In[1]:= Compile[{{x, _Real}}, Sin[x] Cos[x] + x^3, WorkingPrecision -> 40][N[7/5, 40]]  (* opt-in MPFR arithmetic at a fixed 40 digits *)
```

```mathematica
In[1]:= Compile[{{n, _Integer}}, n^20, "BigIntegers" -> True][99]  (* opt-in exact GMP integers inside the compiled program *)
```

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
