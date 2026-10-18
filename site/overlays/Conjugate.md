### Worked examples

```mathematica
In[1]:= Conjugate[3 + 4 I]
```

```mathematica
In[1]:= Conjugate[{1 + I, 2 - 3 I, 5}]  (* Listable; real numbers are fixed *)
```

```mathematica
In[1]:= Conjugate[(2 + I)/(1 - 3 I)]  (* the conjugate of an exact quotient *)
```

```mathematica
In[1]:= z Conjugate[z] /. z -> 3 + 4 I  (* the squared modulus *)
```

```mathematica
In[1]:= Conjugate[Conjugate[z]]  (* an involution, even for symbolic z *)
```

```mathematica
In[1]:= Conjugate[Sqrt[2] Log[3]]  (* an exact real constant is its own conjugate *)
```

```mathematica
In[1]:= Conjugate[Sqrt[-2]]
```

```mathematica
In[1]:= Conjugate[NDArray[{1. + 2. I, 3. - I}]]  (* the NDArray kernel *)
```

```mathematica
In[1]:= Conjugate[z]  (* symbols are not assumed real *)
```

### Notes

- `Conjugate[z]` is `Re[z] - I Im[z]`. Real numbers, exact real constants (tested by evaluating them numerically at machine precision) and calls to the real-valued heads `Re`, `Im`, `Abs`, `Arg` are fixed points, and `Conjugate[Conjugate[z]]` is `z`.
- Symbols are not assumed real, so `Conjugate[z]` stays unevaluated. Calling it with any number of arguments other than one emits `Conjugate::argx`.
- `Listable`, `NumericFunction`, `Protected`. There is a packed/NDArray kernel (`conj`), and `Compile[]` lowers `Conjugate` at scalar and array shapes.
- See also: `Re`, `Im`, `Abs`, `ConjugateTranspose`.
