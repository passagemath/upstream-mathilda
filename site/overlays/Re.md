### Worked examples

```mathematica
In[1]:= Re[3 + 4 I]
```

```mathematica
In[1]:= Re[{1 + I, 2 - 3 I, 5}]  (* Listable *)
```

```mathematica
In[1]:= Re[(1 + I)^10]  (* (1 + I)^10 = 32 I is purely imaginary *)
```

```mathematica
In[1]:= Re[1/(2 + 3 I)]  (* the exact real part of a quotient *)
```

```mathematica
In[1]:= Re[10^30]  (* a bignum is its own real part *)
```

```mathematica
In[1]:= Re[Re[z]]  (* Re, Im, Abs and Arg are real-valued by construction *)
```

```mathematica
In[1]:= Re[Root[#^3 - 2 &, 1]]  (* a real Root object is returned unchanged *)
```

```mathematica
In[1]:= Re[NDArray[{1. + 2. I, 3. - I}]]  (* the NDArray projection kernel *)
```

```mathematica
In[1]:= Re[x + I y]  (* symbols are not assumed real: use ComplexExpand for that *)
```

### Notes

- `Re[z]` is the real part of a numeric `z`. A real number comes back unchanged.
- Symbols are not assumed real, so `Re[x + I y]` stays unevaluated; `ComplexExpand[Re[x + I y]]` gives `x`. A non-real `Root` or `AlgebraicNumber` also stays wrapped, as in Mathematica.
- Limitation: an exact real constant with no complex part written out, such as `Re[Pi]`, is not yet recognised and stays unevaluated (Mathematica gives `Pi`).
- `Listable`, `NumericFunction`, `Protected`. There is a packed/NDArray projection kernel, and `Compile[]` lowers `Re` at scalar and array shapes.
- See also: `Im`, `ReIm`, `Conjugate`, `ComplexExpand`.
