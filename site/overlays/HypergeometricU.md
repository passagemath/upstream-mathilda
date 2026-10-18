### Worked examples

```mathematica
In[1]:= HypergeometricU[1, 2, z]  (* the reduction U(a, a+1, z) = z^-a *)
```

```mathematica
In[1]:= HypergeometricU[1.5, 7, 4.5]  (* machine value *)
```

```mathematica
In[1]:= NIntegrate[Exp[-4.5 t] t^0.5 (1 + t)^4.5, {t, 0, Infinity}]/Gamma[1.5]  (* the Laplace-type integral U = (1/Γ(a))∫ e^(-zt) t^(a-1) (1+t)^(b-a-1) dt gives the same *)
```

```mathematica
In[1]:= N[HypergeometricU[3/2, 7, 9/2], 30]  (* arbitrary precision through FLINT's acb_hypgeom_u *)
```

```mathematica
In[1]:= HypergeometricU[1/2, 1, 2.] - Exp[1.] BesselK[0, 1.]/Sqrt[Pi]  (* integer b: the log series, checked against U(1/2,1,z) = e^(z/2) K_0(z/2)/√π *)
```

```mathematica
In[1]:= HypergeometricU[0.7, 1.3, -2.]  (* real z < 0 lies on the cut: the principal complex value *)
```

```mathematica
In[1]:= HypergeometricU[3 - I, 2 I, 3.0]  (* complex parameters go to arb, demoted to machine precision *)
```

```mathematica
In[1]:= D[HypergeometricU[a, b, z], z]  (* DLMF 13.3.22 *)
```

```mathematica
In[1]:= HypergeometricU[1, 2, {2., 4., 5.}]  (* listable over a packed buffer *)
```

```mathematica
In[1]:= HypergeometricU[-2, b, z]  (* the Laguerre-polynomial case is not expanded for symbolic arguments *)
```

### Notes

- `Listable`, `NumericFunction`, `Protected`. The function is the Tricomi solution of `z w'' + (b − z) w' − a w = 0`. Its branch cut runs along z ∈ (−∞, 0].
- Exact reductions: `U(0, b, z) = 1` and `U(a, a+1, z) = z^(−a)`. All other exact input stays symbolic, including the Laguerre-polynomial case of a non-positive integer a.
- Machine input (real a and b; z real or complex) runs a raw-double kernel. It tries, in order: the asymptotic 2F0 series, accepted only when its least term is below half an ulp; the DLMF 13.2.9 log series for integer b, whose digammas advance by recurrence; the Kummer two-1F1 connection for non-integer b; and an exp-sinh quadrature of the integral form. Packed and `NDArray` buffers share the kernel. A declined point goes to FLINT's `acb_hypgeom_u`, and the result is rounded back to machine precision.
- Complex parameters and precision above machine always use `acb_hypgeom_u`.
- `D` in z follows DLMF 13.3.22. Derivatives in a and b stay inert. There is no origin series and no `FunctionExpand` rule.
- See also: `Hypergeometric1F1`, `WhittakerW`, `MeijerG`.
