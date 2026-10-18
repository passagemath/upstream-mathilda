### Worked examples

```mathematica
In[1]:= Table[LegendreP[n, x], {n, 0, 4}]  (* exact polynomials with rational coefficients *)
```

```mathematica
In[1]:= LegendreP[-3, x]  (* negative order by P_{-n-1} = P_n *)
```

```mathematica
In[1]:= LegendreP[n, 1]  (* every P_n equals 1 at the endpoint *)
```

```mathematica
In[1]:= Integrate[LegendreP[2, x] LegendreP[3, x], {x, -1, 1}]  (* orthogonality on the interval *)
```

```mathematica
In[1]:= LegendreP[200, 0.7]  (* high degree by the machine recurrence *)
```

```mathematica
In[1]:= LegendreP[5, 0.3 + 0.2 I]  (* complex argument *)
```

```mathematica
In[1]:= LegendreP[2, 1, x]  (* the associated function P_2^1 *)
```

```mathematica
In[1]:= LegendreP[30, 4, {0.1, 0.5}]  (* associated Legendre over a machine array, by Horner *)
```

```mathematica
In[1]:= N[LegendreP[1/3, 1/5], 30]  (* non-integer order through the Gauss 2F1 series *)
```

```mathematica
In[1]:= D[LegendreP[n, x], x]  (* the derivative for symbolic order *)
```

```mathematica
In[1]:= LegendreP[0.5, 3.]  (* non-integer order outside |(1-x)/2| < 1: stays unevaluated *)
```

### Notes

- For an exact integer order, `LegendreP[n, x]` is the explicit Legendre
  polynomial `P_n(x)` with exact rational coefficients (orders up to
  `n = 2000`). `LegendreP[n, 1]` is `1` for every order.
- A machine-precision call with an integer order, real or complex `x`, runs the
  three-term recurrence in doubles (degree up to 1000). That kernel is shared
  with packed arrays, `NDArray` arguments and `Compile[]` (scalar, complex and
  rank-1 array), so these all agree bit for bit.
- A non-integer order is evaluated numerically, when an argument is inexact,
  through `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)` at machine or arbitrary precision.
  This needs `|(1-x)/2| < 1`; outside that disc the call stays unevaluated.
- `LegendreP[n, m, x]` is the associated Legendre function, with the
  Condon–Shortley phase `(-1)^m`. `LegendreP[n, m, a, x]` selects type `a`
  (1, 2 or 3). These need integer `n` and integer `m >= 0`. Over a machine array
  the polynomial is evaluated by Horner, which keeps high degree fast.
- Attributes `Listable`, `NumericFunction`, `Protected`.
- See also: `LegendreQ`, `ChebyshevT`, `ChebyshevU`, `HermiteH`,
  `Hypergeometric2F1`.
