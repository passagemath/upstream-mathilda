### Worked examples

```mathematica
In[1]:= LegendreQ[0, x]  (* the second-kind solution carries a logarithm *)
```

```mathematica
In[1]:= LegendreQ[2, x]  (* P_2(x) L(x) plus a polynomial *)
```

```mathematica
In[1]:= N[LegendreQ[3, 1/3], 30]  (* high precision from the closed form *)
```

```mathematica
In[1]:= LegendreQ[3, {0.1, 0.5, 0.9}]  (* Listable, machine kernel *)
```

```mathematica
In[1]:= LegendreQ[2, 0.3 + 0.4 I]  (* complex argument *)
```

```mathematica
In[1]:= LegendreQ[3, 2.]  (* beyond x = 1 the type-1 value is complex *)
```

```mathematica
In[1]:= LegendreQ[2, 1, 0.5]  (* the associated function Q_2^1 *)
```

```mathematica
In[1]:= LegendreQ[1/2, 0]  (* an exact value at the origin for non-integer order *)
```

```mathematica
In[1]:= N[LegendreQ[2.5, 0.4]]  (* non-integer order on the cut, by two Frobenius 2F1 series *)
```

```mathematica
In[1]:= LegendreQ[1.5, 2.]  (* non-integer order off the cut: stays unevaluated *)
```

### Notes

- `LegendreQ[n, x]` is the Legendre function of the second kind `Q_n(x)`
  (type 1). Unlike `P_n` it carries a logarithm even at integer order: the
  result is `P_n(x) L(x) + v_n(x)` with `L(x) = (1/2)(Log[1+x] - Log[1-x])` and
  a polynomial `v_n` (orders up to `n = 2000`). Negative integer orders are
  singular and stay symbolic.
- A machine-precision call with an integer order, real or complex `x`, uses one
  double kernel (`Q_n = P_n Q_0 - W_{n-1}`, with a cancellation-free `Q_0`)
  shared with packed arrays, `NDArray` arguments and `Compile[]`. A real `x`
  with `|x| > 1` gives a complex value, following the branch of the closed form.
  `LegendreQ[n, 1.]` is `Indeterminate`.
- A non-integer order is evaluated through two Frobenius `2F1` series on the cut
  `|x| < 1`. The exact special value `Q_v(0)` is returned for an exact zero
  argument.
- `LegendreQ[n, m, x]` and `LegendreQ[n, m, a, x]` (type `a` in 1, 2, 3) are
  obtained by differentiating `Q_n`, for integer `n, m >= 0`.
- Attributes `Listable`, `NumericFunction`, `Protected`.
- See also: `LegendreP`, `ArcTanh`, `Hypergeometric2F1`.
