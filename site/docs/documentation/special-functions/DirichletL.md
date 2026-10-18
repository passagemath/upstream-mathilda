# DirichletL

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DirichletL[k, j, s]`**

gives the Dirichlet L-function L(chi, s) for the character

<details>
<summary>Notes</summary>

chi = DirichletCharacter\[k, j, .\]; for Re s \> 1 it is Sum\_{n\>=1} chi(n) n^-s. DirichletL\[1, 1, s\] is Zeta\[s\]. Exact values are generated for non-positive integer s (generalized Bernoulli numbers) and for s = 1 with a non-principal character. Evaluates to any precision for real or complex s, and encloses Interval arguments for real characters. Listable.

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= DirichletL[4, 2, 1]
Out[1]= 1/4 Pi

In[2]:= DirichletL[12, 3, 1]
Out[2]= 1/3 Pi

In[3]:= {DirichletL[4, 2, 0], DirichletL[4, 2, -2], DirichletL[4, 2, -3]}
Out[3]= {1/2, -1/2, 0}

In[4]:= N[DirichletL[4, 2, 2], 30]
Out[4]= 0.9159655941772190150546035149325
```

### Applications (12)

The Leibniz series 1 - 1/3 + 1/5 - ... = Pi/4

```mathematica
In[5]:= DirichletL[4, 2, 1]
Out[5]= 1/4 Pi
```

The real character mod 3

```mathematica
In[6]:= DirichletL[3, 2, 1]
Out[6]= Pi/(3 Sqrt[3])
```

An even character: the logarithmic class-number form

```mathematica
In[7]:= DirichletL[8, 2, 1]
Out[7]= -(Log[Sin[1/8 Pi]] - Log[Sin[3/8 Pi]])/Sqrt[2]
```

An imprimitive character: Pi/4 times the Euler factor 1 + 1/3

```mathematica
In[8]:= DirichletL[12, 3, 1]
Out[8]= 1/3 Pi
```

Non-positive integers give generalized Bernoulli numbers; odd characters vanish at negative odd s

```mathematica
In[9]:= {DirichletL[4, 2, 0], DirichletL[4, 2, -2], DirichletL[4, 2, -3]}
Out[9]= {1/2, -1/2, 0}
```

A complex character gives a value in Q(I)

```mathematica
In[10]:= DirichletL[5, 2, -4]
Out[10]= 148/25 + 86/25*I
```

The value at s = 2 for the character mod 4 is Catalan's constant

```mathematica
In[11]:= {N[DirichletL[4, 2, 2], 30], N[Catalan, 30]}
Out[11]= {0.9159655941772190150546035149325, 0.9159655941772190150546035149325}
```

Near the first zero on the critical line

```mathematica
In[12]:= DirichletL[4, 2, 0.5 + 6.0209489 I]
Out[12]= 8.58558e-10 - 6.09043e-09*I
```

A complex character at complex s, 25 digits through arb

```mathematica
In[13]:= N[DirichletL[5, 2, 5/2 + I], 25]
Out[13]= 1.0515319720792428817994336 + 0.12794042693768998717861571*I
```

A packed Range of s runs the machine kernel over the buffer

```mathematica
In[14]:= DirichletL[3, 2, Range[0.5, 2., 0.5]]
Out[14]= {0.480868, 0.6046, 0.703968, 0.781302}
```

A rigorous enclosure for a real character

```mathematica
In[15]:= DirichletL[5, 1, Interval[{1.23, 1.24}]]
Out[15]= Interval[{4.11396, 4.25899}]
```

An exact s >= 2 stays unevaluated, as in Mathematica

```mathematica
In[16]:= DirichletL[4, 2, 2]
Out[16]= DirichletL[4, 2, 2]
```

## Algorithm

dirichlet_l.c -- DirichletL[k, j, s], the Dirichlet L-function L(chi, s) of the character chi = DirichletCharacter[k, j, .].

Attributes: Listable | NumericFunction | Protected. N keeps k and j exact (numericalize_rec in src/numeric.c), as Mathematica's NHoldAll + N rule does.

Evaluation order:

```text
  1. arity (argrx); an NDArray / packed s -> the buffer loop below.
  2. (k, j) validation via the character engine (intp / invl messages).
  3. k = 1: L = Zeta[s] for numeric s (and Interval s).
  4. Interval s, real character: a rigorous arb enclosure.
  5. machine s: the raw-double kernel sf_dirichlet_l_tab_* -- the same one
     the NDArray loop and Compile[] use (src/ndkernels.c ndk_DirichletL_n).
  6. exact s: the closed forms Mathematica gives -- s a non-positive integer
     (generalized Bernoulli numbers, all characters) and s = 1 for a
     non-principal character -- else unevaluated.
  7. inexact s beyond machine (or a machine decline): FLINT acb_dirichlet_l on
     the primitive character plus Euler factors, at the precision of s; the
     fallback build sums HurwitzZeta templates.
```

Exact values live in the cyclotomic field Q(zeta_O), O = order of chi; they are computed as coefficient vectors in Q[x]/(x^O - 1), reduced modulo the cyclotomic polynomial (so the representation is canonical), and rendered as sums of rational multiples of roots of unity in the DirichletCharacter normal form (unit * E^(I Pi r)).

## Implementation notes

**Algorithm.** `builtin_dirichletl` evaluates in this order. (1) Arity is checked, and an `NDArray` or packed `s` goes to `dirichletl_ndarray`, which fetches the character table once and loops the machine kernel. (2) `(k, j)` is validated through the character engine. (3) `k = 1` delegates numbers and intervals to `Zeta[s]`. (4) An `Interval` `s` with a real character gets a rigorous arb enclosure (`flint_dirichlet_l_interval`). That routine bounds monotone pieces by their endpoints and bisects the other pieces up to depth 12. (5) A machine `s` goes to `sf_dirichlet_l_tab_r` (a real character at real `s`) or `sf_dirichlet_l_tab_c`. These reduce `χ` to `χ*` times Euler factors. For `Re s ≥ 1/2` they apply Euler–Maclaurin to `L` itself, with 10 Bernoulli terms, a tail start `M` from the remainder bound, and the integral terms summed through `log1p`/`expm1` so nothing cancels as `s → 1`. For `Re s < 1/2` they use the functional equation with the root number from the table. The principal character at `s = 1` gives `ComplexInfinity`, and any other decline falls back to arb. (6) An exact `s` goes to `dl_exact`. At `s = 1 - m ≤ 0`, `dl_exact_nonpos` computes the generalized Bernoulli numbers of `χ*` times Euler factors. At `s = 1` with a non-principal character, `dl_exact_one` uses the classical odd form (π, Gauss sum) and even form (`Log[Sin[π a/f]]`). Other exact `s` stays unevaluated. (7) An MPFR `s` calls `flint_num_dirichlet_l`, which is arb's `acb_dirichlet_l` on `χ*`, set on FLINT's own generators from Mathilda's character values rather than through Conrey labels, times `nb_dirichlet_euler`, with a precision-doubling ladder. Without FLINT it sums `f^{-s} Σ χ*(a) HurwitzZeta[s, a/f]`.

**Data structures.** Exact values are elements of the cyclotomic field `Q(ζ_O)`, stored as `mpq_t` coefficient vectors (`Cyc`) in `Q[x]/(x^O - 1)` and reduced modulo the cyclotomic polynomial (`cyc_reduce`), so each value has one canonical form. They are rendered in `DirichletCharacter`'s unit-times-`E^(I Pi r)` normal form. The machine kernel reads the dense `χ*` tables of the shared `DchTable`. `N` keeps `k` and `j` exact (`numericalize_rec`). `Compile[]` lowers the head through the n-ary kernel `ndk_DirichletL_n`, at scalar shape and with `s` as an array tile, and it is typed `Real` only for literal labels of a real character. Scalar, packed, `NDArray` and compiled results agree bit for bit.

**Complexity / limits.** A machine evaluation costs about `M f` terms (about 1.5–4 µs per call at small `k`). The machine kernel declines beyond about 2·10⁷ terms or a conductor above 2^22. The exact closed forms are bounded: `f ≤ 2^20` for `s ≤ 0` (`DL_EXACT_FMAX`), and a cyclotomic order of at most 4096 at `s = 1`; beyond these bounds the call stays unevaluated. Exact `s ≥ 2` for `k > 1`, non-integer exact `s`, and an interval containing `s = 1` all stay unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [DirichletCharacter](../../number-theory/DirichletCharacter/), [NDArray](../../linear-algebra/NDArray/)

- DLMF §25.15 — Dirichlet L-functions (functional equation 25.15.5, values at non-positive integers via generalized Bernoulli numbers).
- F. Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Comput. 66 (2017) 1281–1292 — acb_dirichlet_l.
- Source: [`src/special_functions/dirichlet_l.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/dirichlet_l.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_dirichletl.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dirichletl.c)

## Notes & additional examples

### Notes

- `DirichletL[k, j, s]` is `L(χ, s) = Σ χ(n) n^(-s)`, continued analytically to the whole plane, for `χ = DirichletCharacter[k, j, ·]`. It uses the same labelling as `DirichletCharacter`.
- Exact values are given at `s ≤ 0`, where they are generalized Bernoulli numbers, and at `s = 1` for every non-principal character. This includes imprimitive characters, which Mathematica leaves unevaluated at `s = 1`. Exact `s ≥ 2` stays unevaluated, as in Mathematica.
- Numerics: a machine `s` uses an Euler–Maclaurin / functional-equation kernel (max relative error about 6e-14 against arb in testing). Arbitrary precision uses arb's `acb_dirichlet_l`, and an `Interval` `s` with a real character gives a rigorous enclosure. `N` keeps `k` and `j` exact.
- Surfaces: packed and `NDArray` `s` with scalar labels, and `Compile[]` at scalar and array shapes, all share the machine kernel bit for bit.
- Attributes: `Listable`, `NumericFunction`, `Protected`. `DirichletL[1, 1, s]` is `Zeta[s]` for numeric `s`.
- See also: `DirichletCharacter`, `Zeta`, `HurwitzZeta`, `Catalan`.
