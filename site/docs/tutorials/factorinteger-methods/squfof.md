# SQUFOF: Shanks' square forms

`Method -> "ShanksSquareForms"` (alias `"SQUFOF"`) is Daniel Shanks'
square-form factorisation. It walks the continued-fraction expansion of
\(\sqrt{kn}\) until it meets a square, and turns that square into a factor.
Its expected cost is \(O(n^{1/4})\), whatever the size of the factors. Every
number it handles stays below \(2\sqrt{kn}\), so for numbers up to about 62
bits the whole computation runs in single machine words with no bignum
arithmetic. That is why it serves as a backstop in the word engine of every
fast factoring library, Mathilda's included.

```mathematica
In[1]:= n62 = NextPrime[2^30] NextPrime[2^31 + 2^29]; BitLength[n62]
Out[1]= 62

In[2]:= FactorInteger`SQUFOF[n62]
Out[2]= {{1073741827, 1}, {2684354591, 1}}
```

---

## 1. The mathematical idea

Expand \(\sqrt N\), with \(N = kn\) not a square, as a continued fraction. Put
\(s = \lfloor\sqrt N\rfloor\), \(P_0 = s\), \(Q_0 = 1\), \(Q_1 = N - s^2\), and
iterate

\[
b_i = \left\lfloor \frac{s + P_{i-1}}{Q_i} \right\rfloor,\qquad
P_i = b_i Q_i - P_{i-1},\qquad
Q_{i+1} = Q_{i-1} + b_i\,(P_{i-1} - P_i).
\]

The invariant \(P_i^2 + Q_i Q_{i+1} = N\) holds throughout, and
\(0 < P_i, Q_i < 2\sqrt N\). In the language of binary quadratic forms, each
step is one reduction step on a form of discriminant \(4N\).

**Forward phase.** Run until \(Q_i = r^2\) is a perfect square at an even
index \(i\). The form \((Q_i, 2P_{i-1}, \cdot)\) is then the square of a form
in the class group, and \((r, \cdot, \cdot)\) is its square root.

**Reverse phase.** Start again from that square-root form:
\(Q_0' = r\), \(P_0' = P_{i-1} + r\lfloor(s - P_{i-1})/r\rfloor\), and
\(Q_1' = (N - P_0'^2)/r\). Iterate the same recurrence until \(P\) repeats,
\(P_j' = P_{j-1}'\). That symmetry point is an *ambiguous* form, and
\(\gcd(n, P_j')\) is a factor of \(n\) unless it is trivial.

The forward phase needs \(O(N^{1/4})\) steps, since the period of the
expansion is about \(\sqrt N\) and squares occur at density
\(N^{-1/4}\). The reverse phase takes about half as many. Some square forms
give only a trivial factor. When that happens, or a multiplier runs too long,
the method moves on to the next **multiplier** \(k\). Gower and Wagstaff's
analysis suggests racing the 16 squarefree products of \(3, 5, 7, 11\), since
some \(k\) makes the expansion of \(\sqrt{kn}\) much more favourable.

### References

1. D. Shanks, unpublished notes on SQUFOF (c. 1975). Reference 2 gives a full
   account of them.
2. J. E. Gower and S. S. Wagstaff, Jr., "Square form factorization,"
   *Mathematics of Computation* **77** (2008) 551–588.
3. H. Riesel, *Prime Numbers and Computer Methods for Factorization*, 2nd ed.,
   Birkhäuser, 1994 — Ch. 6.
4. H. Cohen, *A Course in Computational Algebraic Number Theory*, Springer,
   1993 — Ch. 8.

---

## 2. How Mathilda realises it

| file | function | range | arithmetic |
|---|---|---|---|
| `fi_u64.c` | `squfof_one`, `fi_u64_squfof` | \(kn < 2^{62}\) | all in `uint64_t` |
| `fi_squfof.c` | `squfof_k`, `fi_split_squfof` | \(kn < 2^{124}\) | GMP for \(\sqrt{kn}\), \(kn - s^2\) and the reverse-phase start; the main loops in words |

Even when \(n\) has two limbs, the forms are bounded by
\(2\sqrt{kn} < 2^{63}\), so the main loops are plain signed word arithmetic.
The bignum code is only needed for the setup and for one subtraction at the
start of the reverse phase.

The multipliers are the 16 Gower–Wagstaff values

\[
k \in \{1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155\},
\]

tried **in turn**. A multiplier for which \(kn\) overflows the range is
skipped. Square candidates are first screened by a quadratic-residue mask
modulo 64. Only values that pass it get an integer square root.

| Suboption | Default | Meaning |
|---|---|---|
| `"Iterations"` | \(4\sqrt{2\sqrt{kn}} + 64\) | forward steps per multiplier, and also the reverse-phase cap |

For a word \(n\) the dispatch tries `fi_u64_squfof` first when \(n < 2^{62}\).
If that fails, it tries the two-word version, which can still use the
multipliers for which \(kn\) exceeds 62 bits but not 124.

**Inside `Automatic`.** SQUFOF is the third splitter in the word engine
(`u64_split` in `fi_u64.c`). The order is: rho for \(n < 2^{40}\), then
micro-ECM with 128 curves, then SQUFOF with its default budget, then a long
rho backstop. It is never used above \(2^{64}\).

---

## 3. Invoking it

```mathematica
In[3]:= FactorInteger[n62, Method -> "ShanksSquareForms"]
Out[3]= {{1073741827, 1}, {2684354591, 1}}

In[4]:= FactorInteger`ShanksSquareForms[n62]
Out[4]= {{1073741827, 1}, {2684354591, 1}}

In[5]:= n80 = NextPrime[2^39] NextPrime[2^40 + 12345]; FactorInteger[n80, Method -> {"SQUFOF", "Iterations" -> 1000}]
Out[5]= {{604462916622637412406697, 1}}
```

`"SQUFOF"` works both as a `Method` string and as a head,
``FactorInteger`SQUFOF``. With only 1000 steps per multiplier, an 80-bit
\(n\), which needs about \(2^{20}\) steps, comes back unfactored.

---

## 4. Where it shines

**The word range.** Up to about 64 bits, SQUFOF is fast and reliable:

```mathematica
In[6]:= FactorInteger`SQUFOF[2^64 + 1]
Out[6]= {{274177, 1}, {67280421310721, 1}}

In[7]:= FactorInteger`SQUFOF[n80]
Out[7]= {{549755813911, 1}, {1099511640127, 1}}
```

Here are 300 semiprimes of about 62 bits, each the product of a 30-bit and a
32-bit prime. All of them are factored, in a fraction of a second in total:

```mathematica
In[8]:= okQ[m_, f_] := Times @@ (Power @@@ f) == m && AllTrue[f[[All, 1]], PrimeQ]

In[9]:= semis = Table[NextPrime[2^29 + 1000 j] NextPrime[2^31 + 777 j], {j, 1, 300}]; AllTrue[semis, okQ[#, FactorInteger`SQUFOF[#]] &]
Out[9]= True
```

**Two-word inputs.** Up to about 110 bits the two-word code still works, with
the cost growing as \(n^{1/4}\). It is about 1 ms at 80 bits (\(2^{20}\)
steps), and tenths of a second at 100–111 bits:

```mathematica
In[10]:= FactorInteger`SQUFOF[NextPrime[2^45] NextPrime[2^46]]
Out[10]= {{35184372088891, 1}, {70368744177679, 1}}

In[11]:= FactorInteger`SQUFOF[NextPrime[2^55] NextPrime[2^56]]
Out[11]= {{36028797018963971, 1}, {72057594037928017, 1}}
```

---

## 5. Where it struggles

**The cost does not depend on the factor.** Rho finds a 21-bit factor in
about \(2^{10}\) steps. SQUFOF on the same 101-bit number still needs its full
\(n^{1/4}\approx2^{25}\) forward steps, about 3 s on the loaded reference
machine against well under a millisecond for rho:

```mathematica
In[12]:= u = NextPrime[2^20] NextPrime[2^80]; FactorInteger`PollardRho[u]
Out[12]= {{1048583, 1}, {1208925819614629174706189, 1}}

In[13]:= FactorInteger`SQUFOF[u]
Out[13]= {{1048583, 1}, {1208925819614629174706189, 1}}
```

**It declines beyond about 124 bits.** The forms must fit a signed word, so
\(kn < 2^{124}\). Above about 113 bits the largest multiplier, 1155, drops
out. At 124
bits only \(k = 1\) is left, and at 125 bits nothing is left. A 126-bit
balanced semiprime is returned unfactored immediately, without a single step:

```mathematica
In[14]:= n126 = NextPrime[2^62] NextPrime[2^63 + 12345]; BitLength[n126]
Out[14]= 126

In[15]:= FactorInteger`SQUFOF[n126]
Out[15]= {{42535295865117366201574668759018862747, 1}}

In[16]:= FactorInteger[n126]
Out[16]= {{4611686018427388039, 1}, {9223372036854788173, 1}}
```

The automatic pipeline factors `n126` in milliseconds with two-limb
micro-ECM and the quadratic sieve. It would not even make sense to run
SQUFOF there: near the top of its range, \(n^{1/4}\approx2^{31}\) steps is
more than ten seconds of work. At 120 bits, ``FactorInteger`SQUFOF`` on a
balanced semiprime took about 7 s on the reference machine.

---

## 6. Limitations

- \(O(n^{1/4})\) regardless of the factors. It cannot take advantage of a small
  factor, which [rho](pollard-rho.md) or [ECM](ecm.md) find far faster.
- Works only while \(kn < 2^{124}\). Beyond that it declines, returning the
  composite with `FactorInteger::nofac`.
- The multipliers are raced *sequentially*, each with a full budget, so an
  input that is unlucky for \(k = 1\) pays for a whole failed multiplier
  before the next one starts.
- In practice it is useful only on machine words. Above about 70 bits,
  two-limb rho, micro-ECM and QS are all faster.

---

## See also

- [FactorInteger methods](index.md): the word engine and where SQUFOF sits in
  it.
- [Fermat and Lehman](fermat-lehman.md): the other word-size methods based on
  squares.
- [CFRAC and Dixon](cfrac-dixon.md): CFRAC uses the same continued-fraction
  expansion of \(\sqrt{kn}\), but collects smooth \(Q_i\) for linear algebra
  instead of waiting for a square.
- `?FactorInteger`ShanksSquareForms` in the REPL.
