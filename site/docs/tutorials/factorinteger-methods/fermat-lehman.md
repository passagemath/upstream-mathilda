# Fermat and Lehman

`Method -> "Fermat"` and `Method -> "Lehman"` look for the factors of \(n\)
from the *middle*, not from below. Every other method in Mathilda prefers a
small factor. Fermat's method prefers two factors that are **close to each
other**, and gets quicker the closer they are. Lehman's method generalises it
to two factors whose **ratio is close to a simple fraction**, and on machine
words it is a complete \(O(n^{1/3})\) algorithm.

```mathematica
In[1]:= p = NextPrime[10^30]; q = NextPrime[p + 10^15]; n = p q; IntegerLength[n]
Out[1]= 61

In[2]:= FactorInteger`Fermat[n]
Out[2]= {{1000000000000000000000000000057, 1}, {1000000000000001000000000000083, 1}}
```

That is a 61-digit semiprime with two 31-digit factors, and it is split in well
under a millisecond. Rho, ECM and the quadratic sieve would all need far
longer, because none of them can make use of the two factors being close.

---

## 1. The mathematical idea

### Fermat (1643)

An odd \(n = pq\) is a difference of two squares:

\[
n = a^2 - b^2 = (a-b)(a+b), \qquad a = \tfrac{p+q}{2},\; b = \tfrac{q-p}{2}.
\]

Fermat's method scans \(a = \lceil\sqrt n\,\rceil, \lceil\sqrt n\,\rceil + 1,
\dots\) and stops when \(a^2 - n\) is a perfect square \(b^2\). Then
\(\gcd(a-b, n)\) is a factor. The number of steps is

\[
a - \sqrt n \;=\; \frac{(\sqrt q - \sqrt p)^2}{2} \;\approx\; \frac{(q-p)^2}{8\sqrt n}.
\]

So everything depends on the gap \(q - p\):

| gap \(q-p\) | steps |
|---|---|
| \(\lesssim n^{1/4}\) | 1 |
| \(n^{1/4}\cdot 2^{k}\) | about \(2^{2k-3}\) |
| comparable with \(\sqrt n\) (an unbalanced \(n\)) | about \(n^{1/2}\), worse than trial division |

### Lehman (1974)

Suppose \(p/q\) is close to a fraction \(u/v\) with small \(uv = k\). Then
\(vp\) and \(uq\) are close, and Fermat's method applied to \(4kn\) finds them:

\[
4kn = (vp + uq)^2 - (vp - uq)^2 .
\]

Lehman's theorem makes this complete. Let \(n\) have no prime factor up to
\(n^{1/3}\). Then for some \(k \le n^{1/3}\), some
\(a \in \bigl[\sqrt{4kn},\; \sqrt{4kn} + n^{1/6}/(4\sqrt k)\,\bigr]\) makes
\(a^2 - 4kn\) a square, and \(\gcd(a + b, n)\) is a proper factor. Trial
division to \(n^{1/3}\), followed by these short Fermat windows, factors any
\(n\) in \(O(n^{1/3})\) steps. Every Dirichlet approximation \(u/v\) of \(p/q\)
is covered by some \(k\).

### References

1. R. S. Lehman, "Factoring large integers," *Mathematics of Computation*
   **28** (1974) 637–646.
2. D. E. Knuth, *The Art of Computer Programming*, Vol. 2, 3rd ed., 1997 —
   §4.5.4, Algorithm C (Fermat's method with residue sieves).
3. J. McKee, "Speeding Fermat's factoring method," *Mathematics of
   Computation* **68** (1999) 1729–1737.
4. W. B. Hart, "A one line factoring algorithm," *Journal of the Australian
   Mathematical Society* **92** (2012) 61–69 — a Lehman variant (not
   implemented here).

---

## 2. How Mathilda realises it

Both methods live in `fi_fermat.c` (bignums) and `fi_u64.c` (machine words).
The dispatch in `fi_core.c` picks the word version when it can be used.

**Fermat on bignums** (`fermat_scan`). Most values of \(a\) can be ruled out
without bignum arithmetic. A square is a quadratic residue modulo every \(m\),
so for nine moduli \(m \in \{64, 63, 65, 11, 17, 19, 23, 29, 31\}\) a table
marks the residues \(a \bmod m\) for which \(a^2 - n\) can be a square. Each
step advances nine word counters and does nine table lookups. Only about 1
candidate in 3000 passes every table and pays for an `mpz_perfect_square_p`
test.

**Fermat on words** (`fi_u64_fermat`, \(n < 2^{62}\)). This version keeps
\(b^2 = a^2 - n\) as a running sum (\(b^2 \mathrel{+}= 2a+1\)) and tests it for
squareness with residue filters modulo 64, 63, 65 and 11.

**Lehman on words** (`fi_u64_lehman`, \(n < 2^{42}\), so that \(4kn < 2^{63}\)).
This is the textbook algorithm. Trial division to \(n^{1/3}\) uses
precomputed inverses, so each prime costs one multiply. Then every
\(k \le n^{1/3}\) gets its window of \(a\) values. Because it runs to the end,
it is **complete**: it always finds a factor of a composite \(n < 2^{42}\).

**Lehman on bignums** (`fi_split_lehman`). Trial division runs to
\(\min(n^{1/3}, 2^{20})\). Then the method runs the Fermat scan on \(4kn\) for
\(k = 1, 2, \dots\), with Lehman's window for each \(k\), until the
`"Iterations"` budget is used up. Above \(2^{42}\) it is therefore a budgeted
"multiplier Fermat", not a complete algorithm.

**Suboptions.**

| Method | Suboption | Default | Meaning |
|---|---|---|---|
| `"Fermat"` | `"Iterations"` | \(2^{26}\) on bignums, \(2^{24}\) on words | values of \(a\) tried |
| `"Lehman"` | `"Iterations"` | \(2^{24}\) | total budget across the \(k\) windows (bignums only; the word version always runs to completion) |

**Inside `Automatic`.** After the rho and micro-ECM tiers, the automatic
pipeline runs a cheap Fermat probe of \(2^{14}\) steps (`AUTO_FERMAT_ITERS`).
That catches any two factors within about \(2^{7} n^{1/4}\) of each other.
Rho, ECM and QS are all at their worst on exactly this construction, and the
probe costs a fraction of a millisecond. Lehman is not used by `Automatic`.

---

## 3. Invoking it

```mathematica
In[3]:= FactorInteger[n, Method -> "Fermat"]
Out[3]= {{1000000000000000000000000000057, 1}, {1000000000000001000000000000083, 1}}

In[4]:= FactorInteger[n, Method -> {"Fermat", "Iterations" -> 10}]
Out[4]= {{1000000000000000000000000000057, 1}, {1000000000000001000000000000083, 1}}

In[5]:= FactorInteger`Lehman[1000003 3000017]
Out[5]= {{1000003, 1}, {3000017, 1}}
```

Ten steps are enough for `n`, because its factors are only \(10^{15}\) apart:
\(10^{30}/(8\cdot10^{30}) < 1\) step.

---

## 4. Where it shines

### Fermat: the gap decides everything

Keep \(p\) and make the gap a thousand times wider, \(10^{18}\) instead of
\(10^{15}\). The step count grows by a factor of a million, to about
\(10^{36}/(8\cdot10^{30}) \approx 1.25\times10^{5}\). That is still a few tens
of milliseconds:

```mathematica
In[6]:= q2 = NextPrime[p + 10^18]; n2 = p q2; FactorInteger`Fermat[n2]
Out[6]= {{1000000000000000000000000000057, 1}, {1000000000001000000000000000203, 1}}
```

The automatic pipeline gets the close-factor case almost for free with its
Fermat probe. Here the two factors are 50-digit primes \(10^{26}\) apart.
That is about \(30\,n^{1/4}\), so Fermat needs about
\(10^{52}/(8\cdot10^{49}) \approx 125\) steps, well inside the probe's
\(2^{14}\):

```mathematica
In[7]:= P = NextPrime[10^49]; Q = NextPrime[P + 10^26]; FactorInteger[P Q]
Out[7]= {{10000000000000000000000000000000000000000000000009, 1}, {10000000000000000000000100000000000000000000000081, 1}}
```

### Lehman: complete on machine words

Below \(2^{42}\), Lehman always finds a factor. As a check, here are the last
3001 integers up to \(2^{41}\). Every one is factored completely, and every
base in the answer is prime:

```mathematica
In[8]:= okQ[m_, f_] := Times @@ (Power @@@ f) == m && AllTrue[f[[All, 1]], PrimeQ]

In[9]:= AllTrue[Range[2^41 - 3000, 2^41], okQ[#, FactorInteger`Lehman[#]] &]
Out[9]= True
```

Lehman does not need the factors to be balanced. Fermat is at its worst on
\(1031 \times 2147483659\), because the factors differ by a factor of two
million. Fermat gives up after its \(2^{24}\) word steps. Lehman finds 1031
during its trial-division stage:

```mathematica
In[10]:= FactorInteger`Fermat[1031 2147483659]
Out[10]= {{2214055652429, 1}}

In[11]:= FactorInteger`Lehman[1031 2147483659]
Out[11]= {{1031, 1}, {2147483659, 1}}
```

### Lehman: a ratio close to a small fraction

Take \(r \approx 3p\). The two factors are very unbalanced for Fermat: the gap
is about \(2p \approx \sqrt n\). For Lehman this is the easy case \(k = 3\),
\(u/v = 1/3\), because \(3p - r\) is tiny:

```mathematica
In[12]:= p3 = NextPrime[10^20]; r = NextPrime[3 p3]; n3 = p3 r; {IntegerLength[n3], 3 p3 - r}
Out[12]= {41, -4}

In[13]:= FactorInteger[n3, Method -> {"Fermat", "Iterations" -> 10^5}]
Out[13]= {{30000000000000000023800000000000000004719, 1}}

In[14]:= FactorInteger`Lehman[n3]
Out[14]= {{100000000000000000039, 1}, {300000000000000000121, 1}}
```

The same works for any small \(u/v\). With \(r \approx 5p/3\), the multiplier
is \(k = 15\):

```mathematica
In[15]:= r5 = NextPrime[Floor[5 p3/3]]; n5 = p3 r5; {5 p3 - 3 r5, FactorInteger`Lehman[n5]}
Out[15]= {-28, {{100000000000000000039, 1}, {166666666666666666741, 1}}}
```

---

## 5. Where it struggles

**Fermat on an ordinary semiprime.** Push the gap to \(10^{20}\) at the same
61-digit size. That needs about \(1.25\times10^{9}\) steps, twenty times the
default budget of \(2^{26}\). Here the budget is cut to \(10^6\) so the failure
comes back at once. With the default budget it fails too, after a couple of
seconds:

```mathematica
In[16]:= q3 = NextPrime[p + 10^20]; FactorInteger[p q3, Method -> {"Fermat", "Iterations" -> 10^6}]
Out[16]= {{1000000000100000000000000000204000000005700000000000000008379, 1}}
```

A random balanced semiprime has a gap comparable with its factors, so Fermat
has essentially no chance on it. The same 61-digit size is a 0.5–1 s job for
the [quadratic sieve](quadratic-sieve.md):

```mathematica
In[17]:= FactorInteger[p NextPrime[3 10^30], Method -> {"Fermat", "Iterations" -> 10^5}]
Out[17]= {{3000000000000000000000000000262000000000000000000000000005187, 1}}
```

**Lehman above \(2^{42}\) is only as good as its budget.** The bignum version
trial-divides only to \(2^{20}\), the size of the shared prime table. A 21-bit
factor of a 61-bit number is above that limit and below \(n^{1/3}\), so the
completeness argument breaks down. With the default \(2^{24}\) budget, the
Fermat windows run for about 3 s and find nothing. Here the budget is cut to
\(10^5\), so the same failure comes back at once:

```mathematica
In[18]:= FactorInteger[NextPrime[2^20] NextPrime[2^41], Method -> {"Lehman", "Iterations" -> 10^5}]
Out[18]= {{2305858402404794557, 1}}
```

!!! warning "Lehman on large bignums never gets past \(k = 1\) (v0.459)"
    In the bignum Lehman, the window for \(k = 1\) is \(n^{1/6}/4\) values of
    \(a\), and it is **not** capped by `"Iterations"`. Beyond about 45 digits,
    that one window is already larger than the \(2^{24}\) budget. The method
    scans all of it, finds the budget spent, and stops before reaching
    \(k = 2, 3, \dots\). So the small-ratio property shown in §4 holds only up
    to about 45 digits. The example below is missed at 51 digits, and
    `"Iterations"` has no effect.

At 51 digits, \(r \approx 3p\) with \(3p - r = -28\) is the same easy
\(k = 3\) case as before, but the bignum Lehman returns it unfactored after
about 0.4 s:

```mathematica
In[19]:= p4 = NextPrime[10^25]; r4 = NextPrime[3 p4]; FactorInteger[p4 r4, Method -> {"Lehman", "Iterations" -> 100}]
Out[19]= {{300000000000000000000001060000000000000000000000871, 1}}
```

Until this is fixed, factor such an \(n\) with the automatic pipeline, which
reaches the quadratic sieve at this size:

```mathematica
In[20]:= FactorInteger[p4 r4]
Out[20]= {{10000000000000000000000013, 1}, {30000000000000000000000067, 1}}
```

---

## 6. Limitations

- Fermat costs about \((q-p)^2/(8\sqrt n)\) steps. It is only useful when the
  two factors agree in roughly their top half of digits. The automatic
  pipeline's \(2^{14}\)-step probe covers the cases worth covering.
- Fermat finds a *pair* \(a \pm b\). When \(n\) has more than two prime factors,
  it finds the split into two divisors closest to \(\sqrt n\), and the driver
  then factors each part again with Fermat.
- Lehman is complete only below \(2^{42}\) (word arithmetic, \(4kn < 2^{63}\)).
  Above that it is a budgeted search. Its trial-division stage stops at
  \(2^{20}\), and on large inputs its budget is spent at \(k = 1\) (see the
  warning above).
- Neither method is used for general factoring at scale. The word engine
  prefers rho, micro-ECM and SQUFOF, and above \(2^{64}\) the sieve and ECM
  dominate.

---

## See also

- [FactorInteger methods](index.md): the automatic pipeline and its Fermat
  probe.
- [SQUFOF](squfof.md): the other \(O(n^{1/4})\)-flavoured method on machine
  words.
- [Quadratic sieve](quadratic-sieve.md): the general-purpose replacement for
  balanced semiprimes.
- [Blake rational base descent](blake-rational-base-descent.md): another
  method for structured factors, here factors close to a power of a rational.
- `?FactorInteger`Fermat` and `?FactorInteger`Lehman` in the REPL.
