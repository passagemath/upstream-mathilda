# CFRAC and Dixon: congruences of squares

`Method -> "CFRAC"` (alias `"ContinuedFraction"`) and `Method -> "Dixon"` are
Mathilda's two older **congruence-of-squares** factoring methods. Neither one
looks for a factor directly. Both collect many small congruences
\(X^2 \equiv (\text{product of small primes}) \pmod N\), and then combine them
with linear algebra over \(\mathrm{GF}(2)\) until the right-hand side is a perfect
square. The quadratic sieve ([quadratic-sieve.md](quadratic-sieve.md)) and the
number field sieve are built on the same idea. What separates the methods is how
the congruences are produced.

```mathematica
In[1]:= n = 84923
Out[1]= 84923

In[2]:= Mod[{513^2, 537^2}, n]
Out[2]= {8400, 33600}

In[3]:= FactorInteger[{8400, 33600}]
Out[3]= {{{2, 4}, {3, 1}, {5, 2}, {7, 1}}, {{2, 6}, {3, 1}, {5, 2}, {7, 1}}}
```

Neither residue is a square. Their product is
\(2^{10}\cdot 3^2\cdot 5^4\cdot 7^2 = (2^5\cdot 3\cdot 5^2\cdot 7)^2\), because each
exponent in the sum is even. So \((513\cdot 537)^2 \equiv 16800^2 \pmod n\), and a
gcd splits \(n\):

```mathematica
In[4]:= Mod[(513 537)^2 - (2^5 3 5^2 7)^2, n]
Out[4]= 0

In[5]:= GCD[513 537 - 2^5 3 5^2 7, n]
Out[5]= 163
```

That is the whole method in miniature (\(84923 = 163 \cdot 521\)). The rest of
this page covers how to find such residues at scale.

---

## 1. The mathematical idea

### Congruence of squares

If \(X^2 \equiv Y^2 \pmod N\) with \(X \not\equiv \pm Y\), then \(N \mid (X-Y)(X+Y)\)
but \(N\) divides neither factor, so \(\gcd(X-Y, N)\) is a proper divisor. For
\(N = pq\), a *random* pair with \(X^2\equiv Y^2\) is useful with probability
\(1/2\). Each further prime factor of \(N\) only raises that probability.

### Factor base and smoothness

Fix a **factor base** \(\mathcal{F} = \{-1, p_1, \dots, p_B\}\) of small primes.
A **relation** is a congruence

\[
X_i^2 \equiv (-1)^{e_{i0}} \prod_{j=1}^{B} p_j^{\,e_{ij}} \pmod N ,
\]

that is, a residue that is \(p_B\)-**smooth**. Write its exponent vector mod 2 as
\(v_i \in \mathrm{GF}(2)^{B+1}\). A subset \(S\) of relations with
\(\sum_{i\in S} v_i = 0\) gives a square on the right:

\[
X = \prod_{i\in S} X_i, \qquad Y = \prod_j p_j^{\frac12\sum_{i\in S} e_{ij}}, \qquad X^2 \equiv Y^2 \pmod N .
\]

Once there are more relations than columns, such a subset must exist; it is a
vector in the left nullspace of the \(R\times(B+1)\) matrix \((v_i)\). In the toy
example, the two vectors \((4,1,2,1)\) and \((6,1,2,1)\) are both
\((0,1,0,1)\) mod 2, so their sum is zero.

### Where the residues come from

How fast a method runs depends almost entirely on how *small* its residues are.
The probability that a random integer \(z\) is \(y\)-smooth is about \(u^{-u}\),
where \(u = \log z/\log y\) (Dickman's \(\rho\)). Halving \(\log z\) is therefore
worth far more than any constant-factor speed-up.

- **Dixon** (1981) takes random \(x\) and keeps \(z = x^2 \bmod N\) when it is smooth.
  In the textbook form \(z\) is as large as \(N\). Dixon proved the first rigorous
  subexponential bound, \(L_N[1/2, 2\sqrt2]\), but the method is slow in practice.
- **CFRAC** (Lehmer–Powers 1931, Morrison–Brillhart 1975) reads the residues off
  the continued fraction of \(\sqrt{kN}\). With convergents \(A_{i}/B_{i}\), the
  recurrence

  \[
  P_i = q_{i-1}Q_{i-1} - P_{i-1},\quad Q_i = \frac{kN - P_i^2}{Q_{i-1}},\quad
  A_{i-1}^2 \equiv (-1)^i\, Q_i \pmod N
  \]

  gives residues with \(0 < Q_i < 2\sqrt{kN}\), which are half the size of
  Dixon's.

The square-root bound is easy to see on the toy number:

```mathematica
In[6]:= cf = ContinuedFraction[Sqrt[n], 10]
Out[6]= {291, 2, 2, 2, 5, 1, 1, 2, 4, 1}

In[7]:= A = Numerator[Table[FromContinuedFraction[Take[cf, i]], {i, 10}]]
Out[7]= {291, 583, 1457, 3497, 18942, 22439, 41381, 105201, 462185, 567386}

In[8]:= Mod[A^2, n, -Floor[n/2]]
Out[8]= {-242, 197, -226, 97, -311, 254, -211, 118, -283, 289}
```

Every least absolute residue is below \(2\sqrt{84923}\approx 583\), against a
modulus of 84923. The signs alternate, as \((-1)^i\) predicts, which is why
\(-1\) is a column of the factor base. The last entry is already a square,
\(289 = 17^2\), so a single relation gives a factor:

```mathematica
In[9]:= GCD[567386 - 17, n]
Out[9]= 521
```

### References

1. M. A. Morrison and J. Brillhart, "A method of factoring and the factorization
   of \(F_7\)," *Mathematics of Computation* **29** (1975) 183–205.
2. J. D. Dixon, "Asymptotically fast factorization of integers," *Mathematics of
   Computation* **36** (1981) 255–260.
3. C. Pomerance, "Analysis and comparison of some integer factoring algorithms,"
   in *Computational Methods in Number Theory*, Math. Centrum Tract 154 (1982)
   89–139 — early abort and the \(L[1/2]\) analyses.
4. D. E. Knuth, *The Art of Computer Programming*, Vol. 2, 3rd ed., §4.5.4 — the
   Knuth–Schroeppel multiplier.
5. R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*,
   2nd ed., Springer 2005, §6.1–6.2.

---

## 2. How Mathilda realises it

Both methods live in `src/numbertheory/factorinteger/` and share all their
back-end code:

| File | Role |
|---|---|
| `fi_cfrac.c` | the continued-fraction walk, the Knuth–Schroeppel multiplier `fi_ks_multiplier` (shared with the QS), and the factor-base size table |
| `fi_dixon.c` | Dixon's candidate generator |
| `fi_smooth.c` | threaded trial division of candidates over the factor base, with early abort |
| `fi_gf2.c` | relation store, large-prime pairing (union–find over a large-prime graph), singleton pruning, dense GF(2) elimination or block Lanczos, square roots |

**CFRAC.**

- **Multiplier.** \(k\) maximises the Knuth–Schroeppel score
  \(-\tfrac12\log k + \sum_p g(p,kN)\log p\) over the 46 square-free
  \(k \le 73\). Here \(g = 2/(p-1)\) when \((kN/p)=1\), and \(g = 1/p\) when
  \(p \mid k\).
- **Factor base.** 2, plus the primes with \((kN/p)\ne -1\) (only these can divide
  some \(Q_i\)). The size is interpolated by digit count: 40 primes at 10 digits,
  260 at 30, 650 at 40, 1500 at 50, 3200 at 60, 5800 at 70.
- **Trial division** runs in word-sized batches of base primes, one bignum
  remainder per batch. It uses **early abort**: once a quarter of the base has
  been tried, a cofactor that still has more than \(0.7\,\log_2 Q_{\max} + 8\) bits
  is dropped.
- **Single large prime.** A leftover cofactor up to \(128\,p_{\max}\) is kept as a
  *partial* relation. Two partials with the same large prime combine into a full
  relation; in `fi_gf2.c` these are 2-cycles of the large-prime graph.
- **Degenerate expansions.** The expansion stops when \(Q_i = 1\) closes the
  period. If three successive solves return only trivial dependencies, CFRAC
  restarts with the next multiplier from a fixed list and a fresh relation store.
- **Threads.** The walk runs in batches of 64 terms, doubling up to 1024 per
  thread. Worker 0 extends the expansion into the next batch while all workers
  trial-divide the current one, and hits are replayed in term order. The relations
  and the divisor are therefore the same for every thread count.

**Dixon.** Textbook Dixon is hopeless beyond toy sizes, so Mathilda's variant draws
each square as \(x = \lfloor\sqrt{kN}\rfloor + j\), with random \(k < 2^{16}\) and
random \(j < 2^{20}\). The residue \(z = x^2 - kN \approx 2j\sqrt{kN}\) is then
square-root sized. That is still about \(\log_{10}(2^{21}\cdot 2^{8}) \approx 9\)
digits larger than a CFRAC residue, which is why Dixon falls behind so quickly.
The choices stay random and nothing is sieved (sieving these residues is the
quadratic sieve). Because \(k\) varies, every prime can divide some \(z\), so the
base is *all* primes, 1.5× CFRAC's base size. Draws are made in batches of 2048
per thread and replayed in draw order. Afterwards the random generator, which is
seeded from \(N\), is left exactly where the serial loop would leave it.

**Linear algebra** (shared with the QS). Singleton columns are pruned repeatedly,
since a relation holding the only odd occurrence of a prime can never be in a
dependency. Below 1000 surviving columns, dense bit-packed Gaussian elimination is
used; above that, block Lanczos (Montgomery 1995). Each dependency yields \(X, Y\)
and a gcd.

### Suboptions

| Method | Suboption | Default | Meaning |
|---|---|---|---|
| `"CFRAC"` | `"Iterations"` | \(2^{32}\) | continued-fraction terms allowed before giving up |
| `"Dixon"` | `"Iterations"` | \(2^{30}\) | squares drawn before giving up |

Neither method has a word-size path. Unlike rho or SQUFOF, a 64-bit input still
runs the full relation machinery, which is useful for teaching but never the fast
choice.

---

## 3. Invoking the methods

There are three equivalent spellings for each method:
`FactorInteger[n, Method -> "CFRAC"]`, `FactorInteger[n, Method -> {"CFRAC",
"Iterations" -> m}]` and ``FactorInteger`CFRAC[n, "Iterations" -> m]``. The alias
string `"ContinuedFraction"` is accepted in `Method`, but no
``FactorInteger`ContinuedFraction`` symbol exists.

The test inputs below are balanced semiprimes built reproducibly. A helper checks
that each answer multiplies back to \(n\) and has only prime bases:

```mathematica
In[10]:= sp[d_] := NextPrime[3 10^(d/2 - 1)] NextPrime[7 10^(d/2 - 1)]
Out[10]= Null

In[11]:= check[n_, f_] := Times @@ (Power @@@ f) == n && AllTrue[f[[All, 1]], PrimeQ]
Out[11]= Null
```

---

## 4. Where they work

**Dixon on small inputs.** At 12 digits (40 bits) Dixon's base has about 80 primes, so it needs
about 130 relations, which take a few tens of milliseconds. At 20
digits it is already well into the sub-second range:

```mathematica
In[12]:= FactorInteger[sp[12], Method -> "Dixon"]
Out[12]= {{300007, 1}, {700001, 1}}

In[13]:= f = FactorInteger`Dixon[sp[20]]
Out[13]= {{3000000019, 1}, {7000000001, 1}}

In[14]:= check[sp[20], f]
Out[14]= True
```

**CFRAC to about 40 digits.** On the same 20-digit input CFRAC needs only
milliseconds, and it reaches 40 digits (133 bits) in about half a second:

```mathematica
In[15]:= FactorInteger`CFRAC[sp[20]]
Out[15]= {{3000000019, 1}, {7000000001, 1}}

In[16]:= f = FactorInteger[sp[40], Method -> "ContinuedFraction"]
Out[16]= {{30000000000000000041, 1}, {70000000000000000013, 1}}

In[17]:= check[sp[40], f]
Out[17]= True
```

Measured wall times on the balanced semiprimes `sp[d]` (16 threads, on a
heavily loaded machine, so treat them as upper bounds):

| digits | `"Dixon"` | `"CFRAC"` | `"QuadraticSieve"` |
|---|---|---|---|
| 12 | ≈ 50 ms | ≈ 2 ms | < 1 ms |
| 16 | ≈ 0.17 s | — | — |
| 20 | ≈ 0.65 s | ≈ 6 ms | ≈ 2 ms |
| 24 | ≈ 2.8 s | ≈ 11 ms | — |
| 28 | > 20 s CPU (aborted) | — | — |
| 30 | — | ≈ 40 ms | ≈ 7 ms |
| 40 | — | ≈ 0.5 s | ≈ 40 ms |
| 50 | — | 2–12 s | ≈ 0.1 s |
| 60 | — | > 10 s CPU (aborted) | ≈ 0.5 s |

---

## 5. Where they struggle

**Budgets.** With a small `"Iterations"` the relation hunt stops early. The
composite is then returned unfactored, with `FactorInteger::nofac` on stderr:

```mathematica
In[18]:= Quiet[FactorInteger`CFRAC[sp[30], "Iterations" -> 500]]
Out[18]= {{210000000000077600000000004539, 1}}

In[19]:= Quiet[FactorInteger`Dixon[sp[20], "Iterations" -> 10000]]
Out[19]= {{21000000136000000019, 1}}
```

**Growth.** Both methods run in \(L_N[1/2, c]\) time,
\(\exp\big(c\sqrt{\log N\log\log N}\big)\), with a larger \(c\) than the quadratic
sieve. The practical gap is bigger still, for two reasons:

- **Dixon's residues are ~9 digits longer than CFRAC's.** Every candidate also
  costs a full trial division, since nothing is sieved. At 40 digits Dixon does not
  finish in any reasonable time.
- **CFRAC's residues are about \(\sqrt{N}\),** the same order as the QS's
  \(M\sqrt{N/2}\). But each \(Q_i\) must be found by trial division. The QS instead
  *sieves*: the polynomial structure lets it add \(\log p\) at every \(p\)-th
  position and trial-divide only the survivors. That changes the cost per candidate
  from about \(B\) operations to roughly \(\sum_{p}1/p \approx \log\log p_B\). By 60
  digits CFRAC cannot finish inside a 10-second CPU budget, while the QS needs
  half a second.

`TimeConstrained` shows this cheaply. Its budget is **process CPU time**, and
these methods use every thread, so 10 seconds of budget is well under a second of
wall time on a 16-thread machine:

```mathematica
In[20]:= TimeConstrained[FactorInteger`Dixon[sp[40]], 10, $Aborted]
Out[20]= $Aborted

In[21]:= TimeConstrained[FactorInteger`CFRAC[sp[60]], 10, $Aborted]
Out[21]= $Aborted

In[22]:= FactorInteger`QuadraticSieve[sp[40]] == FactorInteger`CFRAC[sp[40]]
Out[22]= True

In[23]:= FactorInteger`QuadraticSieve[sp[60]]
Out[23]= {{300000000000000000000000000007, 1}, {700000000000000000000000000033, 1}}
```

**How `Automatic` avoids them.** The `Automatic` cascade never calls CFRAC or
Dixon. Its general-purpose finisher is the quadratic sieve. That is taken once the
ECM effort reaches 2/7 of the composite's digit count, for composites up to 100
digits. CFRAC and Dixon exist so that you can study, compare and benchmark the
classical methods; see [index.md](index.md) for the full cascade.

---

## 6. Limitations

- Neither method suits inputs much beyond 50 digits (CFRAC) or 25 digits (Dixon).
  Use `Automatic` or `"QuadraticSieve"`.
- They use single large primes only. The double-large-prime variation in the
  shared `fi_gf2.c` is fed only by the quadratic sieve. Block Lanczos does take
  over for CFRAC from about 45 digits, where the base passes 1000 primes.
- There is no word-size fast path: a 40-bit input pays for the full relation
  machinery.
- Mathilda's Dixon is a deliberate hybrid (random \(k\) and \(j\) near
  \(\sqrt{kN}\)), not the textbook uniform-\(x\) algorithm, so Dixon's rigorous
  complexity bound does not apply to it as implemented.
- The only tuning knob is `"Iterations"`. The factor-base size, the large-prime
  bound and the early-abort threshold are fixed tables.

---

## See also

- [The quadratic sieve](quadratic-sieve.md): the same relations, found by sieving.
- [Overview of the FactorInteger methods](index.md) and the
  [performance comparison](performance.md).
- Special-purpose methods: [trial division](trial-division.md),
  [Fermat and Lehman](fermat-lehman.md), [Pollard rho](pollard-rho.md),
  [SQUFOF](squfof.md), [p−1 and p+1](pminus1-pplus1.md), [ECM](ecm.md),
  [Blake's rational base descent](blake-rational-base-descent.md).
