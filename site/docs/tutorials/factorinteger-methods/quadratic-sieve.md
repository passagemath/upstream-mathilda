# The quadratic sieve (SIQS)

`Method -> "QuadraticSieve"` (aliases `"QS"`, `"SIQS"`, `"MPQS"`) is Mathilda's
general-purpose factoring engine. It is a native **self-initialising quadratic
sieve**, and it is the fastest method here for composites of roughly 30–100
digits with no small factors. The `Automatic` cascade hands such composites to it
once the cheap special-purpose tiers have failed.
The shoot-out in `benchmarks/98-factorinteger-shootout/REPORT.md` records 60
digits in about 0.6 s, 70 in about 3 s, 80 in about 30 s and 90 in under 5
minutes, on an 8-core i9. That is 5–10× faster than PARI/GP and FLINT from 50
digits up.

```mathematica
In[1]:= sp[d_] := NextPrime[3 10^(d/2 - 1)] NextPrime[7 10^(d/2 - 1)]
Out[1]= Null

In[2]:= FactorInteger`QuadraticSieve[sp[40]]
Out[2]= {{30000000000000000041, 1}, {70000000000000000013, 1}}
```

---

## 1. The mathematical idea

### Relations by sieving

Like [CFRAC and Dixon](cfrac-dixon.md), the quadratic sieve (Pomerance, 1981)
collects congruences \(X^2 \equiv \prod p_j^{e_j} \pmod N\) over a factor base. It
then combines them by \(\mathrm{GF}(2)\) linear algebra into \(X^2\equiv Y^2\), and
\(\gcd(X-Y,N)\) splits \(N\). What is new is where the residues come from: the
values of a *polynomial*

\[
Q(x) = (x + \lfloor\sqrt N\rfloor)^2 - N \approx 2x\sqrt N ,
\]

which are square-root sized near the origin. Because \(Q\) is a polynomial,
\(p \mid Q(x)\) exactly when \(x\) lies in one of two residue classes mod \(p\), given
by the roots \(\pm\sqrt N \bmod p\). Those classes exist only for primes with
\((N/p) = 1\), so only those primes go in the factor base. Instead of
trial-dividing every candidate, the sieve adds \(\log p\) to every \(p\)-th entry
of an array, starting at each root. Positions whose total comes close to
\(\log|Q(x)|\) are almost certainly smooth, and only they are trial-divided. The
cost per position is about \(\sum_{p\le B} 2/p \approx 2\log\log B\) additions,
instead of \(B\) divisions.

Here is a toy instance with \(N = 15347\):

```mathematica
In[3]:= n = 15347
Out[3]= 15347

In[4]:= fb = Select[Prime[Range[10]], # == 2 || JacobiSymbol[n, #] == 1 &]
Out[4]= {2, 17, 23, 29}

In[5]:= smoothQ[v_] := Max[First /@ FactorInteger[v]] <= Max[fb]
Out[5]= Null

In[6]:= rels = Select[Table[{x, x^2 - n}, {x, 124, 244}], smoothQ[Last[#]] &]
Out[6]= {{124, 29}, {126, 529}, {127, 782}, {195, 22678}}

In[7]:= FactorInteger[rels[[All, 2]]]
Out[7]= {{{29, 1}}, {{23, 2}}, {{2, 1}, {17, 1}, {23, 1}}, {{2, 1}, {17, 1}, {23, 1}, {29, 1}}}
```

The relation at \(x = 126\) is already a square, \(23^2\). The product of the
other three is \(2^2\,17^2\,23^2\,29^2\). Either combination factors \(N\):

```mathematica
In[8]:= GCD[126 - 23, n]
Out[8]= 103

In[9]:= X = Mod[124 127 195, n]; Y = 2 17 23 29; GCD[X - Y, n]
Out[9]= 103
```

The sieve finds such positions without division. For \(p = 17\), the square root
of \(N\) mod 17 is \(\pm 8\), so exactly the \(x \equiv \pm 8 \pmod{17}\) are hit:

```mathematica
In[10]:= PowerMod[n, 1/2, 17]
Out[10]= 8

In[11]:= Select[Range[124, 244], Mod[#^2 - n, 17] == 0 &]
Out[11]= {127, 128, 144, 145, 161, 162, 178, 179, 195, 196, 212, 213, 229, 230}
```

These are the classes \(127 \equiv 8\) and \(128 \equiv -8 \pmod{17}\), repeating with
period 17.

### Many polynomials: MPQS and SIQS

With one polynomial, \(|Q(x)|\) grows linearly with the interval, and smoothness
decays quickly. The **multiple-polynomial** QS (Montgomery; Silverman 1987) uses
many polynomials, each sieved over a short interval \([-M, M)\):

\[
g(x) = \frac{(Ax+B)^2 - kN}{A}, \qquad B^2 \equiv kN \pmod A ,
\]

so \((Ax+B)^2 \equiv A\,g(x) \pmod N\). The factor \(A\) is recorded with the
relation. Choosing \(A \approx \sqrt{2kN}/M\) keeps \(|g(x)| \lesssim M\sqrt{kN/2}\)
over the whole interval.

The **self-initialising** QS (Alford–Pomerance; Contini 1997) takes
\(A = q_1\cdots q_s\), a product of factor-base primes. Then \(B\) has \(2^{s}\) square
roots \(B = \pm B_1 \pm \dots \pm B_s\), where each \(B_l\) is a CRT component.
Stepping through the \(2^{s-1}\) sign patterns in **Gray-code** order changes one
sign at a time, so every sieve root moves by a precomputed \(\pm 2B_lA^{-1} \bmod p\).
A new polynomial costs one addition per prime instead of a modular inversion.
The construction, for the 30-digit `sp[30]` with \(M = 2^{15}\) and \(s = 3\):

```mathematica
In[12]:= N30 = sp[30]; M = 32768;
Out[12]= Null

In[13]:= qs = Take[Select[Prime[Range[390, 450]], JacobiSymbol[N30, #] == 1 &], 3]
Out[13]= {2707, 2711, 2719}

In[14]:= A = Times @@ qs; B = ChineseRemainder[PowerMod[N30, 1/2, #] & /@ qs, qs]
Out[14]= 9673830554

In[15]:= Mod[B^2 - N30, A]
Out[15]= 0

In[16]:= g[x_] := ((A x + B)^2 - N30)/A
Out[16]= Null

In[17]:= N[Log10[Max[Abs[g /@ {-M, 0, M}]]]]
Out[17]= 19.0375

In[18]:= N[Log10[(Floor[Sqrt[N30]] + 8 M)^2 - N30]]
Out[18]= 20.3807
```

The four polynomials of this family together sieve \(4\times 2M = 8M\) positions,
with values below \(10^{19}\). A single polynomial over the same \(8M\) positions
reaches \(10^{20.4}\), and every factor of 10 in the size costs smoothness
probability.

### The Knuth–Schroeppel multiplier

The method sieves \(kN\) rather than \(N\), for a small square-free \(k\) chosen so
that many small primes are quadratic residues of \(kN\). Small primes then divide
\(g(x)\) more often. The score maximised is

\[
f(k) = -\tfrac12\log k + \sum_{p} g(p,kN)\log p, \qquad
g = \begin{cases}2/(p-1) & (kN/p)=1\\ 1/p & p\mid k\\ 0 & \text{else}\end{cases}
\]

plus a 2-adic term that depends on \(kN \bmod 8\). Below is a line-for-line
transcription of `fi_ks_multiplier` (`fi_cfrac.c`), with the same 46 candidates
and the odd primes up to 1229:

```mathematica
In[19]:= ks[n_, k_] := Module[{s = -Log[k]/2, m8 = Mod[n k, 8]}, s += Which[m8 == 1, 2 Log[2], m8 == 5, Log[2], True, Log[2]/2]; Do[If[Mod[k, p] == 0, s += Log[p]/p, If[JacobiSymbol[n k, p] == 1, s += 2 Log[p]/(p - 1)]], {p, Prime[Range[2, 201]]}]; N[s]]
Out[19]= Null

In[20]:= K = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73};
Out[20]= Null

In[21]:= Take[SortBy[Table[{k, ks[sp[40], k]}, {k, K}], -Last[#] &], 4]
Out[21]= {{13, 7.83238}, {1, 6.78744}, {10, 6.63672}, {37, 6.38366}}
```

For `sp[40]` the sieve therefore runs on \(13N\). Its score is about 1.04 higher than
that of \(k = 1\), measured in natural-log units of expected small-prime
contribution.

### Large primes and cycles

A position whose cofactor after the factor base is a single prime
\(L \le \texttt{lpmult}\cdot p_{\max}\) is kept as a **partial** relation. Two partials
with the same \(L\) multiply to a full relation, because \(L^2\) is a square. With
**double large primes**, a cofactor \(L_1L_2\) with both primes below the bound is
also kept. Partials then form a graph on the large primes: a one-prime partial is
an edge from a special vertex "1" to \(L\), a two-prime partial an edge \(L_1 - L_2\), and every **cycle**
multiplies out to a full relation. Union–find counts the independent cycles as the
edges arrive (`fi_gf2.c`), so the sieve knows when it has enough relations without
building them.

### Linear algebra

The relation matrix has a few tens of thousands of columns at 80–90 digits, and is
extremely sparse. **Block Lanczos** (Montgomery 1995) finds about 64 nullspace
vectors at once in \(O(wn/64)\) word operations, for \(w\) nonzero entries.
Each dependency splits \(N\) with probability \(\tfrac12\) or better.

### Complexity

The QS runs in \(L_N[\tfrac12, 1] = \exp\big((1+o(1))\sqrt{\ln N\ln\ln N}\big)\).
This is subexponential, and it depends only on the size of \(N\), not on the size
of its factors. That is both its strength (no factor is too large) and its
weakness (no factor is small enough to make it faster): see §4.

### References

1. C. Pomerance, "Analysis and comparison of some integer factoring algorithms,"
   in *Computational Methods in Number Theory*, Math. Centrum Tract 154 (1982)
   89–139; and "The quadratic sieve factoring algorithm," *EUROCRYPT '84*, LNCS
   209 (1985) 169–182.
2. R. D. Silverman, "The multiple polynomial quadratic sieve," *Mathematics of
   Computation* **48** (1987) 329–339.
3. S. P. Contini, *Factoring Integers with the Self-Initializing Quadratic Sieve*,
   M.Sc. thesis, University of Georgia (1997).
4. W. R. Alford and C. Pomerance, "Implementing the self-initializing quadratic
   sieve on a distributed network," in *Number Theoretic and Algebraic Methods in
   Computer Science* (1995) 163–174.
5. A. K. Lenstra and M. S. Manasse, "Factoring with two large primes,"
   *Mathematics of Computation* **63** (1994) 785–798.
6. P. L. Montgomery, "A block Lanczos algorithm for finding dependencies over
   GF(2)," *EUROCRYPT '95*, LNCS 921 (1995) 106–120.
7. D. E. Knuth, *The Art of Computer Programming*, Vol. 2, §4.5.4 (the
   Knuth–Schroeppel function).

---

## 2. How Mathilda realises it

All of the sieve is in `src/numbertheory/factorinteger/fi_qs.c`. The relation
store, cycle counting and dense elimination are in `fi_gf2.c`, block Lanczos in
`fi_lanczos.c`, and thread regions in `fi_par.c`. There is no external dependency.

| Stage | Implementation |
|---|---|
| Multiplier | `fi_ks_multiplier` (above), shared with CFRAC |
| Factor base | 2, then the primes with \((kN/p) \ne -1\), up to the table size. A base prime dividing \(N\) is returned at once |
| Choosing \(A\) | \(s-1\) random distinct base primes of about 11 bits from a window around \((\sqrt{2kN}/M)^{1/s}\), then the last \(q_s\) fitted to the residual. Keys of the used \(A\) are kept so no family repeats |
| Polynomials | \(2^{s-1}\) per \(A\), visited in Gray-code order. Medium-prime roots move in `next_poly`; bucket-prime roots move while their buckets are filled |
| Sieve | byte logarithms in **32 KB blocks** (one L1 data cache). Primes below 32 (below 64 once \(kN\) has 120 bits) are not sieved at all (small-prime variation), and the threshold is lowered to compensate |
| Medium primes | \(p < 2^{15}\): sieved block by block, both roots in one loop |
| Large primes | \(p \ge 2^{15}\): **bucket sieved**. Once per polynomial each hit is dropped into its block's bucket. The buckets double as the trial-division lists for candidates |
| Scan | the array starts at \(128 - T\), so a candidate is any byte with its top bit set, tested 8 bytes at a time. At most 256 candidates per block |
| Trial division | only the unsieved primes, the \(q_l\), the medium primes whose root matches the position, and the bucket primes recorded for it |
| Large primes | single large prime \(\le \texttt{lpmult}\cdot p_{\max}\). From 240-bit \(kN\) (about 72 digits), **double large primes**: a composite cofactor up to \(\texttt{lp}_{\max}^{1.8}\) is split by word micro-ECM (8 curves) or SQUFOF, and the threshold drops by \(0.7\log_2 \texttt{lp}_{\max}\) |
| Stopping rule | solve once the relations reach \(n_{fb}+65\), or the active columns \(+64\). A failed solve asks for 64 more |
| Linear algebra | singleton pruning, then dense Gaussian elimination below 1000 columns and block Lanczos above |
| Threads | from 150-bit \(kN\) (about 45 digits) each worker takes one \(A\) and sieves its whole family into a private relation log. The main thread replays the logs in worker order. Up to 16 threads |
| Below \(2^{64}\) | no sieve: word micro-ECM (256 curves), then SQUFOF |

### The parameter table

The parameters come from `QS_PARAMS`, indexed by the bit size of \(kN\), with
linear interpolation between rows:

- the factor-base size and the large-prime multiplier are interpolated
  linearly;
- the block count snaps to the nearer row;
- the threshold slack is interpolated, and is the number of bits below the full
  log at which a position is accepted.

The table was tuned on an i9-9880H (32 KB L1d, 256 KB L2) against balanced
semiprimes.

| bits of \(kN\) | ≈ digits | factor base | 32 KB blocks per side | interval \(2M\) | LP multiplier | slack (bits) |
|---|---|---|---|---|---|---|
| 64 | 19 | 60 | 1 | 64 KB | 30 | 3 |
| 100 | 30 | 200 | 1 | 64 KB | 40 | 4 |
| 120 | 36 | 400 | 1 | 64 KB | 60 | 5 |
| 133 | 40 | 600 | 1 | 64 KB | 100 | 6 |
| 166 | 50 | 2 500 | 2 | 128 KB | 250 | 8 |
| 200 | 60 | 6 000 | 2 | 128 KB | 250 | 8 |
| 233 | 70 | 18 000 | 3 | 192 KB | 250 | 8 |
| 266 | 80 | 34 000 | 4 | 256 KB | 300 | 9 |
| 300 | 90 | 55 000 | 6 | 384 KB | 350 | 10 |
| 333 | 100 | 80 000 | 8 | 512 KB | 400 | 10 |
| 366 | 110 | 110 000 | 10 | 640 KB | 400 | 11 |

Beyond the last row the parameters are held constant. The method still runs, but
the table stops being tuned there, and above 100 digits the number field sieve
(not in Mathilda) would be the right algorithm. The `-DFI_TUNING` build reads
overrides from environment variables (`FI_QS_FB`, `FI_QS_BLK`, `FI_QS_LP`,
`FI_QS_SLACK`, …) for parameter sweeps. Normal builds compile those hooks out.

---

## 3. Invoking it, and where it shines

The method has **no suboptions**. Passing one is an error (`FactorInteger::mopt`),
and the call stays unevaluated:

```mathematica
In[22]:= Quiet[FactorInteger[sp[20], Method -> {"QuadraticSieve", "Iterations" -> 3}]]
Out[22]= FactorInteger[21000000136000000019, Method -> {"QuadraticSieve", "Iterations" -> 3}]
```

There are four spellings: `Method -> "QuadraticSieve"`, `"QS"`, `"SIQS"` and
`"MPQS"`. The heads ``FactorInteger`QuadraticSieve``, ``FactorInteger`QS`` and
``FactorInteger`SIQS`` also exist; `"MPQS"` is a method string only. Balanced
semiprimes of 20, 30, 40, 50 and 60 digits follow:

```mathematica
In[23]:= FactorInteger`QS[sp[20]]
Out[23]= {{3000000019, 1}, {7000000001, 1}}

In[24]:= FactorInteger`SIQS[sp[30]]
Out[24]= {{300000000000089, 1}, {700000000000051, 1}}

In[25]:= FactorInteger[sp[40], Method -> "MPQS"] == FactorInteger`QuadraticSieve[sp[40]]
Out[25]= True

In[26]:= FactorInteger`QuadraticSieve[sp[50]]
Out[26]= {{3000000000000000000000007, 1}, {7000000000000000000000061, 1}}

In[27]:= f = FactorInteger`QuadraticSieve[sp[60]]
Out[27]= {{300000000000000000000000000007, 1}, {700000000000000000000000000033, 1}}

In[28]:= Times @@ (Power @@@ f) == sp[60] && AllTrue[f[[All, 1]], PrimeQ]
Out[28]= True
```

Approximate wall times in this session, with 16 threads on a loaded machine:

| digits | 20 | 30 | 40 | 50 | 60 |
|---|---|---|---|---|---|
| time | ≈ 2–7 ms | ≈ 7–9 ms | ≈ 40 ms | ≈ 60–100 ms | ≈ 0.3–0.5 s |

The shoot-out report gives the larger sizes: about 3 s at 70 digits, 27–32 s at
80, and about 285 s at 90. In single-threaded development runs, 60 digits took
0.4 s against PARI's 3.2–3.6 s, so the lead is not just parallelism.

**The ≤ 64-bit shortcut.** A word-sized input never reaches the sieve. The
splitter goes straight to micro-ECM (256 curves) and SQUFOF in word arithmetic:

```mathematica
In[29]:= FactorInteger`QuadraticSieve[NextPrime[2^31] NextPrime[2^32]]
Out[29]= {{2147483659, 1}, {4294967311, 1}}
```

Above \(2^{64}\) there is no such shortcut, and even a 150-bit product of two
Mersenne primes is sieved (about 30 ms):

```mathematica
In[30]:= FactorInteger`QS[(2^61 - 1) (2^89 - 1)]
Out[30]= {{2305843009213693951, 1}, {618970019642690137449562111, 1}}
```

---

## 4. Where it is the wrong tool

The sieve's cost depends on \(N\) alone, so the QS gains nothing from a small
factor. Take an 80-digit number with a 10-digit factor and a 71-digit prime
cofactor:

```mathematica
In[31]:= m = NextPrime[10^9] NextPrime[10^70];
Out[31]= Null

In[32]:= FactorInteger`ECM[m]
Out[32]= {{1000000007, 1}, {10000000000000000000000000000000000000000000000000000000000000000000033, 1}}

In[33]:= TimeConstrained[FactorInteger`QuadraticSieve[m], 5, $Aborted]
Out[33]= $Aborted

In[34]:= FactorInteger[m] == FactorInteger`ECM[m]
Out[34]= True
```

ECM finds the 10-digit factor in about 20 ms, because its cost scales with the
size of the factor found ([ecm.md](ecm.md)). The QS would have to do the full
80-digit job, about 30 s of wall time, and it runs out of a 5-second CPU budget at
once. (`TimeConstrained` measures process CPU time, so 16 sieving threads use it
up in a fraction of a second. The abort is taken when the current batch of
polynomial families ends.) `Automatic` gets this right. It runs rho, micro-ECM,
p±1 and rising ECM levels first, and commits to the QS only when

- the composite has at most 100 digits, and
- the ECM effort spent satisfies \(7\cdot\text{effort} \ge 2\cdot\text{digits}\)
  (yafu's crossover rule).

So small factors are peeled off cheaply, and only a genuinely hard balanced
remainder pays for the sieve. On a balanced 40-digit input the cascade reaches the
QS within a few milliseconds:

```mathematica
In[35]:= FactorInteger[sp[40]] == FactorInteger`QuadraticSieve[sp[40]]
Out[35]= True
```

Other weak spots are noted in the shoot-out report:

- At 30–40 digits fixed costs dominate (factor base, multiplier, the first
  polynomial), and PARI is level or ahead.
- 90–126-bit inputs would benefit from a two-limb SIQS, which Mathilda does not
  have; there they go through rho and micro-ECM instead.

---

## 5. Limitations

- It has no suboptions. The parameters come from the fixed table above, and can
  be overridden only in a `-DFI_TUNING` build.
- The factor base is limited to fewer than \(2^{17}\) primes (the bucket entry packs
  the prime index into 17 bits), more than the table ever asks for.
- Its cost depends only on the size of \(N\): it cannot exploit small factors
  (use ECM) or special form (use [p±1](pminus1-pplus1.md),
  [Fermat](fermat-lehman.md) or [rational base descent](blake-rational-base-descent.md)).
- Above about 100 digits the number field sieve is asymptotically faster, and
  Mathilda has none. `Automatic` does not take the QS above 100 digits and keeps
  climbing ECM levels instead.
- Results are deterministic for a given \(N\) and thread count: polynomial
  choices are seeded from \(N\), and logs are replayed in worker order.

---

## See also

- [CFRAC and Dixon](cfrac-dixon.md): the same relations without a sieve.
- [ECM](ecm.md): the right tool when a factor is small.
- [Overview and the Automatic cascade](index.md), [performance](performance.md).
- [Trial division](trial-division.md), [Pollard rho](pollard-rho.md),
  [SQUFOF](squfof.md) (which the ≤ 64-bit shortcut uses).
