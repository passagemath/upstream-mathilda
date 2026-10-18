# Performance: methods head to head

The method pages each describe one algorithm in isolation. This page puts them
side by side, in two ways:

1. **Method against method.** Every `Method` is run on each of twelve *number
   shapes*, so you can see which algorithm suits which number, and how badly
   the others do.
2. **Mathilda against other systems.** `Automatic` is compared with PARI/GP,
   FLINT, Mathematica, Maxima and SymPy on the seeded corpus of the
   FactorInteger shootout.

!!! note "How to read the timings"
    Timings vary from run to run, so they appear only in tables, never in an
    `Out[]`. The method matrix was measured on an Intel i9-9880H (8 cores, 16
    threads) at v0.459 **while other jobs were running** (load average about
    11), so read its numbers as rough upper bounds. The shootout numbers in §2
    were measured on the same machine with nothing else benchmarking. They
    come from `benchmarks/98-factorinteger-shootout/REPORT.md` (v0.453; the
    sieve and ECM code has not changed since).

---

## 1. The method matrix

### The shapes

Each shape is built in a single line of Mathilda, so you can rebuild any of
them yourself. "Smooth p−1" means a 24-digit prime \(p\) whose \(p-1\) has no
prime factor above 61, times a 41-digit prime:

```mathematica
In[1]:= smooth = Product[Prime[i], {i, 2, 18}];
In[2]:= pm1 = 2 NestWhile[# + 1 &, 1, ! PrimeQ[2 # smooth + 1] &] smooth + 1
Out[2]= 234576762718813941966541
In[3]:= FactorInteger[pm1 - 1]
Out[3]= {{2, 2}, {3, 1}, {5, 1}, {7, 1}, {11, 1}, {13, 1}, {17, 1}, {19, 1}, {23, 1}, {29, 1}, {31, 1}, {37, 1}, {41, 1}, {43, 1}, {47, 1}, {53, 1}, {59, 1}, {61, 1}}
```

| Shape | Definition | Digits |
|---|---|--:|
| u64 | `NextPrime[2^31 + 12345] NextPrime[2^32 - 999]` | 19 |
| bal20 | `NextPrime[3 10^9] NextPrime[7 10^9]` | 20 |
| bal30 | `NextPrime[3 10^14] NextPrime[7 10^14]` | 30 |
| bal40 | `NextPrime[3 10^19] NextPrime[7 10^19]` | 40 |
| bal60 | `NextPrime[3 10^29] NextPrime[7 10^29]` | 60 |
| p10·p50 | `NextPrime[10^9] NextPrime[10^50]` | 60 |
| p20·p80 | `NextPrime[10^19] NextPrime[10^79]` | 99 |
| close | `NextPrime[10^30] NextPrime[10^30 + 10^6]` | 61 |
| smooth p−1 | `pm1 NextPrime[10^40]`, with `pm1` as above | 65 |
| smooth p+1 | the same construction with `- 1`, times `NextPrime[10^40]` | 65 |
| RBD | `NextPrime[Floor[3 (22/7)^150]]` × a 40-bit prime (the spec example) | 88 |
| F7 | `2^128 + 1` | 39 |

### The results

Each cell runs `AbsoluteTiming[FactorInteger[n, Method -> m]]` in a fresh
process under a **20-second wall-clock** limit. Legend:

- a time: the answer was checked, its product equals \(n\) and every base
  passes `PrimeQ`;
- **✗**: the method returned with a composite still in the list
  (`FactorInteger::nofac`, often within milliseconds, because it declined or
  ran out of its own iteration budget);
- **—**: still running after 20 s.

TD = TrialDivision, Rho = PollardRho, SQUFOF = ShanksSquareForms,
QS = QuadraticSieve, p−1 = PollardPMinus1, p+1 = WilliamsPPlus1,
RBD = BlakeRationalBaseDescent.

| Shape | Automatic | TD | Rho | Fermat | Lehman | SQUFOF | CFRAC | Dixon | QS | p−1 | p+1 | ECM | RBD |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| u64 | 0.3 ms | ✗ | **0.2 ms** | ✗ | 3.8 ms | 2.9 ms | 7.5 ms | 0.6 s | 0.2 ms | 3.2 ms | 3.9 ms | 4.4 ms | ✗ |
| bal20 | 4.3 ms | ✗ | 6.3 ms | ✗ | 5.4 ms | **2.9 ms** | 11 ms | 0.7 s | 5.5 ms | 139 ms | 398 ms | 6.9 ms | ✗ |
| bal30 | 12 ms | ✗ | 1.4 s | ✗ | **4.8 ms** | 139 ms | 50 ms | — | 11 ms | 154 ms | 407 ms | 76 ms | ✗ |
| bal40 | 51 ms | ✗ | — | ✗ | 50 ms | ✗ | 0.52 s | — | **48 ms** | ✗ | 0.56 s | 2.2 s | ✗ |
| bal60 | 0.80 s | ✗ | — | ✗ | ✗ | ✗ | — | — | **0.49 s** | ✗ | ✗ | 9.5 s | ✗ |
| p10·p50 | 76 ms | ✗ | 25 ms | ✗ | ✗ | ✗ | — | — | 0.58 s | 320 ms | 345 ms | 20 ms | **2.4 ms** |
| p20·p80 | 186 ms | ✗ | — | ✗ | — | ✗ | — | — | — | 290 ms | 0.77 s | 234 ms | **4.4 ms** |
| close | 7.3 ms | ✗ | — | **0.13 ms** | 5.1 ms | ✗ | — | — | 0.72 s | ✗ | ✗ | 5.6 s | ✗ |
| smooth p−1 | 25 ms | ✗ | — | ✗ | — | ✗ | — | — | 1.3 s | 134 ms | 347 ms | 0.48 s | **3.8 ms** |
| smooth p+1 | 93 ms | ✗ | — | ✗ | — | ✗ | — | — | 1.4 s | ✗ | 627 ms | 8.5 s | **2.5 ms** |
| RBD | 47 ms | ✗ | 0.58 s | ✗ | — | ✗ | — | — | — | 183 ms | 0.79 s | 24 ms | **10 ms** |
| F7 | 42 ms | ✗ | — | ✗ | ✗ | ✗ | 0.41 s | — | **41 ms** | ✗ | ✗ | 201 ms | ✗ |

The raw data is `benchmarks/98-factorinteger-shootout/results/method-matrix-v0459.json`.
To regenerate it after `make`:

```
python3 benchmarks/98-factorinteger-shootout/method_matrix.py results.json
```

### What the matrix shows

**Automatic is the only column with no ✗ or —.** It is never the fastest:
it pays for the tiers it tries first. On the balanced rows that overhead stays
under 2× the best pinned method. On bal60, for example, about 0.3 s goes on
rho, the Fermat probe, p−1 and p+1 before the [QS crossover](index.md) sends
the number to the sieve. On the special shapes the right pinned method can be
far ahead: Fermat is 56× faster on "close", and ECM is 4× faster on p10·p50.
Pin a method only when you *know* the shape of your numbers.

**The trap in this matrix is `NextPrime[10^k]`.** RBD "wins" five rows, and
four of them are rows meant to be about something else. The reason is that a
prime just above \(10^{k}\) is \(c \cdot 10^{k}\) with \(c \approx 1\). That
is exactly the structure \(P \approx c\,(a/b)^k\) that
[rational base descent](blake-rational-base-descent.md) looks for, with base
\(a/b = 10\). The shapes built from \(3\cdot10^k\) and \(7\cdot10^k\) are
structured in the same way, but the descent finds a split only when one factor
is much larger than the other, so the balanced rows (bal20–bal60) defeat it.
Benchmark corpora built from round numbers flatter RBD; random primes do not.

```mathematica
In[4]:= n1 = NextPrime[10^9] NextPrime[10^50];
In[5]:= FactorInteger`BlakeRationalBaseDescent[n1]
Out[5]= {{1000000007, 1}, {100000000000000000000000000000000000000000000000151, 1}}
```

**Lehman's bal40 success is also a structured win.** The factors of bal40 are
in the ratio \(7/3\), and Lehman searches \(4kn\) for small \(k\), which
catches a factor ratio close to a small rational. At 60 digits it fails
anyway: for bignums the first \(k = 1\) window uses up the whole iteration
budget (see the warning on the [Fermat and Lehman](fermat-lehman.md) page).
An unstructured method such as p−1 sees nothing in bal40:

```mathematica
In[6]:= bal40 = NextPrime[3 10^19] NextPrime[7 10^19];
In[7]:= FactorInteger`Lehman[bal40]
Out[7]= {{30000000000000000041, 1}, {70000000000000000013, 1}}
In[8]:= FactorInteger`PollardPMinus1[bal40]
Out[8]= {{2100000000000000003260000000000000000533, 1}}
```

**The special-purpose methods are dramatic when they apply, and useless
otherwise.**

- *Fermat* splits the 61-digit "close" number in 0.13 ms, about 5000× faster
  than QS, because \(|p - q| \approx 10^6\) is tiny next to \(n^{1/4}\). On
  every other row it fails.
- *p−1* splits the smooth-p−1 row in 134 ms; on the smooth-p+1 row it fails
  outright.
- *p+1* gets both smooth rows. A Lucas sequence whose discriminant is a
  quadratic residue mod \(p\) works in the group of order \(p-1\), so p+1
  also catches smooth \(p-1\) (see the
  [p−1 and p+1](pminus1-pplus1.md) page).

```mathematica
In[9]:= FactorInteger[NextPrime[10^30] NextPrime[10^30 + 10^6], Method -> "Fermat"]
Out[9]= {{1000000000000000000000000000057, 1}, {1000000000000000000000001000003, 1}}
```

**The cost of rho and ECM depends on the factor's size, not on \(n\).** Rho
finds the 10-digit factor of a 60-digit \(n\) in 25 ms, but fails on every
balanced row from 40 digits up. ECM is the only general method that gets every
row. It is the right tool for p10·p50, where it is 30× faster than QS, and the
wrong tool for balanced semiprimes, where it is 20–45× slower than QS.

**QS is the only method whose cost follows the size of \(n\) and nothing
else.** Its cost is about the same for every number of a given size. It is
the fastest method on bal40, bal60 and F7, and it is pointless when one factor
is small: on p20·p80 it is still sieving a 99-digit number at 20 s, while ECM
finishes in 0.23 s. That is why `Automatic` tries ECM levels before the
crossover, and only up to 100 digits.

**Some methods suit only small numbers.**

- *Trial division* fails every row: its default bound is \(10^7\), and every
  shape here has a factor above that.
- *SQUFOF* is the fastest method at 20 digits and declines above about 124
  bits.
- *Dixon* takes 0.6 s on a 19-digit number, about 3000× slower than rho. It is
  in the library because it is the simplest congruence-of-squares method, a
  stepping stone to [CFRAC and QS](cfrac-dixon.md).
- *CFRAC* is still respectable at 40 digits, about 10× behind QS.

### A rule of thumb

| If you know that… | Pin | Otherwise |
|---|---|---|
| the two factors are very close | `"Fermat"` | Automatic's Fermat probe (\(2^{14}\) steps) finds it after its rho tier |
| some factor has smooth \(p\pm1\) | `"WilliamsPPlus1"` (both) or `"PollardPMinus1"` | Automatic runs one p−1 and one p+1 pass |
| one factor is much smaller than the rest | `"ECM"` | Automatic climbs the same ECM levels |
| \(n\) is a balanced semiprime, 40–100 digits | `"QuadraticSieve"` | Automatic gets there after about 0.3 s of earlier tiers |
| the large factor is near \(c\,(a/b)^k\) | `"BlakeRationalBaseDescent"` | Automatic never tries RBD |
| nothing at all | `Automatic` | — |

---

## 2. Against other systems

From the FactorInteger shootout (`REPORT.md`). These are single-case times for
`Automatic` (or each system's default), and every answer was checked with
SymPy. A timeout is 600 s for PARI and FLINT and 120 s for Mathematica, Maxima
and SymPy.

### Balanced semiprimes (the quadratic sieve)

| Digits | Mathilda | PARI/GP 2.17 | FLINT 3.6 | Mathematica 14 |
|--:|--:|--:|--:|--:|
| 30 | 6–7 ms | **3 ms** | 6–7 ms | 62 ms |
| 40 | 27–30 ms | 26–29 ms | 33–43 ms | 18–331 ms |
| 50 | **76–109 ms** | 0.28–0.31 s | 0.36–0.39 s | 6.2–6.5 s |
| 60 | **0.60 s** | 3.2–3.6 s | 4.7–6.5 s | timeout |
| 70 | **3.1–3.5 s** | 27–31 s | 27–52 s | timeout |
| 80 | **27–32 s** | 322–494 s | 237–346 s | timeout |
| 90 | **≈ 285 s** | timeout | timeout | timeout |

Mathilda's SIQS sieves on up to 16 threads above about 45 digits; PARI and
FLINT run on one thread. The lead does not depend on threads: single-threaded,
Mathilda takes 0.37–0.46 s at 60 digits and 145 s at 80.

### A hidden factor in a 150-digit number (ECM)

| Factor | Mathilda | PARI/GP | FLINT |
|---|--:|--:|--:|
| p20 | 1.2–3.2 s | 20–22 s | 3.0–3.7 s |
| p25 | **3.2–12 s** | 70–188 s | 130–207 s |
| p30 | **26–54 s** | timeout | 451–533 s |
| p35 | **203–370 s** | timeout | timeout |

### Special shapes

| Case | Mathilda | Best other |
|---|--:|---|
| p40 with smooth p−1, × p80 | **32 ms** | every other engine timed out |
| p40 with smooth p+1, × p80 | **≈ 200 ms** | every other engine timed out |
| Fermat-shaped, 60/80/100 digits | 5–6 ms | SymPy 1.4–2.3 ms (PARI 3 s at 60 digits, 265 s at 80) |
| \(F_7 = 2^{128}+1\) | **24.5 ms** | PARI 27 ms |
| \(F_8 = 2^{256}+1\) | **190 ms** | FLINT 337 ms |
| \(M_{149} = 2^{149}-1\) | **30 ms** | PARI 91 ms |
| 2000 random 63-bit integers | 26.2 ms | FLINT 26.1 ms |

### Where Mathilda is behind, and why

- **Below about 40 digits** the cost is fixed overhead, not the sieve.
  `Automatic` spends a few milliseconds on rho and micro-ECM before the
  crossover, and the SIQS setup is a large share of a sub-10 ms run. PARI's 3 ms
  at 30 digits beats Mathilda's 6–7 ms.
- **Fermat-shaped numbers.** SymPy tries a few Fermat steps *first*; Mathilda
  runs its rho tier (about 4 ms) before the Fermat probe. Pinning
  `Method -> "Fermat"` takes 0.13 ms (table above).
- **96-bit semiprimes** (PARI 40 ms against 77 ms for twenty). A minority miss
  the two-limb micro-ECM budget and fall through to the slower mpz tiers.

`REPORT.md` gives the fix proposed for each.

## See also

- [Overview and the Automatic cascade](index.md)
- [Quadratic sieve](quadratic-sieve.md), [ECM](ecm.md),
  [p−1 and p+1](pminus1-pplus1.md), [Fermat and Lehman](fermat-lehman.md),
  [rational base descent](blake-rational-base-descent.md)
- `benchmarks/98-factorinteger-shootout/` — the corpus, the drivers for each
  system, and `verify_methods.py`, the randomised correctness sweep over every
  method.
