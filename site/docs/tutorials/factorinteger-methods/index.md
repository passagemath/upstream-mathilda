# FactorInteger methods

The [number theory tutorial](../05-number-theory.md) introduces
`FactorInteger` and lets it choose an algorithm for you
(`Method -> Automatic`). Underneath, `FactorInteger` is a collection of a
dozen factoring algorithms, from trial division to the self-initialising
quadratic sieve, with a scheduler that runs them cheapest first.

These advanced tutorials open up that collection. Each page covers **one
method**: the mathematics, how Mathilda implements it, its suboptions, and a
set of worked examples. The examples show where the method shines and where it
fails, and the reasons in both cases. Every transcript is reproduced from the
current build (v0.459). The pages are written for readers who know what a
prime factorisation is and want to know *how* one is found. They are also for
anyone who needs to choose a method for a hard number, or to understand why one
number takes a millisecond and its neighbour a minute.

```mathematica
In[1]:= FactorInteger[2^128 + 1]
Out[1]= {{59649589127497217, 1}, {5704689200685129054721, 1}}

In[2]:= FactorInteger[2^128 + 1, Method -> "ECM"]
Out[2]= {{59649589127497217, 1}, {5704689200685129054721, 1}}
```

Both calls return the same answer, since a factorisation is unique. They get
there differently. `Automatic` reaches the quadratic sieve after about 20 ms of
cheap probing. `"ECM"` climbs elliptic-curve levels, which takes about
0.5 s.

---

## How FactorInteger is organised

The subsystem lives in `src/numbertheory/factorinteger/`, one concern per
file:

| file | contents |
|---|---|
| `fi_builtin.c` | `FactorInteger[...]` and every ``FactorInteger`Name[...]`` head: argument and option parsing, rationals, Gaussian integers, messages, output, docstrings |
| `fi_core.c` | the factor list, the per-call context, the \(n\)-seeded RNG, the recursive driver `fi_factor()` and the `Automatic` tier table |
| `fi_u64.c` | the **word engine**: Montgomery arithmetic, deterministic Miller–Rabin, rho, SQUFOF, Lehman, Fermat, micro-ECM and p±1 on `uint64_t` |
| `fi_u128.c` | two-limb Montgomery rho and micro-ECM for \(2^{64} \le n < 2^{126}\) |
| `fi_trial.c` | batched trial division and the shared prime table |
| `fi_rho.c`, `fi_fermat.c`, `fi_squfof.c` | the bignum versions of rho, Fermat/Lehman and SQUFOF |
| `fi_ecm.c` | ECM, Pollard p−1, Williams p+1 (GMP-ECM, with FLINT and native fallbacks) |
| `fi_qs.c`, `fi_lanczos.c` | the self-initialising quadratic sieve and block Lanczos |
| `fi_cfrac.c`, `fi_dixon.c`, `fi_smooth.c`, `fi_gf2.c` | CFRAC, Dixon, threaded smoothness testing, and the GF(2) relation store and solver |
| `fi_rbd*.c` | Blake's rational base descent and its lattice detector |
| `fi_par.c` | thread regions and work claiming |

Every algorithm is a **splitter**. It is given an odd composite that is not a
perfect power, and either returns a proper divisor or gives up when its budget
is spent. A splitter never loops without bound. The driver does the rest. It
strips powers of two, recognises primes and perfect powers, pushes both halves
of every split back onto a work stack, and reports any composite that no
splitter could break. Such a composite is returned with exponent 1 and the
message `FactorInteger::nofac`.

---

## The automatic pipeline

### Below \(2^{64}\): the word engine

A machine integer takes a fast path straight into `fi_u64.c`, with no bignums
and no heap allocation until the result list:

1. Strip powers of 2. Test the odd primes below \(2^{10}\) with a
   multiply-by-inverse divisibility check.
2. A cofactor below \(2^{20}\) with no factor below \(2^{10}\) is prime.
   Otherwise, use deterministic Miller–Rabin, and test whether it is a perfect
   square.
3. Split with `u64_split`. Below \(2^{40}\), run rho with \(c = 1, \dots, 63\)
   at \(2^{22}\) steps each. Then micro-ECM with up to 128 curves, then SQUFOF,
   then a long rho backstop.

```mathematica
In[3]:= FactorInteger[2^64 - 1]
Out[3]= {{3, 1}, {5, 1}, {17, 1}, {257, 1}, {641, 1}, {65537, 1}, {6700417, 1}}
```

### Above \(2^{64}\): the tier table

Larger inputs are trial-divided once to \(2^{16}\). Then every composite on
the work stack is tested for primality (BPSW) and perfect powers, and handed
to `split_auto`. That function walks a **tier table**, cheapest first:

| tier | what runs | budget | catches |
|---|---|---|---|
| rho | Pollard–Brent, one start | \(2^{13}\) steps | factors up to about 8 digits |
| micro-ECM | two-limb ECM (only \(n < 2^{126}\)) | 8 curves, \(B_1 = 150\) | 9–12-digit factors |
| Fermat | residue-sieved Fermat probe | \(2^{14}\) steps | two factors within about \(2^7 n^{1/4}\) |
| p−1 | GMP-ECM P−1 | \(B_1 = \min(2000\,d, 150000)\) | \(p\) with smooth \(p-1\) |
| p+1 | GMP-ECM P+1 | \(B_1 = \min(1500\,d, 100000)\); 1 seed, 3 from 70 digits | \(p\) with smooth \(p+1\) |
| ECM 15, 20, …, 65 | GMP-ECM curves at each level of the table | e.g. 25 curves at 2000, 74 at 11000, … | factors of 15–65 digits |

Here \(d\) is the composite's digit count, estimated as
\(\lfloor 0.30103\cdot\text{bits}\rfloor + 1\).

**The tier is per item.** A work item advances to the next tier only when the
current tier *fails*. After a split, both halves re-enter at the tier that
succeeded. So a number full of small primes keeps peeling them off with cheap
rho, and never pays for ECM.

**The QS crossover.** Before every tier from p−1 on, the scheduler asks
whether the [quadratic sieve](quadratic-sieve.md) would now be cheaper than
searching further. It assigns each finished tier an *effort*, the size of
factor it would have found with good probability: 7 digits after rho and
Fermat, 13 after p−1, 15 after p+1, then 16, 16, 20, 25, … after each ECM
level. QS takes the item, once, when

\[
7 \times \text{effort} \;\ge\; 2 \times d \qquad\text{and}\qquad d \le 100 .
\]

This is yafu's crossover rule: stop looking for a factor once you have spent
about 2/7 of the digit count. The rule works out to:

| composite size \(d\) | QS takes over |
|---|---|
| \(\le 45\) digits | right after the Fermat probe, before p−1 |
| 46–52 | after p−1 |
| 53–56 | after p+1 |
| 57–70 | after ECM level 20 |
| 71–87 | after ECM level 25 |
| 88–100 | after ECM level 30 |
| \(> 100\) | never: ECM levels up to 65 digits, then `nofac` |

There is no trace output. The scheduler is silent, but you can **replay a
tier by hand** by pinning the method with the tier's own budget. Three
examples follow.

**Example 1: a 66-digit number with a smooth \(p-1\).**
\(p = 2\,\mathrm{lcm}(1..60) + 1\), times a 41-digit prime. Here \(d = 66\), so
QS stays out of the way, and the tiers run in order. Rho with its \(2^{13}\)
steps fails, and so does the Fermat probe. Micro-ECM is skipped, because the
number is above \(2^{126}\). Then p−1 at \(B_1 = \min(2000\cdot66, 150000) = 132000\)
succeeds:

```mathematica
In[4]:= L = LCM @@ Range[60]; n1 = (2 L + 1) NextPrime[10^40]; IntegerLength[n1]
Out[4]= 66

In[5]:= FactorInteger[n1, Method -> {"PollardRho", "Iterations" -> 2^13, "Starts" -> 1}]
Out[5]= {{193814243295544634018256010000000000002345152343876090071620897721, 1}}

In[6]:= FactorInteger[n1, Method -> {"Fermat", "Iterations" -> 2^14}]
Out[6]= {{193814243295544634018256010000000000002345152343876090071620897721, 1}}

In[7]:= FactorInteger[n1, Method -> {"PollardPMinus1", "B1" -> 132000}]
Out[7]= {{19381424329554463401825601, 1}, {10000000000000000000000000000000000000121, 1}}

In[8]:= FactorInteger[n1]
Out[8]= {{19381424329554463401825601, 1}, {10000000000000000000000000000000000000121, 1}}
```

The last call takes about 20 ms.

**Example 2: \(F_7 = 2^{128} + 1\), 39 digits.** \(F_7\) has 129 bits, so the
micro-ECM tier is skipped, and the cheap rho tier cannot find a 17-digit
factor. The Fermat probe fails too. Before p−1, the crossover test reads
\(7 \times 13 = 91 \ge 2 \times 39\), so the number goes straight to QS. That
explains the 20 ms of `In[1]`:

```mathematica
In[9]:= FactorInteger[2^128 + 1, Method -> {"Fermat", "Iterations" -> 2^14}]
Out[9]= {{340282366920938463463374607431768211457, 1}}

In[10]:= FactorInteger[2^128 + 1, Method -> "QuadraticSieve"]
Out[10]= {{59649589127497217, 1}, {5704689200685129054721, 1}}
```

**Example 3: a hidden 20-digit factor.** \(p_{20} = 40000000000000000019\)
times an 81-digit prime is 100 digits. The rho, Fermat, p−1 and p+1 tiers
(three seeds at this size) all fail, because \(p_{20} - 1\) and \(p_{20} + 1\) each
have two large prime factors (up to 11 and 10 digits respectively). ECM level 15 has a small chance, and
level 20 (74 curves at \(B_1 = 11000\)) is the one designed for this factor.
QS would only be considered after level 30, since \(7 \times 30 \ge 200\). The
call takes about 0.4–1.2 s on a loaded 8-core i9:

```mathematica
In[11]:= p20 = NextPrime[4 10^19 + 1]; FactorInteger[p20 NextPrime[10^80 + 12345]]
Out[11]= {{40000000000000000019, 1}, {100000000000000000000000000000000000000000000000000000000000000000000000000012397, 1}}
```

---

## Pinning a method

There are three equivalent spellings:

```mathematica
In[12]:= FactorInteger[2^101 - 1, Method -> "PollardRho"]
Out[12]= {{7432339208719, 1}, {341117531003194129, 1}}

In[13]:= FactorInteger[2^101 - 1, Method -> {"PollardRho", "Iterations" -> 10^7, "Starts" -> 2}]
Out[13]= {{7432339208719, 1}, {341117531003194129, 1}}

In[14]:= FactorInteger`PollardRho[2^101 - 1, "Iterations" -> 10^7, "Starts" -> 2]
Out[14]= {{7432339208719, 1}, {341117531003194129, 1}}
```

- `Method -> "Name"` uses the method's defaults.
- `Method -> {"Name", "Sub" -> v, ...}` sets suboptions.
- ``FactorInteger`Name[n, "Sub" -> v, ...]`` calls the method directly. The
  heads are `Listable` and `Protected`, like `FactorInteger` itself, and
  `?FactorInteger`Name` shows each one's help.

**Method names.** `"TrialDivision"`, `"PollardRho"`, `"Fermat"`, `"Lehman"`,
`"ShanksSquareForms"`, `"CFRAC"`, `"Dixon"`, `"QuadraticSieve"`,
`"PollardPMinus1"`, `"WilliamsPPlus1"`, `"ECM"` and
`"BlakeRationalBaseDescent"`. The aliases are `"SQUFOF"`; `"QS"`, `"SIQS"` and
`"MPQS"`; `"ContinuedFraction"`; `"EllipticCurve"`; `"PollardP-1"`,
`"PollardPM1"` and `"P-1"`; `"WilliamsP+1"`, `"WilliamsPP1"` and `"P+1"`. Of
the aliases, only `SQUFOF`, `QS` and `SIQS` also exist as heads
(``FactorInteger`SQUFOF`` and so on).

A method pinned by name is used **at every size**. It keeps its own algorithm
even below \(2^{64}\), where the methods that have a word version (rho,
Fermat, Lehman, SQUFOF, p±1, ECM) run it in word arithmetic. A pinned method
does not fall back to another one. The exceptions are QS, which hands a
word-size composite to micro-ECM and SQUFOF, and the rational base descent,
which retries a stubborn cofactor with a light default descent. A composite it cannot split is returned
unfactored. An unknown method name or suboption leaves the call unevaluated,
with a `FactorInteger::method` or `FactorInteger::mopt` message:

```mathematica
In[15]:= FactorInteger[2^64 + 1, Method -> {"ECM", "Bound" -> 5}]
Out[15]= FactorInteger[18446744073709551617, Method -> {"ECM", "Bound" -> 5}]

In[16]:= FactorInteger[2^64 + 1, Method -> {"TrialDivision", "Bound" -> 1000}]
Out[16]= {{18446744073709551617, 1}}
```

The second call is well formed. Its answer just leaves the composite
\(2^{64}+1 = 274177\times67280421310721\) unfactored, because both factors
are above the bound.

---

## Which method for which number

| the number looks like | best method | why |
|---|---|---|
| anything; you do not know | `Automatic` | It runs the cheap tests first and the expensive ones only when they are needed. |
| below \(2^{64}\) | `Automatic` (the word engine) | Microseconds per number. A pinned method is rarely faster here. |
| many small primes, all below \(10^7\) | [`"TrialDivision"`](trial-division.md) | Batched division. Complete when \(\sqrt n < 10^7\). |
| one or more factors of 6–15 digits | [`"PollardRho"`](pollard-rho.md) | \(\sqrt p\) steps, independent of the size of \(n\). |
| two factors that agree in their leading half of digits | [`"Fermat"`](fermat-lehman.md) | Costs \((q-p)^2/8\sqrt n\) steps: instant when \(\lvert p-q\rvert \lesssim n^{1/4}\). |
| \(p/q\) close to a small fraction \(u/v\); \(n < 2^{42}\) | [`"Lehman"`](fermat-lehman.md) | Fermat on \(4uvn\). Complete \(O(n^{1/3})\) on words. |
| a machine word, as a backstop | [`"SQUFOF"`](squfof.md) | \(O(n^{1/4})\), all in word arithmetic. Declines above about 124 bits. |
| a prime factor with smooth \(p - 1\) | [`"PollardPMinus1"`](pminus1-pplus1.md) | The cost depends on \(B_1\), not on \(p\). |
| a prime factor with smooth \(p + 1\) | [`"WilliamsPPlus1"`](pminus1-pplus1.md) | Same, in the norm-1 group of \(\mathbb F_{p^2}\). Each seed is a coin flip. |
| a 15–40-digit factor of a number of any size | [`"ECM"`](ecm.md) | \(L_p[1/2,\sqrt2]\) in the size of the **factor**. |
| a balanced semiprime of 30–100 digits | [`"QuadraticSieve"`](quadratic-sieve.md) | \(L_n[1/2, 1]\), deterministic, threaded. The fastest general method in this range. |
| teaching, or comparison with QS | [`"CFRAC"` / `"Dixon"`](cfrac-dixon.md) | The historical congruence-of-squares methods that QS improves on. |
| a factor close to \(c\,(a/b)^k\) | [`"BlakeRationalBaseDescent"`](blake-rational-base-descent.md) | A polynomial-time descent for this structure, whatever the size of the factors. |

Timings against PARI/GP, FLINT, Mathematica and SymPy are collected on the
[performance](performance.md) page.

---

## Determinism and threads

**Determinism.** Every random choice is drawn from a generator **seeded by
\(n\) itself** (`fi_ctx_seed` hashes the limbs of \(n\)): rho's starts, the ECM
\(\sigma\) values, the p+1 seeds on words, and the QS polynomials. A given
input therefore always takes the same path. Its running time, which tier
succeeds, and which composite is left over in an incomplete answer all depend
on \(n\) alone, never on what the session factored before. The complete
factorisation is unique anyway. What determinism buys is reproducible timings,
reproducible `nofac` results from pinned methods, and reproducible partial
answers from `FactorInteger[n, k]`.

**Threads.** When Mathilda is built with threads (the default `USE_THREADS=1`),
the expensive methods run on up to \(\min(\text{cores}, 16)\) threads:

| method | parallel unit | from |
|---|---|---|
| QS | one polynomial family (one \(A\)) per worker | \(kn \ge 150\) bits (about 45 digits) |
| ECM | one curve per worker | \(n \ge 160\) bits |
| p+1 | one seed per worker | always (up to 3 seeds) |
| CFRAC, Dixon | chunks of candidates for smoothness testing | always |
| rational base descent | blocks of bases, and lattice depths | always |

Workers keep their results in order. A region's successes are taken in curve,
seed or candidate order, and relations are replayed in worker order. So for
ECM and p+1 the factor found does not depend on the thread count. For QS the
answer is the same for a given thread count. Workers block every signal, and
the main thread waits in `pthread_join`. A `TimeConstrained` timeout is
therefore taken safely between regions. Note that `TimeConstrained` measures
CPU time, which accumulates on all the workers (see the note on the
[ECM](ecm.md) page).

---

## Partial factorisation, rationals and Gaussian integers

**`FactorInteger[n, k]`** returns at most \(k\) entries. The smallest primes
come first, and the last entry may be an unfactored remainder. The driver
always processes the smaller half of a split first, so small factors appear
early, and the search stops as soon as \(k - 1\) primes are known:

```mathematica
In[17]:= FactorInteger[20!, 3]
Out[17]= {{2, 18}, {3, 8}, {1414538125, 1}}

In[18]:= FactorInteger[(2^31 - 1) (2^61 - 1) NextPrime[10^30], 2]
Out[18]= {{2147483647, 1}, {2305843009213693951000000000131433051525180555207, 1}}

In[19]:= FactorInteger[2^64 + 1, 1]
Out[19]= {{18446744073709551617, 1}}
```

`FactorInteger[n, Automatic]` is the same as `FactorInteger[n]`.

**Rationals and signs.** A negative \(n\) gives a leading \(\{-1, 1\}\). The
denominator of a rational gets negative exponents. Both halves are factored
with the chosen method:

```mathematica
In[20]:= FactorInteger[-120/77]
Out[20]= {{-1, 1}, {2, 3}, {3, 1}, {5, 1}, {7, -1}, {11, -1}}
```

**Gaussian integers.** A complex input is factored over \(\mathbb Z[i]\)
automatically. For an ordinary integer, ask with `GaussianIntegers -> True`.
Mathilda factors the norm \(a^2 + b^2\) with the chosen method, then splits
each rational prime:

- \(2 = -i(1+i)^2\);
- a prime \(p \equiv 3 \pmod 4\) stays prime;
- a prime \(p \equiv 1 \pmod 4\) splits as \(p = (c + di)(c - di)\), with
  \(c, d\) from Hermite–Serret (Cornacchia), starting from a square root of
  \(-1\) mod \(p\).

The answer follows Mathematica's normalisation: the unit first, unless it is
1, then the primes rotated into the first quadrant, sorted by real part and
then imaginary part:

```mathematica
In[21]:= FactorInteger[13, GaussianIntegers -> True]
Out[21]= {{-I, 1}, {2 + 3*I, 1}, {3 + 2*I, 1}}

In[22]:= FactorInteger[30, GaussianIntegers -> True]
Out[22]= {{-1, 1}, {1 + I, 2}, {1 + 2*I, 1}, {2 + I, 1}, {3, 1}}

In[23]:= FactorInteger[5 + 3 I]
Out[23]= {{-I, 1}, {1 + I, 1}, {1 + 4*I, 1}}

In[24]:= FactorInteger[(3 + 4 I)/(1 - I)]
Out[24]= {{I, 1}, {1 + I, -1}, {2 + I, 2}}
```

---

## The methods

<div class="grid cards" markdown>

-   :material-division: __[Trial division](trial-division.md)__

    Batched division by every prime up to a bound, with a perfect-power
    rescue. Complete below \(10^{14}\) by default. The first step of every
    pipeline.

-   :material-ray-start-arrow: __[Pollard rho](pollard-rho.md)__

    Brent's cycle-finding variant, with Montgomery arithmetic on one and two
    words. Finds \(p\) in about \(\sqrt p\) steps, whatever the size of \(n\).

-   :material-approximately-equal: __[Fermat and Lehman](fermat-lehman.md)__

    Factors from the middle. Fermat is instant for close factors. Lehman
    handles ratios near a small fraction, and is complete \(O(n^{1/3})\) below
    \(2^{42}\).

-   :material-square-root: __[SQUFOF](squfof.md)__

    Shanks' square forms with 16 multipliers. \(O(n^{1/4})\) in pure word
    arithmetic. The word engine's backstop.

-   :material-plus-minus-variant: __[Pollard p−1 and Williams p+1](pminus1-pplus1.md)__

    Special-purpose methods for primes with smooth \(p \mp 1\). Stage 2, the
    p+1 seed lottery, and the word-size limits.

-   :material-chart-bell-curve-cumulative: __[Elliptic-curve method](ecm.md)__

    Lenstra's ECM through GMP-ECM. The cost follows the factor, not \(n\).
    Effort levels, reproducible \(\sigma\), and threaded curves.

-   :material-fraction-one-half: __[CFRAC and Dixon](cfrac-dixon.md)__

    The congruence-of-squares ancestors of the sieve: continued fractions,
    random squares, large primes, and GF(2) elimination.

-   :material-grid: __[Quadratic sieve](quadratic-sieve.md)__

    The self-initialising quadratic sieve with Gray-code polynomials, bucket
    sieving, double large primes and block Lanczos. The fastest general method
    up to 100 digits.

-   :material-stairs-down: __[Blake rational base descent](blake-rational-base-descent.md)__

    A deterministic descent for factors close to a power of a rational, with
    a Howgrave-Graham lattice detector.

-   :material-speedometer: __[Performance](performance.md)__

    Mathilda's FactorInteger against PARI/GP, FLINT, Mathematica, Maxima and
    SymPy on a seeded corpus, and the reasons behind each row.

</div>

## Following along

Start the REPL with `./Mathilda` and type each `In[...]` line yourself,
without the prompt. Each page is one continuous session, so definitions made
early on a page are used later on the same page. Every `Out[...]` shown was
produced by the actual binary. Timings quoted in the text are approximate
(measured on a loaded 8-core i9), and diagnostic messages such as
`FactorInteger::nofac` go to the message stream rather than into `Out[...]`.
