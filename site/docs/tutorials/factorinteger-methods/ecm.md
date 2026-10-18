# The elliptic-curve method (ECM)

`Method -> "ECM"` (alias `"EllipticCurve"`) is Lenstra's elliptic-curve
method, run through GMP-ECM. Its running time depends on the size of the
**factor it finds**, not on the size of \(n\). A 20-digit prime costs about the
same to extract from a 60-digit number as from a 220-digit one. So ECM is the
method for "medium" factors, from about 15 to 40 digits, of numbers far too
large for the quadratic sieve. It is also the main search engine of Mathilda's
automatic pipeline once the cheap tiers have failed.

```mathematica
In[1]:= FactorInteger`ECM[2^128 + 1]
Out[1]= {{59649589127497217, 1}, {5704689200685129054721, 1}}
```

That is the seventh Fermat number \(F_7\), whose 17-digit factor needed a
mainframe and the continued-fraction method in 1970 (Morrison and Brillhart).
Here it is a fraction of a second of ECM.

---

## 1. The mathematical idea

[Pollard p−1](pminus1-pplus1.md) works in the group \((\mathbb Z/p)^\times\) of
order \(p - 1\). It succeeds only when that one fixed number is smooth. Lenstra
replaced the group by the points of an elliptic curve modulo \(p\):

\[
E:\; By^2 = x^3 + Ax^2 + x \pmod n,\qquad \#E(\mathbb F_p) \in \bigl[p + 1 - 2\sqrt p,\; p + 1 + 2\sqrt p\,\bigr].
\]

The order of \(E(\mathbb F_p)\) varies, roughly uniformly, as the curve
changes. Each new curve is a new draw from the Hasse interval. Stage 1
computes \(Q = E_{B_1}\cdot P\), with \(E_{B_1}\) the product of the prime
powers up to \(B_1\), working modulo \(n\) as if \(n\) were prime. If
\(\#E(\mathbb F_p)\) is \(B_1\)-smooth, then \(Q\) is the identity modulo
\(p\). In projective coordinates that means \(p \mid Z_Q\), and
\(\gcd(Z_Q, n)\) exposes \(p\). Stage 2 allows one further prime up to
\(B_2\), as in p−1.

The probability that one curve succeeds is about the probability that a
random number near \(p\) is smooth. Choosing \(B_1\) optimally gives the
heuristic expected cost

\[
L_p\!\left[\tfrac12, \sqrt2\right] = \exp\!\Bigl((\sqrt2 + o(1))\sqrt{\ln p\,\ln\ln p}\Bigr)
\]

multiplications modulo \(n\). That is subexponential in the size of the
**factor**. The size of \(n\) enters only through the cost of one
multiplication modulo \(n\), which grows polynomially.

**Montgomery curves and Suyama's parametrisation.** Montgomery's form
\(By^2 = x^3 + Ax^2 + x\) allows \(x\)-only arithmetic with a ladder. Suyama's
family, parametrised by an integer \(\sigma\), forces 12 to divide every group
order, which is a free head start on smoothness. GMP-ECM adds an
asymptotically fast stage 2, using polynomial arithmetic over a product tree.
That lets \(B_2\) be hundreds to thousands of times \(B_1\) at modest cost.

### References

1. H. W. Lenstra, Jr., "Factoring integers with elliptic curves," *Annals of
   Mathematics* **126** (1987) 649–673.
2. P. L. Montgomery, "Speeding the Pollard and elliptic curve methods of
   factorization," *Mathematics of Computation* **48** (1987) 243–264.
3. R. P. Brent, "Some integer factorization algorithms using elliptic curves,"
   *Australian Computer Science Communications* **8** (1986) 149–163.
4. P. Zimmermann and B. Dodson, "20 years of ECM," *ANTS VII*, LNCS 4076
   (2006) 525–542.
5. M. A. Morrison and J. Brillhart, "A method of factoring and the
   factorization of \(F_7\)," *Mathematics of Computation* **29** (1975)
   183–205.
6. The GMP-ECM project, [gitlab.inria.fr/zimmerma/ecm](https://gitlab.inria.fr/zimmerma/ecm)
   — the README's table of optimal \(B_1\) and expected curve counts.

---

## 2. How Mathilda realises it

The front end is `fi_ecm.c`. It has three backends, chosen when Mathilda is
built:

| backend | when | curves | stage 2 |
|---|---|---|---|
| **GMP-ECM** | default (this build: ECM 7.0.7) | `ecm_factor`, batch parametrisation, or Suyama when `"Sigma"` is given | GMP-ECM's fast stage 2 |
| FLINT | `USE_ECM=0` with FLINT | `fmpz_factor_ecm` | FLINT's, \(B_2 = 100 B_1\) by default |
| native | neither library | Montgomery curves, Suyama \(\sigma\), prime-power ladder | baby-step giant-step with \(D = 210\), \(B_2 = 50 B_1\) |

Below \(2^{64}\), the method runs the word engine's **micro-ECM** instead
(`fi_u64_ecm`). That is up to 256 curves with \(B_1\) between 47 and 165,
depending on the bit size, and \(B_2 = 25 B_1\), entirely in one-word
Montgomery arithmetic. Primes below \(2^{10}\) are always divided out first,
because GMP-ECM cannot invert a starting point that shares a factor with
\(n\).

**Effort levels.** The level table in `fi_ecm.c` is GMP-ECM's table of
optimal parameters. Running the given number of curves at the given \(B_1\)
finds a factor of that many digits with probability about \(1 - 1/e\):

| digits | \(B_1\) | curves |
|---|---|---|
| 15 | 2000 | 25 |
| 20 | 11000 | 74 |
| 25 | 50000 | 214 |
| 30 | 250000 | 430 |
| 35 | \(10^6\) | 904 |
| 40 | \(3\times10^6\) | 2350 |
| 45 | \(1.1\times10^7\) | 4480 |
| 50 | \(4.3\times10^7\) | 7553 |
| 55 | \(1.1\times10^8\) | 17769 |
| 60 | \(2.6\times10^8\) | 42017 |
| 65 | \(8.5\times10^8\) | 69408 |

**Suboptions.**

| Suboption | Default | Meaning |
|---|---|---|
| `"B1"` | climb the levels | Run at this stage-1 bound. Alone, it uses the curve count of the nearest level. |
| `"Curves"` | that level's count | Number of curves. Alone, without `"B1"`, it runs at \(B_1 = 11000\). |
| `"B2"` | GMP-ECM's default for \(B_1\) | stage-2 bound |
| `"Sigma"` | drawn from the \(n\)-seeded generator | Suyama parameter of the first curve. Curve \(i\) uses \(\sigma + i\). |

With neither `"B1"` nor `"Curves"`, the method **climbs the level table**. It
runs every level up to \(\lceil d/2\rceil + 5\) digits, where \(d\) is the
number of digits of \(n\), and stops at the first success. A curve that catches
every prime of \(n\) at once, so that the gcd is \(n\), is replayed with
\(B_1\) divided by 8 until the primes separate.

**Determinism.** Without `"Sigma"`, each curve's \(\sigma\) is drawn from
Mathilda's generator, which is seeded from \(n\) itself. GMP-ECM's own draw
would be seeded by the clock. So a given \(n\) always runs the same curves.
With `"Sigma" -> s`, the curves are Suyama curves \(s, s+1, s+2, \dots\),
which is the classical convention, so you can reproduce a result exactly and
hand the parameter to another ECM program.

**Threads.** For \(n \ge 160\) bits, curves run one per thread (up to 16)
when Mathilda is built with threads. The \(\sigma\) values are drawn in curve
order, and a batch's successes are taken in curve order. So the factor
reported is always that of the lowest-numbered curve that succeeded, however
many threads there are.

**Inside `Automatic`.** After rho, micro-ECM, Fermat, p−1 and p+1, the
automatic pipeline climbs the same table one level per tier, from 15 to 65
digits. At each level it first checks the quadratic-sieve crossover: QS takes
over a composite of at most 100 digits once \(7 \times\) (effort spent, in
digits) \(\ge 2 \times\) (its digit count). Above 100 digits, ECM is the only
general method the pipeline has.

---

## 3. Invoking it

```mathematica
In[2]:= FactorInteger[2^128 + 1, Method -> "EllipticCurve"]
Out[2]= {{59649589127497217, 1}, {5704689200685129054721, 1}}

In[3]:= FactorInteger[2^128 + 1, Method -> {"ECM", "B1" -> 50000, "Curves" -> 10, "Sigma" -> 1}]
Out[3]= {{59649589127497217, 1}, {5704689200685129054721, 1}}

In[4]:= FactorInteger`ECM[2^128 + 1, "B1" -> 50000, "Curves" -> 10, "Sigma" -> 1]
Out[4]= {{59649589127497217, 1}, {5704689200685129054721, 1}}
```

---

## 4. Where it shines

### A 20-digit factor hidden in a 100-digit number

\(p_{20} = 40000000000000000019\) has no smooth \(p \pm 1\). The largest prime
factors of \(p_{20} - 1\) and \(p_{20} + 1\) have 11 and 10 digits, so p±1
cannot see it. Rho would need about \(6\times10^{9}\) steps. Its cofactor is an
81-digit prime. With no suboptions, `ECM` climbs the levels and finds
\(p_{20}\) in about a second:

```mathematica
In[5]:= p20 = NextPrime[4 10^19 + 1]; c80 = NextPrime[10^80 + 12345]; n = p20 c80; IntegerLength[n]
Out[5]= 100

In[6]:= {FactorInteger[p20 - 1][[-1]], FactorInteger[p20 + 1][[-1]]}
Out[6]= {{50935645921, 1}, {2550536291, 1}}

In[7]:= FactorInteger`ECM[n]
Out[7]= {{40000000000000000019, 1}, {100000000000000000000000000000000000000000000000000000000000000000000000000012397, 1}}
```

### One curve, three sizes of \(n\): the cost follows the factor

Suyama's \(\sigma = 18\) at \(B_1 = 11000\) happens to give a curve whose
order modulo \(p_{20}\) is smooth enough. (It was found by scanning
\(\sigma = 1, \dots, 300\) once. Four values worked: 18, 56, 112 and 254.) The
group \(E(\mathbb F_{p_{20}})\) depends only on \(p_{20}\) and \(\sigma\), not on
the cofactor. So this **single curve** splits every multiple of \(p_{20}\):
60 digits, 100 digits or 220 digits. The time changes only with the cost of
arithmetic modulo \(n\), from about 30 ms to about 120 ms:

```mathematica
In[8]:= ecm1 = {"ECM", "B1" -> 11000, "Curves" -> 1, "Sigma" -> 18};
In[9]:= FactorInteger[p20 NextPrime[10^40 + 12345], Method -> ecm1]
Out[9]= {{40000000000000000019, 1}, {10000000000000000000000000000000000012357, 1}}

In[10]:= FactorInteger[n, Method -> ecm1]
Out[10]= {{40000000000000000019, 1}, {100000000000000000000000000000000000000000000000000000000000000000000000000012397, 1}}

In[11]:= Length[FactorInteger[p20 NextPrime[10^200], Method -> ecm1]]
Out[11]= 2
```

The same holds at the 15-digit level. Twenty-five default curves at
\(B_1 = 2000\) pull a 15-digit prime out of a 75-digit and a 215-digit number
alike:

```mathematica
In[12]:= p15 = NextPrime[3 10^14 + 7]; FactorInteger[p15 NextPrime[10^60], Method -> {"ECM", "B1" -> 2000, "Curves" -> 25}][[1]]
Out[12]= {300000000000089, 1}

In[13]:= FactorInteger[p15 NextPrime[10^200], Method -> {"ECM", "B1" -> 2000, "Curves" -> 25}][[1]]
Out[13]= {300000000000089, 1}
```

### Reproducible runs

The neighbouring \(\sigma = 17\) fails. A run of 17 curves starting at
\(\sigma = 2\) covers \(\sigma = 2, \dots, 18\), so it includes the lucky
curve and succeeds:

```mathematica
In[14]:= FactorInteger[n, Method -> {"ECM", "B1" -> 11000, "Curves" -> 1, "Sigma" -> 17}]
Out[14]= {{4000000000000000001900000000000000000000000000000000000000000000000000000000495880000000000000235543, 1}}

In[15]:= FactorInteger[n, Method -> {"ECM", "B1" -> 11000, "Curves" -> 17, "Sigma" -> 2}]
Out[15]= {{40000000000000000019, 1}, {100000000000000000000000000000000000000000000000000000000000000000000000000012397, 1}}
```

### Beyond the reach of the quadratic sieve

At 121 digits, \(p_{20}\cdot\mathrm{NextPrime}(10^{100})\) is past the
automatic pipeline's 100-digit QS limit, so its ECM tiers do the work. It
takes about a second on the loaded reference machine:

```mathematica
In[16]:= FactorInteger[p20 NextPrime[10^100]]
Out[16]= {{40000000000000000019, 1}, {10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000267, 1}}
```

The FactorInteger shootout (`benchmarks/98-factorinteger-shootout/REPORT.md`)
puts the automatic pipeline at about 1–3 s for a hidden 20-digit factor, 3–12 s
for 25 digits, 26–54 s for 30 digits, and 3–6 minutes for 35 digits, on an
8-core i9 with 16 threads. That is the \(L_p[1/2,\sqrt2]\) growth, against
PARI and FLINT timing out at 30–35 digits.

---

## 5. Where it struggles

### \(B_1\) too low for the factor

The 15-digit level (\(B_1 = 2000\), whose 25 curves `"B1" -> 2000` selects)
is very unlikely to find a 20-digit factor. The probability that one curve's
group order is 2000-smooth near \(4\times10^{19}\) is tiny:

```mathematica
In[17]:= FactorInteger[n, Method -> {"ECM", "B1" -> 2000}]
Out[17]= {{4000000000000000001900000000000000000000000000000000000000000000000000000000495880000000000000235543, 1}}
```

### Balanced semiprimes

ECM is the wrong tool when the smallest factor has half the digits of \(n\). A
59-digit product of two 30-digit primes would need the 30-digit level, about
430 curves at \(B_1 = 250000\), and minutes of work, with no guarantee of
success. A 20-curve sample at the 20-digit level shows the shape of the
search. The [quadratic sieve](quadratic-sieve.md) factors the same number
deterministically in about half a second:

```mathematica
In[18]:= bal = NextPrime[10^29] NextPrime[3 10^29]; FactorInteger[bal, Method -> {"ECM", "B1" -> 11000, "Curves" -> 20}]
Out[18]= {{30000000000000000000000000096400000000000000000000000002233, 1}}

In[19]:= FactorInteger[bal, Method -> "QS"]
Out[19]= {{100000000000000000000000000319, 1}, {300000000000000000000000000007, 1}}
```

The automatic pipeline encodes exactly this trade-off. For a 59-digit
composite, the crossover \(7 \times 20 \ge 2 \times 59\) is met after the
20-digit level, and the number goes to QS.

!!! note "`TimeConstrained` and threaded ECM"
    `TimeConstrained` counts the process's CPU time (`ITIMER_PROF`, as
    documented in `src/core.c`), which accumulates on every worker thread. With threaded curves (\(n \ge 160\) bits), a limit of
    \(t\) seconds can therefore expire after about \(t/16\) seconds of wall
    time on a 16-thread machine. For example,
    `TimeConstrained[FactorInteger`ECM[bal], 3]` returns `$Aborted` after
    about 0.3 s. Prefer limiting the work with `"Curves"` when you want a
    bounded ECM run.

### Small numbers

On machine words, micro-ECM is excellent:

```mathematica
In[20]:= FactorInteger`ECM[NextPrime[2^31] NextPrime[2^31 + 2^30]]
Out[20]= {{2147483659, 1}, {3221225473, 1}}
```

Between \(2^{64}\) and about \(2^{100}\), GMP-ECM's per-curve setup dominates.
There the automatic pipeline's two-limb rho and two-limb micro-ECM
(`fi_u128.c`), followed by QS, are much cheaper than pinned `ECM`.

---

## 6. Limitations

- Probabilistic. A level's curve count finds a factor of that size with
  probability about \(1 - 1/e\), not certainty. A factor can survive a level
  meant for it, and turn up in the next one.
- \(L_p[1/2,\sqrt2]\) is still exponential in \(\sqrt{\log p}\). Factors past
  about 40 digits take hours. The automatic pipeline's table stops at 65
  digits. A composite that survives every level is returned unfactored with
  `FactorInteger::nofac`.
- `"Sigma"` selects Suyama curves. Without it, GMP-ECM's faster batch
  parametrisation is used, so a run with `"Sigma"` is reproducible but not
  identical to the default run.
- `"B2"` is ignored by micro-ECM on machine words.

---

## See also

- [FactorInteger methods](index.md): the ECM levels inside the automatic
  pipeline, and the QS crossover.
- [Pollard p−1 and Williams p+1](pminus1-pplus1.md): the fixed-group
  ancestors of ECM.
- [Quadratic sieve](quadratic-sieve.md): the method for balanced composites
  up to 100 digits.
- [Performance](performance.md): Mathilda's FactorInteger against PARI, FLINT
  and Mathematica.
- `?FactorInteger`ECM` in the REPL.
