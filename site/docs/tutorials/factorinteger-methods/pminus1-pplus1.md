# Pollard p−1 and Williams p+1

`Method -> "PollardPMinus1"` and `Method -> "WilliamsPPlus1"` are
**special-purpose** methods. They find a prime \(p \mid n\) quickly whenever
\(p - 1\) (respectively \(p + 1\)) is *smooth*, meaning that all its prime
factors are small. The size of \(p\) hardly matters. A 26-digit prime with a
smooth \(p-1\) falls out in a millisecond, when rho would need about
\(10^{13}\) steps and ECM at its cheapest levels would very likely miss it.

```mathematica
In[1]:= L = LCM @@ Range[60]; p = 2 L + 1; {PrimeQ[p], IntegerLength[p]}
Out[1]= {True, 26}

In[2]:= FactorInteger[p - 1]
Out[2]= {{2, 6}, {3, 3}, {5, 2}, {7, 2}, {11, 1}, {13, 1}, {17, 1}, {19, 1}, {23, 1}, {29, 1}, {31, 1}, {37, 1}, {41, 1}, {43, 1}, {47, 1}, {53, 1}, {59, 1}}

In[3]:= q = NextPrime[10^40]; n = p q; IntegerLength[n]
Out[3]= 66

In[4]:= FactorInteger`PollardPMinus1[n, "B1" -> 100]
Out[4]= {{19381424329554463401825601, 1}, {10000000000000000000000000000000000000121, 1}}
```

---

## 1. The mathematical idea

### Pollard p−1 (1974)

By Fermat's little theorem, \(a^{p-1} \equiv 1 \pmod p\) for \(p \nmid a\).
Suppose \(p - 1\) divides a number \(E\) that we can write down without knowing
\(p\). Then \(p \mid a^E - 1\), and

\[
g = \gcd(a^E - 1,\; n)
\]

reveals \(p\). Stage 1 takes \(E = \prod_{\ell \le B_1} \ell^{\lfloor \log_\ell B_1\rfloor}\),
the product of all prime powers up to \(B_1\). It succeeds when \(p-1\) is
\(B_1\)-*powersmooth*. **Stage 2** then allows one more prime \(\ell\) with
\(B_1 < \ell \le B_2\). It tests \(a^{E\ell} - 1\) for every such \(\ell\) at
far less than the cost of raising \(B_1\) to \(B_2\). Overall, the method
finds \(p\) when \(p - 1 = (\text{$B_1$-powersmooth}) \times (\text{at most one prime} \le B_2)\).

It **fails** in two ways. If no prime of \(n\) has a smooth \(p-1\), then
\(g = 1\). If *every* prime of \(n\) is caught at once, then \(g = n\). The
second failure is cured by backing off: replay with a smaller bound, so that
the primes separate.

### Williams p+1 (1982)

The same idea in the group of norm-1 elements of \(\mathbb F_{p^2}\), whose
order is \(p+1\), works through Lucas sequences:

\[
V_0 = 2,\quad V_1 = x_0,\quad V_{2k} = V_k^2 - 2,\quad V_{2k+1} = V_k V_{k+1} - x_0 .
\]

After stage 1, \(\gcd(V_E - 2, n)\) exposes \(p\) when the group order divides
\(E\). The catch is which group you are in. With discriminant
\(\Delta = x_0^2 - 4\), the order is

\[
p - \left(\frac{\Delta}{p}\right) =
\begin{cases} p + 1 & \text{if } \Delta \text{ is a non-residue mod } p,\\
              p - 1 & \text{if } \Delta \text{ is a residue mod } p. \end{cases}
\]

Since \(p\) is unknown, you cannot choose \(x_0\) to land in the
\(p+1\) group. Each seed is a coin flip between "test \(p+1\)" and
"test \(p-1\)". That is why the method runs several seeds.

Both methods are *smoothness lotteries* with a fixed group order. For a random
26-digit prime, \(p \pm 1\) is very unlikely to be \(10^6\)-smooth. [ECM](ecm.md)
removes that limitation by drawing a fresh group order for every curve.

### References

1. J. M. Pollard, "Theorems on factorization and primality testing,"
   *Proceedings of the Cambridge Philosophical Society* **76** (1974) 521–528.
2. H. C. Williams, "A \(p+1\) method of factoring," *Mathematics of
   Computation* **39** (1982) 225–234.
3. P. L. Montgomery, "Speeding the Pollard and elliptic curve methods of
   factorization," *Mathematics of Computation* **48** (1987) 243–264 —
   stage 2.
4. P. L. Montgomery and A. Kruppa, "Improved stage 2 to \(P\pm1\) factoring
   algorithms," *ANTS VIII*, LNCS 5011 (2008) 180–195 — the fast stage 2
   in GMP-ECM.
5. P. Zimmermann and B. Dodson, "20 years of ECM," *ANTS VII*, LNCS 4076
   (2006) 525–542 — the GMP-ECM implementation, including P±1.

---

## 2. How Mathilda realises it

Both methods are in `fi_ecm.c`. The word versions are in `fi_u64.c`. There
are three backends, chosen when Mathilda is built:

| backend | p−1 | p+1 |
|---|---|---|
| **GMP-ECM** (default; this build has ECM 7.0.7) | `ecm_factor` with `ECM_PM1`, start \(a = 3\), fast stage 2 | `ECM_PP1`, fixed seeds, fast stage 2 |
| FLINT (`USE_ECM=0`) | native stage 1 + prime-gap stage 2 | `fmpz_factor_pp1` |
| native (neither) | native | native Lucas ladder, stage 1 only |

**Seeds.** The p+1 seeds are fixed, so a run is reproducible:

| seed | \(x_0\) | \(\Delta = x_0^2 - 4\) up to squares | tests \(p+1\) when |
|---|---|---|---|
| 0 | \(2/7\) | \(-3\) | \(\left(\frac{-3}{p}\right) = -1\), i.e. \(p \equiv 2 \pmod 3\) |
| 1 | \(6/5\) | \(-1\) | \(\left(\frac{-1}{p}\right) = -1\), i.e. \(p \equiv 3 \pmod 4\) |
| 2 | \(3\) | \(5\) | \(\left(\frac{5}{p}\right) = -1\), i.e. \(p \equiv \pm2 \pmod 5\) |

The first two are GMP-ECM's recommended seeds. They make 6 and 4 divide the
group order, so those factors come for free. Seeds run one per thread when
threads are available. Successes are taken in seed order, so the answer does
not depend on the thread count.

**Overshoot.** Before any run, every prime below \(2^{10}\) is split off by
division, because GMP-ECM cannot invert a starting value that shares a factor
with \(n\). On a gcd of \(n\), the run is repeated with \(B_1\) divided by 16
(without stage 2) until the primes separate.

**Machine words.** For \(n < 2^{64}\), both methods use an allocation-free word
walk (`smooth_walk` in `fi_u64.c`). It takes the prime powers up to \(B_1\) in
batches of 32 with one gcd per batch, and replays an overshooting batch one
prime at a time. The base is \(a = 2\) for p−1. The p+1 seeds are three
pseudo-random values \(3 + (r \bmod 1000)\), drawn from the generator seeded
by \(n\). On words there is **no stage 2**, and \(B_1\) is capped at
\(2^{20}\).

**Suboptions.**

| Suboption | Default | Meaning |
|---|---|---|
| `"B1"` | \(10^6\) | stage-1 bound |
| `"B2"` | GMP-ECM's default for \(B_1\) | stage-2 bound. Ignored on words. With the FLINT or native backend, it defaults to \(100\,B_1\). |

**Inside `Automatic`.** After rho, micro-ECM and the Fermat probe, the
automatic pipeline runs one p−1 pass and then one p+1 pass, scaled to the
size \(d\) (in digits) of the composite:

- p−1 at \(B_1 = \min(2000\,d,\ 150000)\);
- p+1 at \(B_1 = \min(1500\,d,\ 100000)\), with **one** seed below 70 digits
  and all three from 70 digits up.

Both use GMP-ECM's default \(B_2\), so a small composite is not charged a
large-number search. Composites up to about 45 digits never reach these tiers,
because the quadratic sieve takes over first (see the
[crossover rule](index.md#the-automatic-pipeline)).

---

## 3. Invoking it

```mathematica
In[5]:= FactorInteger[n, Method -> {"PollardPMinus1", "B1" -> 100}]
Out[5]= {{19381424329554463401825601, 1}, {10000000000000000000000000000000000000121, 1}}

In[6]:= FactorInteger[n, Method -> "P-1"]
Out[6]= {{19381424329554463401825601, 1}, {10000000000000000000000000000000000000121, 1}}
```

The p−1 aliases are `"PollardP-1"`, `"PollardPM1"` and `"P-1"`. The p+1 aliases
are `"WilliamsP+1"`, `"WilliamsPP1"` and `"P+1"`. The aliases are accepted only
as `Method` strings. The heads are ``FactorInteger`PollardPMinus1`` and
``FactorInteger`WilliamsPPlus1``.

---

## 4. Where they shine

### p−1 where nothing else works

The 26-digit factor \(p = 2\,\mathrm{lcm}(1..60) + 1\) of the 66-digit `n` is
far beyond rho. It needs about \(\sqrt p \approx 4\times10^{12}\) steps, and
\(2\times10^{5}\) is nowhere near. It is also unlikely for ECM at a 15-digit
level. Each of those 25 curves has only a small chance of a smooth group order
modulo a 26-digit prime. Both come back with `n` unfactored, while p−1 at
\(B_1 = 100\) succeeded above in about 3 ms:

```mathematica
In[7]:= FactorInteger[n, Method -> {"PollardRho", "Iterations" -> 2 10^5, "Starts" -> 1}]
Out[7]= {{193814243295544634018256010000000000002345152343876090071620897721, 1}}

In[8]:= FactorInteger[n, Method -> {"ECM", "B1" -> 2000, "Curves" -> 25}]
Out[8]= {{193814243295544634018256010000000000002345152343876090071620897721, 1}}
```

### p+1 where p−1 fails

Now take \(p' = 5\,\mathrm{lcm}(1..60) - 1\). Its \(p'+1\) is 59-smooth, but
\(p'-1\) has a 13-digit prime factor. p−1 at \(B_1 = 10^5\) fails. p+1 at
\(B_1 = 100\) succeeds:

```mathematica
In[9]:= pp = 5 L - 1; np = pp q; {IntegerLength[pp], FactorInteger[pp + 1][[-1]], FactorInteger[pp - 1][[-1]]}
Out[9]= {26, {59, 1}, {2842500523517, 1}}

In[10]:= FactorInteger[np, Method -> {"PollardPMinus1", "B1" -> 10^5}]
Out[10]= {{484535608238861585045639990000000000005862880859690225179052243879, 1}}

In[11]:= FactorInteger[np, Method -> {"WilliamsPPlus1", "B1" -> 100}]
Out[11]= {{48453560823886158504563999, 1}, {10000000000000000000000000000000000000121, 1}}
```

The seed table explains it. For \(p'\), the discriminants of seeds 0 and 1 are
non-residues, so those two seeds work in the \(p'+1\) group:

```mathematica
In[12]:= {JacobiSymbol[-3, pp], JacobiSymbol[-1, pp], JacobiSymbol[5, pp]}
Out[12]= {-1, -1, 1}
```

The automatic pipeline finds \(p'\) in its p+1 tier. Its p−1 tier runs first,
at \(B_1 = 2000 \times 66 = 132000\), and fails because of the 13-digit prime.
Then the single p+1 seed, \(x_0 = 2/7\), lands in the \(p'+1\) group because
\(p' \equiv 2 \pmod 3\). The whole call takes a few tens of milliseconds:

```mathematica
In[13]:= FactorInteger[np, Method -> {"PollardPMinus1", "B1" -> 132000}]
Out[13]= {{484535608238861585045639990000000000005862880859690225179052243879, 1}}

In[14]:= FactorInteger[np]
Out[14]= {{48453560823886158504563999, 1}, {10000000000000000000000000000000000000121, 1}}
```

### Stage 2 rescues a single large prime

Give \(p_2 - 1\) one prime \(Q = 1000003\) beyond a small \(B_1\). Stage 1 at
\(B_1 = 200\) cannot reach it. GMP-ECM's default \(B_2\) for such a tiny
\(B_1\) is also too small. Setting `"B2" -> 10^7` puts \(Q\) inside stage 2,
and the factor appears in about 15 ms. Raising \(B_1\) to \(2\times10^6\)
would also work, at many times the cost:

```mathematica
In[15]:= Q = NextPrime[10^6]; k = SelectFirst[Range[1000], PrimeQ[# L Q + 1] &]; p2 = k L Q + 1; n2 = p2 q; FactorInteger[p2 - 1][[-1]]
Out[15]= {1000003, 1}

In[16]:= FactorInteger[n2, Method -> {"PollardPMinus1", "B1" -> 200}]
Out[16]= {{387629649476549041304316109536010000004690318758666243399782224925385721, 1}}

In[17]:= FactorInteger[n2, Method -> {"PollardPMinus1", "B1" -> 200, "B2" -> 10^7}]
Out[17]= {{38762964947654904130431610953601, 1}, {10000000000000000000000000000000000000121, 1}}
```

---

## 5. Where they struggle

### \(B_1\) too small

With \(B_1 = 40\), stage 1 leaves out the primes 41, 43, 47, 53 and 59 of
\(p - 1\), and one factor of 2 (\(2^6 > 40\)). Stage 2 can supply at most one
of them, so the run fails:

```mathematica
In[18]:= FactorInteger[n, Method -> {"PollardPMinus1", "B1" -> 40}]
Out[18]= {{193814243295544634018256010000000000002345152343876090071620897721, 1}}
```

### The p+1 seed lottery

Here is a prime whose \(p+1\) is 71-smooth, but for which **all three**
discriminants \(-3, -1, 5\) are quadratic residues. Every seed then computes
in the \(p-1\) group instead. This \(p_b\) has \(p_b - 1 = 2^6\cdot3\cdot(\text{26-digit prime})\),
so p+1 fails even at \(B_1 = 10^5\), a thousand times the bound its
\(p_b+1\) needs. The construction: \(p_b + 1 = 2 c m\) with \(m = 5\cdot\prod_{7\le\ell\le71}\ell\),
and \(c\) chosen so that \(p_b \equiv 1 \pmod{12}\) and \(p_b \equiv -1 \pmod 5\):

```mathematica
In[19]:= m0 = 5 Times @@ Prime[Range[4, 20]]; c = SelectFirst[Range[1, 3000, 2], Mod[#, 3] != 0 && Mod[# m0, 3] == 1 && PrimeQ[2 # m0 - 1] && Max[FactorInteger[#][[All, 1]]] < 60 &]; pb = 2 c m0 - 1; nb = pb q;
In[20]:= {IntegerLength[pb], FactorInteger[pb + 1][[-1]], {JacobiSymbol[-3, pb], JacobiSymbol[-1, pb], JacobiSymbol[5, pb]}}
Out[20]= {28, {71, 1}, {1, 1, 1}}

In[21]:= FactorInteger[nb, Method -> {"WilliamsPPlus1", "B1" -> 10^5}]
Out[21]= {{24177435972156954975254666890000000000292546975263099155200581469369, 1}}

In[22]:= FactorInteger[pb - 1]
Out[22]= {{2, 6}, {3, 1}, {12592414568831747382945139, 1}}
```

The coin can also land well. If a discriminant is a residue, p+1 becomes a
p−1 run with the same \(B_1\). For the original \(p = 2\,\mathrm{lcm}(1..60)+1\),
all three discriminants are residues, so ``WilliamsPPlus1`` finds it through
\(p - 1\):

```mathematica
In[23]:= {JacobiSymbol[-3, p], JacobiSymbol[-1, p], JacobiSymbol[5, p]}
Out[23]= {1, 1, 1}

In[24]:= FactorInteger[n, Method -> {"WilliamsPPlus1", "B1" -> 100}]
Out[24]= {{19381424329554463401825601, 1}, {10000000000000000000000000000000000000121, 1}}
```

### No stage 2 on machine words

Below \(2^{64}\), the word walk is stage 1 only, and `"B2"` is silently
ignored. Here \(p_w - 1 = k\cdot\mathrm{lcm}(1..20)\cdot10007\). On a bignum,
`"B1" -> 20, "B2" -> 10^5` would find \(p_w\). On a word, \(B_1\) has to cover
10007 itself:

```mathematica
In[25]:= L20 = LCM @@ Range[20]; kw = SelectFirst[Range[1000], PrimeQ[# L20 10007 + 1] &]; pw = kw L20 10007 + 1; qw = 700907; nw = pw qw; {FactorInteger[qw - 1], nw < 2^63}
Out[25]= {{{2, 1}, {350453, 1}}, True}

In[26]:= FactorInteger[nw, Method -> {"PollardPMinus1", "B1" -> 20, "B2" -> 10^5}]
Out[26]= {{1632801510063864347, 1}}

In[27]:= FactorInteger[nw, Method -> {"PollardPMinus1", "B1" -> 10^4 + 10}]
Out[27]= {{700907, 1}, {2329555147921, 1}}
```

### A random prime

Most primes have no smooth \(p \pm 1\). A 26-digit prime \(p\) has \(p-1\)
\(10^6\)-smooth with probability \(\rho(26/6)\) (Dickman's
function, \(\rho(4.3) \approx 2\times10^{-3}\)). Stage 2 helps only modestly. So p±1 is a cheap side bet, one tier
of the automatic pipeline. It is never the main search. That is what
[ECM](ecm.md) is for.

---

## 6. Limitations

- They succeed only when \(p - 1\) (or \(p + 1\)) is \(B_1\)-powersmooth apart
  from one prime up to \(B_2\). For a random \(p\) that is rare.
- p+1 depends on a Legendre symbol per seed. Three fixed seeds can all land in
  the \(p-1\) group, as with \(p_b\) above.
- On machine words there is no stage 2. `"B2"` is ignored, and \(B_1\) is
  capped at \(2^{20}\).
- When every prime of \(n\) is smooth at \(B_1\), the gcd is \(n\). The
  back-off divides \(B_1\) by 16 at a time and drops stage 2, which can hide a
  prime that needed it.

---

## See also

- [FactorInteger methods](index.md): the p±1 tiers in the automatic pipeline.
- [ECM](ecm.md): the generalisation to random group orders.
- [Pollard rho](pollard-rho.md): the general-purpose method for small factors.
- `?FactorInteger`PollardPMinus1` and `?FactorInteger`WilliamsPPlus1` in the
  REPL.
