# Trial division

`Method -> "TrialDivision"` factors \(n\) by dividing it by every prime up to a
bound \(B\), and by nothing else. It is the oldest factoring algorithm. It is
also part of every other method in Mathilda, because it is the fastest way to
remove the small primes that the clever algorithms are worst at.

```mathematica
In[1]:= FactorInteger`TrialDivision[2^64 - 1]
Out[1]= {{3, 1}, {5, 1}, {17, 1}, {257, 1}, {641, 1}, {65537, 1}, {6700417, 1}}
```

Every prime factor of \(2^{64}-1\) is below \(10^7\), so a plain scan up to the
default bound finds all of them.

---

## 1. The mathematical idea

If \(n\) is composite then it has a prime factor \(p \le \sqrt n\). Trial
division tests the primes \(p = 2, 3, 5, \dots\) in order. Each time \(p \mid n\),
it divides \(p\) out as many times as it goes, and records the exponent. The scan
stops at the first of these:

- the bound \(B\) is passed;
- the cofactor \(m\) left over drops to 1;
- \(p^2 > m\). Then \(m\) has no prime factor \(\le p\), so it is prime.

Finding a prime factor \(p\) costs about \(\pi(p) \approx p/\ln p\) divisions. That
is exponential in the *size of the factor*, so the method is perfect for small
factors and hopeless for large ones. A 12-digit factor already needs about
\(3.8\times10^{10}\) divisions. Pollard rho finds the same factor in about
\(10^6\) steps (see [Pollard rho](pollard-rho.md)).

When \(B \ge \sqrt n\), trial division is a *complete* algorithm. It is the only
method here that can prove its answer without a separate primality test.

### References

1. D. E. Knuth, *The Art of Computer Programming*, Vol. 2, *Seminumerical
   Algorithms*, 3rd ed., Addison-Wesley, 1997 — §4.5.4, Algorithm A.
2. R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*,
   2nd ed., Springer, 2005 — §3.2 (sieving) and §5.1 (trial division).
3. S. S. Wagstaff, Jr., *The Joy of Factoring*, AMS Student Mathematical
   Library 68, 2013 — Ch. 5.

---

## 2. How Mathilda realises it

The code is in `src/numbertheory/factorinteger/fi_trial.c`. The main routine is
`fi_trial_divide(n, bound, mult, out)`, and the method entry is the
`FI_TRIAL` branch of `fi_factor()` in `fi_core.c`.

**Batched remainders.** Dividing a bignum by one word costs a full pass over
its limbs. Mathilda does not make that pass once per prime. It groups the
primes into products \(M = p_1 p_2 \cdots p_k < 2^{64}\), with \(k \le 8\), and
computes \(r = n \bmod M\) once per group. After that, each prime costs only
the single-word operation \(r \bmod p_i\). Primes below \(2^{16}\) go four to a
pass, primes below \(2^{21}\) three, and primes up to \(2^{32}\) two. When a
prime divides \(n\), it is removed with `mpz_divexact_ui` until it no longer
divides.

**Prime supply.** All odd primes below \(2^{20}\) come from a table that is
sieved once and shared with the rest of the subsystem (`fi_small_primes`).
Primes above \(2^{20}\) come from a segmented sieve with \(2^{18}\)-wide
segments, which are themselves sieved by the table primes.

**Early stop.** After each batch, if the cofactor fits a word and is smaller
than the square of the largest prime tried, it is prime. It is recorded and the
scan ends at once.

**What happens to the leftover.** Once the scan is over, the cofactor \(m > 1\)
is classified:

| Cofactor | Result |
|---|---|
| prime (BPSW; deterministic Miller–Rabin below \(2^{64}\)) | recorded as \(\{m, 1\}\) |
| a perfect power \(b^e\) with \(b\) prime | recorded as \(\{b, e\}\) |
| anything else | returned unfactored as \(\{m, 1\}\), with `FactorInteger::nofac` |

**Suboption.**

| Suboption | Default | Meaning |
|---|---|---|
| `"Bound"` | \(\lfloor\sqrt n\rfloor\) if that is below \(10^7\), else \(10^7\) | Largest prime tried. Capped at \(2^{40}\). |

With the default bound, every \(n < 10^{14}\) is factored completely. For
larger \(n\), the default amounts to "take every prime below ten million".

**Inside `Automatic`.** Trial division is the first step of the automatic
pipeline (see the [overview](index.md#the-automatic-pipeline)). Below \(2^{64}\)
the word engine tests the odd primes below \(2^{10}\) with a multiply-by-inverse
divisibility test, which is cheaper than a division. Above \(2^{64}\),
`fi_trial_divide` runs once to \(2^{16}\) before any splitter is called.

---

## 3. Invoking it

The three spellings are equivalent:

```mathematica
In[2]:= FactorInteger[30!, Method -> "TrialDivision"]
Out[2]= {{2, 26}, {3, 14}, {5, 7}, {7, 4}, {11, 2}, {13, 2}, {17, 1}, {19, 1}, {23, 1}, {29, 1}}

In[3]:= FactorInteger[2^64 - 1, Method -> {"TrialDivision", "Bound" -> 100}]
Out[3]= {{3, 1}, {5, 1}, {17, 1}, {72340172838076673, 1}}

In[4]:= FactorInteger`TrialDivision[2^64 - 1, "Bound" -> 100]
Out[4]= {{3, 1}, {5, 1}, {17, 1}, {72340172838076673, 1}}
```

With `"Bound" -> 100`, only 3, 5 and 17 are found. The cofactor
\(72340172838076673 = 257\cdot641\cdot65537\cdot6700417\) is composite, so it
is returned unfactored and `FactorInteger::nofac` is printed on the message
channel. Read a `TrialDivision` answer with this in mind: its last entry may be
composite.

Rational and negative inputs work as they do for every method. The
denominator's primes get negative exponents:

```mathematica
In[5]:= FactorInteger[-1001/1024, Method -> "TrialDivision"]
Out[5]= {{-1, 1}, {2, -10}, {7, 1}, {11, 1}, {13, 1}}
```

---

## 4. Where it shines

**Highly composite numbers.** \(30! + 1\) has four prime factors below
\(1.1\times10^6\). The scan divides them out and runs on to the default bound
of \(10^7\). The 16-digit cofactor left at the end passes the primality test,
so it is recorded as a prime:

```mathematica
In[6]:= FactorInteger`TrialDivision[30! + 1]
Out[6]= {{31, 1}, {12421, 1}, {82561, 1}, {1080941, 1}, {7719068319927551, 1}}
```

**A big prime cofactor is free.** The scan never has to reach a large prime
factor, because the primality test settles the cofactor. Here two 7-digit
primes are found, and a 31-digit prime is left over and certified:

```mathematica
In[7]:= n = Prime[10^5] Prime[10^5 + 1] NextPrime[10^30]
Out[7]= 1689259081189000000000000000096287767627773

In[8]:= FactorInteger`TrialDivision[n]
Out[8]= {{1299709, 1}, {1299721, 1}, {1000000000000000000000000000057, 1}}
```

**Perfect-power cofactors.** \(p^3\) with \(p \approx 10^{12}\) is far beyond
the default bound of \(10^7\). The cofactor is not prime, but the perfect-power
test sees that it is the cube of a prime:

```mathematica
In[9]:= FactorInteger`TrialDivision[NextPrime[10^12]^3]
Out[9]= {{1000000000039, 3}}

In[10]:= FactorInteger`TrialDivision[NextPrime[10^6]^5 NextPrime[10^50]]
Out[10]= {{1000003, 5}, {100000000000000000000000000000000000000000000000151, 1}}
```

**Early stop on a large bound.** The cap of \(2^{40}\) is not a cost you pay in
full. \(2^{67}-1\) (Cole's 1903 factorisation) has the 9-digit factor
193707721. Once it is removed, the 12-digit cofactor is below \(193707721^2\),
so the scan stops there and never goes near \(2^{40}\):

```mathematica
In[11]:= FactorInteger[2^67 - 1, Method -> {"TrialDivision", "Bound" -> 2^40}]
Out[11]= {{193707721, 1}, {761838257287, 1}}
```

---

## 5. Where it struggles

**A factor just above the bound.** The default bound is \(10^7\), and
`NextPrime[10^7]` \(= 10000019\) is just above it. That factor and everything
after it stay in one composite cofactor:

```mathematica
In[12]:= FactorInteger`TrialDivision[NextPrime[10^6] NextPrime[10^7] NextPrime[10^8]]
Out[12]= {{1000003, 1}, {1000001970000133, 1}}
```

**Two medium factors.** A semiprime of two 9-digit primes is out of reach of
the default bound. To find the smaller factor, the scan must reach it: that is
about 5.8 million primes, and roughly a second on a loaded 8-core i9. The
automatic pipeline factors the same number in well under a millisecond, using
the word engine's micro-ECM:

```mathematica
In[13]:= s = NextPrime[10^8] NextPrime[2 10^8]
Out[13]= 20000004700000231

In[14]:= FactorInteger`TrialDivision[s]
Out[14]= {{20000004700000231, 1}}

In[15]:= FactorInteger`TrialDivision[s, "Bound" -> 10^8 + 10]
Out[15]= {{100000007, 1}, {200000033, 1}}

In[16]:= FactorInteger[s]
Out[16]= {{100000007, 1}, {200000033, 1}}
```

**A square of a medium prime times a large prime.** The perfect-power rescue
only applies when the *whole* cofactor is a prime power. \(p^2 q\) with
\(p \approx 10^9\) and \(q \approx 10^{40}\) is neither prime nor a perfect
power, so trial division returns it whole:

```mathematica
In[17]:= FactorInteger[NextPrime[10^9]^2 NextPrime[10^40], Method -> "TrialDivision"]
Out[17]= {{10000000140000000490000000000000000000121000001694000005929, 1}}
```

The cost is linear in the bound, give or take the logarithm. Going from
\(B = 10^7\) to \(B = 10^8\) costs ten times as much, and \(B = 10^9\) about
ten times more again (several seconds). The \(2^{40}\) cap is there to
stop a runaway request, not as a recommended setting. For any factor beyond
about \(10^7\), use [Pollard rho](pollard-rho.md) for factors up to about 15
digits, and [ECM](ecm.md) beyond that.

---

## 6. Limitations

- Exponential in the size of the factor: about \(p/\ln p\) divisions to find
  \(p\). In practice it is only useful below about \(10^8\).
- `"Bound"` is capped at \(2^{40}\).
- A cofactor that is not prime and not a prime power is returned unfactored
  with `FactorInteger::nofac`. There is no fallback to another method. Use
  `Automatic` when you want a complete answer.
- `TrialDivision` has no word-size fast path. Small inputs go through the
  general driver, which is still microseconds, but slower than the word engine
  that `Automatic` uses.

---

## See also

- [FactorInteger methods](index.md): the overview, the automatic pipeline, and
  which method suits which number.
- [Pollard rho](pollard-rho.md): the next step up, \(O(\sqrt p)\) instead of
  \(O(p)\).
- [Fermat and Lehman](fermat-lehman.md): the opposite case, where the two
  factors are close together.
- `?FactorInteger`TrialDivision` in the REPL gives the built-in help string.
