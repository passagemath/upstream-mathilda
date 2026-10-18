# Pollard rho

`Method -> "PollardRho"` finds a prime factor \(p\) of \(n\) in about
\(\sqrt p\) steps. Its cost depends on the size of the factor, not the size of
\(n\), and it needs no memory. That makes it the standard tool for factors from
about 6 to 15 digits, the range just beyond trial division. Mathilda runs
Brent's improved variant with Montgomery arithmetic on one and two machine
words, and GMP above that.

```mathematica
In[1]:= p12 = NextPrime[10^11 + 12345]; q49 = NextPrime[3 10^48]; n = p12 q49; IntegerLength[n]
Out[1]= 60

In[2]:= FactorInteger`PollardRho[n]
Out[2]= {{100000012349, 1}, {3000000000000000000000000000000000000000000000103, 1}}
```

A 12-digit factor of a 60-digit number comes out in a few tens of
milliseconds. Trial division would need about \(4\times10^{9}\) divisions to
reach it.

---

## 1. The mathematical idea

Iterate a pseudo-random polynomial map modulo \(n\):

\[
x_{i+1} = f(x_i) = x_i^2 + c \pmod n .
\]

Reduced modulo an unknown prime \(p \mid n\), the sequence lives in a set of
size \(p\), so it must repeat. By the birthday paradox, the tail and the
cycle each have expected length about \(\sqrt{\pi p/8}\), so the first repeat
\(x_i \equiv x_j \pmod p\) comes after about
\(\sqrt{\pi p/2} \approx 1.25\sqrt p\) steps. From then on
the sequence mod \(p\) runs round a cycle, and its tail-plus-cycle shape is
the Greek letter \(\rho\). A collision mod \(p\) shows up as

\[
1 < \gcd(x_i - x_j,\; n),
\]

which is usually exactly \(p\), and almost never \(n\) itself.

**Cycle finding.** Pollard used Floyd's tortoise and hare, comparing
\(x_i\) with \(x_{2i}\). Brent (1980) compares \(x_{r}\) with
\(x_{r+1}, \dots, x_{2r}\) for \(r = 1, 2, 4, \dots\), which is about 24%
faster on average. He also **batches the gcds**: he multiplies \(m\) values
of \(|x - y|\) together modulo \(n\) and takes one gcd of the product. If a
batch overshoots, so that the product's gcd is \(n\), the batch is replayed one
step at a time from a saved point (*backtracking*).

The cost is \(O(\sqrt p)\) multiplications modulo \(n\). That is
exponential in the size of the factor, but only half as steep as trial
division:

| factor \(p\) | steps \(\approx\sqrt p\) | trial divisions \(\approx p/\ln p\) |
|---|---|---|
| \(10^{8}\) | \(10^{4}\) | \(5\times10^{6}\) |
| \(10^{12}\) | \(10^{6}\) | \(4\times10^{10}\) |
| \(10^{16}\) | \(10^{8}\) | \(3\times10^{14}\) |
| \(10^{20}\) | \(10^{10}\) | \(2\times10^{18}\) |

### References

1. J. M. Pollard, "A Monte Carlo method for factorization," *BIT* **15** (1975)
   331–334.
2. R. P. Brent, "An improved Monte Carlo factorization algorithm," *BIT*
   **20** (1980) 176–184.
3. R. P. Brent and J. M. Pollard, "Factorization of the eighth Fermat
   number," *Mathematics of Computation* **36** (1981) 627–630.
4. P. L. Montgomery, "Modular multiplication without trial division,"
   *Mathematics of Computation* **44** (1985) 519–521.

---

## 2. How Mathilda realises it

There are three engines, chosen by the size of the composite being split:

| size of \(n\) | engine | file | gcd batch | arithmetic |
|---|---|---|---|---|
| \(< 2^{64}\) | `fi_u64_rho` | `fi_u64.c` | 128 | one-word Montgomery (\(R = 2^{64}\)) |
| \(2^{64}\) to \(2^{126}\) | `fi_u128_rho` | `fi_u128.c` | 128 | two-limb Montgomery (\(R = 2^{128}\)) |
| \(\ge 2^{126}\) | `rho_mpz` | `fi_rho.c` | 256 | GMP `mpz_mul` + `mpz_tdiv_r` |

All three run Brent's cycle finder with batched gcds and backtracking, start at
\(y_0 = 2\), and use \(f(y) = y^2 + c\). Start number \(s\) uses \(c = s\), so
start 1 is \(y^2+1\), start 2 is \(y^2+2\), and so on. The two Montgomery
engines are roughly an order of magnitude faster per step than GMP, because
GMP makes two calls, a normalisation and a quotient for every modular
multiplication on a two-limb number.

**Suboptions.**

| Suboption | Default | Meaning |
|---|---|---|
| `"Iterations"` | \(2^{24}\) on bignums, \(2^{26}\) on words | steps per start |
| `"Starts"` | 4 on bignums, 16 on words | number of polynomials \(y^2 + c\), \(c = 1, 2, \dots\) |

The default budget, \(4\times2^{24}\approx6.7\times10^{7}\) steps on bignums,
finds factors up to about 15 digits with high probability. When every start
fails, the composite is returned unfactored with `FactorInteger::nofac`.

**Inside `Automatic`.** Rho appears three times:

- In the word engine, numbers below \(2^{40}\) are split by rho with
  \(c = 1, \dots, 63\) at \(2^{22}\) steps each.
- In the word engine, the last-resort backstop (after micro-ECM and SQUFOF)
  is rho with \(c = 1, \dots, 1023\) at \(2^{26}\) steps each.
- Above \(2^{64}\), the first tier is one cheap start of \(2^{13}\) steps
  (`AUTO_RHO_ITERS`). It peels off factors up to about 8 digits (\(2^{24}\))
  before anything more expensive runs. After a split, both parts re-enter at
  the same tier, so a number full of small factors keeps going through cheap
  rho.

---

## 3. Invoking it

```mathematica
In[3]:= FactorInteger[n, Method -> "PollardRho"]
Out[3]= {{100000012349, 1}, {3000000000000000000000000000000000000000000000103, 1}}

In[4]:= FactorInteger[n, Method -> {"PollardRho", "Iterations" -> 2 10^6, "Starts" -> 1}]
Out[4]= {{100000012349, 1}, {3000000000000000000000000000000000000000000000103, 1}}

In[5]:= FactorInteger`PollardRho[n, "Iterations" -> 10^5, "Starts" -> 1]
Out[5]= {{300000037047000000000000000000000000000000000010300001271947, 1}}
```

The last call allows only \(10^5\) steps. The expected collision for
\(p \approx 10^{11}\) comes at about \(3\times10^{5}\) steps, so this budget
runs out first, and `n` comes back whole.

---

## 4. Where it shines

**Scaling with the factor, not with \(n\).** Here is the same cofactor with a
growing hidden factor, timed on a loaded 8-core i9 with GMP arithmetic (a
48-digit cofactor) and with the two-limb engine (a 21-digit cofactor):

| hidden factor | 60-digit \(n\) (GMP) | about 30-digit \(n\) (two-limb) |
|---|---|---|
| 6 digits | 0.5 ms | 0.1 ms |
| 8 digits | 0.6 ms | 0.1 ms |
| 10 digits | 15 ms | 3 ms |
| 12 digits | 30–110 ms | 11 ms |
| 14 digits | about 14 s | about 2 s |

Each two extra digits multiply the work by about \(\sqrt{100} = 10\), and the
run-to-run spread is large, because the collision time is random. The two-limb
engine's advantage is the constant factor in the table above.

**Two-limb Montgomery.** \(2^{101}-1\) is a 31-digit Mersenne number with a
13-digit factor. It is below \(2^{126}\), so every step runs on `fi_u128_rho`:

```mathematica
In[6]:= FactorInteger`PollardRho[2^101 - 1]
Out[6]= {{7432339208719, 1}, {341117531003194129, 1}}

In[7]:= FactorInteger`PollardRho[(2^31 - 1) (2^61 - 1)]
Out[7]= {{2147483647, 1}, {2305843009213693951, 1}}
```

**Machine words.** On a 63-bit semiprime with two 32-bit factors, the word
engine needs about \(2^{16}\) steps:

```mathematica
In[8]:= FactorInteger`PollardRho[NextPrime[2^31] NextPrime[2^31 + 2^30]]
Out[8]= {{2147483659, 1}, {3221225473, 1}}
```

**Many small factors.** The driver re-runs the method on every part it splits
off, so rho takes a number apart one factor at a time:

```mathematica
In[9]:= FactorInteger`PollardRho[Product[NextPrime[10^k], {k, 3, 9}]]
Out[9]= {{1009, 1}, {10007, 1}, {100003, 1}, {1000003, 1}, {10000019, 1}, {100000007, 1}, {1000000007, 1}}
```

Brent and Pollard factored the eighth Fermat number this way in 1980.
\(F_8 = 2^{256}+1\) has the 16-digit factor 1238926361552897.
``FactorInteger`PollardRho[2^256 + 1]`` finds it in about 6 s on the loaded
reference machine. That run is not shown here because it is slow. The
automatic pipeline factors \(F_8\) in about 0.2 s with ECM
(`benchmarks/98-factorinteger-shootout/REPORT.md`, row F8).

---

## 5. Where it struggles

**Balanced semiprimes.** A 39-digit product of two 20-digit primes needs about
\(10^{10}\) rho steps. That is 150 times the whole default budget, and an hour
or more of GMP arithmetic. With a small budget the failure comes back at once:

```mathematica
In[10]:= b1 = NextPrime[10^19]; b2 = NextPrime[3 10^19]; m = b1 b2; IntegerLength[m]
Out[10]= 39

In[11]:= FactorInteger[m, Method -> {"PollardRho", "Iterations" -> 10^5, "Starts" -> 2}]
Out[11]= {{300000000000000001940000000000000002091, 1}}
```

With the default budget it simply runs until its \(6.7\times10^{7}\) steps
are spent. `TimeConstrained` is the cheap way to show that:

```mathematica
In[12]:= TimeConstrained[FactorInteger`PollardRho[m], 2]
Out[12]= $Aborted
```

The automatic pipeline does not wait for rho here. Its cheap rho tier fails
in a few milliseconds, and once the effort spent reaches \(2/7\) of the digit
count it hands the number to the quadratic sieve, which takes about 20 ms:

```mathematica
In[13]:= FactorInteger[m]
Out[13]= {{10000000000000000051, 1}, {30000000000000000041, 1}}
```

**Medium factors of big numbers.** Above about 15 digits, the elliptic-curve
method wins clearly. Its cost per curve is fixed, and the number of curves
needed grows only subexponentially with the factor
(\(L_p[1/2, \sqrt2]\) against rho's \(p^{1/2}\)). A 20-digit factor costs
rho about \(10^{10}\) steps and ECM about 74 curves at \(B_1 = 11000\) (see
[ECM](ecm.md)).

---

## 6. Limitations

- The expected work is \(\sqrt p\) for the smallest prime \(p\), with large
  variance. A budget that is "enough on average" fails some of the time. That
  is why there are several starts.
- Every step is a full multiplication modulo \(n\), so the constant grows with
  \(n\). The two Montgomery engines cover only \(n < 2^{126}\).
- Rho cannot exploit structure: smooth \(p \pm 1\) or close factors. Other
  methods handle those (see [p−1 / p+1](pminus1-pplus1.md) and
  [Fermat](fermat-lehman.md)).
- `"Iterations"` counts steps per start. The total work is
  `"Iterations"` \(\times\) `"Starts"`.

---

## See also

- [FactorInteger methods](index.md): where rho sits in the automatic pipeline.
- [Trial division](trial-division.md): the method rho replaces once factors
  pass about \(10^7\).
- [SQUFOF](squfof.md): the word-size backstop that rho shares the word engine
  with.
- [ECM](ecm.md): the method that takes over from rho beyond about 15 digits.
- `?FactorInteger`PollardRho` in the REPL.
