# FactorInteger

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FactorInteger[n]`**

gives a list of the prime factors of the integer n, together with their exponents: {{p1, e1}, {p2, e2}, ...}.  For a negative n the first entry is {-1, 1}; for a rational n the denominator's primes have negative exponents; a Gaussian integer or Gaussian rational is factored over Z\[i\] (unit first, primes in the first quadrant).

**`FactorInteger[n, k]`**

does a partial factorisation: at most k entries, the last of which may be an unfactored remainder.

**`FactorInteger[n, Automatic]`**

is FactorInteger\[n\].

<details>
<summary>Notes</summary>

Options: GaussianIntegers -\> True factors over the Gaussian integers. Method -\> m chooses the algorithm; m is a name below, Automatic, or {name, "Sub" -\> v, ...} with suboptions.  Each method is also callable directly as FactorInteger\`name\[n, "Sub" -\> v, ...\] (see ?FactorInteger\`name). Methods:

```text
Automatic                  word engine below 2^64 (micro-ECM, rho, SQUFOF);
                           above: trial division to 2^16, rho, micro-ECM,
                           Fermat probe, p-1, p+1, then ECM levels, switching
                           to the quadratic sieve for composites up to 100 digits
"TrialDivision"            batched trial division        ("Bound")
"PollardRho"               Pollard-Brent rho             ("Iterations", "Starts")
"Fermat"                   Fermat, residue-sieved        ("Iterations")
"Lehman"                   Lehman, O(n^(1/3))            ("Iterations")
"ShanksSquareForms"        SQUFOF; alias "SQUFOF"        ("Iterations")
"CFRAC"                    continued fractions           ("Iterations")
"Dixon"                    random squares                ("Iterations")
"QuadraticSieve"           SIQS; aliases "QS", "SIQS"
"PollardPMinus1"           p-1; alias "PollardP-1"       ("B1", "B2")
"WilliamsPPlus1"           p+1; alias "WilliamsP+1"      ("B1", "B2")
"ECM"                      elliptic curves (GMP-ECM)     ("B1", "B2", "Curves", "Sigma")
"BlakeRationalBaseDescent" rational base descent         ("Base", "Variant", ...)
```

A composite no factor of which is found within a method's bounds is returned with exponent 1 and FactorInteger::nofac.

</details>

## Examples (25)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= FactorInteger`BlakeRationalBaseDescent[1001, "Base" -> 22/7]
Out[1]= {{7, 1}, {11, 1}, {13, 1}}

In[2]:= FactorInteger`BlakeRationalBaseDescent[742245361114170698080188160556986382238609541182118959142360194921632978579427957105431]
Out[2]= {{623347347967, 1}, {1190741187132579502456731851112421085210615890827643793135496647841845583593, 1}}

In[3]:= Select[Quiet[FactorInteger`BlakeRationalBaseDescentAllHits[2^1041 + 1, "Base" -> 2]][[All, 1]], # < 10^120 &]
Out[3]= {3, 2083, 95562442332919646317117537304253622533190207882011713489066201641121786503686867002917439712921903606443}

In[4]:= First /@ FactorInteger`BlakeRationalBaseDescentLattice[195539142822724777401163665279552940381006352574199957151137, "Base" -> 22/7]
Out[4]= {296454907933147663952202071141, 659591518271709676788151702157}

In[5]:= First /@ FactorInteger`BlakeRationalBaseDescentCoppersmith[195539142822724777401163665279552940381006352574199957151137, "Base" -> 22/7, "Coefficient" -> 3, "Exponent" -> 59]
Out[5]= {296454907933147663952202071141, 659591518271709676788151702157}
```

### Scope (5)

```mathematica
In[6]:= FactorInteger[-12]
Out[6]= {{-1, 1}, {2, 2}, {3, 1}}

In[7]:= FactorInteger[3/4]
Out[7]= {{2, -2}, {3, 1}}

In[8]:= FactorInteger[100, 2]
Out[8]= {{2, 2}, {25, 1}}

In[9]:= FactorInteger`ECM[2^128 + 1]
Out[9]= {{59649589127497217, 1}, {5704689200685129054721, 1}}

In[10]:= FactorInteger`QuadraticSieve[NextPrime[10^29] NextPrime[3 10^29]]
Out[10]= {{100000000000000000000000000319, 1}, {300000000000000000000000000007, 1}}
```

### Options (1)

```mathematica
In[11]:= FactorInteger[5, GaussianIntegers -> True]
Out[11]= {{-I, 1}, {1 + 2*I, 1}, {2 + I, 1}}
```

### Applications (14)

The Mersenne number Cole factored by hand in 1903

```mathematica
In[12]:= FactorInteger[2^67 - 1]
Out[12]= {{193707721, 1}, {761838257287, 1}}
```

A Gaussian rational is factored over Z[i]: 3 + 4 I = (2 + I)^2

```mathematica
In[13]:= FactorInteger[(3 + 4 I)/(1 + I)]
Out[13]= {{1 + I, -1}, {2 + I, 2}}
```

At most three entries; the last is the unfactored remainder

```mathematica
In[14]:= FactorInteger[10^30 + 1, 3]
Out[14]= {{61, 1}, {101, 1}, {162311313098522967050803441, 1}}
```

A balanced 60-digit semiprime: Automatic's ECM tiers hand over to the quadratic sieve

```mathematica
In[15]:= FactorInteger[300000000000000000000000000966100000000000000000000000029029]
Out[15]= {{100000000000000000000000000319, 1}, {3000000000000000000000000000091, 1}}
```

Rho needs about Sqrt[p] steps for the 13-digit factor p

```mathematica
In[16]:= FactorInteger`PollardRho[2^101 - 1]
Out[16]= {{7432339208719, 1}, {341117531003194129, 1}}
```

Two 21-digit primes about 10^7 apart: Fermat finds them at once

```mathematica
In[17]:= FactorInteger`Fermat[10000000000001000028200000000000390009477]
Out[17]= {{100000000000000000039, 1}, {100000000000010000243, 1}}
```

Square forms, O(n^(1/4)) work

```mathematica
In[18]:= FactorInteger`ShanksSquareForms[3000000000130000000000507]
Out[18]= {{1000000000039, 1}, {3000000000013, 1}}
```

The same semiprime by continued fractions, Morrison and Brillhart

```mathematica
In[19]:= FactorInteger`CFRAC[3000000000130000000000507]
Out[19]= {{1000000000039, 1}, {3000000000013, 1}}
```

A 40-digit semiprime by the self-initialising quadratic sieve

```mathematica
In[20]:= FactorInteger[3000000000000000015830000000000000002703, Method -> "QuadraticSieve"]
Out[20]= {{10000000000000000051, 1}, {300000000000000000053, 1}}
```

P - 1 = 2^4 3 5 7 11 13 17^2 19 ... 59 for the 21-digit factor, so a tiny B1 suffices

```mathematica
In[21]:= FactorInteger[1049760769253203207210000000000000000012702105307963758807241, Method -> {"PollardPMinus1", "B1" -> 300}]
Out[21]= {{104976076925320320721, 1}, {10000000000000000000000000000000000000121, 1}}
```

Here p + 1 is smooth instead, and p - 1 is not

```mathematica
In[22]:= FactorInteger`WilliamsPPlus1[1046620906533602399990000000000000000012664112969056589039879]
Out[22]= {{104662090653360239999, 1}, {10000000000000000000000000000000000000121, 1}}
```

The Fermat number F7, which Morrison and Brillhart first split with CFRAC

```mathematica
In[23]:= FactorInteger`ECM[2^128 + 1]
Out[23]= {{59649589127497217, 1}, {5704689200685129054721, 1}}
```

The larger factor is near 3 (22/7)^150; the base 22/7 is found by search

```mathematica
In[24]:= FactorInteger`BlakeRationalBaseDescent[742245361114170698080188160556986382238609541182118959142360194921632978579427957105431]
Out[24]= {{623347347967, 1}, {1190741187132579502456731851112421085210615890827643793135496647841845583593, 1}}
```

No prime below the bound divides it: the composite comes back whole, with FactorInteger::nofac

```mathematica
In[25]:= FactorInteger`TrialDivision[2^67 - 1, "Bound" -> 1000]
Out[25]= {{147573952589676412927, 1}}
```

## Algorithm

fi_builtin.c — FactorInteger[...] and the FactorInteger`Method[...] entry points: argument and option parsing, rational and Gaussian inputs, partial factorisation, messages, output, registration and docstrings.

```text
  FactorInteger[n]                       complete factorisation
  FactorInteger[n, k]                    at most k entries; the last may be
                                         an unfactored remainder
  FactorInteger[n, Automatic]            same as FactorInteger[n]
  FactorInteger[n, Method -> m]          m a method name, Automatic, or
                                         {"Name", "Sub" -> v, ...}
  FactorInteger[n, GaussianIntegers -> True]   factor over Z[i]
  FactorInteger`Name[n, ...]             = FactorInteger[n, Method -> {"Name", ...}]
```

n may be an integer, a rational (the denominator's primes get negative exponents), or a Gaussian integer / Gaussian rational (factored over Z[i] automatically, unit first, primes in the first quadrant sorted by real then imaginary part — Mathematica's normalisation).

## Implementation notes

**Algorithm.** `builtin_factorinteger` (`fi_builtin.c`) has a heap-free fast path for a machine integer `n ≥ 2`: `fi_u64_factor_array` factors it entirely in the **word engine** (`fi_u64.c`) — inverse-multiply trial division by the odd primes below 2^10, a deterministic Miller–Rabin test, then `u64_split`: Pollard–Brent rho below 2^40, **micro-ECM** above (Montgomery curves with Suyama's parametrisation, prime-power ladder stage 1, baby-step/giant-step stage 2), with SQUFOF and long rho runs as backstops — all on single-word Montgomery arithmetic. Everything else goes through `factorinteger_core` (rationals, Gaussian inputs, `k`, options) into the driver `fi_factor` (`fi_core.c`), a work stack of `(composite, multiplicity, stage)` items: powers of two (and, for `Automatic`, every prime up to 2^16, by batched trial division) come off once; each item is then a word (→ word engine), a prime (→ output), a perfect power (→ `primitive_root`, multiplicity scaled), or is handed to the method's splitter and both parts are pushed back. **`Automatic`'s splitter** `split_auto` is a tier table resumed per item: cheap rho (2^13 steps; two-limb Montgomery `fi_u128_rho` below 2^126) → two-limb micro-ECM (8 curves at B1 = 150, below 2^126) → a Fermat probe (2^14 steps, for two factors near √n) → one Pollard p−1 (B1 = min(2000·digits, 1.5·10^5)) → one Williams p+1 (B1 = min(1500·digits, 10^5); three seeds from 70 digits, one below) → ECM at the 15, 20, …, 65-digit levels of GMP-ECM's table. Before each tier from p−1 on, the **quadratic sieve** takes over once the effort already spent (in digits) reaches 2/7 of the composite's digit count — yafu's crossover — for composites of at most 100 digits. A tier advances only on failure; after a split both parts restart at the tier that succeeded, so small factors keep peeling off with cheap rho. **Named methods** (`Method -> "Name"` or the head ``FactorInteger`Name``) keep their own algorithm at every size, including word inputs (`fi_word_method_array`, allocation-free): `TrialDivision` (word-sized products of primes, one bignum remainder per 2–4 primes, segmented sieve above 2^20; default bound √n when below 10^7, else 10^7); `PollardRho` (Brent cycle finding, 256 products per gcd, backtracking); `Fermat` (residue tables modulo nine moduli; ~1 candidate in 3000 pays a bignum square test); `Lehman` (complete below 2^42, budgeted above); `ShanksSquareForms` (16 Gower–Wagstaff multipliers; word arithmetic while kN < 2^124, declined beyond); `CFRAC` and `Dixon` (Knuth–Schroeppel multiplier, early-abort batched trial division, single large primes; Dixon draws squares near √(kN) so residues are square-root sized); `QuadraticSieve` (native SIQS: Gray-code self-initialisation, 32 KB block sieve with bucket sieving for primes ≥ 2^15, no sieving below 32/64, single large primes plus double large primes once kN ≥ 2^240); `PollardPMinus1` / `WilliamsPPlus1` / `ECM` (`fi_ecm.c`); and `BlakeRationalBaseDescent` (`fi_rbd*.c`): descents `Q ← ⌊bQ/a⌋` over supplied or Stern–Brocot-enumerated bases a/b in six variants — `Exact` (default: scaled exact descent on S·n with a four-value window, deferred testing above √(cmax n) reached by one division, batched gcds), `Classic` (the paper's R1–R4), `AllHits` (every hit refined into coprime divisors; truncation divisibility for Cunningham numbers), `Lattice` (Howgrave-Graham approximate-divisor detector via LLL, integer roots by Hensel lifting), `TwoPass`, and `Coppersmith` (one detector call at known a/b, c, k).

**Data structures.** `FiFactors` is a growable array of `{mpz_t p, int64 e}` merged on insert and sorted by `mpz_cmp` at the end (negative `e` for denominator primes of a rational). `FiCtx` carries the method, `k` limit, suboptions (`Bound`, `Iterations`, `Starts`, `B1`, `B2`, `Curves`, `Sigma`, and `FiRbdOpts` for the RBD family) and an xorshift64* generator that `fi_ctx_seed` seeds from the limbs of `n`; every random choice (rho constants, ECM sigmas — GMP-ECM's own RNG is reseeded from it — p+1 seeds, QS polynomials, Dixon draws) comes from it, so the answer, the time and any `nofac` depend on `n` alone. The congruence-of-squares methods share `fi_gf2.c`: a relation store whose one- and two-large-prime partials are edges of a union-find graph on the large primes (every cycle is a full relation), singleton-column pruning, then dense bit-packed elimination below 1000 columns and **block Lanczos** (`fi_lanczos.c`, 64-bit blocks, Montgomery 1995) above. Worker threads write relations to `FiRelLog`s replayed by the main thread in a fixed order. ECM backends, best first: system **GMP-ECM** (batch-parametrised curves, fast stage 2; a user `"Sigma"` is a Suyama parameter), **FLINT** (`fmpz_factor_ecm`, `fmpz_factor_pp1`) when built `USE_ECM=0`, and a native mpz implementation otherwise; QS and every other method are native. Gaussian inputs factor the norm, then split each p ≡ 1 (mod 4) by Cornacchia two-squares and rotate primes into the first quadrant, unit first.

**Complexity / limits.** Below 2^64 the word engine finishes in microseconds. Rho costs ~√p for the smallest prime p; ECM's cost is subexponential in the size of the *factor*, not of `n`; SIQS is `L_n[1/2, 1]` and is the general-purpose finisher for 30–100 digits; p∓1 succeed exactly when p∓1 is B1-smooth up to one prime ≤ B2 (prime powers must also be ≤ B1); Fermat/Lehman are instant only for factors close to √n (or to a small rational ratio); RBD succeeds only when the larger factor is near c (a/b)^k and is no threat to unstructured moduli. Every splitter carries a budget, so nothing loops unboundedly: a composite that survives its method's bounds (or Automatic's whole table, which stops after the 65-digit ECM level on a composite above 100 digits) is returned with exponent 1 and `FactorInteger::nofac`. Threads (`fi_par.c`, builds with `MATHILDA_THREADS`, up to 16): QS polynomial families from kN ≥ 2^150, ECM curves from n ≥ 2^160, p+1 seeds, CFRAC/Dixon candidate trial division and RBD base search run in parallel regions bracketed by `tc_async_defer`, with results taken in index order so the answer does not depend on the thread count; a `TimeConstrained` abort is taken on the main thread at the end of the current short region (and budgets, being CPU time, drain faster). `FactorInteger` is `Listable` and `Protected`; it has no packed-array kernel or `Compile[]` lowering — it returns a ragged list of exact pairs, not a machine array.

- `Listable`, `Protected`; `Options[FactorInteger]` is `{GaussianIntegers -> False, Method -> Automatic}`.
- Perfect powers are reduced to a primitive base before any search, so `(2^61 - 1)^2` factors correctly.
- `EulerPhi` and `facint_factor_complete` decline (stay unevaluated or report incompleteness) rather than
  build an answer on a composite that survived the search.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Quiet](../../control-flow/Quiet/), [Check](../../control-flow/Check/), [TimeConstrained](../../time-and-date/TimeConstrained/), [EulerPhi](../../number-theory/EulerPhi/)

- J. M. Pollard, "A Monte Carlo method for factorization", BIT 15 (1975).
- R. P. Brent, "An improved Monte Carlo factorization algorithm", BIT 20 (1980).
- J. M. Pollard, "Theorems on factorization and primality testing", Proc. Cambridge Phil. Soc. 76 (1974) (p−1 method).
- H. C. Williams, "A p+1 method of factoring", Math. Comp. 39 (1982).
- H. W. Lenstra Jr., "Factoring integers with elliptic curves", Annals of Mathematics 126 (1987).
- P. L. Montgomery, "Speeding the Pollard and elliptic curve methods of factorization", Math. Comp. 48 (1987).
- R. S. Lehman, "Factoring large integers", Math. Comp. 28 (1974).
- J. E. Gower and S. S. Wagstaff Jr., "Square form factorization", Math. Comp. 77 (2008).
- M. A. Morrison and J. Brillhart, "A method of factoring and the factorization of F_7", Math. Comp. 29 (1975) (CFRAC).
- J. D. Dixon, "Asymptotically fast factorization of integers", Math. Comp. 36 (1981).
- S. Contini, "Factoring integers with the self-initializing quadratic sieve", M.Sc. thesis, University of Georgia (1997).
- P. L. Montgomery, "A block Lanczos algorithm for finding dependencies over GF(2)", EUROCRYPT '95, LNCS 921 (1995).
- N. Howgrave-Graham, "Approximate integer common divisors", CaLC 2001, LNCS 2146 (2001).
- H. Cohen, A Course in Computational Algebraic Number Theory, Springer (1993) (integral LLL, Algorithm 2.6.7).
- S. Blake, "Rational Base Descent: A Deterministic Algorithm for Factoring Structured Semiprimes" and "Lattice-Augmented Rational Base Descent" (manuscripts).
- Source: [`src/numbertheory/factorinteger/fi_builtin.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/factorinteger/fi_builtin.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_factorinteger.c`](https://github.com/stblake/mathilda/blob/main/tests/test_factorinteger.c)
- Tests: [`tests/test_partitionsp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_partitionsp.c)

## Notes & additional examples

### Notes

- `FactorInteger[n]` gives `{prime, exponent}` pairs in ascending order; a negative `n` starts with `{-1, 1}`, a rational's denominator primes get negative exponents, `0` gives `{{0, 1}}` and `1` gives `{{1, 1}}`. `FactorInteger[n, k]` returns at most `k` entries, the last of which may be composite (Mathematica's contract).
- Gaussian integers and Gaussian rationals factor over Z[i] automatically, as does any integer under `GaussianIntegers -> True`; the unit comes first when it is not 1, and primes are put in the first quadrant, sorted by real then imaginary part, as in Mathematica.
- `Method -> m` takes `Automatic`, a method name, or `{name, "Sub" -> v, ...}`. The names are `"TrialDivision"`, `"PollardRho"`, `"Fermat"`, `"Lehman"`, `"ShanksSquareForms"` (alias `"SQUFOF"`), `"CFRAC"` (alias `"ContinuedFraction"`), `"Dixon"`, `"QuadraticSieve"` (aliases `"QS"`, `"SIQS"`, `"MPQS"`), `"PollardPMinus1"` (aliases `"PollardP-1"`, `"PollardPM1"`, `"P-1"`), `"WilliamsPPlus1"` (aliases `"WilliamsP+1"`, `"WilliamsPP1"`, `"P+1"`), `"ECM"` (alias `"EllipticCurve"`) and `"BlakeRationalBaseDescent"`.
- Context-qualified entry points have no pages of their own, so they are listed here. ``FactorInteger`Name[n, "Sub" -> v, ...]`` is `FactorInteger[n, Method -> {"Name", "Sub" -> v, ...}]` and takes the same `k` and `GaussianIntegers` arguments. The heads are ``FactorInteger`TrialDivision``, ``FactorInteger`PollardRho``, ``FactorInteger`Fermat``, ``FactorInteger`Lehman``, ``FactorInteger`ShanksSquareForms`` (and ``FactorInteger`SQUFOF``), ``FactorInteger`CFRAC``, ``FactorInteger`Dixon``, ``FactorInteger`QuadraticSieve`` (and ``FactorInteger`QS``, ``FactorInteger`SIQS``; there is no ``FactorInteger`MPQS`` head), ``FactorInteger`PollardPMinus1``, ``FactorInteger`WilliamsPPlus1``, ``FactorInteger`ECM`` and ``FactorInteger`BlakeRationalBaseDescent``. Each RBD variant also has a head, ``FactorInteger`BlakeRationalBaseDescent<Variant>`` for `Exact`, `Classic`, `AllHits`, `Lattice`, `TwoPass` and `Coppersmith`. `?FactorInteger`Name` prints each method's docstring.
- Suboptions: `"Bound"` (TrialDivision, at most 2^40); `"Iterations"` (PollardRho, Fermat, Lehman, ShanksSquareForms, CFRAC, Dixon); `"Starts"` (PollardRho); `"B1"` and `"B2"` (PollardPMinus1, WilliamsPPlus1, ECM; B1 defaults to 10^6 for p∓1, and B2 to GMP-ECM's choice); `"Curves"` and `"Sigma"` (ECM, where a given Sigma is a Suyama parameter). BlakeRationalBaseDescent takes `"Base"`, `"Variant"`, `"Order"`, `"Height"`, `"Bases"`, `"MaxCoefficient"`, `"Multiplier"`, `"Tolerance"`, `"Beta"`, `"Multiplicity"`, `"Dimension"`, `"BaseRange"`, `"Coefficient"` and `"Exponent"`. The Coppersmith variant needs one Base, a Coefficient and an Exponent.
- With `Automatic`, a machine integer is factored entirely by the word engine. Above 2^64 the tiers are: trial division to 2^16, then rho, then micro-ECM, then a Fermat probe, then one p−1 and one p+1 pass, then ECM at levels of 15 to 65 digits. The quadratic sieve takes over for composites of up to 100 digits once the effort spent reaches 2/7 of the digit count. A named method keeps its own algorithm at every size, so ``FactorInteger`Fermat`` really runs Fermat even on a word.
- Results are deterministic. Every random choice is seeded from `n`, and threaded methods take their results in a fixed order, so the answer, the time and any message depend on `n` alone, never on earlier calls or on the thread count.
- Messages (all respect `Quiet` and `Check`): `nofac` (a composite survived the search bounds and is returned with exponent 1), `method` (unknown method), `mopt` (bad suboption), `optx` (unknown option), `exact` (inexact or symbolic argument), and `rbdcop` (incomplete Coppersmith setup). `FactorInteger` is `Listable` and `Protected`.
- See also the advanced tutorial [FactorInteger methods](../../tutorials/factorinteger-methods/index.md).
