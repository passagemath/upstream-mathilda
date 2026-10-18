### Worked examples

```mathematica
In[1]:= FactorInteger[2^67 - 1]  (* the Mersenne number Cole factored by hand in 1903 *)
```

```mathematica
In[1]:= FactorInteger[(3 + 4 I)/(1 + I)]  (* a Gaussian rational is factored over Z[i]: 3 + 4 I = (2 + I)^2 *)
```

```mathematica
In[1]:= FactorInteger[10^30 + 1, 3]  (* at most three entries; the last is the unfactored remainder *)
```

```mathematica
In[1]:= FactorInteger[300000000000000000000000000966100000000000000000000000029029]  (* a balanced 60-digit semiprime: Automatic's ECM tiers hand over to the quadratic sieve *)
```

```mathematica
In[1]:= FactorInteger`PollardRho[2^101 - 1]  (* rho needs about Sqrt[p] steps for the 13-digit factor p *)
```

```mathematica
In[1]:= FactorInteger`Fermat[10000000000001000028200000000000390009477]  (* two 21-digit primes about 10^7 apart: Fermat finds them at once *)
```

```mathematica
In[1]:= FactorInteger`ShanksSquareForms[3000000000130000000000507]  (* square forms, O(n^(1/4)) work *)
```

```mathematica
In[1]:= FactorInteger`CFRAC[3000000000130000000000507]  (* the same semiprime by continued fractions, Morrison and Brillhart *)
```

```mathematica
In[1]:= FactorInteger[3000000000000000015830000000000000002703, Method -> "QuadraticSieve"]  (* a 40-digit semiprime by the self-initialising quadratic sieve *)
```

```mathematica
In[1]:= FactorInteger[1049760769253203207210000000000000000012702105307963758807241, Method -> {"PollardPMinus1", "B1" -> 300}]  (* p - 1 = 2^4 3 5 7 11 13 17^2 19 ... 59 for the 21-digit factor, so a tiny B1 suffices *)
```

```mathematica
In[1]:= FactorInteger`WilliamsPPlus1[1046620906533602399990000000000000000012664112969056589039879]  (* here p + 1 is smooth instead, and p - 1 is not *)
```

```mathematica
In[1]:= FactorInteger`ECM[2^128 + 1]  (* the Fermat number F7, which Morrison and Brillhart first split with CFRAC *)
```

```mathematica
In[1]:= FactorInteger`BlakeRationalBaseDescent[742245361114170698080188160556986382238609541182118959142360194921632978579427957105431]  (* the larger factor is near 3 (22/7)^150; the base 22/7 is found by search *)
```

```mathematica
In[1]:= FactorInteger`TrialDivision[2^67 - 1, "Bound" -> 1000]  (* no prime below the bound divides it: the composite comes back whole, with FactorInteger::nofac *)
```

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
