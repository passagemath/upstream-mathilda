# Number Theory

Integer and number-theoretic functions: greatest common divisor and least common multiple (`GCD`, `LCM`, `ExtendedGCD`), modular arithmetic (`PowerMod`, `MultiplicativeOrder`, `PrimitiveRoot`, `ChineseRemainder`), primality and factorization (`PrimeQ`, `Prime`, `PrimePi`, `NextPrime`, `FactorInteger`, `SquareFreeQ`), arithmetic functions (`EulerPhi`), partitions (`IntegerPartitions`, `PartitionsP`, `PartitionsQ`), and continued fractions (`ContinuedFraction`, `FromContinuedFraction`).

## GCD, LCM

- `GCD[n1, n2, ...]`: Greatest common divisor.
- `LCM[n1, n2, ...]`: Least common multiple.

Both fold through GMP when any argument is a bigint, so results that exceed `int64` (e.g. `LCM[20!, 10^100 + 3]`) are computed exactly rather than left symbolic.

> **Packed arrays.** Over an `int64` buffer both run on the machine words —
> `GCD[list, n]`, `LCM[list, n]` and the two-array forms — and answer with a
> packed result. A product that would overflow `int64` abandons the buffer so
> the ordinary path answers exactly in GMP; never a wrapped value. See
> [Packed arrays](packed-arrays.md).

## Divisible

- `Divisible[n, m]`: yields `True` if `n` is divisible by `m`, and `False` otherwise. Attributes: `Listable`, `Protected`.

`n` is divisible by `m` when `n` is an integer multiple of `m`; this is effectively `Mod[n, m] == 0`. `Divisible[n, m]` returns `False` unless `n` and `m` are manifestly divisible.

**Features**:
- Machine integers and GMP bigints: tested directly with `mpz_divisible_p`, so large cases such as `Divisible[10^3000 + 1, 16001]` → `True` are exact. By the GMP convention, divisibility by `0` holds iff `n == 0` (`Divisible[0, 0]` → `True`, `Divisible[6, 0]` → `False`); sign is ignored (`Divisible[10, -2]` → `True`).
- Gaussian integers, rationals, and exact numeric quantities: the quotient `n/m` is formed and evaluated; the result is `True` iff it reduces to an integer or a Gaussian integer. So `Divisible[3 + I, 1 - I]` → `True`, `Divisible[3/2, 1/2]` → `True`, `Divisible[2 Pi, Pi/2]` → `True`, while `Divisible[Sqrt[6], Sqrt[2]]` → `False`.
- `Listable`: threads element-wise over lists, e.g. `Divisible[{1, 2, 3, 4, 5, 6}, 2]` → `{False, True, False, True, False, True}`.
- Symbolic, non-numeric arguments leave the call unevaluated (e.g. `Divisible[x, 2]`).
- Diagnostics: too few arguments emit `Divisible::argm`, too many emit `Divisible::argt`; both leave the call unevaluated.

```
In[1]:= Divisible[10, 2]
Out[1]= True

In[2]:= Divisible[5, 2]
Out[2]= False

In[3]:= Divisible[3 + I, 1 - I]
Out[3]= True

In[4]:= Divisible[2 Pi, Pi/2]
Out[4]= True

In[5]:= Divisible[Sqrt[6], Sqrt[2]]
Out[5]= False

In[6]:= Divisible[{1, 2, 3, 4, 5, 6}, 2]
Out[6]= {False, True, False, True, False, True}
```

## CoprimeQ

- `CoprimeQ[n1, n2, ...]`: yields `True` if the arguments are pairwise relatively prime, and `False` otherwise. Attributes: `Listable`, `Orderless`, `Protected`.

Integers are relatively prime when their greatest common divisor is `1`. `CoprimeQ` returns `False` unless the arguments are manifestly coprime integers or Gaussian integers.

**Features**:
- Machine integers and GMP bigints, handled uniformly through `mpz_gcd`, so large cases are exact: `CoprimeQ[2^100 - 1, 3^100 - 1]` → `False` (both even), `CoprimeQ[2^127 - 1, 2^61 - 1]` → `True`. Sign is ignored; `GCD(0, n) = |n|`, so `CoprimeQ[0, 1]` → `True` but `CoprimeQ[0, 5]` → `False`.
- More than two arguments are tested *pairwise*: `CoprimeQ[6, 35, 143]` → `True`, while `CoprimeQ[2, 3, 4]` → `False` (2 and 4 share a factor).
- Gaussian integers: with `GaussianIntegers -> True`, or when any argument is an exact Gaussian integer, coprimality is tested over `Z[i]` via the Gaussian Euclidean algorithm (round-to-nearest division). `CoprimeQ[5 + I, 1 - I]` → `False` (both divisible by `1 + I`); `CoprimeQ[2, 5, GaussianIntegers -> True]` → `True`, while `CoprimeQ[2, 10, GaussianIntegers -> True]` → `False`.
- `Orderless`: argument order is irrelevant, and the `GaussianIntegers` option may appear at any position.
- `Listable`: threads element-wise over lists, e.g. `CoprimeQ[{1, 2, 3, 4, 5}, 6]` → `{True, False, False, False, True}`.
- As a `*Q` predicate it always returns a Boolean: `CoprimeQ[]` → `False`, `CoprimeQ[n]` → `True` (no pairs), and anything not a manifestly coprime integer or Gaussian integer — rationals, reals, symbols, malformed options — yields `False` (e.g. `CoprimeQ[a, b]` → `False`).

```
In[1]:= CoprimeQ[8, 11]
Out[1]= True

In[2]:= CoprimeQ[2, 4]
Out[2]= False

In[3]:= CoprimeQ[2, 3, -5, 7]
Out[3]= True

In[4]:= CoprimeQ[5 + I, 1 - I]
Out[4]= False

In[5]:= CoprimeQ[{1, 2, 3, 4, 5}, 6]
Out[5]= {True, False, False, False, True}
```

## ExtendedGCD

- `ExtendedGCD[n1, n2, ...]`: returns `{g, {r1, r2, ...}}` where `g == GCD[n1, ...]` and `g == r1 n1 + r2 n2 + ...` (a multi-argument Bézout / extended GCD). Attributes: `Listable`, `Protected` (deliberately *not* `Flat`/`Orderless`/`OneIdentity`, since the cofactor list is positional).

**Features**:
- Integer-only; both machine integers and GMP bigints are accepted, and cofactors demote back to machine integers when they fit.
- Computed by folding GMP's `mpz_gcdext` pairwise: `gcd(running_g, n_i) = s·running_g + t·n_i`, scaling the accumulated cofactors by `s` and appending `t` each step. The running gcd stays non-negative, so `g` is exactly `GCD` and the sign convention is the conventional one.
- Threads element-wise over lists, e.g. `ExtendedGCD[3, {5, 15}]` → `{{1, {2, -1}}, {3, {1, 0}}}`.
- Edge cases: `ExtendedGCD[]` → `{0, {}}`; `ExtendedGCD[n]` → `{|n|, {±1}}`; zeros and negatives handled (`g` non-negative).
- Diagnostics: an inexact (Real) argument emits `ExtendedGCD::exact`; an exact non-integer (Rational) argument emits `ExtendedGCD::egcd`; symbolic arguments leave the call unevaluated silently.

## ChineseRemainder

- `ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}]`: the smallest `x ≥ 0` satisfying every integer congruence `x mod mi == ri mod mi`, lying in `0 ≤ x < LCM[m1, m2, ...]`. Attributes: `Protected` (deliberately *not* `Listable` — the residue and modulus lists are atomic arguments).
- `ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}, d]`: the smallest such `x ≥ d`, lying in `d ≤ x < d + LCM[m1, m2, ...]`.

**Features**:
- Integer-only; machine integers and GMP bigints share one path, so a residue-number-system recovery over large coprime moduli returns the full bignum (`ChineseRemainder[data, keys]` with prime keys is the CRT-encryption idiom).
- The moduli need **not** be pairwise coprime: a solution exists iff every pair of congruences agrees modulo the gcd of their moduli. When the system is inconsistent, `ChineseRemainder` is left unevaluated, e.g. `ChineseRemainder[{1, 2}, {6, 10}]` (since `1 ≢ 2 (mod gcd(6, 10) = 2)`).
- Computed by a streaming pairwise CRT fold: an accumulator congruence `x ≡ x_acc (mod m_acc)` is merged with each `(ri, |mi|)` via `mpz_gcdext` (`x_acc += ((ri − x_acc)/g)·u·m_acc`, `m_acc = lcm`), reducing modulo the running lcm each step. The final `x = d + ((x_acc − d) mod L)` gives the smallest `x ≥ d` and handles a negative `d` and empty lists (`L = 1`).
- If all `0 ≤ ri < mi` the result satisfies `x mod mi == ri`; unreduced or negative residues are handled (the result gives `x mod mi == ri mod mi`).
- Edge cases: `ChineseRemainder[{}, {}]` → `0`; a zero modulus, a length mismatch between the two lists, and non-integer or symbolic arguments leave the call unevaluated.
- Diagnostic: `ChineseRemainder::argt` when called with other than 2 or 3 arguments.

```
In[1]:= ChineseRemainder[{3, 4}, {4, 5}]
Out[1]= 19

In[2]:= ChineseRemainder[{2, 3, 5}, {3, 5, 7}]
Out[2]= 68

In[3]:= ChineseRemainder[{2, 3, 5}, {3, 5, 7}, 100]
Out[3]= 173

In[4]:= ChineseRemainder[{1, 2}, {6, 10}]
Out[4]= ChineseRemainder[{1, 2}, {6, 10}]

In[5]:= Mod[ChineseRemainder[{123, 456, 789}, {1000000007, 1000000009, 1000000021}], {1000000007, 1000000009, 1000000021}]
Out[5]= {123, 456, 789}
```

## PowerMod

Gives modular exponentiations, inverses, and roots.
- `PowerMod[a, b, m]`: Gives $a^b \pmod m$.
- `PowerMod[a, -1, m]`: Finds the modular inverse of `a` modulo `m`.
- `PowerMod[a, 1/r, m]`: Finds a modular `r`-th root of `a`.

**Features**:
- `Protected`, `Listable`.
- Evaluates much more efficiently than `Mod[a^b, m]`.
- Integer-exponent path uses GMP `mpz_powm` / `mpz_invert`; `a`, `b`, `m`
  may all be arbitrary-precision bignums.
- Rational-exponent path `PowerMod[a, p/q, m]` solves the modular `q`-th
  root of `a^p` via Tonelli–Shanks, the coprime closed form, Hensel
  lifting, and CRT over the factorisation of `m`. The
  Adleman–Manders–Miller branch (`gcd(q, p-1) > 1` for `q > 2` over a
  large prime) is not yet supported.
- Returns unevaluated if the corresponding inverse or root does not exist
  (or cannot be computed by the implemented algorithms).
- Allows threading over lists natively.

```mathematica
In[1]:= PowerMod[2, 10, 3]
Out[1]= 1

In[2]:= PowerMod[3, -2, 7]
Out[2]= 4

In[3]:= PowerMod[3, 1/2, 2]
Out[3]= 1

In[4]:= PowerMod[2, {10, 11, 12, 13, 14}, 5]
Out[4]= {4, 3, 1, 2, 4}

In[5]:= PowerMod[100, 1/2, 17 * 19 * 23]
Out[5]= 10

In[6]:= PowerMod[2, 1/2, 10^18 + 9]
Out[6]= 742174169206529574
```

> **Packed arrays.** `PowerMod[list, e, m]` over an `int64` buffer is
> square-and-multiply on the machine words, with no GMP allocation per
> element. Taken for a non-negative `e` and a modulus small enough that the
> squaring step cannot overflow; a negative exponent is a modular inverse,
> which exists only for some elements, and takes the ordinary path.

## MultiplicativeOrder

Gives the multiplicative order of `k` modulo `n` — the smallest positive
integer `m` such that $k^m \equiv 1 \pmod{n}$.
- `MultiplicativeOrder[k, n]`: order of `k` modulo `n`.
- `MultiplicativeOrder[k, n, {r1, r2, ...}]`: smallest positive `m` with
  $k^m \equiv r_i \pmod{n}$ for some `i` (a multi-target discrete log).

The order is the size of the cyclic subgroup generated by `k` in the units
modulo `n`. Reversing the map — given `k` and a target `b`, find `x` with
$k^x \equiv b \pmod n$ — is the **discrete logarithm problem**, and the
three-argument form solves exactly that. Computing the order is easy (it divides
`EulerPhi[n]`); the discrete logarithm it inverts is believed hard, and that
presumed hardness underlies Diffie–Hellman key exchange and the ElGamal and DSA
schemes. See the *Handbook of Applied Cryptography* (Menezes, van Oorschot,
Vanstone, 1996), §3.6, for the standard algorithms (baby-step giant-step,
Pohlig–Hellman, index calculus).

**Features**:
- `Protected`.
- All arithmetic uses GMP `mpz_t`, so `k`, `n`, and any `r_i` may be
  arbitrary-precision bignums.
- Negative `n` is treated as `|n|`; `k` (and each `r_i`) is reduced modulo
  `n` before the search, so negative or out-of-range inputs work
  transparently.
- Returns unevaluated when `gcd(k, n) != 1`, when `n` is zero, or when
  no power of `k` lands in the residue set (3-arg form).
- The 2-arg form factors `phi(n)` and successively strips prime factors
  from `phi(n)` whose corresponding exponent still maps `k` to 1 — so
  the work is dominated by factoring `phi(n)`, not by walking the orbit.
- The 3-arg form walks `k^m mod n` for `m = 1, ..., order(k, n)`. To
  guard against pathological group sizes, the call returns unevaluated
  if the order exceeds `10^8` or does not fit in an unsigned long.
- Non-integer numeric inputs (`Real`, `Complex`, `Rational`) and symbolic
  arguments flow through unevaluated with no diagnostic.
- Diagnostic: `MultiplicativeOrder::argt` when called with anything other
  than 2 or 3 arguments.

```mathematica
In[1]:= MultiplicativeOrder[5, 8]
Out[1]= 2

In[2]:= MultiplicativeOrder[5, 7]
Out[2]= 6

In[3]:= MultiplicativeOrder[-5, 7]
Out[3]= 3

In[4]:= MultiplicativeOrder[5, 7, {3, 11}]
Out[4]= 2

In[5]:= MultiplicativeOrder[10^10000, 7919]
Out[5]= 3959

In[6]:= Select[Range[43], MultiplicativeOrder[#, 43] == EulerPhi[43] &]
Out[6]= {3, 5, 12, 18, 19, 20, 26, 28, 29, 30, 33, 34}

In[7]:= MultiplicativeOrder[10, 22]
Out[7]= MultiplicativeOrder[10, 22]
```

## JacobiSymbol

Gives the Jacobi symbol $(n/m)$.
- `JacobiSymbol[n, m]`: the Jacobi symbol of `n` with respect to `m`,
  returned as `-1`, `0`, or `1`.

**Features**:
- `Protected`, `Listable` — threads element-wise over lists and arrays.
- For prime `m` the Jacobi symbol reduces to the Legendre symbol, equal to
  `±1` according to whether `n` is a quadratic residue modulo `m`, and `0`
  when `m` divides `n`.
- Following the Wolfram Language, this is the full Kronecker-symbol
  generalisation: the second argument `m` may be even or non-positive, and
  the first argument `n` may be negative.
- Computed with GMP `mpz_kronecker` in $O((\log m)^2)$ time, so `n` and `m`
  may be arbitrary-precision bignums.
- Non-integer numeric inputs and symbolic arguments flow through
  unevaluated with no diagnostic.
- Diagnostic: `JacobiSymbol::argrx` when called with other than 2 arguments.

```mathematica
In[1]:= JacobiSymbol[10, 5]
Out[1]= 0

In[2]:= Table[JacobiSymbol[n, m], {n, 0, 10}, {m, 1, n, 2}]
Out[2]= {{}, {1}, {1}, {1, 0}, {1, 1}, {1, -1, 0}, {1, 0, 1},
         {1, 1, -1, 0}, {1, -1, -1, 1}, {1, 0, 1, 1, 0}, {1, 1, 0, -1, 1}}

In[3]:= JacobiSymbol[10^10 + 1, Prime[1000]]
Out[3]= 1

In[4]:= JacobiSymbol[7, 6]
Out[4]= 1

In[5]:= JacobiSymbol[{2, 3, 5, 7, 11}, 3]
Out[5]= {-1, 0, -1, 1, -1}

In[6]:= JacobiSymbol[-3, {1, 3, 5, 7}]
Out[6]= {1, 0, -1, 1}
```

## DirichletCharacter

- `DirichletCharacter[k, j, n]` — the Dirichlet character χ_{k,j}(n), modulus
  k, index j = 1..EulerPhi[k] (j = 1 principal).

Implemented in `src/special_functions/dirichlet_character.c` on the shared
engine `src/special_functions/dirichlet.{c,h}`, registered via
`dirichletcharacter_init()`. Attributes `Protected` (as in Mathematica); it
threads over lists in every argument itself.

- **Labelling (Mathematica's, verified on 148,138 values for k ≤ 100 and large
  prime powers / 2-adic / composite moduli).** Factor k with primes ascending;
  (Z/kZ)* splits into cyclic components — one per odd p^a, generated by the
  smallest primitive root mod p^a; one for 4 (generator 3); two for 2^a, a ≥ 3:
  ⟨−1⟩ then ⟨5⟩. j − 1 is the mixed-radix number over the components, the first
  component most significant, and χ(g_c) = exp(2πi d_c/ord_c). So
  `DirichletCharacter[k, j, PrimitiveRoot[k]] = Exp[2 Pi I (j-1)/EulerPhi[k]]`
  for prime k, and `DirichletCharacter[45, 8, n] =
  DirichletCharacter[9, 2, n] DirichletCharacter[5, 4, n]`.
- **Values** are exact: 0 when GCD[n, k] > 1, otherwise a product of the
  component values in Mathematica's normal form — the components that are units
  (±1, ±I) multiply into a unit, the rest combine into one E^(I π r),
  r ∈ (−1, 1), folded back into the unit when it lands on one:
  `DirichletCharacter[45, 8, 7] = -I E^(-2 I Pi/3)`,
  `DirichletCharacter[32, 10, 3] = -E^(3 I Pi/4)`,
  `DirichletCharacter[20!, 300, 23] = E^(I Pi/9)`.
- **Engine.** Per (k, j) a table (cached, immutable once published, safe from
  threaded Compile loops): components, digits, order, parity, conductor and the
  primitive χ*, the Euler-factor primes, and — on demand, for f ≤ 2^22 — dense
  χ* tables (exact exponents, machine values, root number) built by walking the
  component generators. Single values use Pohlig–Hellman discrete logs with
  baby-step giant-step in each prime-order subgroup (k < 2^63; a component
  order whose largest prime factor exceeds 2^44 leaves the call unevaluated).
- **Messages.** `argrx`, `intp` (position 1 or 2), `invl` (j > EulerPhi[k]),
  `int` (a non-integer n).
- **Surfaces.** An int64 packed / `NDArray` n: a real character returns an
  int64 buffer of 0, ±1; a complex one a List of exact roots of unity (memoised
  per residue). Exact-valued, so no `Compile[]` lowering (compile-coverage
  `EXEMPT`).

```mathematica
In[1]:= Table[DirichletCharacter[5, 2, n], {n, 1, 5}]
Out[1]= {1, I, -I, -1, 0}

In[2]:= Table[DirichletCharacter[8, j, n], {j, 1, 4}, {n, 1, 8}]
Out[2]= {{1, 0, 1, 0, 1, 0, 1, 0}, {1, 0, -1, 0, -1, 0, 1, 0}, {1, 0, -1, 0, 1, 0, -1, 0}, {1, 0, 1, 0, -1, 0, -1, 0}}

In[3]:= DirichletCharacter[45, 8, 7]
Out[3]= -I E^((-2/3*I) Pi)
```

## PrimitiveRoot

Gives a primitive root of `n`, i.e. a generator of the multiplicative group
of integers modulo `n` that are coprime to `n`.
- `PrimitiveRoot[n]`: a primitive root of `n`.
- `PrimitiveRoot[n, k]`: the smallest primitive root of `n` greater than or
  equal to `k`.

**Features**:
- `Protected`, `Listable`.
- Returns unevaluated unless `n` is 2, 4, an odd prime power $p^k$, or
  twice an odd prime power $2 p^k$ (the moduli for which $(\mathbb{Z}/n\mathbb{Z})^*$
  is cyclic). For all other `n`, the call is left unevaluated.
- The 1-argument form returns a canonical primitive root: smallest for
  $n \in \{2, 4\}$ and odd prime powers; for $n = 2 p^k$ the formula
  $g$ if $g$ is odd else $g + p^k$ is applied, where $g$ is the smallest
  primitive root of $p^k$. This matches Mathematica's convention so that,
  e.g. `PrimitiveRoot[10] == 7` while `PrimitiveRoot[10, 1] == 3`.
- The 2-argument form walks forward from `k`; if `k > n - 1` the call is
  left unevaluated.
- All arithmetic uses GMP `mpz_t`, so machine integers, bignums, and
  symbolic bignum products like `Prime[1000000]^1000000` are handled
  uniformly. The prime-power detection iteratively strips prime exponents
  via `mpz_root`, which runs in $O(\omega(k))$ root extractions.
- Diagnostics:
  - `PrimitiveRoot::argt` if not called with 1 or 2 arguments.
  - `PrimitiveRoot::intg` if `n` (or the 2nd-arg `k` when numeric) is not
    an integer greater than 1.

```mathematica
In[1]:= PrimitiveRoot[9]
Out[1]= 2

In[2]:= PrimitiveRoot[10]
Out[2]= 7

In[3]:= PrimitiveRoot[10, 1]
Out[3]= 3

In[4]:= PrimitiveRoot[10, 4]
Out[4]= 7

In[5]:= PrimitiveRoot[{9, 7, 19}]
Out[5]= {2, 3, 2}

In[6]:= PrimitiveRoot[12]
Out[6]= PrimitiveRoot[12]
```

## PrimitiveRootList

Gives the sorted list of all primitive roots of `n` in the canonical
residues $\{1, \ldots, n-1\}$.
- `PrimitiveRootList[n]`

**Features**:
- `Protected`, `Listable`.
- Returns `{}` if `n` is not 2, 4, an odd prime power, or twice an odd
  prime power.
- Enumerates the $\varphi(\varphi(n))$ primitive roots as $g^i \bmod n$
  for $i \in \{1, \ldots, \varphi(n)\}$ with $\gcd(i, \varphi(n)) = 1$,
  where $g$ is the smallest primitive root of `n`. The list is sorted
  ascending.
- Falls back to unevaluated when $\varphi(n)$ exceeds `unsigned long`,
  since the enumeration cannot be represented.
- Non-integer numeric inputs (e.g. `11.0`, `11 + I`) flow through
  unevaluated with no diagnostic, matching Mathematica.
- Diagnostic: `PrimitiveRootList::argx` if not called with exactly 1
  argument.

```mathematica
In[1]:= PrimitiveRootList[9]
Out[1]= {2, 5}

In[2]:= PrimitiveRootList[19]
Out[2]= {2, 3, 10, 13, 14, 15}

In[3]:= PrimitiveRootList[12]
Out[3]= {}

In[4]:= Union[Table[PowerMod[2, i, 9], {i, 6}]]
Out[4]= {1, 2, 4, 5, 7, 8}
```

## PrimeQ

- `PrimeQ[n]`: Returns `True` if `n` is a prime number, `False` otherwise.
- `PrimeQ[z]`: For a Gaussian integer `z = a + b I`, returns `True` if `z` is a Gaussian prime.
- `PrimeQ[n, GaussianIntegers -> True]`: Tests primality of `n` in the
  Gaussian integers `Z[i]`. A rational integer `n` is a Gaussian prime
  iff `|n|` is prime in `Z` AND `n ≡ 3 (mod 4)` (primes `≡ 1 (mod 4)`
  split as `(a+bi)(a-bi)`, and `2 = -i(1+i)^2` is reducible).

**Features**:
- `Listable`, `Protected`.
- Always returns `True` or `False`. For non-integer / non-Gaussian
  inputs (symbols, `Sqrt[2]`, `Exp[2 Pi I/3]`, strings, etc.) returns
  `False` — `*Q` predicates never remain symbolic.
- A Gaussian integer `a + b I` is a Gaussian prime if:
  - Both `a` and `b` are nonzero and `a^2 + b^2` is an ordinary prime, or
  - One of `a`, `b` is zero and the absolute value of the other is a prime congruent to 3 mod 4.

```mathematica
In[1]:= PrimeQ[7]
Out[1]= True

In[2]:= PrimeQ[1 + I]
Out[2]= True

In[3]:= PrimeQ[1 + 2 I]
Out[3]= True

In[4]:= PrimeQ[3 I]
Out[4]= True

In[5]:= PrimeQ[5 I]
Out[5]= False

In[6]:= PrimeQ[2 + 2 I]
Out[6]= False

In[7]:= PrimeQ[5, GaussianIntegers -> True]
Out[7]= False

In[8]:= PrimeQ[3, GaussianIntegers -> True]
Out[8]= True

In[9]:= PrimeQ[Exp[2 Pi I/3]]
Out[9]= False
```

## Prime

- `Prime[n]`: the `n`-th prime $p_n$ (`Prime[1] = 2`, `Prime[2] = 3`, ...).

**Features**:
- `Listable`, `Protected`.
- Small `n` is read directly from a sieve table of the primes below $10^6$.
- Large `n` inverts `PrimePi`: a Cipolla asymptotic estimate for $p_n$ is refined
  by a Newton step against the exact prime counter, then `NextPrime`/`PrevPrime`
  steps land exactly on $p_n$.  Exact for `n` up to about $1.4 \times 10^{12}$
  ($p_n \le 5 \times 10^{13}$); beyond that the call is left unevaluated.
- A non-positive-integer argument emits `Prime::intpp`; a wrong argument count
  emits `Prime::argx`; both leave the call unevaluated.

```mathematica
In[1]:= Prime[100]
Out[1]= 541

In[2]:= Prime[{1, 3, 4, 10}]
Out[2]= {2, 5, 7, 29}

In[3]:= Prime[10^10]
Out[3]= 252097800623
```

> **Packed arrays.** `Prime[list]` over an `int64` buffer is one **sieve** to
> the largest requested index followed by a gather, rather than a separate
> prime count per element — so the whole array costs about what its largest
> single index used to.

## PrimePi

- `PrimePi[x]`: the number of primes less than or equal to `x`.
- `PrimePi[x, Method -> m]`: count with a specific algorithm.

**Features**:
- `Listable`, `Protected`.  Default `Method -> Automatic`.
- Several recognised prime-counting algorithms are available via `Method`:
  `"Sieve"` (segmented sieve of Eratosthenes), `"Legendre"`, `"Meissel"`,
  `"Lehmer"`, `"LMO"` (Lagarias-Miller-Odlyzko, segmented special leaves),
  `"DelegliseRivat"`, and `"LucyHedgehog"`.  `"DelegliseRivat"` refines `"LMO"`
  by classifying special leaves as trivial / easy / hard and answering the
  trivial and easy ones with $O(1)$ prime-count lookups, sending only the hard
  leaves through the segmented sieve; it is the fastest combinatorial method here
  (≈ 1.6 s at $10^{12}$, 9.9 s at $10^{13}$).  `Automatic` selects an efficient
  method for the magnitude of `x` (a direct sieve table below $10^6$,
  `"LucyHedgehog"` up to $10^9$, then `"DelegliseRivat"`).  The combinatorial
  methods exact-count up to $5 \times 10^{13}$; `"Sieve"`/`"Legendre"` are
  intended for smaller `x`.
- An unrecognised `Method` setting emits `PrimePi::method` and leaves the call
  unevaluated.

```mathematica
In[1]:= PrimePi[10^9]
Out[1]= 50847534

In[2]:= PrimePi[10^9, Method -> "LMO"]
Out[2]= 50847534

In[3]:= PrimePi[{10, 100}]
Out[3]= {4, 25}
```

## NextPrime

- `NextPrime[x]`: Smallest prime above `x`.
- `NextPrime[x, k]`: $k$-th next prime above `x` (or previous if $k$ is negative).

**Features**:
- `Protected`.
- Supports negative $k$ for finding previous primes.
- Remains unevaluated if no such prime exists (e.g., `NextPrime[2, -1]`).

## FactorInteger

- `FactorInteger[n]`: the prime factorisation `{{p1, e1}, {p2, e2}, ...}`, primes ascending.
- `FactorInteger[n, k]`: partial factorisation with **at most `k` entries**; the last entry may be an
  unfactored remainder (Mathematica's contract: `FactorInteger[100, 1]` is `{{100, 1}}`,
  `FactorInteger[100, 2]` is `{{2, 2}, {25, 1}}`). `k < 1` leaves the call unevaluated.
- `FactorInteger[n, Automatic]`: the same as `FactorInteger[n]`.
- `FactorInteger[n, GaussianIntegers -> True]`: factorisation over Z[i]. A Gaussian integer or Gaussian
  rational argument (`3 + 4 I`, `(3 + 4 I)/(1 + I)`) is factored over Z[i] automatically. The unit comes
  first, when it is not 1. Primes are given in the first quadrant (Re > 0, Im >= 0), sorted by real part
  and then imaginary part. Every output matches Mathematica 14.
- `FactorInteger[n, Method -> m]`: the algorithm is chosen by `m`, which is a method name, `Automatic`,
  or `{name, "Sub" -> v, ...}`. Every method is also a builtin of its own,
  `FactorInteger`name[n, "Sub" -> v, ...]`, taking the same `k` and `GaussianIntegers` arguments.
- Rationals: the denominator's primes get negative exponents. Negative `n` starts with `{-1, 1}`;
  `0` gives `{{0, 1}}`, `1` gives `{{1, 1}}`.
- Messages, all routed through the `Quiet`/`Check` funnel:
  - `FactorInteger::nofac`: a composite was left unfactored, once per such composite;
  - `FactorInteger::method`: unknown method;
  - `FactorInteger::mopt`: bad suboption;
  - `FactorInteger::optx`: unknown option;
  - `FactorInteger::exact`: inexact or symbolic argument, as in Mathematica.

The implementation is in `src/numbertheory/factorinteger/` (one file per algorithm; see `fi_internal.h`).
Every randomised choice (rho constants, ECM sigmas, p+1 seeds, QS polynomials) is drawn from a generator
seeded by `n`, so the answer, the time taken and any `nofac` depend on `n` alone.

| Method (`FactorInteger`name`) | Algorithm | Suboptions |
|---|---|---|
| `Automatic` | Below 2^64 the word engine handles everything: Montgomery arithmetic, a deterministic Miller–Rabin test, micro-ECM with Suyama curves and a baby-step/giant-step stage 2, and rho/SQUFOF as backstops. Above 2^64 the tiers are: trial division to 2^16 → Pollard rho (two-limb Montgomery below 2^126) → two-limb micro-ECM (below 2^126) → a Fermat probe → one p-1 → one p+1 (three seeds from 70 digits) → ECM levels 15, 20, ... digits, switching to the quadratic sieve once the ECM effort reaches 2/7 of the composite's digit count, for composites up to 100 digits. After a split, both parts restart at the tier that succeeded. | — |
| `"TrialDivision"` | Primes are tested in word-sized products, one bignum remainder per 2–4 primes; a segmented sieve takes over above 2^20. The default bound is sqrt(n) when that is below 10^7, else 10^7. | `"Bound"` (≤ 2^40) |
| `"PollardRho"` | Brent's variant with batched gcds and backtracking; Montgomery arithmetic on one and two words, GMP above that. | `"Iterations"`, `"Starts"` |
| `"Fermat"` | Fermat's method with quadratic-residue filters modulo nine small moduli. | `"Iterations"` |
| `"Lehman"` | Lehman's O(n^(1/3)) method; complete below 2^42, budgeted above. | `"Iterations"` |
| `"ShanksSquareForms"` (alias `"SQUFOF"`) | SQUFOF over the 16 Gower–Wagstaff multipliers, in word arithmetic for kN < 2^124. | `"Iterations"` |
| `"CFRAC"` | Morrison–Brillhart: Knuth–Schroeppel multiplier, early-abort trial division, single large primes, GF(2) elimination. Relation collection is pipelined across threads: one thread extends the continued fraction into the next batch while every thread trial-divides the current one, and the hits are replayed in term order. The relations, restarts and divisor are therefore those of the serial method for any thread count. | `"Iterations"` |
| `"Dixon"` | Random squares drawn near sqrt(kN), so the residues are square-root sized, with the same relation machinery as CFRAC. The draws are taken in batches; the candidates are built and trial-divided in parallel, and the hits are replayed in draw order. Results and the random stream are those of the serial method. | `"Iterations"` |
| `"QuadraticSieve"` (aliases `"QS"`, `"SIQS"`, `"MPQS"`; the heads `FactorInteger`QS` and `FactorInteger`SIQS` exist, `"MPQS"` is a `Method` string only) | Native self-initialising QS: Gray-code polynomial switching, 32 KB block sieve with bucket sieving for primes ≥ 2^15, batched candidate trial division. Single large primes, plus double large primes from 72 digits, combined by cycle-finding in the large-prime graph. Polynomial families are sieved in parallel from 45 digits (one A per thread). Linear algebra is block Lanczos (Montgomery 1995) above 1000 columns, dense elimination below. | — |
| `"PollardPMinus1"` (alias `"PollardP-1"`) | GMP-ECM p-1, with a B1 back-off when every prime is caught at once; a word walk with prime-by-prime backtracking below 2^64. | `"B1"` (10^6), `"B2"` |
| `"WilliamsPPlus1"` (alias `"WilliamsP+1"`) | GMP-ECM p+1 with seeds 2/7, 6/5, 3, and the same back-off and word walk. The seeds run in parallel, one per thread, and the first success in seed order is taken. | `"B1"` (10^6), `"B2"` |
| `"ECM"` | GMP-ECM, climbing its optimal levels (B1 = 2000, 11000, 5e4, ...) with the expected curve count at each. Curves run in parallel, one per thread, for n ≥ 2^160. They use the batch parametrisation with sigmas drawn deterministically, and successes are taken in curve order, so the result does not depend on scheduling. A user `"Sigma"` is a Suyama parameter. There is a B1 back-off on gcd = n. | `"B1"`, `"B2"`, `"Curves"`, `"Sigma"` |
| `"BlakeRationalBaseDescent"` | Blake's rational base descent, for a factor close to c (a/b)^k: six variants (exact descent, the reference descent, all-hits, lattice, two-pass, Coppersmith), threaded base search. See [below](#blakerationalbasedescent). | `"Base"`, `"Variant"`, `"Order"`, `"Height"`, `"Bases"`, `"MaxCoefficient"`, `"Multiplier"`, `"Tolerance"`, `"Beta"`, `"Multiplicity"`, `"Dimension"`, `"BaseRange"`, `"Coefficient"`, `"Exponent"` |

Fallbacks: builds without GMP-ECM use FLINT's ECM and p+1 (and a native p-1); builds without either use
native mpz ECM, p-1 and p+1. The quadratic sieve and every other method are native, with no external
dependency.

Threads: the quadratic sieve, ECM, p+1, CFRAC, Dixon and BlakeRationalBaseDescent use every online CPU (up to 16) in builds with
`MATHILDA_THREADS`. `TimeConstrained` budgets are process CPU time, so they are consumed correspondingly
faster while they run. An abort is taken on the main thread when the current region ends: one polynomial
family per QS thread, one curve per ECM thread, one seed per p+1 thread, or a few milliseconds of
CFRAC/Dixon candidates.

### BlakeRationalBaseDescent

`FactorInteger`BlakeRationalBaseDescent[n]` (the same as `FactorInteger[n, Method -> "BlakeRationalBaseDescent"]`)
is S. Blake's *rational base descent* (RBD), from two papers: "Rational Base Descent: A Deterministic
Algorithm for Factoring Structured Semiprimes" and "Lattice-Augmented Rational Base Descent". It is a
special-purpose method. It succeeds when the **larger** factor P of `n = P M` is close to a power of a
rational number, `P ≈ c (a/b)^k`, and is no threat to unstructured (e.g. RSA) moduli.

**Mechanism.** Starting from `Q = n`, repeat `Q <- floor(b Q / a)`. Each step divides by the base. After
`k` steps `Q` has divided `P` down to about `c`, so `Q` lies within a small constant of the integer `c M`,
and `gcd(Q - j, n)` for a small `j` exposes `M` (or `M gcd(c, P)`). The paper's theorem: if
`(a/b)^k > M` and `|P - floor(c (a/b)^k)| < (a/b)^k / M`, a factor is found at step `k`. One descent tests
every exponent and every coefficient at once, in time polynomial in `log n` for a fixed base.

The base is supplied with `"Base"` or found by enumeration. The bases that succeed form a comb of isolated
teeth of relative width ~1/(k M) around the representations of P (plus harmonic teeth at `a/d`, `d | a`).
No signal points towards a tooth from a distance, so bases are enumerated breadth-first over the
Stern–Brocot tree (every rational once, coarse to fine), skipping perfect powers `(x/y)^m`. For an integer
base the descent truncates the base-`a` digits of `n`, and truncations inherit algebraic divisors. This
*truncation divisibility* is what recovers the cyclotomic factors of Cunningham numbers and repunits.

| Variant (`"Variant"` -> name; also the head `FactorInteger`BlakeRationalBaseDescent<name>`) | Algorithm |
|---|---|
| `"Exact"` (default) | The scaled exact descent. The recurrence runs on `S n`, with `S = 2^s ≥ 2a/(a-b)`, so the target always lies in a base-independent window of four consecutive values. Nothing is tested above `sqrt(cmax n)`, since no target `c M` can lie there; that depth is reached by one exact division. The window products are folded into an accumulator mod `n`, with one gcd per block of steps; a block is replayed candidate by candidate only when its gcd is non-trivial. |
| `"Classic"` | The paper's reference algorithm R1–R4. It tests `gcd(Q - j, n)` for `|j| ≤ ceil(2 + a/(a-b))` at every depth, searching before the termination test. |
| `"AllHits"` | The reference descent run to the end. Every distinct non-trivial gcd is recorded, and the record is refined into pairwise-coprime divisors of `n`; later calls on the cofactors split from that cache. Composite pieces are re-run with AllHits, which recovers nested algebraic divisors. |
| `"Lattice"` | Lattice-augmented descent. Howgrave-Graham's approximate-divisor detector (LLL on `N^(t-i) (Q+x)^i`, `x^j (Q+x)^t`) replaces the gcd window. It finds `c M` from `Q` within `X ≤ n^(β²)` when `M ≥ n^β`, widening the per-step tolerance from O(a/(a-b)) to `X`. The detector runs only in the depth band `X ≤ Q ≤ sqrt(cmax n) + X`, in parallel over the depths. Integer roots are found p-adically on the square-free part of each short polynomial. |
| `"TwoPass"` | The paper's practical algorithm. Pass 1 is Exact over the bases. Pass 2 is Lattice over a coarse grid of `"BaseRange" -> {lo, hi}`, with relative spacing `X/(2 c k n^β)` and each grid point realised by the simplest rational near it; without a range, pass 2 uses the Stern–Brocot bases of height 16 (or `"Height"`). |
| `"Coppersmith"` | A known base, coefficient and exponent: one detector call at `floor(c (a/b)^k)`. This is factoring with high bits known. It needs `"Base"`, `"Coefficient"` and `"Exponent"` (else `FactorInteger::rbdcop`). |

Suboptions (defaults in parentheses):
- `"Base"` -> `a/b` or a list of bases. A base is reduced, and inverted if below 1; `1` is rejected with
  `FactorInteger::mopt`. (Default: enumerate.)
- `"Order"` -> `"SternBrocot"` (default) or `"Sum"` (the papers' `a + b` order).
- `"Height"` -> `H`: largest numerator enumerated (1000). `"Bases"` -> `B`: enumeration budget (200000).
- `"MaxCoefficient"` -> `cmax` (2^64): sets the deferred-testing ceiling `sqrt(cmax n)`; `Infinity`
  tests every depth.
- `"Multiplier"` -> `m`: descend from `m n`, with gcds still taken with `n` (the paper's Open Question 2).
- `"Tolerance"` -> `X`, `"Beta"` -> `β` (1/2), `"Dimension"` -> `d` (12), `"Multiplicity"` -> `t`: the
  lattice detector. Given a tolerance, the least dimension whose LLL guarantee covers it is chosen.
  Otherwise `X` is the largest tolerance the guarantee covers at dimension `d`, less two bits.
- `"BaseRange"` -> `{lo, hi}` (TwoPass grid), `"Coefficient"` -> `c`, `"Exponent"` -> `k` (Coppersmith).

Behaviour:
- The suboptions describe the input. A composite cofactor that the configured search cannot split gets a
  light default search (Exact, Stern–Brocot height 64, at most 2000 bases); AllHits re-runs AllHits
  instead. Anything still composite is returned with `FactorInteger::nofac`.
- Threads: a many-base search runs blocks of 2048 bases on every thread. A lattice base runs its band two
  depths per thread per region. The lowest-indexed success wins (a hit caps the work claim), so the answer
  does not depend on the thread count, and a `TimeConstrained` abort is taken within milliseconds.
- Lattice reduction is FLINT's `fmpz_lll` when available, else a native fraction-free integral LLL
  (de Weger / Cohen 2.6.7, δ = 0.99): exact, but slower.
- The core (`fi_rbd.h`; `fi_rbd_descent.c`, `fi_rbd_bases.c`, `fi_rbd_lattice.c`, `fi_rbd_lll.c`) uses only
  GMP (and optionally FLINT), with no global state; the Mathilda glue is `fi_rbd.c`.

Performance, measured on 500 Stern–Brocot bases against a random semiprime, single thread:

| n (bits) | Exact (µs per base) | the v0.454 implementation | speedup |
|---|---|---|---|
| 256 | 61 | 2240 | 37× |
| 512 | 151 | 11400 | 75× |
| 1024 | 546 | 54800 | 100× |

The paper's Example 1 (1140 bits, base 22/7) takes 2.8 ms with Exact and 27 ms with Classic. The paper's
2048-bit challenge modulus is factored from scratch, with no base given, in 106 s on 16 threads.

```mathematica
In[1]:= FactorInteger`BlakeRationalBaseDescent[1001, "Base" -> 22/7]
Out[1]= {{7, 1}, {11, 1}, {13, 1}}

(* P = NextPrime[Floor[3 (22/7)^150]] times a 40-bit prime: the base is found by search *)
In[2]:= FactorInteger`BlakeRationalBaseDescent[742245361114170698080188160556986382238609541182118959142360194921632978579427957105431]
Out[2]= {{623347347967, 1}, {1190741187132579502456731851112421085210615890827643793135496647841845583593, 1}}

(* truncation divisibility: the algebraic factors of 2^1041 + 1 *)
In[3]:= Select[Quiet[FactorInteger`BlakeRationalBaseDescentAllHits[2^1041 + 1, "Base" -> 2]][[All, 1]], # < 10^120 &]
Out[3]= {3, 2083, 95562442332919646317117537304253622533190207882011713489066201641121786503686867002917439712921903606443}

(* a planted miss of 2^28 at depth 59: invisible to the gcd window, found by the lattice detector *)
In[4]:= First /@ FactorInteger`BlakeRationalBaseDescentLattice[195539142822724777401163665279552940381006352574199957151137, "Base" -> 22/7]
Out[4]= {296454907933147663952202071141, 659591518271709676788151702157}

In[5]:= First /@ FactorInteger`BlakeRationalBaseDescentCoppersmith[195539142822724777401163665279552940381006352574199957151137, "Base" -> 22/7, "Coefficient" -> 3, "Exponent" -> 59]
Out[5]= {296454907933147663952202071141, 659591518271709676788151702157}
```

**Features**:
- `Listable`, `Protected`; `Options[FactorInteger]` is `{GaussianIntegers -> False, Method -> Automatic}`.
- Perfect powers are reduced to a primitive base before any search, so `(2^61 - 1)^2` factors correctly.
- `EulerPhi` and `facint_factor_complete` decline (stay unevaluated or report incompleteness) rather than
  build an answer on a composite that survived the search.

```mathematica
In[1]:= FactorInteger[-12]
Out[1]= {{-1, 1}, {2, 2}, {3, 1}}

In[2]:= FactorInteger[3/4]
Out[2]= {{2, -2}, {3, 1}}

In[3]:= FactorInteger[100, 2]
Out[3]= {{2, 2}, {25, 1}}

In[4]:= FactorInteger[5, GaussianIntegers -> True]
Out[4]= {{-I, 1}, {1 + 2 I, 1}, {2 + I, 1}}

In[5]:= FactorInteger`ECM[2^128 + 1]
Out[5]= {{59649589127497217, 1}, {5704689200685129054721, 1}}

In[6]:= FactorInteger`QuadraticSieve[NextPrime[10^29] NextPrime[3 10^29]]
Out[6]= {{100000000000000000000000000319, 1}, {300000000000000000000000000007, 1}}
```

## Divisors

- `Divisors[n]`: the ascending list of positive integers that divide `n`.
- `Divisors[n, GaussianIntegers -> True]`: includes Gaussian-integer divisors.

**Options**:
- `GaussianIntegers -> False` (default): ordinary integer divisors. Set to
  `True` to compute divisors in the Gaussian integers `Z[i]`. A non-real
  Gaussian-integer input (e.g. `6 + 4 I`) auto-enables Gaussian mode.

**Features**:
- `Listable`, `Protected`.
- Machine integers and GMP bigints are handled uniformly; the result promotes to
  a big-integer list when needed.
- The sign of `n` is ignored (`Divisors[-12] == Divisors[12]`).
- Divisors are computed from the prime factorization (the divisor lattice), so
  cost scales with the number of divisors rather than `Sqrt[n]`.
- In Gaussian mode each divisor is the canonical first-quadrant representative of
  its associate class (`Re > 0`, `Im >= 0`), and the list is sorted by
  `(Re, Im)`. Rational primes are lifted to `Z[i]`: `2` ramifies as `1 + I`,
  primes `p ≡ 1 (mod 4)` split via a sum-of-two-squares (Cornacchia)
  decomposition, and primes `p ≡ 3 (mod 4)` stay inert.
- `Divisors[0]`, non-integer arguments, and calls whose divisor count would
  overflow (e.g. `Divisors[100!]`) are left unevaluated; `Divisors[]` issues a
  `Divisors::argx` message.

```mathematica
In[1]:= Divisors[1729]
Out[1]= {1, 7, 13, 19, 91, 133, 247, 1729}

In[2]:= Divisors[6]
Out[2]= {1, 2, 3, 6}

In[3]:= Divisors[{605, 871, 824}]
Out[3]= {{1, 5, 11, 55, 121, 605}, {1, 13, 67, 871}, {1, 2, 4, 8, 103, 206, 412, 824}}

In[4]:= Divisors[6 + 4 I]
Out[4]= {1, 1 + I, 1 + 5 I, 2, 3 + 2 I, 6 + 4 I}

In[5]:= Divisors[2, GaussianIntegers -> True]
Out[5]= {1, 1 + I, 2}

In[6]:= Divisors[3, GaussianIntegers -> True]
Out[6]= {1, 3}
```

## DivisorSigma

- `DivisorSigma[k, n]`: the divisor function `sigma_k(n)`, the sum of the
  `k`-th powers of the (positive) divisors of `n`.
- `DivisorSigma[k, n, GaussianIntegers -> True]`: the sum over Gaussian-integer
  divisors.

**Options**:
- `GaussianIntegers -> False` (default): ordinary integer divisors. Set to
  `True` to sum over the divisors in `Z[i]`. A non-real Gaussian-integer input
  (e.g. `3 + I`) auto-enables Gaussian mode.

**Features**:
- `Listable`, `NHoldAll`, `Protected`.
- Computed from the multiplicative formula
  `sigma_k(n) = Product_i (p_i^((e_i+1) k) - 1) / (p_i^k - 1)` for
  `n = Product_i p_i^e_i`, so a single path serves every exponent type: exact
  integers and rationals for integer `k`, and symbolic / radical forms for
  symbolic or rational `k`. `k == 0` returns the divisor count `sigma_0(n)`.
- The sign of `n` is ignored; machine integers and GMP bigints are handled
  uniformly.
- In Gaussian mode the product runs over the first-quadrant associates
  (`Re > 0`, `Im >= 0`) of the Gaussian prime factors of `n`. This is the
  multiplicative definition — note it differs from naively summing
  `d^k` over `Divisors[n, GaussianIntegers -> True]`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `DivisorSigma::argrx` message.

```mathematica
In[1]:= DivisorSigma[1, 20]
Out[1]= 42

In[2]:= DivisorSigma[2, 20]
Out[2]= 546

In[3]:= DivisorSigma[0, 12]
Out[3]= 6

In[4]:= DivisorSigma[-2, 10]
Out[4]= 13/10

In[5]:= DivisorSigma[1/2, 12]
Out[5]= (2 (-1 + 2 Sqrt[2]))/((-1 + Sqrt[2]) (-1 + Sqrt[3]))

In[6]:= DivisorSigma[k, {2, 3, 6}]
Out[6]= {(-1 + 2^(2 k))/(-1 + 2^k), (-1 + 3^(2 k))/(-1 + 3^k), ((-1 + 2^(2 k)) (-1 + 3^(2 k)))/((-1 + 2^k) (-1 + 3^k))}

In[7]:= DivisorSigma[2, {1, 2, 3, 4, 5}]
Out[7]= {1, 5, 10, 21, 26}

In[8]:= DivisorSigma[1, 3 + I]
Out[8]= 2 + 6 I

In[9]:= DivisorSigma[2, 6, GaussianIntegers -> True]
Out[9]= -30 + 20 I
```

> **Packed arrays.** `DivisorSigma[k, list]` over an `int64` buffer factors
> each element by trial division in `int64`, with no GMP allocation per
> element. A non-negative integer `k` only: `DivisorSigma[-1, n]` is a
> `Rational`, which no buffer holds.

## MoebiusMu

- `MoebiusMu[n]`: the Möbius function `mu(n)`. For
  `n = u p_1^k_1 ... p_m^k_m` with `u` a unit and `p_i` primes, `mu(n)` is `0`
  unless every `k_i` equals `1`, in which case it is `(-1)^m`. Equivalently
  `mu(1) = 1`, `mu(n) = (-1)^m` for a product of `m` distinct primes, and
  `mu(n) = 0` when `n` has a squared prime factor.

**Features**:
- `Listable`, `Protected`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `0`, `1`, or `-1`.
- The sign of `n` is ignored (`mu(-n) = mu(n)`).
- A non-real Gaussian-integer argument `Complex[a, b]` is handled over `Z[i]`:
  the input is factored into Gaussian primes (the unit factor does not count),
  giving `0` for a repeated Gaussian prime factor and `(-1)^m` otherwise.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `MoebiusMu::argx` message.

```mathematica
In[1]:= MoebiusMu[11]
Out[1]= -1

In[2]:= MoebiusMu[10]
Out[2]= 1

In[3]:= MoebiusMu[1440]
Out[3]= 0

In[4]:= MoebiusMu[{4, 10, 17, 20}]
Out[4]= {0, 1, -1, 0}

In[5]:= MoebiusMu[10^50 + 1]
Out[5]= -1

In[6]:= MoebiusMu[5 + 6 I]
Out[6]= -1
```

> **Packed arrays.** Runs on an `int64` buffer. `MoebiusMu[0]` is undefined,
> so an array containing `0` takes the ordinary path and leaves that element
> unevaluated exactly as the scalar does.

## LiouvilleLambda

- `LiouvilleLambda[n]`: the Liouville function `lambda(n) = (-1)^Omega(n)`, where
  `Omega(n)` is the number of prime factors of `n` counted **with multiplicity**.
  For `n = u p_1^k_1 ... p_m^k_m` with `u` a unit and `p_i` primes,
  `lambda(n) = (-1)^(k_1 + ... + k_m)`. Equivalently `lambda(1) = 1`, and
  `lambda` is `1` when `n` is a product of an even number of primes and `-1`
  otherwise.

**Features**:
- `Listable`, `Protected`.
- Completely multiplicative: `lambda(m n) = lambda(m) lambda(n)`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `1` or `-1`.
- The sign of `n` is ignored (`lambda(-n) = lambda(n)`).
- Gaussian integers: `LiouvilleLambda[n, GaussianIntegers -> True]`, or a
  non-real Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]`
  and counts the Gaussian prime factors with multiplicity. Because `2` factors
  as `-i (1 + i)^2` in `Z[i]`, e.g. `LiouvilleLambda[2, GaussianIntegers -> True]`
  is `1` while `LiouvilleLambda[2]` is `-1`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `LiouvilleLambda::argt` message.

```mathematica
In[1]:= LiouvilleLambda[8]
Out[1]= -1

In[2]:= LiouvilleLambda[9]
Out[2]= 1

In[3]:= LiouvilleLambda[{1, 2, 3, 4, 5, 6}]
Out[3]= {1, -1, -1, 1, -1, 1}

In[4]:= LiouvilleLambda[10^30 + 1]
Out[4]= -1

In[5]:= LiouvilleLambda[2 + I]
Out[5]= -1

In[6]:= LiouvilleLambda[8, GaussianIntegers -> True]
Out[6]= 1
```

## PrimeOmega

- `PrimeOmega[n]`: the number of prime factors of `n` counted **with
  multiplicity**, `Omega(n)`. For `n = u p_1^k_1 ... p_m^k_m` with `u` a unit
  and `p_i` primes, `PrimeOmega[n] = k_1 + ... + k_m`. This is the exponent sum
  that `LiouvilleLambda` reduces via `lambda(n) = (-1)^Omega(n)`.

**Features**:
- `Listable`, `Protected`.
- Completely additive: `Omega(m n) = Omega(m) + Omega(n)`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly).
- `PrimeOmega[1]` (and `PrimeOmega[-1]`) is `0`; the sign of `n` is ignored
  (`Omega(-n) = Omega(n)`).
- Gaussian integers: `PrimeOmega[n, GaussianIntegers -> True]`, or a non-real
  Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]` and counts
  the Gaussian prime factors with multiplicity. Because `2` factors as
  `-i (1 + i)^2` in `Z[i]`, `PrimeOmega[12, GaussianIntegers -> True]` is `5`
  (from `(1 + i)^4 3`) while `PrimeOmega[12]` is `3`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `PrimeOmega::argt` message.

```mathematica
In[1]:= PrimeOmega[30]
Out[1]= 3

In[2]:= PrimeOmega[12]
Out[2]= 3

In[3]:= PrimeOmega[{4, 12, 24}]
Out[3]= {2, 3, 4}

In[4]:= PrimeOmega[30!]
Out[4]= 59

In[5]:= PrimeOmega[5 + 9 I]
Out[5]= 2

In[6]:= PrimeOmega[12, GaussianIntegers -> True]
Out[6]= 5
```

## PrimeNu

- `PrimeNu[n]`: the number of **distinct** prime factors of `n`, `nu(n)`. For
  `n = u p_1^k_1 ... p_m^k_m` with `u` a unit and `p_i` distinct primes,
  `PrimeNu[n] = m` (independent of the exponents). `nu` and `Omega`
  (`PrimeOmega`) coincide exactly when `n` is square-free; for a prime power
  `PrimeNu[n]` is `1`.

**Features**:
- `Listable`, `Protected`.
- Additive on coprime arguments: `nu(m n) = nu(m) + nu(n)` when
  `GCD[m, n] == 1`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly).
- `PrimeNu[1]` (and `PrimeNu[-1]`) is `0`; the sign of `n` is ignored
  (`nu(-n) = nu(n)`).
- Gaussian integers: `PrimeNu[n, GaussianIntegers -> True]`, or a non-real
  Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]` and counts
  the distinct Gaussian prime factors. Because a rational prime `p ≡ 1 (mod 4)`
  splits into two conjugate Gaussian primes, e.g.
  `PrimeNu[105, GaussianIntegers -> True]` is `4` (from `3`, the split pair over
  `5`, and `7`) while `PrimeNu[105]` is `3`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `PrimeNu::argt` message.
- Relations: for a square-free `n`, `MoebiusMu[n] == (-1)^PrimeNu[n]` and
  `LiouvilleLambda[n] == (-1)^PrimeNu[n]`.

```mathematica
In[1]:= PrimeNu[24]
Out[1]= 2

In[2]:= PrimeNu[105]
Out[2]= 3

In[3]:= PrimeNu[{4, 28, 180}]
Out[3]= {1, 2, 3}

In[4]:= PrimeNu[50!]
Out[4]= 15

In[5]:= PrimeNu[3 + I]
Out[5]= 2

In[6]:= PrimeNu[105, GaussianIntegers -> True]
Out[6]= 4
```

## EulerPhi

Gives the Euler totient function $\phi(n)$.
- `EulerPhi[n]`

**Features**:
- `Listable`, `Protected`.
- Counts the number of positive integers less than or equal to $n$ that are relatively prime to $n$.
- Returns 0 for $n = 0$, and handles negative integers via $\phi(-n) = \phi(n)$.
- Accepts arbitrary-precision integers (`BigInt`). Factorization runs in GMP
  through the same trial-division / Pollard-rho / ECM cascade used by
  `FactorInteger`, so inputs of cryptographic size are tractable.
- For a prime decomposition $n = \prod p_i^{k_i}$, computes
  $\phi(n) = n \prod (1 - 1/p_i)$ as $(n / \prod p_i) \prod (p_i - 1)$
  with exact integer arithmetic.

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4

In[2]:= EulerPhi[2^89 - 1]
Out[2]= 618970019642690137449562110
```

> **Packed arrays.** Runs on an `int64` buffer, factoring by trial division
> in `int64`. `EulerPhi[3.]` is not `EulerPhi[3]`, so a real buffer takes the
> ordinary path.

## SquareFreeQ

- `SquareFreeQ[expr]`: Returns `True` if `expr` is a square-free polynomial or
  number, and `False` otherwise.
- `SquareFreeQ[expr, vars]`: Returns `True` if `expr` is square-free with
  respect to the variables `vars` (a `Symbol` or a `List` of `Symbol`s).
- `SquareFreeQ[..., GaussianIntegers -> True | False | Automatic]`: Tests
  square-freeness in `Z[i]` when `True`, in `Z` when `False`. `Automatic`
  (the default) enables Gaussian mode iff the input is a `Complex` literal
  with integer real and imaginary parts.

**Features**:
- `Protected`. Not `Listable` -- passing a list of inputs treats the list as
  the expression (`SquareFreeQ[{1, 2, 3}]` returns `False`).
- Always returns `True` or `False` on a structurally valid call. For inputs
  that are neither a recognised number nor a manifest polynomial -- reals,
  `Sqrt[2]`, `Sin[x]`, `Pi`, strings -- the result is `False`, never
  symbolic.
- An integer `n` is square-free iff `|n|` has no rational prime factor of
  multiplicity `>= 2`. `0` is not square-free; `+/-1` and `+/-p` (for any
  prime `p`) are.
- A rational `p/q` is square-free iff both `p` and `q` are square-free
  integers.
- A polynomial in `vars` is square-free iff for every variable `x_i` in
  `vars` that the polynomial actually depends on, `PolynomialGCD(p, dp/dx_i)`
  has degree `0` in `x_i`. Implementation routes the derivative through
  the `D` builtin and the gcd through `PolynomialGCD`.
- For `GaussianIntegers -> True` (or `Automatic` on a `Complex[Integer,
  Integer]` input), the test factors `N(z) = a^2 + b^2` over `Z` and
  classifies each rational prime by residue `mod 4`:
  - `p == 2`: the Gaussian prime above `2` is `1 + I`; its multiplicity in
    `z` equals the multiplicity of `2` in `N(z)`.
  - `p ≡ 3 (mod 4)`: `p` itself is a Gaussian prime with `N(p) = p^2`;
    its multiplicity in `z` is half the multiplicity of `p` in `N(z)`.
  - `p ≡ 1 (mod 4)`: `p` splits as `pi * conj(pi)`; ambiguous cases (the
    `e_in_norm == 2` slice) are resolved by computing `pi` via a
    Cornacchia search and stripping the multiplicity directly.
- The `Modulus -> p` option is parsed but only `Modulus -> 0` (the default
  no-modulus path) is wired in; any other value -- non-zero integer,
  non-integer, or symbolic -- emits `SquareFreeQ::modnotimpl` and leaves
  the call unevaluated until a polynomial sqfree-mod-`p` test is added.

Diagnostics:
- `SquareFreeQ[]` emits `SquareFreeQ::argb` and stays unevaluated.
- A non-`Rule` past the optional vars slot (e.g. `SquareFreeQ[1, 2, 3]`)
  emits `SquareFreeQ::nonopt` and stays unevaluated.
- `Modulus -> n` with `n != 0` emits `SquareFreeQ::modnotimpl` and stays
  unevaluated (e.g. `SquareFreeQ[x^2 + 1, Modulus -> 2]`).

```mathematica
In[1]:= SquareFreeQ[10]
Out[1]= True

In[2]:= SquareFreeQ[4]
Out[2]= False

In[3]:= SquareFreeQ[20]
Out[3]= False

In[4]:= SquareFreeQ[3 + 2 I]
Out[4]= True

In[5]:= SquareFreeQ[2, GaussianIntegers -> True]
Out[5]= False

In[6]:= SquareFreeQ[2/3]
Out[6]= True

In[7]:= SquareFreeQ[6 + 6 x + x^2]
Out[7]= True

In[8]:= SquareFreeQ[x^3 - x^2 y]
Out[8]= False

In[9]:= SquareFreeQ[x y^2, x]
Out[9]= True

In[10]:= SquareFreeQ[x y^2, y]
Out[10]= False

In[11]:= SquareFreeQ[10^70 + 3]
Out[11]= True
```

## IntegerPartitions

- `IntegerPartitions[n]`: all partitions of `n` in reverse-lexicographic order. `Length[IntegerPartitions[n]] == PartitionsP[n]`.
- `IntegerPartitions[n, k]`: partitions into at most `k` parts.
- `IntegerPartitions[n, {k}]`: exactly `k` parts.
- `IntegerPartitions[n, {kmin, kmax}]`: between `kmin` and `kmax` parts.
- `IntegerPartitions[n, {kmin, kmax, dk}]`: into `kmin, kmin+dk, …` parts.
- `IntegerPartitions[n, kspec, sspec]`: parts drawn only from `sspec` (a list of allowed values). `kspec` of `All` is `{0, Infinity}`; `sspec` of `All` is `Range[n]`.
- `IntegerPartitions[n, kspec, sspec, m]`: the first `m` partitions (`m > 0`) or the last `|m|` (`m < 0`); `m` of `All` is `Infinity`.

`n` and the `s_i` may be rational and/or negative. Attributes: `Protected`.

**Implementation**: a single count-vector enumerator (`src/partitions.c`) over `P = reverse(sspec)`, assigning a count to each allowed part in descending order with exact GMP-rational (`mpq_t`) arithmetic. The default `Range[n]` path is the classic descending-part recursion (bounded by `remaining/part`); explicit `sspec` becomes a coin-problem enumeration bounded by `kmax`. Within each partition, parts appear in `P` order (descending for the default form).

**Edge cases & messages**:
- `IntegerPartitions[1/2]` → `{}` (no integer parts); `IntegerPartitions[1/2, All, {1/6, 1/3}]` → `{{1/3, 1/6}, {1/6, 1/6, 1/6}}`.
- `IntegerPartitions[0]` → `{{}}` (the empty partition).
- `IntegerPartitions::undef` — an unbounded part set (mixed signs or a zero part) with no finite length cap and no positive `m` would be infinitely large; the call is left unevaluated.
- `IntegerPartitions::take` — fewer partitions exist than the requested `m`; the available ones are returned.
- `IntegerPartitions::argb` — fewer than 1 or more than 4 arguments.

```
In[1]:= IntegerPartitions[5]
Out[1]= {{5}, {4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}}

In[2]:= IntegerPartitions[50, All, {6, 9, 20}]
Out[2]= {{20, 9, 9, 6, 6}, {20, 6, 6, 6, 6, 6}}

In[3]:= IntegerPartitions[5, 10, {1, -1}]
Out[3]= {{-1, -1, 1, 1, 1, 1, 1, 1, 1}, {-1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}}
```

## PartitionsP

Gives the number `p(n)` of unrestricted partitions of the integer `n` — the
number of ways to write `n` as a sum of positive integers, where order is
ignored. `PartitionsP[n] == Length[IntegerPartitions[n]]`.

- `PartitionsP[n]`: the partition count `p(n)` as an exact (arbitrary-precision)
  integer.

**Features**:
- `Protected`, `Listable` — `PartitionsP[{2, 4, 6}]` → `{2, 5, 11}`.
- Two engines, dispatched by the size of `n` (threshold `n = 1000`):
  - **Small `n`** — Euler's pentagonal-number-theorem recurrence
    `p(m) = Σ_k (-1)^(k-1) [p(m - g_k) + p(m - g_k')]` over the generalized
    pentagonal numbers `g_k = k(3k∓1)/2`, using exact GMP integers only.
    O(n^1.5) integer additions; O(n·√n)-bit table.
  - **Large `n`** — the non-recursive **Hardy–Ramanujan–Rademacher** exact
    formula evaluated in MPFR (Johansson, arXiv:1205.5991): a convergent
    series `p(n) = Σ_{k=1}^N √(3/k)·(4/(24n−1))·A_k(n)·U(C/k)` with
    `U(x) = cosh x − sinh x / x`, `C = (π/6)√(24n−1)`. The exponential sums
    `A_k(n)` use Selberg's cosine-sum identity (no Dedekind sums). The number
    of terms `N` is chosen from Rademacher's rigorous remainder bound so the
    rounded sum is provably exact; working precision is `≈ C/ln 2` bits plus
    guard bits. Uses only O(√n)-bit precision, so memory stays tiny while the
    recurrence's table would be prohibitive.
- Both engines are cross-validated to agree exactly (see
  `tests/test_partitionsp.c`). The result auto-promotes to a big integer
  (`p(4096)` has 67 digits; `p(10000)` has 106).
- `p(0) = 1`; `p(n) = 0` for `n < 0`.
- Non-integer, symbolic, or astronomically large (big-integer) arguments are
  left unevaluated; `PartitionsP` called with other than one argument emits
  `PartitionsP::argx` and is left unevaluated.

```
In[1]:= Table[PartitionsP[k], {k, 0, 12}]
Out[1]= {1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77}

In[2]:= PartitionsP[100]
Out[2]= 190569292

In[3]:= PartitionsP[4096]
Out[3]= 6927233917602120527467409170319882882996950147283323368445315320451

In[4]:= Table[Times @@ PartitionsP[Last /@ FactorInteger[n]], {n, 12}]
Out[4]= {1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2}
```

## PartitionsQ

Gives the number `q(n)` of partitions of the integer `n` into **distinct**
parts — equivalently (Euler), into **odd** parts (OEIS A000009). It is the
counting companion to `IntegerPartitions` with non-repeating parts.

- `PartitionsQ[n]`: the distinct-parts partition count `q(n)` as an exact
  (arbitrary-precision) integer.

**Features**:
- `Protected`, `Listable` — `PartitionsQ[{2, 4, 6}]` → `{1, 2, 4}`.
- Two engines, dispatched by the size of `n` (threshold `n = 1000`):
  - **Small `n`** — an exact GMP recurrence derived from the Euler identity
    `∏(1−x^k)·∏(1+x^k) = ∏(1−x^{2k})`:
    `q(m) = r(m) + Σ_k (−1)^(k−1) [q(m − g_k) + q(m − g_k')]` over the same
    generalized pentagonal numbers `g_k = k(3k∓1)/2` as `PartitionsP`, with an
    inhomogeneous term `r(m) = (−1)^k` when `m/2` is the generalized pentagonal
    number `g_k` (else `0`). O(n^1.5) integer additions; O(n·√n)-bit table.
  - **Large `n`** — the **Hardy–Ramanujan–Rademacher / Hagis** exact convergent
    series evaluated in MPFR (Hagis, *Amer. J. Math.* 85 (1963) 213–222):
    `q(n) = (π/√(24n+1)) · Σ_{k odd} (1/k)·A_k(n)·I₁(π√(48n+2)/(12k))`, summing
    over **odd `k` only**, with the modified Bessel function `I₁` (computed by
    its everywhere-positive power series, as MPFR has no `I₁`) and the character
    sum `A_k(n) = Σ_{gcd(h,k)=1} cos(π(s(h,k) − s(2h,k)) − 2πnh/k)` built from
    exact Dedekind sums `s(h,k)`. Working precision `≈ π√(n/3)/ln 2` bits plus
    guard bits; memory stays tiny where the recurrence's table would be
    prohibitive.
- Unlike `PartitionsP`, the q-series has no clean closed-form Rademacher
  remainder bound, so the HRR engine's term count is convergence-driven (sum
  `≈√n` odd terms; accept only when the sum lies within `1/4` of an integer and
  the last term is below `1/8`, growing terms/precision otherwise). Correctness
  is pinned by an exhaustive HRR-equals-recurrence cross-check (every `n` up to
  1500 plus a strided sample and large spot values to `n = 100000`; see
  `tests/test_partitionsq.c`).
- `q(0) = 1`; `q(n) = 0` for `n < 0`. The result auto-promotes to a big integer
  (`q(1000)` has 22 digits; `q(100000)` has 245).
- Non-integer, symbolic, or astronomically large (big-integer) arguments are
  left unevaluated; `PartitionsQ` called with other than one argument emits
  `PartitionsQ::argx` and is left unevaluated.

```
In[1]:= Table[PartitionsQ[k], {k, 0, 20}]
Out[1]= {1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 32, 38, 46, 54, 64}

In[2]:= PartitionsQ[100]
Out[2]= 444793

In[3]:= PartitionsQ[{2, 4, 6}]
Out[3]= {1, 2, 4}
```

## ContinuedFraction

Gives the simple continued-fraction expansion of a number.
- `ContinuedFraction[x, n]`: a list of the first `n` terms.
- `ContinuedFraction[x]`: all terms determinable from `x`.

The list `{a1, a2, a3, ...}` corresponds to `a1 + 1/(a2 + 1/(a3 + ...))`.

**Features**:
- `Protected`, `Listable`.
- **Exact rationals** (Integer / BigInt / Rational) use the Euclidean
  algorithm and return the canonical terminating form (last term `>= 2`,
  never `{..., k-1, 1}`). A finite rational may yield fewer than `n` terms.
- **Quadratic surds** `Sqrt[d]` with `d` a non-square integer use the
  periodic surd recurrence. Without a count the result is
  `{a1, ..., {b1, ...}}`, where the bracketed block repeats cyclically; with
  a count the periodic sequence is unrolled to exactly `n` terms. (General
  quadratic irrationals are not recognised symbolically — pass an explicit
  `n` to use the numeric path.) The no-count form is declined for a `d` whose
  period would be astronomically long.
- **Inexact reals** (machine `Real` or arbitrary-precision MPFR) yield terms
  only as far as the input precision determines them, tracking the value's
  uncertainty and stopping when the next term is no longer pinned down.
- **Exact symbolic reals with an explicit `n`** (e.g. `Pi`, `Sqrt[E]`,
  `Exp[Pi Sqrt[163]]`) are numericised at adaptively increasing precision
  until `n` terms are confirmed by two consecutive evaluations.
- Left unevaluated for an exact non-rational, non-quadratic value with no
  count, or a non-real numeric value.

```mathematica
In[1]:= ContinuedFraction[47/17]
Out[1]= {2, 1, 3, 4}

In[2]:= ContinuedFraction[Sqrt[13]]
Out[2]= {3, {1, 1, 1, 1, 6}}

In[3]:= ContinuedFraction[Pi, 20]
Out[3]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2}

In[4]:= ContinuedFraction[N[Pi]]
Out[4]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}

In[5]:= ContinuedFraction[Exp[Pi Sqrt[163]], 10]
Out[5]= {262537412640768743, 1, 1333462407511, 1, 8, 1, 1, 5, 1, 4}
```

## FromContinuedFraction

The inverse of `ContinuedFraction`: reconstructs a value from its continued-fraction terms.
- `FromContinuedFraction[{a1, a2, ..., an}]`: gives `a1 + 1/(a2 + 1/(a3 + ... + 1/an))`.
- `FromContinuedFraction[{a1, ..., am, {b1, ..., bk}}]`: gives the exact
  quadratic irrational whose terms begin with the `ai` then cycle through the
  `bi` forever.

**Features**:
- `Protected` (not `Listable` — the argument is the whole term list).
- The `ai` of the finite form may be **symbolic**; the result is the convergent
  `h_n / k_n` built from the fundamental recurrence
  `h_i = a_i h_{i-1} + h_{i-2}`, `k_i = a_i k_{i-1} + k_{i-2}`, kept in nested
  (un-expanded) form — `Together` collapses it to a flat rational.
- The **periodic** form requires all `ai` and `bi` to be integers. The purely
  periodic tail solves the quadratic
  `k_{k-1} x^2 + (k_{k-2} - h_{k-1}) x - h_{k-2} = 0` (h, k the period's
  convergents); its positive root is then pushed through the leading terms by a
  Möbius transform and rationalised to a single `(P + Q Sqrt[R]) / S` in lowest
  terms.
- `FromContinuedFraction[{}]` is `0`; `FromContinuedFraction[{x}]` is `x`.
- Left unevaluated for a non-list argument, a sub-list anywhere but last, an
  empty period block, or non-integer terms in a periodic form.

```mathematica
In[1]:= FromContinuedFraction[{2, 1, 3, 4}]
Out[1]= 47/17

In[2]:= FromContinuedFraction[{a, b, c, d}]
Out[2]= (1 + a b + (a + (1 + a b) c) d)/(b + (1 + b c) d)

In[3]:= FromContinuedFraction[{8, {2, 2, 1, 7, 1, 2, 2, 16}}]
Out[3]= Sqrt[71]

In[4]:= FromContinuedFraction[{{1, 2, 3, 4}}]
Out[4]= 1/15 (9 + 2 Sqrt[39])

In[5]:= FromContinuedFraction[ContinuedFraction[Pi, 3]]
Out[5]= 333/106
```

