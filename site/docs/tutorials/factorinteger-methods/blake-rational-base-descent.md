# Blake's rational base descent

`Method -> "BlakeRationalBaseDescent"` is a **special-purpose** factoring method
introduced by Sam Blake. It factors \(N = P\,M\) quickly when the larger factor is
close to a power of a rational number,

\[
P \;\approx\; c\left(\frac{a}{b}\right)^{k}, \qquad a > b \ge 1,\ \gcd(a,b)=1,
\]

with \(c\) moderate. The cofactor \(M\) can be anything of a suitable size. Like
Pollard \(p-1\) (smooth \(p-1\)) or Fermat (close factors), the method exploits a
property of the factor rather than the size of \(N\). Given the right base it is
polynomial in \(\log N\). Against an unstructured semiprime it is hopeless, and the
analysis shows why.

```mathematica
In[1]:= FactorInteger`BlakeRationalBaseDescent[1001, "Base" -> 22/7]
Out[1]= {{7, 1}, {11, 1}, {13, 1}}
```

---

## 1. The mathematical idea

### The descent

Start from \(Q_0 = N\) and repeatedly divide by the base, rounding down:

\[
Q_k = \left\lfloor \frac{b\,Q_{k-1}}{a} \right\rfloor .
\]

After \(k\) steps \(Q_k \approx (b/a)^k N = (b/a)^k P M\). If \(P\approx c(a/b)^k\),
the powers cancel and \(Q_k\approx c\,M\), an integer multiple of the hidden
factor. A small gcd search then finds it: some \(\gcd(Q_k - j, N)\) with
\(|j|\le J\) equals \(M\gcd(c,P)\).

**Theorem** (Blake, Theorem 1). Let \(\Delta = P - \lfloor c(a/b)^k\rfloor\).
Suppose

\[
\left(\frac{a}{b}\right)^{k} > M \qquad\text{and}\qquad |\Delta| < \frac{(a/b)^k}{M}.
\]

Then \(|Q_k - cM| \le J = \lceil 2 + a/(a-b)\rceil\), so the window
\(\gcd(Q_k - j, N)\), \(|j|\le J\), returns \(M\gcd(c,P)\), a non-trivial factor.

The proof unrolls the recurrence: \(Q_k = (b/a)^k N - E_k\), with truncation error
\(0\le E_k < a/(a-b)\). Three terms are each bounded:

- the \(\Delta\) term, \(|\Delta(b/a)^k M| < 1\);
- the fractional part of \(c(a/b)^k\), times \((b/a)^kM < 1\);
- \(E_k\).

Some consequences:

- **One descent tests every exponent \(k\) and every coefficient \(c\) at once.**
  For a fixed base the cost is \(O(\log_{a/b}N)\) steps of \(2J+1\) gcds, which is
  \(O(\log^3 N)\) when \(a/b\) is bounded away from 1.
- **The structure must be on the larger factor.** The hypotheses force
  \(P > cM \ge M\). For a balanced semiprime only \(c = 1\) can work.
- **The tolerance is generous.** \(|\Delta|\) may be as large as \((a/b)^k/M\),
  which for a small cofactor is astronomically larger than 1.

### The search over bases, and why it is a comb

If the base is unknown, it must be searched for. The companion paper proves a
*miss-distance law*: if the trial base has relative error \(u\), the descent
misses its target by \(\eta \approx -c\,k\,M\,u\). So the bases that succeed form a
**comb** of isolated teeth of relative width \(\tilde O(1/(kM))\). There are also
harmonics at \(a/d\) for divisors \(d\mid a\), with an inflated coefficient
\(c\,d^k\). Between the teeth \(\eta\) is enormous and carries no usable signal.
In fact, any efficient monotone "closer/further" signal would factor every
semiprime by bisection. Bases are therefore *enumerated*, not optimised.
Mathilda walks the Stern–Brocot tree breadth-first, which visits every rational
once in lowest terms, coarse to fine, and skips perfect powers \((x/y)^m\), since
they duplicate \(x/y\).

### Lattices widen the teeth

The detection step asks whether \(Q_k\) is near a multiple of an unknown divisor
of \(N\). That is a partial approximate-common-divisor problem.
Howgrave-Graham's lattice method solves it for a tolerance up to
\(X \le N^{\beta^2}\) when \(M \ge N^\beta\), against \(J = O(a/(a-b))\) for the
gcd window. Every tooth then widens by a factor of about \(X\), and the expected
base search shrinks from \(\tilde\Theta(M)\) to \(\tilde O(N^{\beta(1-\beta)})\),
which is at most \(\tilde O(N^{1/4})\). With a known base, coefficient and exponent,
this is exactly Coppersmith-style factoring with high bits known.

### Truncation divisibility

For an integer base \(B\) the descent simply truncates base-\(B\) digits:
\(Q_k = \lfloor N/B^k\rfloor\). If \(d\mid N\), \(\gcd(d,B)=1\) and
\(d \mid N \bmod B^k\), then \(d \mid \lfloor N/B^k\rfloor\). For Cunningham numbers
\(B^m\pm1\) and repunits the remainders are themselves structured. So the descent
meets their algebraic (cyclotomic) divisors at aligned depths, with no knowledge
of the cyclotomic structure.

### References

1. S. Blake, "Rational Base Descent: A Deterministic Algorithm for Factoring
   Structured Semiprimes," preprint (2026).
2. S. Blake, "Lattice-Augmented Rational Base Descent," preprint (2026).
3. N. Howgrave-Graham, "Approximate integer common divisors," *CaLC 2001*, LNCS
   2146, Springer (2001) 51–66.
4. D. Coppersmith, "Small solutions to polynomial equations, and low exponent RSA
   vulnerabilities," *Journal of Cryptology* **10** (1997) 233–260.
5. A. May, "Using LLL-reduction for solving RSA and factorization problems," in
   *The LLL Algorithm*, Springer (2010) 315–348.
6. A. K. Lenstra, H. W. Lenstra Jr., L. Lovász, "Factoring polynomials with
   rational coefficients," *Mathematische Annalen* **261** (1982) 515–534.
7. R. L. Graham, D. E. Knuth, O. Patashnik, *Concrete Mathematics*, 2nd ed.,
   §4.5 (the Stern–Brocot tree).
8. F. Lawrence, "Factorisation of numbers," *Messenger of Mathematics* **24**
   (1895) 100–109, and W. Hart, "A one line factoring algorithm," *J. Aust. Math.
   Soc.* **92** (2012) 61–69 — earlier methods for factors near known ratios or
   powers.

---

## 2. How Mathilda realises it

The implementation is in `src/numbertheory/factorinteger/`:

- `fi_rbd.c` holds the Mathilda glue: options, variants, the threaded base search
  and the hit cache.
- A reentrant GMP/FLINT core (`fi_rbd.h`) provides the rest:
  - `fi_rbd_descent.c`: the iterated and exact descents, and the batched gcds;
  - `fi_rbd_bases.c`: the Stern–Brocot and \(a+b\) enumerations, and the
    perfect-power filter;
  - `fi_rbd_lattice.c`: the Howgrave-Graham detector;
  - `fi_rbd_lll.c`: FLINT's `fmpz_lll`, or a native fraction-free integral LLL.

### The six variants

Each variant is selected with `"Variant" -> name`, and each also has its own head,
``FactorInteger`BlakeRationalBaseDescent<name>``.

| Variant | What it does | Use it when |
|---|---|---|
| `"Exact"` (default) | Scaled exact descent. The recurrence runs on \(S N\), \(S = 2^s \ge 2a/(a-b)\), so the target is always one of **four** consecutive values, whatever the base. Testing is deferred while \(Q > \sqrt{c_{\max}N}\), since no \(cM\) can lie there; that depth is reached with one division. Window products are folded mod \(N\), with one gcd per block, replayed only on a hit. | The general case: structured larger factor, base known or searched. |
| `"Classic"` | The paper's reference algorithm R1–R4. It tests \(\gcd(Q-j,N)\), \(\lvert j\rvert\le\lceil 2+a/(a-b)\rceil\), at every depth. | Reproducing the paper, or as a cross-check. About 5–10× slower than Exact. |
| `"AllHits"` | Runs the reference descent to the end, records **every** distinct non-trivial gcd, and refines the record into pairwise-coprime divisors. Composite pieces are re-run with AllHits. | Cunningham numbers and repunits with an integer base (truncation divisibility). |
| `"Lattice"` | Replaces the window by the Howgrave-Graham detector. LLL runs on \(N^{t-i}(Q+x)^i\) and \(x^j(Q+x)^t\), only in the depth band \(X \le Q \le \sqrt{c_{\max}N}+X\), in parallel over depths. | The factor is near \(c(a/b)^k\), but the miss is far beyond the window. |
| `"TwoPass"` | Pass 1 is Exact over the bases. Pass 2 is Lattice over a coarse grid of `"BaseRange"` (relative spacing \(X/(2c\,k\,N^\beta)\), each point realised by its simplest rational), or over the Stern–Brocot bases of height 16. | The paper's practical search. Expensive: see §5. |
| `"Coppersmith"` | Base, coefficient and exponent all known: one detector call at \(\lfloor c(a/b)^k\rfloor\). | Factoring with high bits known. |

### Suboptions

| Suboption | Default | Meaning |
|---|---|---|
| `"Variant"` | `"Exact"` | one of the six names above |
| `"Base"` | enumerate | `a/b` or a list of bases. Each is reduced, and inverted if below 1; `1` is rejected (`FactorInteger::mopt`) |
| `"Order"` | `"SternBrocot"` | base enumeration order: breadth-first Stern–Brocot, or `"Sum"` (the papers' \(a+b\) order) |
| `"Height"` | 1000 (16 for the TwoPass lattice pass) | largest numerator enumerated |
| `"Bases"` | 200000 | enumeration budget |
| `"MaxCoefficient"` | \(2^{64}\) | largest \(c\) searched for. Sets the deferred-testing ceiling \(\sqrt{c_{\max}N}\); `Infinity` tests every depth |
| `"Multiplier"` | none | descend from \(mN\), with gcds still taken with \(N\) |
| `"Tolerance"` | from the LLL bound | lattice tolerance \(X\). Given \(X\), the least dimension whose guarantee covers it is chosen |
| `"Beta"` | 1/2 (Coppersmith: from the size of the guess) | the lattice assumes a divisor \(M \ge N^\beta\), \(0<\beta\le 1\) |
| `"Dimension"` | 12 | lattice dimension \(d\) (2–48). Without `"Tolerance"`, \(X\) is the largest the guarantee allows at \(d\), less two bits |
| `"Multiplicity"` | chosen with \(d\) | Howgrave-Graham multiplicity \(t\) (at most 47) |
| `"BaseRange"` | none | `{lo, hi}` with \(1 < lo < hi\): the TwoPass lattice grid |
| `"Coefficient"` | — | \(c\), required by `"Coppersmith"` |
| `"Exponent"` | — | \(k\), required by `"Coppersmith"` |

**Behaviour.**

- The suboptions describe the *input*. If a composite cofactor cannot be split by
  the configured search, it gets a light default search: Exact, Stern–Brocot
  height 64, at most 2000 bases. AllHits instead re-runs AllHits.
- Anything still composite is returned with `FactorInteger::nofac`.
- Many-base searches run blocks of 2048 bases on every thread. The
  lowest-indexed hit wins, so the answer does not depend on the thread count.
- Unlike `Automatic`, the method does **no** trial division. Only powers of 2
  and perfect powers are removed before the splitter runs.

---

## 3. Invoking it

The method can be called in three ways:

- `FactorInteger[n, Method -> "BlakeRationalBaseDescent"]`;
- `FactorInteger[n, Method -> {"BlakeRationalBaseDescent", "Variant" -> "Lattice", "Base" -> 22/7}]`;
- ``FactorInteger`BlakeRationalBaseDescentLattice[n, "Base" -> 22/7]``.

The examples below construct their inputs in the page.

---

## 4. Worked examples

### 4.1 Watching the descent by hand

Plant \(P = \text{nextprime}(\lfloor 3\,(22/7)^{150}\rfloor)\), so \(c=3\) and
\(k=150\), and multiply by a 40-bit prime. Then run 150 descent steps directly:

```mathematica
In[1]:= p = NextPrime[Floor[3 (22/7)^150]]; q = NextPrime[10^12]; n = p q;
Out[1]= Null

In[2]:= Q = Nest[Floor[7 #/22] &, n, 150]; Q - 3 q
Out[2]= -1

In[3]:= Ceiling[2 + 22/(22 - 7)]
Out[3]= 4

In[4]:= GCD[Q + 1, n] == q
Out[4]= True
```

After 150 steps \(Q_{150}\) sits one below \(3q\), well inside the window
\(J = 4\). The miss is non-positive, as the paper's one-sided remark predicts. The
builtin does the same thing, with the base given or found by search:

```mathematica
In[5]:= FactorInteger`BlakeRationalBaseDescent[n, "Base" -> 22/7]
Out[5]= {{1000000000039, 1}, {1190741187132579502456731851112421085210615890827643793135496647841845583593, 1}}

In[6]:= FactorInteger`BlakeRationalBaseDescent[n]
Out[6]= {{1000000000039, 1}, {1190741187132579502456731851112421085210615890827643793135496647841845583593, 1}}
```

With the base given, the run takes about 0.5 ms. The search over Stern–Brocot
bases takes about 13 ms on 16 threads. The search has to reach the base, so the
height bound matters. Below the numerator 22 nothing fires, and the composite is
returned (`nofac` is printed on stderr):

```mathematica
In[7]:= Quiet[FactorInteger`BlakeRationalBaseDescent[n, "Height" -> 10]]
Out[7]= {{1190741187179018408754902451708233627404000313150857812877774755774129952858831977760127, 1}}

In[8]:= FactorInteger`BlakeRationalBaseDescent[n, "Height" -> 22, "Order" -> "Sum"] == FactorInteger[n]
Out[8]= True
```

### 4.2 The harmonic tooth and `"MaxCoefficient"`

Because \(22 = 2\cdot 11\), the base \(11/7\) represents the same factor with
coefficient \(c' = 3\cdot 2^{150}\). This is far above the default \(c_{\max} = 2^{64}\),
so Exact's deferred testing skips the depth where the target lies. Raising the
ceiling (or removing it) finds the harmonic:

```mathematica
In[9]:= Quiet[FactorInteger`BlakeRationalBaseDescent[n, "Base" -> 11/7]] == {{n, 1}}
Out[9]= True

In[10]:= FactorInteger`BlakeRationalBaseDescent[n, "Base" -> 11/7, "MaxCoefficient" -> Infinity] == FactorInteger[n]
Out[10]= True
```

`Classic` never defers, so it finds the harmonic as well, at about 4× the cost.

### 4.3 The paper's Example 1: a 1140-bit modulus at the boundary

Example 1 of the first paper has an 886-bit structured prime \(p_{\min}\) (with
\(c = 1002003004005\), \(k = 512\)) and a 255-bit cofactor, the first 77 digits of
\(\pi\). The defect is set *exactly at* the theorem's boundary,
\(\Delta = -\lfloor (22/7)^{512}/q\rfloor\):

```mathematica
In[11]:= n1 = 13434917067328449643383271289062122492729438008563207396013386022436439625786814860310635092648849808707283885079425559380060603741711887741983827702536077124606526390026531236424033896417798802339681855054916007583273591721147874872971014173899334436432277372280112100845364029259802969568473368054409653616559493278947418005745364854816116209;
Out[11]= Null

In[12]:= q255 = 31415926535897932384626433832795028841971693993751058209749445923078164119021;
Out[12]= Null

In[13]:= {BitLength[n1], PrimeQ[n1/q255], BitLength[n1/q255]}
Out[13]= {1140, True, 886}

In[14]:= N[n1/q255 - Floor[1002003004005 (22/7)^512]]
Out[14]= -1.35852e+178

In[15]:= MemberQ[FactorInteger`BlakeRationalBaseDescent[n1, "Base" -> 22/7][[All, 1]], q255]
Out[15]= True

In[16]:= MemberQ[FactorInteger`BlakeRationalBaseDescentClassic[n1, "Base" -> 22/7][[All, 1]], q255]
Out[16]= True
```

The defect is \(10^{178}\) in absolute terms, yet the descent finds the 255-bit
factor at step 512. The scaled exact descent takes about 4 ms (1536 gcds in the
paper's count) and Classic about 22 ms (4604 gcds). These are this session's
timings on a loaded machine; `docs/spec/builtins/number-theory.md` records 2.8 ms
and 27 ms.

### 4.4 A visible signature in base 3

For an integer base the descent truncates digits. A structured factor therefore
leaves a long run of zeros in the base-\(B\) expansion of \(N\). Its length is the
base-\(B\) logarithm of the slack in the theorem. The first paper's Open Question
3 gives an example:

```mathematica
In[17]:= n3 = 14046133535443172227430326880457961558090460914355488360756493;
Out[17]= Null

In[18]:= Max[Length /@ Select[Split[IntegerDigits[n3, 3]], First[#] == 0 &]]
Out[18]= 26

In[19]:= FactorInteger`BlakeRationalBaseDescent[n3, "Base" -> 3]
Out[19]= {{59604644762897191, 1}, {235655016338368235499067731947072878028222123, 1}}
```

### 4.5 `"AllHits"`: Cunningham numbers and repunits

The first paper's Example 2 is \(2^{1041}+1\) with base 2. The full descent finds
3, 2083, and the 105-digit algebraic divisor \(2^{347}+1\), whose odd part is the
Wagstaff prime \((2^{347}+1)/3\):

```mathematica
In[20]:= Select[Quiet[FactorInteger`BlakeRationalBaseDescentAllHits[2^1041 + 1, "Base" -> 2]][[All, 1]], # < 10^120 &] == {3, 2083, (2^347 + 1)/3}
Out[20]= True
```

The fourth piece is a 206-digit cofactor, which AllHits leaves unsplit (with
`nofac`). Repunits \(R_m = (10^m-1)/9\) show the mechanism most transparently.
With base 10, \(Q_k = R_{m-k}\) exactly, and \(\gcd(R_i,R_j) = R_{\gcd(i,j)}\):

```mathematica
In[21]:= R[m_] := (10^m - 1)/9
Out[21]= Null

In[22]:= f = Quiet[FactorInteger`BlakeRationalBaseDescentAllHits[R[60], "Base" -> 10]]; Times @@ (Power @@@ f) == R[60]
Out[22]= True

In[23]:= Select[f[[All, 1]], ! PrimeQ[#] &]
Out[23]= {520801, 99009901, 165573604901641}
```

The hit record splits \(R_{60}\) into 17 coprime pieces in under a millisecond.
Three of them are composite: \(520801 = 241\cdot2161\), and so on. No window
of the base-10 descent separates their prime factors, and re-running AllHits on
the pieces does not either. The right workflow is to hand the pieces to `Automatic`:

```mathematica
In[24]:= Union @@ (FactorInteger[#][[All, 1]] & /@ f[[All, 1]]) == FactorInteger[R[60]][[All, 1]]
Out[24]= True
```

The paper's \(R_{2224}\) (2224 digits) takes about 0.3 s. Its record yields 13
pieces, the largest the 1104-digit cofactor of \(R_{1112}\):

```mathematica
In[25]:= IntegerLength /@ Quiet[FactorInteger`BlakeRationalBaseDescentAllHits[R[2224], "Base" -> 10]][[All, 1]]
Out[25]= {2, 2, 2, 3, 3, 3, 4, 7, 136, 139, 273, 552, 1104}
```

### 4.6 `"Multiplier"`: the paper's Open Question 2

Take \(p = \text{nextprime}(2^{304})\) and \(N = 4779846103\,p\). Here \(p\) is close to
a power of 2, but base 3 sees nothing. The multiplier
\(m = \text{round}(3^{211}/2^{304})\) moves \(mp\) within \(1.6\times10^{-10}\) of
\(3^{211}\), and the descent on \(mN\) then finds the cofactor. The cofactor
fallback splits \(4779846103 = 24473\cdot195311\) further:

```mathematica
In[26]:= nm = 4779846103 NextPrime[2^304]; m = Round[3^211/2^304]
Out[26]= 1443663791

In[27]:= Quiet[FactorInteger`BlakeRationalBaseDescent[nm, "Base" -> 3]] == {{nm, 1}}
Out[27]= True

In[28]:= FactorInteger`BlakeRationalBaseDescent[nm, "Base" -> 3, "Multiplier" -> m]
Out[28]= {{24473, 1}, {195311, 1}, {32592575621351777380295131014550050576823494298654980010178247189670100796213387298934358053, 1}}
```

Finding such an \(m\) without knowing \(p\) costs about as much as the base search
it replaces. This is a demonstration, not a shortcut.

### 4.7 `"Lattice"`: a miss of 2²⁸

The companion paper's experiment is a balanced 197-bit semiprime. Its 98-bit \(p\)
is pushed away from \(\lfloor 3(22/7)^{59}\rfloor\), so the descent at depth 59
misses \(3q\) by about \(2^{28}\). The gcd window (\(J=4\)) cannot see that, and the
lattice detector can:

```mathematica
In[29]:= nl = 195539142822724777401163665279552940381006352574199957151137;
Out[29]= Null

In[30]:= Quiet[FactorInteger`BlakeRationalBaseDescent[nl, "Base" -> 22/7]] == {{nl, 1}}
Out[30]= True

In[31]:= First /@ FactorInteger`BlakeRationalBaseDescentLattice[nl, "Base" -> 22/7]
Out[31]= {296454907933147663952202071141, 659591518271709676788151702157}

In[32]:= First /@ FactorInteger`BlakeRationalBaseDescentLattice[nl, "Base" -> 22/7, "Tolerance" -> 2^30, "Multiplicity" -> 6, "Dimension" -> 12]
Out[32]= {296454907933147663952202071141, 659591518271709676788151702157}
```

The first Lattice call, with automatic \((t,d)\) and \(X\), takes about 60 ms. The
second uses the paper's \((6, 12)\) with \(X = 2^{30}\). This instance sits *outside*
the plain theorem's regime: \((22/7)^{59} < q\), so \(cM = 3q\) exceeds \(P\). The
depth band must therefore reach above \(\sqrt{4N}\). A `"MaxCoefficient"` of 4 cuts
the band just below the target and fails; 5 succeeds, about 5× faster than the
default \(2^{64}\):

```mathematica
In[33]:= Quiet[FactorInteger`BlakeRationalBaseDescentLattice[nl, "Base" -> 22/7, "MaxCoefficient" -> 4]] == {{nl, 1}}
Out[33]= True

In[34]:= First /@ FactorInteger`BlakeRationalBaseDescentLattice[nl, "Base" -> 22/7, "MaxCoefficient" -> 5]
Out[34]= {296454907933147663952202071141, 659591518271709676788151702157}
```

### 4.8 `"Coppersmith"` and `"TwoPass"`

When the base, coefficient and exponent are all known, a single detector call at
\(\lfloor 3(22/7)^{59}\rfloor\) recovers the factor in a few milliseconds. The wrong
exponent gives nothing, and missing parameters leave the call unevaluated with
`FactorInteger::rbdcop`:

```mathematica
In[35]:= First /@ FactorInteger`BlakeRationalBaseDescentCoppersmith[nl, "Base" -> 22/7, "Coefficient" -> 3, "Exponent" -> 59]
Out[35]= {296454907933147663952202071141, 659591518271709676788151702157}

In[36]:= Quiet[FactorInteger`BlakeRationalBaseDescentCoppersmith[nl, "Base" -> 22/7, "Coefficient" -> 3, "Exponent" -> 58]] == {{nl, 1}}
Out[36]= True

In[37]:= Quiet[FactorInteger`BlakeRationalBaseDescentCoppersmith[1001, "Base" -> 2]]
Out[37]= FactorInteger`BlakeRationalBaseDescentCoppersmith[1001, "Base" -> 2]

In[38]:= First /@ FactorInteger`BlakeRationalBaseDescentTwoPass[nl, "Base" -> {3/2, 22/7}]
Out[38]= {296454907933147663952202071141, 659591518271709676788151702157}
```

In the TwoPass call, pass 1 (Exact) fails on both bases. Pass 2 then runs the
lattice detector over both, about 0.3 s in all, most of it in the many depths of
the slow base \(3/2\).

---

## 5. Where it fails: generic semiprimes

Take a 40-digit semiprime whose factors are not close to a rational power:

```mathematica
In[39]:= g = NextPrime[Floor[Sqrt[2] 10^19]] NextPrime[Floor[Sqrt[3] 10^19]]
Out[39]= 244948974278317811047525616383344641411

In[40]:= Quiet[FactorInteger`BlakeRationalBaseDescent[g]]
Out[40]= {{244948974278317811047525616383344641411, 1}}

In[41]:= FactorInteger[g]
Out[41]= {{14142135623730950533, 1}, {17320508075688772967, 1}}
```

The full default search tries 200 000 Stern–Brocot bases in about a second and
finds nothing. `Automatic` (rho, then the QS) factors the number in about 30 ms.
This is no tuning failure. It is the method's complexity:

- A typical real \((p/c)^{1/k}\) has rational approximations \(a/b\) with error
  about \(1/b^2\). The theorem's tolerance then requires \(c\,k\,q < ab\).
- For a balanced \(N\) that means \(ab \gtrsim \sqrt N \log N\). The number of
  bases to try is therefore \(\tilde\Theta(\sqrt N)\), worse than trial division.
- The lattice detector improves this to \(\tilde O(N^{1/4})\), the classical
  deterministic baseline, and still exponential.
- The papers bound the chance that an honestly generated RSA-2048 prime is
  vulnerable to any base with \(ab \le 10^{12}\): below \(2^{-960}\).

**RSA moduli from a sound key generator are not at risk.** The realistic use is
the defender's: screening keys for a planted rational-base backdoor (a *weak*
backdoor in the kleptographic sense, since anyone can run the same search), and
recovering the algebraic structure of Cunningham-type numbers.

TwoPass's lattice grid is the expensive part of a blind search. Its spacing
\(X/(2c_{\max}kN^\beta)\) is tiny with the default \(c_{\max}=2^{64}\). Even a range as
narrow as `"BaseRange" -> {314285/100000, 314287/100000}` around \(22/7\) on `nl`
does not finish within a 30-second CPU budget. Prefer explicit `"Base"` lists, a
small `"MaxCoefficient"`, or plain Lattice when the base is known.

---

## 6. Limitations

- It is special-purpose: it succeeds only when the **larger** factor is close to
  \(c(a/b)^k\) with a base low in the Stern–Brocot tree (or one you supply).
- There is no trial division. A small prime factor is found only if some descent
  window happens to catch it. Combine with `Automatic` for cleanup.
- AllHits returns pairwise-coprime pieces, some composite (with `nofac`). It does
  not finish the factorisation.
- The default \(c_{\max} = 2^{64}\) hides harmonics with huge coefficients. Use
  `"MaxCoefficient" -> Infinity` (or `"Classic"`) to see them.
- The lattice detector's guarantee follows the LLL bound, which is pessimistic.
  A `"Tolerance"` beyond it is still attempted at \((t,d)=(6,12)\).
- `TimeConstrained` budgets are process CPU time and the search is threaded, so
  a budget runs out up to 16× faster than wall-clock time.

---

## See also

- [Overview and the Automatic cascade](index.md), [performance](performance.md).
- Other special-purpose methods: [Fermat and Lehman](fermat-lehman.md) (close
  factors), [p−1 and p+1](pminus1-pplus1.md) (smooth group orders),
  [ECM](ecm.md) (small factors).
- General-purpose methods: [CFRAC and Dixon](cfrac-dixon.md),
  [the quadratic sieve](quadratic-sieve.md).
- `?FactorInteger`BlakeRationalBaseDescent` in the REPL for the full docstring.
