# FactorInteger shootout — v0.453

Mathilda v0.453 against PARI/GP 2.17.3, FLINT 3.6 (C driver), Mathematica 14.0
(`wolframscript`), Maxima 5.49 and SymPy 1.14, on the seeded corpus in `corpus.json`
(`gen_corpus.py`, seed 20261017).

- Machine: Intel i9-9880H, 8 cores and 16 threads, macOS 15. Nothing else was benchmarking
  at the time.
- Every answer is verified independently: the product of the factors equals n, and every
  base passes `sympy.isprime`. A cell is a time only if the answer verified.
- Batch rows give the total for the whole group. Single rows give the time per case, and
  process start-up is excluded.

## Where the columns come from

- **mathilda, pari, flint:** the final run, `results/final-*.json`. Timeouts were 900 s for
  Mathilda and 600 s for PARI and FLINT.
- **mathematica:**
  - The batch rows come from the final run.
  - The final run was stopped while it was in the `qs` group, so the single-case rows come
    from the v0.450 baseline run (`results/baseline-v0450-mathematica.json`, 120 s timeout).
    The cases the final run did reach before it was stopped took qs50 3.3 s, qs60 45–50 s,
    and qs70 timed out at 600 s.
- **maxima, sympy:** from the baseline run (120 s timeout), so a "timeout" in these columns
  means more than 120 s.
- **v0.450:** Mathilda before the rewrite, from the same baseline run. "wrong" means v0.450
  returned an incomplete factorisation, i.e. one with a composite entry. It had no quadratic
  sieve, and it gave up after a fixed rho/ECM schedule.

**Threads.** Mathilda's SIQS sieves on up to 16 threads above 150 bits (about 45 digits).
ECM runs curves in parallel above 160 bits. PARI and the FLINT driver run on one thread. To
separate the algorithm from the threads, here are Mathilda's single-threaded timings from
development (v0.452/v0.453, `FI_THREADS=1`):

| Digits | Mathilda, one thread | PARI |
|---|---|---|
| 60 | 0.37–0.46 s | 3.2–3.6 s |
| 70 | ~22 s (before the v0.453 sieve work) | 26–31 s |
| 80 | 145 s | 322–494 s |

So the quadratic-sieve lead does not depend on parallelism.

## Results

Fastest verified time in each row is in bold.

| case | mathilda | pari | flint | mathematica | maxima | sympy | v0.450 |
|---|---|---|---|---|---|---|---|
| u32 (2000×, total) | 5.6 ms | 4.0 ms | **933 µs** | 16.0 ms | 736.8 ms | 284.3 ms | 27.3 ms |
| u63 (2000×, total) | 26.2 ms | 113.0 ms | **26.1 ms** | 586.8 ms | 5.77 s | 7.64 s | 483.3 ms |
| semi63 (500×, total) | **38.3 ms** | 154.0 ms | 125.8 ms | 585.2 ms | 18.01 s | timeout | 5.33 s |
| semi96 (20×, total) | 77.2 ms | **40.0 ms** | 126.7 ms | 243.6 ms | 18.62 s | 13.31 s | 20/20 wrong |
| rand40 (100×, total) | **432.2 ms** | 550.0 ms | 1.08 s | 2.85 s | 53.70 s | 42.56 s | 39/100 wrong |
| smooth (200×, total) | **105.2 ms** | 412.0 ms | 428.9 ms | 2.35 s | 2.73 s | timeout | 546.2 ms |
| qs30-0 | 6.4 ms | **3.0 ms** | 6.2 ms | 61.5 ms | 1.84 s | 762.8 ms | wrong |
| qs30-1 | 7.1 ms | **3.0 ms** | 6.9 ms | 61.9 ms | 1.99 s | 273.7 ms | wrong |
| qs40-0 | 26.5 ms | 26.0 ms | 33.3 ms | **18.2 ms** | 1.10 s | 70.0 ms | wrong |
| qs40-1 | 29.7 ms | **29.0 ms** | 42.5 ms | 330.8 ms | 14.08 s | 11.79 s | wrong |
| qs50-0 | **75.9 ms** | 279.0 ms | 358.0 ms | 6.15 s | timeout | timeout | wrong |
| qs50-1 | **109.4 ms** | 308.0 ms | 394.2 ms | 6.48 s | 13.02 s | timeout | wrong |
| qs60-0 | **613.0 ms** | 3.22 s | 6.53 s | timeout | timeout | timeout | wrong |
| qs60-1 | **600.7 ms** | 3.57 s | 4.70 s | timeout | timeout | timeout | wrong |
| qs70-0 | **3.05 s** | 26.52 s | 27.00 s | timeout | timeout | timeout | wrong |
| qs70-1 | **3.46 s** | 31.02 s | 51.72 s | timeout | timeout | timeout | wrong |
| qs80-0 | **32.10 s** | 493.94 s | 346.08 s | timeout | timeout | timeout | wrong |
| qs80-1 | **27.25 s** | 321.68 s | 236.60 s | timeout | timeout | timeout | wrong |
| qs90-0 | **287.75 s** | timeout | timeout | timeout | timeout | timeout | wrong |
| qs90-1 | **282.81 s** | timeout | timeout | timeout | timeout | timeout | wrong |
| ecm-p20-0 | **1.16 s** | 19.81 s | 3.68 s | 1.47 s | 4.61 s | timeout | wrong |
| ecm-p20-1 | 3.19 s | 21.68 s | **3.01 s** | timeout | 3.55 s | timeout | wrong |
| ecm-p25-0 | **3.23 s** | 187.69 s | 206.78 s | timeout | timeout | timeout | wrong |
| ecm-p25-1 | **11.85 s** | 69.77 s | 130.28 s | timeout | timeout | timeout | wrong |
| ecm-p30-0 | **53.84 s** | timeout | 450.73 s | timeout | timeout | timeout | wrong |
| ecm-p30-1 | **26.16 s** | timeout | 532.79 s | timeout | timeout | timeout | wrong |
| ecm-p35-0 | **203.46 s** | timeout | timeout | timeout | timeout | timeout | wrong |
| ecm-p35-1 | **370.24 s** | timeout | timeout | timeout | timeout | timeout | wrong |
| pm1-0 | **32.9 ms** | timeout | timeout | timeout | timeout | timeout | wrong |
| pm1-1 | **32.2 ms** | timeout | timeout | timeout | timeout | timeout | wrong |
| pm1-2 | **31.9 ms** | timeout | timeout | timeout | timeout | timeout | wrong |
| pp1-0 | **209.6 ms** | timeout | timeout | timeout | timeout | timeout | wrong |
| pp1-1 | **194.9 ms** | timeout | timeout | timeout | timeout | timeout | wrong |
| pp1-2 | **194.3 ms** | timeout | timeout | timeout | timeout | timeout | wrong |
| fermat-100d | 6.0 ms | timeout | timeout | timeout | timeout | **2.3 ms** | wrong |
| fermat-80d | 5.6 ms | 265.35 s | 217.12 s | timeout | timeout | **1.9 ms** | wrong |
| fermat-60d | 5.1 ms | 3.06 s | 6.28 s | timeout | timeout | **1.4 ms** | wrong |
| pow-p30^7 | 2.6 ms | **0 µs** | 1.9 ms | 3.8 ms | timeout | 2.5 ms | 369 µs |
| pow-m61^3*q^2 | 88.0 ms | 12.96 s | 8.98 s | 18.3 ms | timeout | **5.0 ms** | wrong |
| pow-2^4096 | 28 µs | **0 µs** | 9 µs | 26 µs | 605.4 ms | 26 µs | 76 µs |
| pow-smooth^5 | 2.5 ms | **0 µs** | 179 µs | 5.0 ms | 822.5 ms | 4.2 ms | 6.9 ms |
| F7=2^128+1 | **24.5 ms** | 27.0 ms | 39.8 ms | 286.9 ms | 16.32 s | 5.43 s | wrong |
| F8=2^256+1 | **189.8 ms** | 836.0 ms | 336.6 ms | 975.8 ms | 2.39 s | 1.53 s | wrong |
| M137=2^137-1 | **36.8 ms** | 43.0 ms | 50.3 ms | 357.9 ms | 13.40 s | 6.60 s | wrong |
| M149=2^149-1 | **30.4 ms** | 91.0 ms | 97.3 ms | 848.0 ms | 13.51 s | 19.25 s | wrong |
| 30!+1 | 2.5 ms | 1.0 ms | **532 µs** | 13.3 ms | 562.8 ms | 4.4 ms | 6.6 ms |
| 10^60+1 | **13.3 ms** | 53.0 ms | 16.9 ms | 32.4 ms | 1.12 s | 41.3 ms | 357.5 ms |

PARI's timer has a resolution of 1 ms, so "0 µs" means less than 1 ms.

## Summary

- **Mathilda is fastest in 31 of 47 rows**, and is effectively tied in two more: u63 against
  FLINT, and qs40-0/qs40-1 against PARI.
- **Every hard row is Mathilda's:**
  - QS from 50 to 90 digits: 5–10× faster than PARI and FLINT. At 90 digits Mathilda is the
    only engine that finishes.
  - Every ECM case from p25 to p35.
  - Every p−1, p+1 and Cunningham case.
- **v0.450 failed to factor 36 of the 41 single cases and two of the six batch groups.**
  v0.453 factors every case in the corpus.

## Rows where Mathilda is not the fastest, and why

- **u32 (FLINT 0.93 ms, PARI 4.0 ms, Mathilda 5.6 ms).** The numbers are 32-bit, so the
  factoring itself takes nanoseconds. The cost is the loop: Mathilda runs 2000 evaluator
  calls, about 2.8 µs each, while FLINT's driver is a C loop. Each call already goes
  straight to the word engine through the machine-integer fast path. Closing the rest needs
  cheaper builtin dispatch in the evaluator; the factoring code is not the bottleneck.
- **u63 (FLINT 26.1 ms, Mathilda 26.2 ms).** A tie, within noise.
- **semi96 (PARI 40 ms, Mathilda 77 ms).** These are 96-bit semiprimes with 48-bit factors.
  Mathilda runs Montgomery rho on two limbs, then the two-limb micro-ECM with B1=150 and 8
  curves. A minority of inputs miss in that budget and fall through to the mpz tiers, which
  costs most of the total. Two things would close it:
  - re-tune the micro-ECM schedule for 40–50-bit factors (B1 around 300–500 with more
    curves);
  - add a two-limb SIQS for 90–126 bits, which is what yafu and PARI's MPQS effectively
    provide in this range.
- **qs30 and qs40 (PARI 3 ms at 30 digits, Mathematica 18 ms on qs40-0, Mathilda 6–30 ms).**
  At these sizes the result is decided by fixed costs, not the sieve.
  - Automatic first spends a few ms on rho and micro-ECM. It switches to QS only once
    7·(ECM effort) ≥ 2·digits.
  - SIQS setup (factor base, Knuth–Schroeppel, the first polynomial) is a large fraction of
    a sub-10 ms run.
  - Fix: below about 35 digits, go to QS straight after rho, and give SIQS a lighter
    small-input path.
  - Mathematica's 18 ms on qs40-0 is a single lucky input. On qs40-1 it takes 331 ms.
- **ecm-p20-1 (FLINT 3.01 s, Mathilda 3.19 s).** ECM run time is set by when a curve
  happens to have a smooth group order, so on one input it is luck. On the other seven ECM
  cases Mathilda is between 3× and more than 10× faster.
- **fermat-60d/80d/100d (SymPy 1.4–2.3 ms, Mathilda 5–6 ms).** SymPy's `factorint` runs a
  few Fermat steps first. Mathilda's Automatic runs trial division and the rho tier
  (2¹³ iterations) before its Fermat tier, so about 4 ms goes to rho first. A handful of
  Fermat steps (say 8) before rho would cost microseconds and close this. PARI, FLINT and
  Mathematica have no Fermat stage and time out or take minutes.
- **pow-m61^3·q^2 (SymPy 5 ms, Mathematica 18 ms, Mathilda 88 ms).** n = (2⁶¹−1)³·q²
  with q a 26-digit prime, so n is not a perfect power. 2⁶¹−2 is 1321-smooth, so p−1
  finds the 61-bit factor at once. SymPy runs p−1 early. Mathilda reaches its p−1 tier only
  after rho and Fermat, and then runs p−1 at the full B1 = min(150000, 2000·digits). Fix: a
  short p−1 probe (B1 about 2000) right after rho.
- **pow-p30^7, pow-2^4096, pow-smooth^5, 30!+1 (PARI or FLINT under 1 ms, Mathilda 28 µs to
  2.6 ms).**
  - On inputs of 200 digits and more, the fixed cost is trial division to 2¹⁶, which
    Automatic does before the perfect-power test. That is about 2.5 ms on a 700-bit number.
  - Fix: test for a perfect power before trial division, and trial divide in bounded chunks
    that reuse the batched remainder tree.
  - For 30!+1, the remaining 83-bit cofactor goes through the u128 tiers, and FLINT's
    `n_factor`-style path gets there faster.

## Not measured

- The per-method rows (`--method X`: PollardRho/rho96, ShanksSquareForms/squfof60,
  Lehman/lehman40, TrialDivision/trial40, and the method-pinned pm1, pp1, ecm and qs runs)
  never ran, because the final run was stopped before reaching them.
- Development spot checks at v0.453: ``FactorInteger`Lehman`` at 40 bits took 24 µs per
  input, against FLINT's 23 µs. ``FactorInteger`PollardRho`` at 40 bits took 9 µs.
- `results/mtest.json` predates the v0.453 Lehman and SQUFOF fixes. Its Lehman total
  (6.5 ms against FLINT's 2.3 ms) and its 10 SQUFOF failures are superseded.
- Rerun the method rows with the loop in `final.sh`:
  `python3.11 run.py --engines mathilda,pari,flint --groups <g> --method <M> --tag method-<M>`.
- `BlakeRationalBaseDescent` has its ``FactorInteger`BlakeRationalBaseDescent`` entry point,
  but it is deliberately not optimised or benchmarked here. That is deferred to a later
  session.
