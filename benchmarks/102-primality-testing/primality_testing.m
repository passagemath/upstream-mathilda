(* ==========================================================================
   Experiment 102 -- primality testing (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  Probable-prime testing at cryptographic sizes, as each
   system's builtin ("prim:") and as the two textbook tests written in the
   language ("app:"): Miller-Rabin with the first 12 prime bases, and
   Baillie-PSW (a strong base-2 Miller-Rabin round, then a strong Lucas test
   with Selfridge's method-A parameters D = 5, -7, 9, ..., P = 1,
   Q = (1 - D)/4, built on JacobiSymbol).

   Inputs, all deterministic:
     - 30 primes, 10 each of 512 / 1024 / 2048 bits (primeOf, shared LCG);
     - 2000 random odd 1024-bit integers (the early-reject path);
     - 47 composites built to fool weak tests: 8 Chernick Carmichael numbers
       (6k+1)(12k+1)(18k+1) found by searching k upward from an LCG start;
       7 literal strong pseudoprimes to the first 4..13 prime bases (verified:
       each is composite and fools exactly the stated bases); the 27
       composite Mersenne numbers 2^p - 1, 1000 < p < 1200, and the Fermat
       numbers F8..F11 -- every composite Mersenne / Fermat number is a
       strong pseudoprime to base 2, so they exercise BPSW's Lucas half;
     - the Mersenne prime 2^4423 - 1.

   Every check is a COUNT of "prime" verdicts, so the systems agree only if
   every individual verdict agrees in total.  Expected: the 12-base
   Miller-Rabin is fooled by exactly two inputs (the spsp's to bases 2..37
   and 2..41) -- that is the algorithm being honest, not a bug -- and BPSW
   and every builtin by none.

   The proven-primality row uses ProvablePrimeQ (PrimalityProving` in WL);
   Mathilda lacks it today, so that case is benchIf'd (ABSENT).
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

Quiet[Needs["PrimalityProving`"]];

require[{"PrimeQ", "PowerMod", "JacobiSymbol", "IntegerExponent", "ProvablePrimeQ"}];

(* ---- inputs (untimed) ------------------------------------------------------ *)

primes = Join[Table[primeOf[s, 512], {s, 401, 410}],
              Table[primeOf[s, 1024], {s, 411, 420}],
              Table[primeOf[s, 2048], {s, 421, 430}]];
odds1024 = Table[bigOf[s, 1024] + 1 - Mod[bigOf[s, 1024], 2], {s, 5001, 7000}];
odds500 = odds1024[[;; 500]];

(* First k >= start with 6k+1, 12k+1, 18k+1 all prime; their product. *)
chernick[start_] := Module[{k = start},
  While[! (PrimeQ[6 k + 1] && PrimeQ[12 k + 1] && PrimeQ[18 k + 1]), k++];
  (6 k + 1) (12 k + 1) (18 k + 1)];

carmichaels = Table[chernick[bigOf[s, 64]], {s, 3001, 3008}];
spsps = {3215031751, 2152302898747, 3474749660383, 341550071728321,
         3825123056546413051, 318665857834031151167461,
         3317044064679887385961981};
mersenneC = Table[2^p - 1, {p, Select[Range[1000, 1199], PrimeQ]}];
fermatC = Table[2^(2^k) + 1, {k, 8, 11}];
pseudos = Join[carmichaels, spsps, mersenneC, fermatC];

(* ---- app: the tests written in the language ------------------------------- *)

$MRBases = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

(* One strong-probable-prime round to base a, n - 1 = d 2^s. *)
mrRound[n_, a_, d_, s_] := Module[{x = PowerMod[a, d, n], r = 1},
  If[x == 1 || x == n - 1, Return[True]];
  While[r < s,
    x = Mod[x^2, n];
    If[x == n - 1, Return[True]];
    r++];
  False];

millerRabin[n_] := Module[{s, d},
  Which[n < 2, False, EvenQ[n], n == 2,
    True,
    s = IntegerExponent[n - 1, 2]; d = (n - 1)/2^s;
    AllTrue[$MRBases, Mod[#, n] == 0 || mrRound[n, #, d, s] &]]];

half[x_, n_] := If[OddQ[x], (x + n)/2, x/2];

(* Strong Lucas probable-prime test, Selfridge method A. *)
strongLucas[n_] := Module[{dd = 5, j, p = 1, q, d, s, bits, u, v, qk, r},
  While[True,
    j = JacobiSymbol[dd, n];
    If[j == -1, Break[]];
    If[j == 0 && Abs[dd] != n, Return[False]];
    dd = If[dd > 0, -dd - 2, -dd + 2]];
  q = (1 - dd)/4;
  s = IntegerExponent[n + 1, 2]; d = (n + 1)/2^s;
  bits = Rest[IntegerDigits[d, 2]];
  u = 1; v = p; qk = Mod[q, n];
  Do[
    {u, v} = {Mod[u v, n], Mod[v^2 - 2 qk, n]};
    qk = Mod[qk^2, n];
    If[b == 1,
      {u, v} = {half[Mod[p u + v, n], n], half[Mod[dd u + p v, n], n]};
      qk = Mod[qk q, n]],
    {b, bits}];
  If[u == 0 || v == 0, Return[True]];
  r = 1;
  While[r < s,
    v = Mod[v^2 - 2 qk, n];
    If[v == 0, Return[True]];
    qk = Mod[qk^2, n];
    r++];
  False];

bpsw[n_] := Module[{s, d},
  Which[n < 2, False, EvenQ[n], n == 2,
    True,
    s = IntegerExponent[n - 1, 2]; d = (n - 1)/2^s;
    Which[! mrRound[n, 2, d, s], False,
          IntegerQ[Sqrt[n]], False,
          True, strongLucas[n]]]];

(* The same two tests with every early exit expressed as a loop condition
   instead of Return[] -- see README "Findings": Mathilda's Return[] inside a
   While returns from the While, not from the enclosing function, so the
   natural spellings above answer wrongly there.  Same algorithm. *)
mrRound2[n_, a_, d_, s_] := Module[{x = PowerMod[a, d, n], r = 1},
  If[x == 1 || x == n - 1, True,
    While[r < s && x != n - 1, x = Mod[x^2, n]; r++];
    x == n - 1]];

millerRabin2[n_] := Module[{s, d},
  Which[n < 2, False, EvenQ[n], n == 2,
    True,
    s = IntegerExponent[n - 1, 2]; d = (n - 1)/2^s;
    AllTrue[$MRBases, Mod[#, n] == 0 || mrRound2[n, #, d, s] &]]];

strongLucas2[n_] := Module[{dd = 5, p = 1, q, d, s, bits, u, v, qk, r = 1},
  While[JacobiSymbol[dd, n] == 1, dd = If[dd > 0, -dd - 2, -dd + 2]];
  If[JacobiSymbol[dd, n] == 0 && Abs[dd] != n, False,
    q = (1 - dd)/4;
    s = IntegerExponent[n + 1, 2]; d = (n + 1)/2^s;
    bits = Rest[IntegerDigits[d, 2]];
    u = 1; v = p; qk = Mod[q, n];
    Do[
      {u, v} = {Mod[u v, n], Mod[v^2 - 2 qk, n]};
      qk = Mod[qk^2, n];
      If[b == 1,
        {u, v} = {half[Mod[p u + v, n], n], half[Mod[dd u + p v, n], n]};
        qk = Mod[qk q, n]],
      {b, bits}];
    If[u == 0 || v == 0, True,
      While[r < s && v != 0,
        v = Mod[v^2 - 2 qk, n];
        qk = Mod[qk^2, n];
        r++];
      v == 0]]];

bpsw2[n_] := Module[{s, d},
  Which[n < 2, False, EvenQ[n], n == 2,
    True,
    s = IntegerExponent[n - 1, 2]; d = (n - 1)/2^s;
    Which[! mrRound2[n, 2, d, s], False,
          IntegerQ[Sqrt[n]], False,
          True, strongLucas2[n]]]];

(* ---- prim -------------------------------------------------------------------- *)

bench["prim: PrimeQ 512/1024/2048-bit primes x30", c1 = Count[PrimeQ /@ primes, True];];
check["prim: PrimeQ 512/1024/2048-bit primes x30", c1];

bench["prim: PrimeQ 1024-bit odd x2000", c2 = Count[PrimeQ /@ odds1024, True];];
check["prim: PrimeQ 1024-bit odd x2000", c2];

bench["prim: PrimeQ pseudoprimes (Carmichael/spsp/Mersenne/Fermat) x47",
      c3 = Count[PrimeQ /@ pseudos, True];];
check["prim: PrimeQ pseudoprimes (Carmichael/spsp/Mersenne/Fermat) x47", c3];

bench["prim: PrimeQ Mersenne 2^4423-1", c4 = Boole[PrimeQ[2^4423 - 1]];];
check["prim: PrimeQ Mersenne 2^4423-1", c4];

proofPrimes = Table[primeOf[s, 256], {s, 431, 434}];
If[haveQ["ProvablePrimeQ"],
   benchOnce["prim: proven primality 256-bit prime x4",
             c5 = Count[ProvablePrimeQ /@ proofPrimes, True];];
   check["prim: proven primality 256-bit prime x4", c5],
   Print["SKIP\tprim: proven primality 256-bit prime x4\tProvablePrimeQ"]];

(* ---- app ---------------------------------------------------------------------- *)

bench["app: Miller-Rabin 12 bases primes x30", a1 = Count[millerRabin /@ primes, True];];
check["app: Miller-Rabin 12 bases primes x30", a1];

bench["app: Miller-Rabin 12 bases pseudoprimes x47", a2 = Count[millerRabin /@ pseudos, True];];
check["app: Miller-Rabin 12 bases pseudoprimes x47", a2];

bench["app: Baillie-PSW primes x30", a3 = Count[bpsw /@ primes, True];];
check["app: Baillie-PSW primes x30", a3];

bench["app: Baillie-PSW pseudoprimes x47", a4 = Count[bpsw /@ pseudos, True];];
check["app: Baillie-PSW pseudoprimes x47", a4];

bench["app: Baillie-PSW 1024-bit odd x500", a5 = Count[bpsw /@ odds500, True];];
check["app: Baillie-PSW 1024-bit odd x500", a5];

bench["app: Miller-Rabin 12 bases (no-Return) primes x30", b1 = Count[millerRabin2 /@ primes, True];];
check["app: Miller-Rabin 12 bases (no-Return) primes x30", b1];

bench["app: Miller-Rabin 12 bases (no-Return) pseudoprimes x47", b2 = Count[millerRabin2 /@ pseudos, True];];
check["app: Miller-Rabin 12 bases (no-Return) pseudoprimes x47", b2];

bench["app: Baillie-PSW (no-Return) primes x30", b3 = Count[bpsw2 /@ primes, True];];
check["app: Baillie-PSW (no-Return) primes x30", b3];

bench["app: Baillie-PSW (no-Return) 1024-bit odd x500", b4 = Count[bpsw2 /@ odds500, True];];
check["app: Baillie-PSW (no-Return) 1024-bit odd x500", b4];
