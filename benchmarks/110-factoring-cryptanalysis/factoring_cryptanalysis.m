(* ==========================================================================
   Experiment 110 -- Factoring cryptanalysis (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  Three classical integer-factoring algorithms written in
   the language ("app:"), each on the family of semiprimes it is designed to
   break, then the native FactorInteger primitive ("prim:") on balanced
   semiprimes (both factors ~half the modulus -- the honest, hard case).

   app: POLLARD RHO (Brent's variant) on semiprimes with one small (~28-bit)
        factor; the cycle-finding loop is deterministic (start 2, c = 1) so
        all systems find the same factor.
   app: POLLARD p-1 (stage 1) on semiprimes whose p-1 is B-smooth: a = 2 is
        raised to pr^e for each prime pr <= B1 (exponents from a precomputed
        prime list), then gcd(a - 1, n) reveals p.
   app: FERMAT on semiprimes with |p - q| tiny: a from ceil(sqrt n) upward
        until a^2 - n is a perfect square.
   prim: FactorInteger / factorint / factor on balanced semiprimes at
        80/100/120/160 bits.

   Every check is a SINGLE INTEGER: the sum of the recovered smaller factors
   (or, for FactorInteger, the smaller prime of each balanced semiprime),
   folded mod $CheckMod, plus a "recovered (count)" guard so a silent miss
   cannot hide behind a matching sum.  Inputs come from crypto_data.m's LCG.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"FactorInteger", "PowerMod", "GCD", "NextPrime", "PrimeQ", "Prime", "PrimePi"}];

smaller[n_, f_] := Min[f, n/f];

(* ---- app: Pollard rho (Brent) ------------------------------------------ *)

rhoBrent[n_] := Module[{y = 2, c = 1, m = 128, g = 1, r = 1, q = 1, x, ys, k},
  If[EvenQ[n], Return[2]];
  While[g == 1,
    x = y; Do[y = Mod[y^2 + c, n], {r}];
    k = 0;
    While[k < r && g == 1,
      ys = y;
      Do[y = Mod[y^2 + c, n]; q = Mod[q Abs[x - y], n], {Min[m, r - k]}];
      g = GCD[q, n]; k += m];
    r *= 2];
  If[g == n, g = 1;
    While[g == 1, ys = Mod[ys^2 + c, n]; g = GCD[Abs[x - ys], n]]];
  g];

(* semiprime with a `small`-bit factor and `bits` total. *)
weakRho[s_, small_, bits_] :=
  NextPrime[bigOf[s, small]] NextPrime[bigOf[s + 1, bits - small]];

rhoNs = {weakRho[2501, 28, 64], weakRho[2511, 28, 80],
         weakRho[2521, 28, 96], weakRho[2531, 28, 112]};

bench["app: Pollard rho (Brent) 28-bit factor x4",
      rhoF = rhoBrent /@ rhoNs;];
check["app: Pollard rho (Brent) 28-bit factor x4",
      Mod[Total[MapThread[smaller, {rhoNs, rhoF}]], $CheckMod] + $CheckMod * (Count[MapThread[Mod[#1, #2] == 0 &, {rhoNs, rhoF}], True])];

(* ---- app: Pollard p-1 (stage 1) ---------------------------------------- *)

$Primes = Prime[Range[PrimePi[20000]]];

pm1[n_, B1_] := Module[{a = 2, e, pr, i, g},
  Do[pr = $Primes[[i]]; If[pr > B1, Break[]];
     e = Floor[Log[B1]/Log[pr]];
     a = PowerMod[a, pr^e, n], {i, Length[$Primes]}];
  GCD[a - 1, n]];

(* p with p-1 B-smooth: p = 2*(prod primes <= b)*t + 1 prime. *)
smoothPrime[b_, t0_] := Module[{base = 2 Apply[Times, Prime[Range[PrimePi[b]]]], t = t0, p},
  p = base t + 1;
  While[! PrimeQ[p], t++; p = base t + 1];
  p];

pm1Keys = {{smoothPrime[40, 1] NextPrime[bigOf[2601, 45]], 60},
           {smoothPrime[50, 3] NextPrime[bigOf[2602, 55]], 80},
           {smoothPrime[60, 2] NextPrime[bigOf[2603, 60]], 100},
           {smoothPrime[70, 5] NextPrime[bigOf[2604, 64]], 120}};

bench["app: Pollard p-1 stage 1 x120",
      Do[pm1F = pm1[#[[1]], #[[2]]] & /@ pm1Keys, {30}];];
check["app: Pollard p-1 stage 1 x120",
      Mod[Total[MapThread[smaller, {pm1Keys[[All, 1]], pm1F}]], $CheckMod] + $CheckMod * (Count[MapThread[1 < #2 < #1 && Mod[#1, #2] == 0 &, {pm1Keys[[All, 1]], pm1F}], True])];

(* ---- app: Fermat (close primes) ---------------------------------------- *)

fermat[n_] := Module[{a = Ceiling[Sqrt[n]], b2},
  b2 = a^2 - n;
  While[! IntegerQ[Sqrt[b2]], a++; b2 = a^2 - n];
  a - Sqrt[b2]];

closeSemiprime[s_, bits_, gap_] := Module[{p = NextPrime[bigOf[s, bits]]},
  p NextPrime[p + gap]];

fermatNs = {closeSemiprime[2701, 40, 2^25], closeSemiprime[2702, 40, 2^25],
            closeSemiprime[2703, 40, 2^26], closeSemiprime[2704, 40, 2^24]};

bench["app: Fermat (close primes) x4",
      ferF = fermat /@ fermatNs;];
check["app: Fermat (close primes) x4",
      Mod[Total[MapThread[smaller, {fermatNs, ferF}]], $CheckMod] + $CheckMod * (Count[MapThread[1 < #2 < #1 && Mod[#1, #2] == 0 &, {fermatNs, ferF}], True])];

(* ---- prim: FactorInteger on balanced semiprimes ------------------------ *)

bal[s_, bits_] := NextPrime[bigOf[s, bits/2]] NextPrime[bigOf[s + 1, bits/2]];
smallestFactor[n_] := FactorInteger[n][[1, 1]];

bal80 = Table[bal[2800 + 2 k, 80], {k, 0, 24}];     (* 25 semiprimes *)
benchOnce["prim: FactorInteger balanced 80-bit x25",
          sf80 = smallestFactor /@ bal80;];
check["prim: FactorInteger balanced 80-bit x25", Mod[Total[sf80], $CheckMod]];

bal100 = Table[bal[2850 + 2 k, 100], {k, 0, 9}];    (* 10 semiprimes *)
benchOnce["prim: FactorInteger balanced 100-bit x10",
          sf100 = smallestFactor /@ bal100;];
check["prim: FactorInteger balanced 100-bit x10", Mod[Total[sf100], $CheckMod]];

bal120 = Table[bal[2880 + 2 k, 120], {k, 0, 2}];    (* 3 semiprimes *)
benchOnce["prim: FactorInteger balanced 120-bit x3",
          sf120 = smallestFactor /@ bal120;];
check["prim: FactorInteger balanced 120-bit x3", Mod[Total[sf120], $CheckMod]];

n160 = bal[2900, 160];
benchOnce["prim: FactorInteger balanced 160-bit x1",
          sf160 = smallestFactor[n160];];
check["prim: FactorInteger balanced 160-bit x1", Mod[sf160, $CheckMod]];
