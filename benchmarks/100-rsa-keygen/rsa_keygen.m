(* ==========================================================================
   Experiment 100 -- RSA key generation (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  The integer kernels an RSA key generator leans on --
   NextPrime / PrimeQ on 1024-2048-bit operands, the private-exponent inverse
   d = e^-1 mod lambda(n), and the CRT parameters -- first as each system's
   builtin ("prim:" cases), then as a whole keygen written in the language
   ("app:" cases: an incremental odd-candidate search driven by the system's
   own primality test, exactly as a textbook generator is written).

   Inputs come from crypto_data.m's shared LCG, so every system searches from
   the same starting integers and must land on the same primes; each check
   folds the result to a residue mod $CheckMod.

   ModularInverse is the WL spelling of the inverse; Mathilda lacks it today,
   so that case is benchIf'd (ABSENT).  PowerMod[e, -1, m] is the portable
   spelling and is timed alongside it.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"NextPrime", "PrimeQ", "PowerMod", "ModularInverse", "LCM", "GCD"}];

$E = 65537;

(* ---- prim: the builtins, on precomputed inputs ------------------------- *)

starts1024 = Table[bigOf[s, 1024], {s, 101, 104}];
bench["prim: NextPrime 1024-bit x4", nps = NextPrime /@ starts1024;];
check["prim: NextPrime 1024-bit x4", Mod[Total[nps], $CheckMod]];

primes2048 = Table[primeOf[s, 2048], {s, 201, 208}];
bench["prim: PrimeQ 2048-bit prime x8", cnt = Count[PrimeQ /@ primes2048, True];];
check["prim: PrimeQ 2048-bit prime x8", cnt];

odds2048 = Table[bigOf[s, 2048] + 1 - Mod[bigOf[s, 2048], 2], {s, 301, 700}];
bench["prim: PrimeQ 2048-bit odd x400", cnt = Count[PrimeQ /@ odds2048, True];];
check["prim: PrimeQ 2048-bit odd x400", cnt];

(* lambda values of the shape a keygen inverts against: lcm(p-1, q-1) of
   1024-bit primes is ~2047 bits; a product of LCG words coprime to e stands
   in: 2000 of them, inverted 10 times over, so the case clears the PARI
   timer's 1 ms resolution. *)
lams = Select[Table[bigOf[s, 2048], {s, 1001, 3100}], GCD[#, $E] == 1 &][[;; 2000]];
bench["prim: PowerMod[e,-1,lambda] 2048-bit x20000",
      Do[ds = PowerMod[$E, -1, #] & /@ lams, {10}];];
check["prim: PowerMod[e,-1,lambda] 2048-bit x20000", Mod[Total[ds], $CheckMod]];

benchIf["prim: ModularInverse 2048-bit x20000", "ModularInverse",
        Do[ds2 = ModularInverse[$E, #] & /@ lams, {10}];];
checkIf["prim: ModularInverse 2048-bit x20000", "ModularInverse",
        Mod[Total[ModularInverse[$E, #] & /@ lams], $CheckMod]];

(* ---- app: a whole key generator, written in the language --------------- *)

(* The first prime >= start (start forced odd) with gcd(p - 1, e) = 1. *)
genPrime[start_] := Module[{c = start + 1 - Mod[start, 2]},
  While[! (PrimeQ[c] && GCD[c - 1, $E] == 1), c += 2];
  c];

(* {n, d, dp, dq, qinv} from two LCG seeds. *)
rsaKey[bits_, s_] := Module[{p, q, lam, d},
  p = genPrime[bigOf[s, bits/2]];
  q = genPrime[bigOf[s + 1, bits/2]];
  lam = LCM[p - 1, q - 1];
  d = PowerMod[$E, -1, lam];
  {p q, d, Mod[d, p - 1], Mod[d, q - 1], PowerMod[q, -1, p]}];

keyDigest[keys_] := Mod[Total[Flatten[keys]], $CheckMod];

bench["app: RSA-1024 keygen x8", k1 = Table[rsaKey[1024, 2 s], {s, 1, 8}];];
check["app: RSA-1024 keygen x8", keyDigest[k1]];

benchOnce["app: RSA-2048 keygen x2", k2 = Table[rsaKey[2048, 2 s + 50], {s, 1, 2}];];
check["app: RSA-2048 keygen x2", keyDigest[k2]];
