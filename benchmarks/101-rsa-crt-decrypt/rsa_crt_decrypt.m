(* ==========================================================================
   Experiment 101 -- RSA-2048 encrypt / decrypt, plain vs CRT (group F)
   ==========================================================================
   WHAT IT MEASURES.  The modular-exponentiation workload of textbook RSA once
   a key exists: encryption with e = 65537 (17 squarings), decryption with the
   full 2047-bit private exponent d, and the CRT (Garner) decryption every real
   implementation uses instead -- two half-size exponentiations mod p and mod q
   recombined as m = m2 + q (qinv (m1 - m2) mod p), about 4x less work.

   "prim:" cases time each system's builtin: PowerMod (scalar, and over a
   whole List -- PowerMod is Listable in both CAS), and ChineseRemainder for
   the recombination step.  "app:" cases write the CRT decryption and the
   Garner recombination in the language from PowerMod / Mod.

   The key is built exactly as 100-rsa-keygen builds one (an incremental
   odd-candidate search from shared LCG starts), OUTSIDE every timed region;
   the messages are shared LCG integers below n.  Every decryption check folds
   the recovered plaintexts mod $CheckMod, so a wrong decryption cannot pass.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod", "ChineseRemainder", "PrimeQ", "LCM", "GCD"}];

$E = 65537;

(* ---- the key (untimed) --------------------------------------------------- *)

genPrime[start_] := Module[{c = start + 1 - Mod[start, 2]},
  While[! (PrimeQ[c] && GCD[c - 1, $E] == 1), c += 2];
  c];

rsaP = genPrime[bigOf[1101, 1024]];
rsaQ = genPrime[bigOf[1102, 1024]];
rsaN = rsaP rsaQ;
rsaD = PowerMod[$E, -1, LCM[rsaP - 1, rsaQ - 1]];
rsaDp = Mod[rsaD, rsaP - 1];
rsaDq = Mod[rsaD, rsaQ - 1];
rsaQinv = PowerMod[rsaQ, -1, rsaP];

(* messages: 2000-bit LCG integers, all < n.  The recombination cases loop
   the 4000 pairs 5 times so the PARI timer clears 20 ms. *)
msgs = Table[bigOf[s, 2000], {s, 2001, 6000}];
msgs200 = msgs[[;; 200]];
cts = PowerMod[msgs, $E, rsaN];
cts50 = cts[[;; 50]];
cts200 = cts[[;; 200]];

(* ---- prim: encryption ---------------------------------------------------- *)

bench["prim: PowerMod encrypt e=65537 2048-bit x4000",
      enc = PowerMod[#, $E, rsaN] & /@ msgs;];
check["prim: PowerMod encrypt e=65537 2048-bit x4000", Mod[Total[enc], $CheckMod]];

bench["prim: PowerMod[list] encrypt e=65537 2048-bit x4000",
      encL = PowerMod[msgs, $E, rsaN];];
check["prim: PowerMod[list] encrypt e=65537 2048-bit x4000", Mod[Total[encL], $CheckMod]];

(* ---- prim: plain decryption --------------------------------------------- *)

bench["prim: PowerMod decrypt c^d mod n 2048-bit x50",
      dec = PowerMod[#, rsaD, rsaN] & /@ cts50;];
check["prim: PowerMod decrypt c^d mod n 2048-bit x50", Mod[Total[dec], $CheckMod]];

(* ---- app: CRT (Garner) decryption --------------------------------------- *)

crtDecrypt[c_] := Module[{m1, m2, h},
  m1 = PowerMod[c, rsaDp, rsaP];
  m2 = PowerMod[c, rsaDq, rsaQ];
  h = Mod[rsaQinv (m1 - m2), rsaP];
  m2 + h rsaQ];

bench["app: RSA-2048 CRT decrypt x200", decC = crtDecrypt /@ cts200;];
check["app: RSA-2048 CRT decrypt x200", Mod[Total[decC], $CheckMod]];

(* ---- recombination only: the two half-results precomputed --------------- *)

m1s = PowerMod[cts, rsaDp, rsaP];
m2s = PowerMod[cts, rsaDq, rsaQ];

bench["prim: ChineseRemainder recombine 2048-bit x20000",
      Do[recC = MapThread[ChineseRemainder[{#1, #2}, {rsaP, rsaQ}] &, {m1s, m2s}], {5}];];
check["prim: ChineseRemainder recombine 2048-bit x20000", Mod[Total[recC], $CheckMod]];

bench["app: Garner recombine 2048-bit x20000",
      Do[recG = MapThread[#2 + Mod[rsaQinv (#1 - #2), rsaP] rsaQ &, {m1s, m2s}], {5}];];
check["app: Garner recombine 2048-bit x20000", Mod[Total[recG], $CheckMod]];

(* ---- app: full round trip ------------------------------------------------ *)

bench["app: RSA-2048 encrypt+CRT decrypt round trip x200",
      rt = Count[crtDecrypt[PowerMod[#, $E, rsaN]] - # & /@ msgs200, 0];];
check["app: RSA-2048 encrypt+CRT decrypt round trip x200", rt];
