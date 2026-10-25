(* ==========================================================================
   Experiment 118 -- quadratic-residuosity cryptography (group F)
   ==========================================================================
   WHAT IT MEASURES.  Two classic constructions whose security rests on the
   quadratic-residuosity problem modulo a Blum integer n = p q
   (p = q = 3 mod 4, 1024-bit n):

     * Blum-Blum-Shub, x -> x^2 mod n, written in the language ("app:"):
       one parity bit per step, and the variant that extracts
       log2 log2 n = 10 low bits per step;
     * Goldwasser-Micali: encrypt bit b as x^2 y^b mod n with the pseudo-
       square y = n - 1, decrypt by the Legendre symbol of c mod p --
       timed through the native symbol (prim: JacobiSymbol) and through
       Euler's criterion PowerMod[c, (p - 1)/2, p] written in the language;

   plus the raw symbol kernels on 1024-bit operands:  JacobiSymbol on odd
   moduli, and KroneckerSymbol on signed / even arguments.

   p and q are found by the same deterministic search in every system: from
   bigOf[seed, 512] step to the first candidate = 3 mod 4, then by 4 until
   PrimeQ.  Every other input comes from crypto_data.m's LCG; checks fold
   bit streams with FromDigits[..., 2] (or a running 2^j-ary fold) mod
   $CheckMod.

   KroneckerSymbol is absent in Mathilda today, so that case is benchIf'd
   (ABSENT); JacobiSymbol on the same signed/even inputs is not a substitute
   (Jacobi is undefined for even moduli), so no portable twin is timed.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod", "JacobiSymbol", "KroneckerSymbol", "PrimeQ", "FromDigits"}];

(* the first prime = 3 mod 4 at or above start *)
blumPrime[start_] := Module[{c = start - Mod[start, 4] + 3},
  While[! PrimeQ[c], c += 4];
  c];

bbsP = blumPrime[bigOf[1181, 512]];
bbsQ = blumPrime[bigOf[1182, 512]];
bbsN = bbsP bbsQ;
bbsX0 = Mod[bigOf[1183, 1000]^2, bbsN];         (* a quadratic residue seed *)

(* ---- app: Blum-Blum-Shub ----------------------------------------------- *)

(* parity of each successive square, folded as a binary number mod $CheckMod *)
bbsParity[x0_, n_, steps_] := Module[{x = x0, h = 0},
  Do[x = Mod[x^2, n]; h = Mod[2 h + Mod[x, 2], $CheckMod], {steps}];
  h];

(* the low 10 = log2 log2 n bits of each square *)
bbsLow10[x0_, n_, steps_] := Module[{x = x0, h = 0},
  Do[x = Mod[x^2, n]; h = Mod[1024 h + Mod[x, 1024], $CheckMod], {steps}];
  h];

bench["app: BBS-1024 parity x20000 bits", bbs1 = bbsParity[bbsX0, bbsN, 20000];];
check["app: BBS-1024 parity x20000 bits", bbs1];

bench["app: BBS-1024 10 bits/step x20000 steps", bbs2 = bbsLow10[bbsX0, bbsN, 20000];];
check["app: BBS-1024 10 bits/step x20000 steps", bbs2];

(* ---- app / prim: Goldwasser-Micali ------------------------------------- *)

gmY = bbsN - 1;                                  (* Jacobi +1, non-residue mod p, q *)
gmMsg = smallList[1184, 2000, 2];
gmR = Table[bigOf[2000 + i, 1000], {i, 2000}];  (* encryption randomness *)

bench["app: GM-1024 encrypt x2000 bits x10",
      Do[gmC = MapThread[Mod[#1^2 gmY^#2, bbsN] &, {gmR, gmMsg}], {10}];];
check["app: GM-1024 encrypt x2000 bits x10", Mod[Total[gmC], $CheckMod]];

(* decryption: b = 0 iff c is a square mod p *)
bench["prim: JacobiSymbol GM-1024 decrypt x2000 bits x4",
      Do[gmD1 = (1 - JacobiSymbol[#, bbsP])/2 & /@ gmC, {4}];];
check["prim: JacobiSymbol GM-1024 decrypt x2000 bits x4",
      Mod[FromDigits[gmD1, 2], $CheckMod] + If[gmD1 === gmMsg, 10^10, 0]];

gmE = (bbsP - 1)/2;
bench["app: GM-1024 decrypt Euler criterion x2000 bits",
      gmD2 = If[PowerMod[#, gmE, bbsP] == 1, 0, 1] & /@ gmC;];
check["app: GM-1024 decrypt Euler criterion x2000 bits",
      Mod[FromDigits[gmD2, 2], $CheckMod] + If[gmD2 === gmMsg, 10^10, 0]];

(* ---- prim: the symbol kernels on 1024-bit operands --------------------- *)

jacA = Table[bigOf[3000 + i, 1024], {i, 5000}];
jacM = Table[With[{b = bigOf[9000 + i, 1024]}, b + 1 - Mod[b, 2]], {i, 5000}];
bench["prim: JacobiSymbol 1024-bit x5000x4",
      Do[jac = MapThread[JacobiSymbol, {jacA, jacM}], {4}];];
check["prim: JacobiSymbol 1024-bit x5000x4", Mod[Total[jac], $CheckMod]];

(* signed numerators, even (and odd) denominators: the Kronecker extension *)
kroA = Table[(-1)^i jacA[[i]], {i, 5000}];
kroM = Table[jacM[[i]] - Mod[i, 2], {i, 5000}];   (* odd i -> even modulus *)
benchIf["prim: KroneckerSymbol signed/even 1024-bit x5000x4", "KroneckerSymbol",
        Do[kro = MapThread[KroneckerSymbol, {kroA, kroM}], {4}];];
checkIf["prim: KroneckerSymbol signed/even 1024-bit x5000x4", "KroneckerSymbol",
        Mod[Total[kro], $CheckMod]];
