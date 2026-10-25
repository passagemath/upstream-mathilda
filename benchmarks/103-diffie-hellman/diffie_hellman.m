(* ==========================================================================
   Experiment 103 -- finite-field Diffie-Hellman (group F)
   ==========================================================================
   WHAT IT MEASURES.  The number theory behind classic (mod-p) Diffie-Hellman:
     - safe-prime generation: the first q >= an LCG start with q and 2q + 1
       both prime, written in the language around the system's own PrimeQ;
     - a primitive root of a safe prime p = 2q + 1, as the builtin
       (PrimitiveRoot) and as the textbook test (the least g >= 2 with
       g^2 != 1 and g^q != 1 mod p, since p - 1 = 2q has only those two
       prime factors);
     - the RFC 3526 MODP-2048 exchange (g = 2): A = g^a, B = g^b, then each
       side's shared secret B^a, A^b, with 256-bit private exponents;
     - public-key validation: a peer key y must satisfy 1 < y < p - 1 and lie
       in the prime-order-q subgroup, i.e. be a quadratic residue -- tested as
       y^q == 1 mod p (Euler's criterion via PowerMod, the "app" form) and as
       JacobiSymbol[y, p] == 1 (the "prim" form).  Half the candidate keys are
       honest 2^a (2 is a QR since p == 7 mod 8), half arbitrary LCG integers.

   All inputs from crypto_data.m; safe primes for the PrimitiveRoot cases are
   found untimed.  Checks: sums of primes / roots / secrets mod $CheckMod,
   and counts.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PrimeQ", "PrimitiveRoot", "PowerMod", "JacobiSymbol"}];

p2048 = $MODP2048;
q2048 = (p2048 - 1)/2;

(* ---- app: safe-prime generation ------------------------------------------ *)

(* p = 2q + 1 for the first odd q >= start with q and 2q + 1 both prime. *)
safePrime[start_] := Module[{q = start + 1 - Mod[start, 2]},
  While[! (PrimeQ[q] && PrimeQ[2 q + 1]), q += 2];
  2 q + 1];

bench["app: safe-prime search 256-bit x4",
      sp256 = Table[safePrime[bigOf[s, 255]], {s, 1, 4}];];
check["app: safe-prime search 256-bit x4", Mod[Total[sp256], $CheckMod]];

benchOnce["app: safe-prime search 512-bit x2",
          sp512 = Table[safePrime[bigOf[s, 511]], {s, 9, 10}];];
check["app: safe-prime search 512-bit x2", Mod[Total[sp512], $CheckMod]];

(* ---- primitive roots of safe primes -------------------------------------- *)

safes = Table[safePrime[bigOf[s, 255]], {s, 1, 8}];   (* untimed *)

bench["prim: PrimitiveRoot 256-bit safe prime x800",
      Do[roots = PrimitiveRoot /@ safes, {100}];];
check["prim: PrimitiveRoot 256-bit safe prime x800", Total[roots]];

(* least g >= 2 with g^2 != 1 and g^q != 1 (mod p), p = 2q + 1 *)
safeRoot[p_] := Module[{q = (p - 1)/2, g = 2},
  While[Mod[g^2, p] == 1 || PowerMod[g, q, p] == 1, g++];
  g];

bench["app: primitive root of safe prime 256-bit x800",
      Do[roots2 = safeRoot /@ safes, {100}];];
check["app: primitive root of safe prime 256-bit x800", Total[roots2]];

(* ---- MODP-2048 key exchange ---------------------------------------------- *)

privA = Table[bigOf[s, 256], {s, 6001, 6050}];
privB = Table[bigOf[s, 256], {s, 7001, 7050}];

bench["prim: PowerMod 2^a mod MODP-2048 (256-bit a) x100",
      pubs = PowerMod[2, #, p2048] & /@ Join[privA, privB];];
check["prim: PowerMod 2^a mod MODP-2048 (256-bit a) x100", Mod[Total[pubs], $CheckMod]];

dhExchange[a_, b_] := Module[{pa, pb, sa, sb},
  pa = PowerMod[2, a, p2048];
  pb = PowerMod[2, b, p2048];
  sa = PowerMod[pb, a, p2048];
  sb = PowerMod[pa, b, p2048];
  If[sa == sb, sa, 0]];

bench["app: MODP-2048 DH exchange x50",
      secrets = MapThread[dhExchange, {privA, privB}];];
check["app: MODP-2048 DH exchange x50", Mod[Total[secrets], $CheckMod]];

(* ---- peer-key validation ------------------------------------------------- *)

cands = Join[PowerMod[2, privA, p2048],
             Table[Mod[bigOf[s, 2048], p2048], {s, 8001, 8050}]];

validQ[y_] := 1 < y < p2048 - 1 && PowerMod[y, q2048, p2048] == 1;

bench["app: subgroup validation y^q mod p MODP-2048 x100",
      nv = Count[validQ /@ cands, True];];
check["app: subgroup validation y^q mod p MODP-2048 x100", nv];

bench["prim: JacobiSymbol validation MODP-2048 x10000",
      Do[nj = Count[cands, y_ /; 1 < y < p2048 - 1 && JacobiSymbol[y, p2048] == 1], {100}];];
check["prim: JacobiSymbol validation MODP-2048 x10000", nj];
