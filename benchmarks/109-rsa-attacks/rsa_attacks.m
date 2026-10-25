(* ==========================================================================
   Experiment 109 -- RSA attacks (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  Three textbook attacks on misused RSA, each written in
   the language ("app:") over primitives that crypto number theory leans on,
   plus the native continued-fraction primitives ("prim:") behind Wiener.

   (a) WIENER.  A key with a small private exponent d < N^(1/4)/3 (here d is
       120-bit, N is 1024-bit, bound ~254-bit) leaks d: d is the denominator
       of a convergent of the continued fraction of e/N.  The recovery is an
       in-language convergent recurrence on ContinuedFraction[e/N], testing
       each convergent (k_num, d) by reconstructing phi and checking that
       x^2 - (N - phi + 1) x + N has integer roots (a perfect-square
       discriminant).  "prim: ContinuedFraction" times the builtin that
       drives it; "prim: Convergents" is the WL head that would yield the
       convergents directly -- ABSENT in Mathilda (benchIf -> SKIP).

   (b) COMMON MODULUS.  Two coprime public exponents on the same modulus:
       given c1 = m^e1, c2 = m^e2, ExtendedGCD[e1, e2] = {1, {a, b}} recovers
       m = c1^a c2^b mod N (negative powers via modular inverse).

   (c) HASTAD BROADCAST.  e = 3, one message m sent under three coprime
       moduli: ChineseRemainder lifts {c_i} to m^3 mod (n1 n2 n3) = m^3
       exactly (m^3 < the product), then an in-language integer cube root
       (Newton on integers) recovers m.

   Inputs come from crypto_data.m's shared LCG; each check folds a recovered
   secret (or a count) to a residue mod $CheckMod.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"ContinuedFraction", "Convergents", "ExtendedGCD", "ChineseRemainder",
         "PowerMod", "NextPrime"}];

(* ---- (a) Wiener -------------------------------------------------------- *)

(* A weak RSA key {e, N, d} with d 120-bit (<< N^(1/4)/3), e = d^-1 mod phi. *)
wienerKey[s_] := Module[{p, q, n, phi, d, e},
  p = NextPrime[bigOf[s, 512]]; q = NextPrime[bigOf[s + 1, 512]];
  n = p q; phi = (p - 1) (q - 1);
  d = bigOf[s + 2, 120]; If[EvenQ[d], d++]; While[GCD[d, phi] != 1, d += 2];
  e = PowerMod[d, -1, phi]; {e, n, d}];

(* Recover d from the convergents of e/N (in-language recurrence). *)
wiener[e_, n_] := Module[
  {cf, hm1 = 1, hm2 = 0, km1 = 0, km2 = 1, h, k, i, phi, b, disc, found = 0},
  cf = ContinuedFraction[e/n];
  Do[
    h = cf[[i]] hm1 + hm2; k = cf[[i]] km1 + km2;   (* convergent h/k *)
    If[h > 0 && k > 0 && Mod[e k - 1, h] == 0,
      phi = (e k - 1)/h; b = n - phi + 1; disc = b^2 - 4 n;
      If[disc >= 0 && IntegerQ[Sqrt[disc]], found = k]];
    hm2 = hm1; hm1 = h; km2 = km1; km1 = k,
    {i, Length[cf]}];
  found];

wkeys = Table[wienerKey[1700 + 3 s], {s, 0, 19}];   (* 20 weak keys *)

bench["prim: ContinuedFraction e/N 1024-bit x600",
      Do[cfs = ContinuedFraction[#[[1]]/#[[2]]] & /@ wkeys, {30}];];
check["prim: ContinuedFraction e/N 1024-bit x600",
      Mod[Total[Total /@ cfs] + Total[Length /@ cfs], $CheckMod]];

(* Convergents returns the convergent list directly; from it the same disc
   test recovers d.  ABSENT in Mathilda (SKIP); Mathematica has Convergents. *)
wienerFromConv[cv_, e_, n_] := Module[{found = 0, h, k, phi, b, disc},
  Do[h = Numerator[c]; k = Denominator[c];
     If[h > 0 && k > 0 && Mod[e k - 1, h] == 0,
        phi = (e k - 1)/h; b = n - phi + 1; disc = b^2 - 4 n;
        If[disc >= 0 && IntegerQ[Sqrt[disc]], found = k]], {c, cv}];
  found];
benchIf["prim: Convergents e/N 1024-bit x40", "Convergents",
        Do[cvrec = wienerFromConv[Convergents[#[[1]]/#[[2]]], #[[1]], #[[2]]] & /@ wkeys, {2}];];
checkIf["prim: Convergents e/N 1024-bit x40", "Convergents",
        Mod[Total[cvrec], $CheckMod]];

bench["app: Wiener attack 1024-bit x20",
      wrec = wiener[#[[1]], #[[2]]] & /@ wkeys;];
check["app: Wiener attack 1024-bit x20", Mod[Total[wrec], $CheckMod] + $CheckMod * (Count[MapThread[Equal, {wrec, wkeys[[All, 3]]}], True])];

(* ---- (b) common modulus ------------------------------------------------ *)

cmN = NextPrime[bigOf[1801, 512]] NextPrime[bigOf[1802, 512]];
$E1 = 65537; $E2 = 65539;

commonModulus[m_] := Module[{c1, c2, g, a, b},
  c1 = PowerMod[m, $E1, cmN]; c2 = PowerMod[m, $E2, cmN];
  {g, {a, b}} = ExtendedGCD[$E1, $E2];
  Mod[PowerMod[c1, a, cmN] PowerMod[c2, b, cmN], cmN]];

cmMsgs = Table[Mod[bigOf[s, 900], cmN], {s, 1810, 1909}];   (* 100 messages *)

bench["app: common-modulus attack 1024-bit x600",
      Do[cmRec = commonModulus /@ cmMsgs, {6}];];
check["app: common-modulus attack 1024-bit x600", Mod[Total[cmRec], $CheckMod] + $CheckMod * (Count[MapThread[Equal, {cmRec, cmMsgs}], True])];

(* ---- (c) Hastad broadcast (e = 3) -------------------------------------- *)

(* Three fixed coprime 1024-bit moduli. *)
hMods = {NextPrime[bigOf[1900, 512]] NextPrime[bigOf[1901, 512]],
         NextPrime[bigOf[1902, 512]] NextPrime[bigOf[1903, 512]],
         NextPrime[bigOf[1904, 512]] NextPrime[bigOf[1905, 512]]};

(* Integer cube root by integer Newton iteration (robust for any size). *)
icbrt[x_] := Module[{r, rn},
  If[x == 0, Return[0]];
  r = 2^Ceiling[IntegerLength[x, 2]/3];
  rn = Quotient[2 r + Quotient[x, r^2], 3];
  While[rn < r, r = rn; rn = Quotient[2 r + Quotient[x, r^2], 3]];
  While[r^3 > x, r--];
  r];

hastad[m_] := Module[{cs},
  cs = Mod[m^3, #] & /@ hMods;               (* ciphertexts under e = 3 *)
  icbrt[ChineseRemainder[cs, hMods]]];        (* CRT lift, then cube root *)

hMsgs = Table[bigOf[s, 300], {s, 1950, 1999}];   (* 50 messages, m^3 < prod *)

bench["app: Hastad broadcast e=3 x1500", Do[hRec = hastad /@ hMsgs, {30}];];
check["app: Hastad broadcast e=3 x1500", Mod[Total[hRec], $CheckMod] + $CheckMod * (Count[MapThread[Equal, {hRec, hMsgs}], True])];
