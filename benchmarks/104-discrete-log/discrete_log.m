(* ==========================================================================
   Experiment 104 -- discrete logarithms (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  Generic discrete-log algorithms in (Z/pZ)^*, written in
   the language ("app:"), against each system's builtin ("prim:"):
     - baby-step giant-step in a prime-order-r subgroup, with the baby-step
       table built two ways: inserted one key at a time (tab[y] = i -- the
       textbook loop), and built in bulk from the list of powers
       (AssociationThread / dict(zip) / Map(matrix));
     - Pollard rho for logarithms: the standard 3-partition walk on x mod 3
       (square / times g / times h), Floyd cycle finding, restart from g^t on
       a degenerate collision;
     - Pohlig-Hellman over a 131-bit prime p with p - 1 = 2^3 3 q1...q8 (eight
       primes just above 40000, 45000, ..., 75000), digit-by-digit per prime
       power with the bulk BSGS for each prime, recombined by
       ChineseRemainder;
     - prim: the WL discrete-log form MultiplicativeOrder[g, p, {h}] /
       sympy.ntheory.discrete_log / PARI znlog.

   Instances: r = primeOf[s, bits] (24- and 32-bit), p = 2 k r + 1 the first
   prime of that shape, g = 2^((p-1)/r) (next base if that is 1), and
   h = g^x with x = bigOf[t, 64] mod r planted; for the smooth p, g is the
   least primitive root and x = bigOf[t, 128] mod (p - 1).  Every check is the
   sum of recovered exponents mod $CheckMod, which must equal the sum of the
   planted ones (24-bit: 39973122, 32-bit: 595305666, smooth: 250044235).
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"MultiplicativeOrder", "ChineseRemainder", "PowerMod", "AssociationThread",
         "KeyExistsQ", "FactorInteger"}];

(* ---- instances (untimed) ------------------------------------------------- *)

subgroup[bits_, s_] := Module[{r = primeOf[s, bits], k = 1, p, b = 2},
  While[! PrimeQ[2 k r + 1], k++];
  p = 2 k r + 1;
  While[PowerMod[b, (p - 1)/r, p] == 1, b++];
  {p, r, PowerMod[b, (p - 1)/r, p]}];

instance[bits_, s_, t_] := Module[{p, r, g, x},
  {p, r, g} = subgroup[bits, s];
  x = Mod[bigOf[t, 64], r];
  {p, r, g, PowerMod[g, x, p], x}];

inst24 = Table[instance[24, 9000 + i, 9100 + i], {i, 0, 3}];
inst32 = Table[instance[32, 9200 + i, 9300 + i], {i, 0, 1}];

qs = Table[NextPrime[40000 + 5000 i], {i, 0, 7}];
smoothL = Times @@ qs;
smoothK = 1;
While[! PrimeQ[2 smoothK smoothL + 1], smoothK++];
pS = 2 smoothK smoothL + 1;
fS = FactorInteger[pS - 1];
gS = 2;
While[AnyTrue[fS[[All, 1]], PowerMod[gS, (pS - 1)/#, pS] == 1 &], gS++];
xsSmooth = Table[Mod[bigOf[s, 128], pS - 1], {s, 9401, 9410}];
hsSmooth = PowerMod[gS, xsSmooth, pS];

(* ---- app: baby-step giant-step -------------------------------------------- *)

(* Baby steps inserted one at a time into an Association. *)
bsgsInc[g_, h_, p_, r_] := Module[{m = Floor[Sqrt[r]] + 1, tab = <||>, y = 1, c, j = 0},
  Do[tab[y] = i; y = Mod[y g, p], {i, 0, m - 1}];
  c = PowerMod[g, -m, p];
  y = h;
  While[! KeyExistsQ[tab, y], y = Mod[y c, p]; j++];
  j m + tab[y]];

(* Baby-step table built in bulk from the list of powers. *)
bsgs[g_, h_, p_, r_] := Module[{m = Floor[Sqrt[r]] + 1, tab, c, y = h, j = 0},
  tab = AssociationThread[NestList[Mod[# g, p] &, 1, m - 1] -> Range[0, m - 1]];
  c = PowerMod[g, -m, p];
  While[! KeyExistsQ[tab, y], y = Mod[y c, p]; j++];
  j m + tab[y]];

bench["app: BSGS 24-bit subgroup (incremental table) x4",
      b1 = bsgsInc[#[[3]], #[[4]], #[[1]], #[[2]]] & /@ inst24;];
check["app: BSGS 24-bit subgroup (incremental table) x4", Mod[Total[b1], $CheckMod]];

bench["app: BSGS 24-bit subgroup (bulk table) x4",
      b2 = bsgs[#[[3]], #[[4]], #[[1]], #[[2]]] & /@ inst24;];
check["app: BSGS 24-bit subgroup (bulk table) x4", Mod[Total[b2], $CheckMod]];

bench["app: BSGS 32-bit subgroup (bulk table) x2",
      b3 = bsgs[#[[3]], #[[4]], #[[1]], #[[2]]] & /@ inst32;];
check["app: BSGS 32-bit subgroup (bulk table) x2", Mod[Total[b3], $CheckMod]];

(* ---- app: Pollard rho for logarithms -------------------------------------- *)

rhoStep[{x_, a_, b_}, g_, h_, p_, r_] := Switch[Mod[x, 3],
  0, {Mod[x^2, p], Mod[2 a, r], Mod[2 b, r]},
  1, {Mod[x g, p], Mod[a + 1, r], b},
  _, {Mod[x h, p], a, Mod[b + 1, r]}];

pollardRho[g_, h_, p_, r_] := Module[{t = 0, ans = -1, tort, hare, den},
  While[ans < 0,
    tort = {PowerMod[g, t, p], t, 0};
    hare = tort;
    tort = rhoStep[tort, g, h, p, r];
    hare = rhoStep[rhoStep[hare, g, h, p, r], g, h, p, r];
    While[tort[[1]] != hare[[1]],
      tort = rhoStep[tort, g, h, p, r];
      hare = rhoStep[rhoStep[hare, g, h, p, r], g, h, p, r]];
    den = Mod[tort[[3]] - hare[[3]], r];
    If[den != 0, ans = Mod[(hare[[2]] - tort[[2]]) PowerMod[den, -1, r], r], t++]];
  ans];

benchOnce["app: Pollard rho 32-bit subgroup x2",
          r1 = pollardRho[#[[3]], #[[4]], #[[1]], #[[2]]] & /@ inst32;];
check["app: Pollard rho 32-bit subgroup x2", Mod[Total[r1], $CheckMod]];

(* ---- app: Pohlig-Hellman over a smooth p ----------------------------------- *)

pohligHellman[g_, h_, p_, fs_] := Module[{n = p - 1, rs = {}, ms = {}, q, e, gq, xk, hk},
  Do[
    {q, e} = qe;
    gq = PowerMod[g, n/q, p];
    xk = 0;
    Do[
      hk = PowerMod[Mod[PowerMod[g, -xk, p] h, p], n/q^(j + 1), p];
      xk += bsgs[gq, hk, p, q] q^j,
      {j, 0, e - 1}];
    AppendTo[rs, xk]; AppendTo[ms, q^e],
    {qe, fs}];
  ChineseRemainder[rs, ms]];

bench["app: Pohlig-Hellman 131-bit smooth p x10",
      ph = pohligHellman[gS, #, pS, fS] & /@ hsSmooth;];
check["app: Pohlig-Hellman 131-bit smooth p x10", Mod[Total[ph], $CheckMod]];

(* ---- prim: MultiplicativeOrder[g, p, {h}] ----------------------------------- *)

(* Repeated 10x so PARI's znlog clears its 1-ms timer; at ~350 ms per call
   in Mathilda that makes this the one case over 5 s, hence benchOnce. *)
benchOnce["prim: discrete log 24-bit subgroup x40",
          Do[d1 = MultiplicativeOrder[#[[3]], #[[1]], {#[[4]]}] & /@ inst24, {10}];];
check["prim: discrete log 24-bit subgroup x40", Mod[Total[d1], $CheckMod]];

(* The 32-bit and smooth-p rows are capped with TimeConstrained (20 s) and run
   once: Mathematica's MultiplicativeOrder[g, p, {h}] is an exhaustive search
   (~350 s per 32-bit call, unbounded for the 131-bit p -- no Pohlig-Hellman),
   and Mathilda returns these unevaluated; see README "Findings".  A capped
   run prints $Aborted, an unevaluated one the held form -- either is a
   CHECK-FAIL against the planted exponents, which is the honest result. *)
benchOnce["prim: discrete log 32-bit subgroup x10",
          TimeConstrained[Do[d2 = MultiplicativeOrder[#[[3]], #[[1]], {#[[4]]}] & /@ inst32, {5}], 20, d2 = $Aborted];];
check["prim: discrete log 32-bit subgroup x10", Mod[Total[d2], $CheckMod]];

benchOnce["prim: discrete log 131-bit smooth p x40",
          TimeConstrained[Do[d3 = MultiplicativeOrder[gS, pS, {#}] & /@ hsSmooth, {4}], 20, d3 = $Aborted];];
check["prim: discrete log 131-bit smooth p x40", Mod[Total[d3], $CheckMod]];
