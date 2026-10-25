(* ==========================================================================
   Experiment 112 -- Rabin cryptosystem & modular square roots (group F)
   ==========================================================================
   WHAT IT MEASURES.  Square roots modulo a prime -- the kernel behind the
   Rabin cryptosystem and the final step of index-calculus / point
   decompression.  Two classical algorithms written in the language
   ("app:"): Tonelli-Shanks and Cipolla, both exercised on the p = 1 mod 8
   case (S >= 3 in p - 1 = Q 2^S, the branch the easy (p+1)/4 shortcut cannot
   take).  Then the WL primitive route ("prim:"): PowerMod[a, 1/2, p], which
   is how a square root mod p is spelled in the Wolfram Language -- there is
   no SqrtMod head (it is absent in Mathilda AND does not exist in WL), so the
   "prim: SqrtMod" label is benchIf'd and reported ABSENT, while
   PowerMod[a, 1/2, p] is timed as the real WL primitive.

   Every root is canonicalised to min(r, p - r) before checking, so the four
   systems -- which may return either of the two roots -- agree on one
   integer.  The Rabin case decrypts a ciphertext to its four square roots
   mod n = p q (p = q = 3 mod 4) via CRT and sums their canonical forms.

   Inputs come from crypto_data.m's shared LCG.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod", "JacobiSymbol", "PrimeQ", "NextPrime"}];

(* First prime >= primeOf[seed, bits] congruent to r mod m. *)
primeMod[seed_, bits_, r_, m_] := Module[{p = primeOf[seed, bits]},
  While[Mod[p, m] != r, p = NextPrime[p]]; p];

canon[r_, p_] := Min[r, p - r];

(* ---- Tonelli-Shanks ---------------------------------------------------- *)

tonelli[n_, p_] := Module[{Q, S0, z, M, c, t, R, i, b, tt},
  Q = p - 1; S0 = 0; While[EvenQ[Q], Q /= 2; S0++];
  If[S0 == 1, Return[PowerMod[n, (p + 1)/4, p]]];
  z = 2; While[PowerMod[z, (p - 1)/2, p] != p - 1, z++];
  M = S0; c = PowerMod[z, Q, p]; t = PowerMod[n, Q, p];
  R = PowerMod[n, (Q + 1)/2, p];
  While[t != 1,
    i = 0; tt = t; While[tt != 1, tt = Mod[tt^2, p]; i++];
    b = PowerMod[c, 2^(M - i - 1), p]; M = i; c = Mod[b^2, p];
    t = Mod[t c, p]; R = Mod[R b, p]];
  R];


(* ---- Cipolla ----------------------------------------------------------- *)

qrQ[x_, p_] := x == 0 || PowerMod[x, (p - 1)/2, p] == 1;

(* Multiply in F_p(sqrt w): (x1 + y1 s)(x2 + y2 s) = (x1 x2 + y1 y2 w) + (...) s *)
fmul[{x1_, y1_}, {x2_, y2_}, w_, p_] :=
  {Mod[x1 x2 + y1 y2 w, p], Mod[x1 y2 + y1 x2, p]};

fpow[base_, e_, w_, p_] := Module[{result = {1, 0}, b = base, ee = e},
  While[ee > 0,
    If[OddQ[ee], result = fmul[result, b, w, p]];
    b = fmul[b, b, w, p]; ee = Quotient[ee, 2]];
  result];

cipolla[n_, p_] := Module[{a0 = 1, w},
  While[qrQ[Mod[a0^2 - n, p], p], a0++];
  w = Mod[a0^2 - n, p];
  fpow[{a0, 1}, (p + 1)/2, w, p][[1]]];

(* ---- inputs: p = 1 mod 8 (hard case), several sizes ------------------- *)

p256 = primeMod[1201, 256, 1, 8];
as256 = Table[Mod[bigOf[s, 256]^2, p256], {s, 10, 109}];   (* 100 QRs *)

p1024 = primeMod[1251, 1024, 1, 8];
as1024 = Table[Mod[bigOf[s, 1024]^2, p1024], {s, 10, 19}];  (* 10 QRs *)

(* ---- app: Tonelli-Shanks ---------------------------------------------- *)

bench["app: Tonelli-Shanks 256-bit p=1mod8 x100",
      ts256 = canon[tonelli[#, p256], p256] & /@ as256;];
check["app: Tonelli-Shanks 256-bit p=1mod8 x100", Mod[Total[ts256], $CheckMod]];

benchOnce["app: Tonelli-Shanks 1024-bit p=1mod8 x10",
          ts1024 = canon[tonelli[#, p1024], p1024] & /@ as1024;];
check["app: Tonelli-Shanks 1024-bit p=1mod8 x10", Mod[Total[ts1024], $CheckMod]];

(* ---- app: Cipolla ------------------------------------------------------ *)

bench["app: Cipolla 256-bit p=1mod8 x100",
      cp256 = canon[cipolla[#, p256], p256] & /@ as256;];
check["app: Cipolla 256-bit p=1mod8 x100", Mod[Total[cp256], $CheckMod]];

benchOnce["app: Cipolla 1024-bit p=1mod8 x10",
          cp1024 = canon[cipolla[#, p1024], p1024] & /@ as1024;];
check["app: Cipolla 1024-bit p=1mod8 x10", Mod[Total[cp1024], $CheckMod]];

(* ---- prim: PowerMod[a, 1/2, p] (the WL spelling) ----------------------- *)

bench["prim: PowerMod[a,1/2,p] 256-bit x500",
      Do[pm256 = canon[PowerMod[#, 1/2, p256], p256] & /@ as256, {5}];];
check["prim: PowerMod[a,1/2,p] 256-bit x500", Mod[Total[pm256], $CheckMod]];

bench["prim: PowerMod[a,1/2,p] 1024-bit x40",
      Do[pm1024 = canon[PowerMod[#, 1/2, p1024], p1024] & /@ as1024, {4}];];
check["prim: PowerMod[a,1/2,p] 1024-bit x40", Mod[Total[pm1024], $CheckMod]];

(* SqrtMod is the non-WL spelling; absent in Mathilda and in WL itself. *)
benchIf["prim: SqrtMod 256-bit x500", "SqrtMod",
        Do[sm256 = canon[SqrtMod[#, p256], p256] & /@ as256, {5}];];
checkIf["prim: SqrtMod 256-bit x500", "SqrtMod",
        Mod[Total[sm256], $CheckMod]];

(* ---- app: Rabin decrypt (four roots via CRT) --------------------------- *)

pr = primeMod[1301, 256, 3, 4];
qr = primeMod[1302, 256, 3, 4];
nr = pr qr;
(* yp = p^-1 mod q, yq = q^-1 mod p (fixed for this n). *)
ypr = PowerMod[pr, -1, qr]; yqr = PowerMod[qr, -1, pr];

(* Four roots of c mod n, summed by canonical form. *)
rabinRoots[c_] := Module[{mp, mq, u, v},
  mp = PowerMod[c, (pr + 1)/4, pr]; mq = PowerMod[c, (qr + 1)/4, qr];
  u = Mod[mp qr yqr + mq pr ypr, nr]; v = Mod[mp qr yqr - mq pr ypr, nr];
  Mod[canon[u, nr] + canon[v, nr], nr]];

msgs = Table[Mod[bigOf[s, 256], nr], {s, 500, 999}];       (* 500 messages *)
cts = Mod[#^2, nr] & /@ msgs;

bench["app: Rabin decrypt 512-bit n x500", rr = rabinRoots /@ cts;];
check["app: Rabin decrypt 512-bit n x500", Mod[Total[rr], $CheckMod]];
