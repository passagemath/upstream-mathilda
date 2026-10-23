(* BesselHyper.m -- indefinite integrals of
     c x^m Z1(w) [Z2(w)],   w = b x^p,   Z in {BesselJ, BesselY, BesselI, BesselK}
   in closed form as HypergeometricPFQ / MeijerG.

   Loaded lazily by the indefinite-integration stage "BesselHypergeometric"
   (src/calculus/integrate.c), which runs after "Holonomic": an antiderivative
   that closes in Bessel functions (Lommel's x Z_nu^2) is found there first,
   and this stage covers the integrals that have none (Z_nu^2, Z_nu, J_mu J_nu
   with mu - nu non-integer, ...).  All symbols are qualified in BesselHyper`.

   DERIVATION (every identity below is proved, not tabled).
   D1  J_nu(w) = (w/2)^nu / Gamma(nu+1) 0F1(; nu+1; -w^2/4)   (defining series;
       I_nu with +w^2/4).
   D2  J_mu J_nu: the Cauchy product of the two series has, after scaling by
       k! Gamma(mu+k+1) Gamma(nu+k+1), inner sum
         Sum_j C(k,j) (mu+j+1)_(k-j) (nu+k-j+1)_j = (mu+nu+k+1)_k
       (Chu-Vandermonde), so
         J_mu J_nu = (w/2)^(mu+nu) / (Gamma(mu+1) Gamma(nu+1))
                     2F3((mu+nu+1)/2, (mu+nu+2)/2; mu+1, nu+1, mu+nu+1; -w^2)
       (I_mu I_nu with +w^2).
   D3  For a non-integer order: Y_nu = (J_nu Cos[nu Pi] - J_-nu)/Sin[nu Pi],
       K_nu = Pi (I_-nu - I_nu) / (2 Sin[nu Pi]); products expand bilinearly.
   D4  Term-by-term: d/dx [x^(m+1) w^(lam+2k)] = (m+1+p(lam+2k)) x^m w^(lam+2k)
       holds on every branch, so with alpha = (m+1+p lam)/(2p)
         Int x^m w^lam pFq(A; B; kap w^2) dx
           = x^(m+1) w^lam/(m+1+p lam) p+1Fq+1(A, alpha; B, alpha+1; kap w^2).
       If alpha = -j (j = 0, 1, ...) the k = j term is a Log[x] and the tail is
       re-indexed into a shifted p+2Fq+2.
   D5  Integer or symbolic orders of Y / K (where D3 divides by zero): the
       Mellin-Barnes representations, Gamma factors chosen so the residue sum
       reproduces the D1/D2 series and fixed by the behaviour at infinity
       (t = w^2, verified numerically to 1e-40 at several points and orders):
         Y_nu      = G^{2,0}_{1,3}(w^2/4 | ; -(nu+1)/2 | nu/2, -nu/2 ; -(nu+1)/2)
         K_nu      = 1/2 G^{2,0}_{0,2}(w^2/4 | nu/2, -nu/2)
         K_mu K_nu = Sqrt[Pi]/2 G^{4,0}_{2,4}(t | ; 0, 1/2 | (+-mu+-nu)/2)
         I_mu K_nu = 1/(2 Sqrt[Pi]) G^{2,2}_{2,4}(t | 0, 1/2 ; |
                       (mu+nu)/2, (mu-nu)/2 ; (nu-mu)/2, -(mu+nu)/2)
         J_nu Y_nu = -1/Sqrt[Pi] G^{2,0}_{1,3}(t | ; 1/2 | nu, 0 ; -nu)
         J_mu Y_nu = G^{2,1}_{2,4}(t | 1/2 ; 0 | (mu+nu)/2, (mu-nu)/2 ;
                       (nu-mu)/2, -(mu+nu)/2) / (Sqrt[Pi] Sin[Pi (mu-nu)/2])
                     - Cot[Pi (mu-nu)/2] J_mu J_nu          (mu - nu not even)
         J_nu^2 + Y_nu^2 = 2 Cos[nu Pi]/Pi^(5/2) G^{3,1}_{1,3}(t | 1/2 ; | nu, 0, -nu)
       and Int t^(beta-1) G(t) dt = t^beta G^{m,n+1}_{p+1,q+1}(t | 1-beta, a ; b, -beta)
       (integrate the contour integrand under the integral sign).
   D6  Orders differing by a non-zero integer: the three-term recurrence writes
       Z_(nu0+i) in the basis {Z_nu0, Z_(nu0+1)}, and with
       X = Z_(nu0+1) W_nu0, Y = Z_nu0 W_(nu0+1), D = Z_nu0 W_nu0,
         sZ X + sW Y = 2 nu0/w D - D'      (Z' = nu/w Z - sZ Z_(nu+1))
         sZ X - sW Y = Wronskian[Z, W]     (2/(Pi w), -1/w, or 0)
       so every such product is a Laurent combination of the diagonal products
       D(nu0), D(nu0+1) and D'(nu0); the D' terms are integrated by parts.  This
       is what removes the empty Mellin strips of J_0 Y_1, I_0 K_1, Y_0 Y_2, ...

   BOUNDARY (declines with $Failed): three or more Bessel factors; different
   arguments J(a x) J(b x) (an Appell F4 double series); mixed ordinary /
   modified families (J I, J K, Y K: the D2 inner sum does not close);
   parameter sets where a series denominator or a Mellin-Barnes Gamma pair
   degenerates (half-integer pairs such as J_(-1/2) J_(-1/2), which are
   elementary and closed by other stages). *)

BesselHyper`fail = BesselHyper`fail;

BesselHyper`Heads = BesselJ | BesselY | BesselI | BesselK;

(* w = b x^p with b, p free of x, p != 0: returns {b, p} or fails *)
BesselHyper`ArgForm[w_, x_] := Module[{p, b},
  If[w === x, Return[{1, 1}]];
  If[FreeQ[w, x], Throw[$Failed, BesselHyper`fail]];
  p = Simplify[x D[w, x]/w];
  If[!FreeQ[p, x] || p === 0, Throw[$Failed, BesselHyper`fail]];
  b = Simplify[w/x^p];
  If[!FreeQ[b, x], b = Simplify[PowerExpand[w/x^p]]];
  If[!FreeQ[b, x], Throw[$Failed, BesselHyper`fail]];
  {b, p}];

BesselHyper`IntQ[e_] := IntegerQ[Simplify[e]];
BesselHyper`NonIntNumQ[nu_] := NumericQ[nu] && !BesselHyper`IntQ[nu];
BesselHyper`NonPosIntQ[e_] := With[{s = Simplify[e]}, IntegerQ[s] && s <= 0];
BesselHyper`PosIntQ[e_] := With[{s = Simplify[e]}, IntegerQ[s] && s > 0];

BesselHyper`Ordinary[h_] := MemberQ[{BesselJ, BesselY}, h];
BesselHyper`Rank[h_] := Switch[h, BesselJ, 1, BesselY, 2, BesselI, 1, BesselK, 2];

(* ---------- D4: integrate c x^m w^lam pFq(A; B; kap w^2) ---------- *)
BesselHyper`IntSeries[c_, m_, lam_, A0_List, B0_List, kap_, w_, p_, x_] := Module[
  {A = Expand /@ A0, B = Expand /@ B0}, BesselHyper`IntSeries1[c, m, lam, A, B, kap, w, p, x]];
BesselHyper`IntSeries1[c_, m_, lam_, A_List, B_List, kap_, w_, p_, x_] := Module[
  {alpha, j, tk, z = kap w^2, head, tail, logt},
  If[AnyTrue[B, BesselHyper`NonPosIntQ], Throw[$Failed, BesselHyper`fail]];
  alpha = Expand[Simplify[(m + 1 + p lam)/(2 p)]];
  If[!BesselHyper`NonPosIntQ[alpha],
    Return[c x^(m + 1) w^lam/(2 p alpha) HypergeometricPFQ[Append[A, alpha], Append[B, alpha + 1], z]]];
  j = -alpha;
  tk[k_] := (Times @@ (Pochhammer[#, k] & /@ A))/(Times @@ (Pochhammer[#, k] & /@ B))/k! kap^k;
  head = Sum[c tk[k] x^(m + 1) w^(lam + 2 k)/(m + 1 + p (lam + 2 k)), {k, 0, j - 1}];
  (* x^(m+1) w^(lam+2j) is locally constant: its derivative vanishes *)
  logt = c tk[j] x^(m + 1) w^(lam + 2 j) Log[x];
  tail = c tk[j + 1] x^(m + 1) w^(lam + 2 j + 2)/(2 p) *
    HypergeometricPFQ[Join[A + j + 1, {1, 1}], Join[B + j + 1, {j + 2, 2}], z];
  head + logt + tail];

(* ---------- D5: integrate c x^m G(kap w^2) ---------- *)
BesselHyper`Degenerate[an_List, bm_List] :=
  AnyTrue[Flatten[Outer[Subtract, an, bm]], BesselHyper`PosIntQ];

BesselHyper`IntMeijer[c_, m_, {an0_List, ap0_List}, {bm0_List, bq0_List}, kap_, w_, p_, x_] := Module[
  {beta = Expand[Simplify[(m + 1)/(2 p)]], an = Expand /@ an0, ap = Expand /@ ap0,
   bm = Expand /@ bm0, bq = Expand /@ bq0, an2, bq2},
  If[BesselHyper`Degenerate[an, bm], Throw[$Failed, BesselHyper`fail]];
  an2 = Prepend[an, 1 - beta]; bq2 = Append[bq, -beta];
  If[BesselHyper`Degenerate[an2, bm], Throw[$Failed, BesselHyper`fail]];
  c x^(m + 1)/(2 p) MeijerG[{an2, ap}, {bm, bq2}, kap w^2]];

(* ---------- D3: replace a Y / K factor of non-integer order ---------- *)
BesselHyper`Split[{BesselY, nu_}] :=
  {{Cot[nu Pi], {BesselJ, nu}}, {-Csc[nu Pi], {BesselJ, -nu}}};
BesselHyper`Split[{BesselK, nu_}] :=
  {{Pi Csc[nu Pi]/2, {BesselI, -nu}}, {-Pi Csc[nu Pi]/2, {BesselI, nu}}};

(* Negative integer orders: J_-n = (-1)^n J_n, Y_-n = (-1)^n Y_n, I_-n = I_n, K_-n = K_n *)
BesselHyper`Normalize[{h_, nu_}] := With[{s = Simplify[nu]},
  If[IntegerQ[s] && s < 0,
    {If[BesselHyper`Ordinary[h], (-1)^s, 1], {h, -s}},
    {1, {h, nu}}]];

(* ---------- the dispatcher: Int c x^m Prod Z(w) dx ---------- *)
BesselHyper`Int[c_, m_, fs_List, w_, b_, p_, x_] := Module[{nf, sp, pos, h1, n1, h2, n2, d, bs},
  nf = BesselHyper`Normalize /@ fs;
  If[!AllTrue[nf[[All, 1]], # === 1 &],
    Return[BesselHyper`Int[c Times @@ nf[[All, 1]], m, nf[[All, 2]], w, b, p, x]]];
  (* D3 first for every Y / K of numeric non-integer order: pFq output *)
  pos = Position[fs, {BesselY | BesselK, _?BesselHyper`NonIntNumQ}, {1}, 1];
  If[pos =!= {},
    pos = pos[[1, 1]];
    Return[Total[(BesselHyper`Int[c #[[1]], m, ReplacePart[fs, pos -> #[[2]]], w, b, p, x] &) /@
      BesselHyper`Split[fs[[pos]]]]]];
  Switch[Length[fs],
   1,
    {h1, n1} = fs[[1]];
    Switch[h1,
      BesselJ, BesselHyper`IntSeries[c 2^-n1/Gamma[n1 + 1], m, n1, {}, {n1 + 1}, -1/4, w, p, x],
      BesselI, BesselHyper`IntSeries[c 2^-n1/Gamma[n1 + 1], m, n1, {}, {n1 + 1}, 1/4, w, p, x],
      BesselY, BesselHyper`IntMeijer[c, m, {{}, {-(n1 + 1)/2}}, {{n1/2, -n1/2}, {-(n1 + 1)/2}}, 1/4, w, p, x],
      BesselK, BesselHyper`IntMeijer[c/2, m, {{}, {}}, {{n1/2, -n1/2}, {}}, 1/4, w, p, x]],
   2,
    {{h1, n1}, {h2, n2}} = SortBy[fs, BesselHyper`Rank[#[[1]]] &];
    If[BesselHyper`Ordinary[h1] =!= BesselHyper`Ordinary[h2], Throw[$Failed, BesselHyper`fail]];
    If[h1 === h2 && (h1 === BesselJ || h1 === BesselI),
      Return[BesselHyper`IntSeries[c 2^(-n1 - n2)/(Gamma[n1 + 1] Gamma[n2 + 1]), m, n1 + n2,
        {(n1 + n2 + 1)/2, (n1 + n2 + 2)/2}, {n1 + 1, n2 + 1, n1 + n2 + 1},
        If[h1 === BesselJ, -1, 1], w, p, x]]];
    bs = {(n1 + n2)/2, (n1 - n2)/2, (n2 - n1)/2, -(n1 + n2)/2};
    d = Simplify[n1 - n2];
    If[d === 0 && Simplify[m + 1] === 0 && IntegerQ[Simplify[n1]] && Simplify[n1] >= 1,
      Return[BesselHyper`Telescope[c, {h1, h2}, Simplify[n1], w, p, x]]];
    If[h1 === BesselK,
      Return[BesselHyper`IntMeijer[c Sqrt[Pi]/2, m, {{}, {0, 1/2}}, {bs, {}}, 1, w, p, x]]];
    If[IntegerQ[d] && d != 0, Return[BesselHyper`Reduce[c, m, {h1, n1}, {h2, n2}, w, b, p, x]]];
    If[d === 0,
      Return[Switch[{h1, h2},
        {BesselI, BesselK}, BesselHyper`IntMeijer[c/(2 Sqrt[Pi]), m, {{1/2}, {}}, {{n1, 0}, {-n1}}, 1, w, p, x],
        {BesselJ, BesselY}, BesselHyper`IntMeijer[-c/Sqrt[Pi], m, {{}, {1/2}}, {{n1, 0}, {-n1}}, 1, w, p, x],
        {BesselY, BesselY}, BesselHyper`IntMeijer[2 c Cos[n1 Pi]/Pi^(5/2), m, {{1/2}, {}}, {{n1, 0, -n1}, {}}, 1, w, p, x] -
                            BesselHyper`Int[c, m, {{BesselJ, n1}, {BesselJ, n1}}, w, b, p, x]]]];
    (* orders differing by a non-integer *)
    Switch[{h1, h2},
      {BesselI, BesselK}, BesselHyper`IntMeijer[c/(2 Sqrt[Pi]), m, {{0, 1/2}, {}}, {bs[[{1, 2}]], bs[[{3, 4}]]}, 1, w, p, x],
      {BesselJ, BesselY}, BesselHyper`IntMeijer[c/(Sqrt[Pi] Sin[Pi d/2]), m, {{1/2}, {0}}, {bs[[{1, 2}]], bs[[{3, 4}]]}, 1, w, p, x] -
                          BesselHyper`Int[c Cot[Pi d/2], m, {{BesselJ, n1}, {BesselJ, n2}}, w, b, p, x],
      (* Y_mu Y_nu with a symbolic order: D3 on it (generic in that order) *)
      {BesselY, BesselY}, With[{k = If[NumericQ[n1], 2, 1]},
        Total[(BesselHyper`Int[c #[[1]], m, ReplacePart[{{h1, n1}, {h2, n2}}, k -> #[[2]]], w, b, p, x] &) /@
          BesselHyper`Split[{{h1, n1}, {h2, n2}}[[k]]]]]],
   _, Throw[$Failed, BesselHyper`fail]]];

(* ---------- D7: Int D_nu(w) dx/x for an integer order nu >= 1 ----------
   With D_nu = Z_nu W_nu the derivative rules give
     (D_nu + s D_(nu+1))' = 2 nu/w D_nu - 2 s (nu+1)/w D_(nu+1)
   (s = -1 for I^2 and K^2, +1 otherwise), so
     Int D_nu/w dw = (2 (nu-1) Int D_(nu-1)/w dw - (D_(nu-1) + s D_nu))/(2 s nu),
   which telescopes to nu = 0, where the coefficient vanishes.  This is the
   log-degenerate Mellin-Barnes case (beta = 0) left by D5/D6 (J_n Y_(n+1)). *)
BesselHyper`Telescope[c_, {h1_, h2_}, nu_, w_, p_, x_] := Module[
  {s = If[h1 === h2 && (h1 === BesselI || h1 === BesselK), -1, 1], dm, dn},
  dm = h1[nu - 1, w] h2[nu - 1, w]; dn = h1[nu, w] h2[nu, w];
  c/(2 s nu) (If[nu == 1, 0, 2 (nu - 1) BesselHyper`Telescope[1, {h1, h2}, nu - 1, w, p, x]] -
    (dm + s dn)/p)];

(* ---------- D6: orders differing by a non-zero integer ---------- *)
BesselHyper`Rec[h_] := Switch[h, BesselJ | BesselY, {1, -1}, BesselI, {-1, 1}, BesselK, {1, 1}];
BesselHyper`Sigma[h_] := If[h === BesselI, -1, 1];
(* Z_(nu0+i) = a Z_nu0 + c Z_(nu0+1): {a, c} *)
BesselHyper`Basis[h_, nu0_, i_, w_] := Module[{s0 = {1, 0}, s1 = {0, 1}, s2, r = BesselHyper`Rec[h]},
  If[i == 0, Return[s0]];
  Do[s2 = r[[1]] 2 (nu0 + k + 1)/w s1 + r[[2]] s0; s0 = s1; s1 = s2, {k, 0, i - 2}];
  s1];
BesselHyper`Wronskian[h1_, h2_, w_] := Switch[{h1, h2},
  {BesselJ, BesselY}, 2/(Pi w), {BesselY, BesselJ}, -2/(Pi w),
  {BesselI, BesselK}, -1/w, {BesselK, BesselI}, 1/w, _, 0];

BesselHyper`Reduce[c_, m_, {h1_, n1_}, {h2_, n2_}, w_, b_, p_, x_] := Module[
  {d = Simplify[n1 - n2], nu0, za, wa, X, Yv, expr, ws = BesselHyper`w, dd, dd1, ddp, res = 0,
   sZ = BesselHyper`Sigma[h1], sW = BesselHyper`Sigma[h2], wr, terms, k, cf, e, prod},
  nu0 = If[d > 0, n2, n1];
  za = BesselHyper`Basis[h1, nu0, Simplify[n1 - nu0], ws];
  wa = BesselHyper`Basis[h2, nu0, Simplify[n2 - nu0], ws];
  wr = BesselHyper`Wronskian[h1, h2, ws];
  X  = ((2 nu0/ws dd - ddp) + wr)/(2 sZ);
  Yv = ((2 nu0/ws dd - ddp) - wr)/(2 sW);
  expr = Expand[za[[1]] wa[[1]] dd + za[[2]] wa[[1]] X + za[[1]] wa[[2]] Yv + za[[2]] wa[[2]] dd1];
  prod = h1[nu0, w] h2[nu0, w];
  Do[
    cf = Expand[If[sym === 1, expr /. {dd -> 0, dd1 -> 0, ddp -> 0}, Coefficient[expr, sym]]];
    terms = If[Head[cf] === Plus, List @@ cf, {cf}];
    Do[
      If[t =!= 0,
        k = Exponent[t, ws];
        If[!IntegerQ[k], k = -Exponent[1/t, ws]];
        If[!IntegerQ[k] || !FreeQ[t/ws^k, ws], Throw[$Failed, BesselHyper`fail]];
        (* c t = (c t/ws^k) b^k x^(p k) *)
        cf = Simplify[c t/ws^k] b^k; e = m + p k;
        res += Switch[sym,
          dd,  BesselHyper`Int[cf, e, {{h1, nu0}, {h2, nu0}}, w, b, p, x],
          dd1, BesselHyper`Int[cf, e, {{h1, nu0 + 1}, {h2, nu0 + 1}}, w, b, p, x],
          (* Int cf x^e dD/dw dx = cf/(p b) Int x^(e+1-p) dD/dx dx, by parts *)
          ddp, cf/(p b) (x^(e + 1 - p) prod -
                 If[Simplify[e + 1 - p] === 0, 0,
                   (e + 1 - p) BesselHyper`Int[1, e - p, {{h1, nu0}, {h2, nu0}}, w, b, p, x]]),
          1,   If[Simplify[e + 1] === 0, cf Log[x], cf x^(e + 1)/(e + 1)]]],
      {t, terms}],
    {sym, {dd, dd1, ddp, 1}}];
  res];

(* ---------- entry point ---------- *)
BesselHyper`Term[f_, x_] := Module[{fs = {}, c = 1, m = 0, args = {}, bp, w},
  Scan[Function[t,
    Which[
      FreeQ[t, x], c *= t,
      t === x, m += 1,
      MatchQ[t, Power[x, k_ /; FreeQ[k, x]]], m += t[[2]],
      MatchQ[t, (BesselHyper`Heads)[nu_ /; FreeQ[nu, x], _]],
        AppendTo[fs, {Head[t], t[[1]]}]; AppendTo[args, t[[2]]],
      MatchQ[t, Power[(BesselHyper`Heads)[nu_ /; FreeQ[nu, x], _], 2]],
        fs = Join[fs, {{Head[t[[1]]], t[[1, 1]]}, {Head[t[[1]]], t[[1, 1]]}}];
        args = Join[args, {t[[1, 2]], t[[1, 2]]}],
      True, Throw[$Failed, BesselHyper`fail]]],
    If[Head[f] === Times, List @@ f, {f}]];
  If[fs === {} || Length[fs] > 2, Throw[$Failed, BesselHyper`fail]];
  w = args[[1]];
  If[Length[args] == 2 && Simplify[args[[2]] - w] =!= 0, Throw[$Failed, BesselHyper`fail]];
  bp = BesselHyper`ArgForm[w, x];
  BesselHyper`Int[c, m, fs, w, bp[[1]], bp[[2]], x]];

BesselHyper`Integrate[f_, x_Symbol] := Catch[
  If[Head[f] === Plus, Total[BesselHyper`Term[#, x] & /@ List @@ f], BesselHyper`Term[f, x]],
  BesselHyper`fail];
