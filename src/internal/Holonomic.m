(* Holonomic.m -- indefinite integrals of
     x^alpha L(x) * (monomials in Bessel / Airy / trig kernels)
   by a rational ansatz in the kernels' finite-dimensional differential module.

   Loaded lazily by the indefinite-integration stage "Holonomic"
   (src/calculus/integrate.c).  All symbols are qualified in Holonomic`.

   MODULE.  Each kernel family spans a 2-dimensional space closed under d/dx
   with Laurent-polynomial coefficients (z = b x):
     Bessel J, Y:  p = Z_nu(z), q = Z_(nu+1)(z)
                   p' = (nu/x) p - b q,        q' = b p - ((nu+1)/x) q
     Bessel I:     p' = (nu/x) p + b q,        q' = b p - ((nu+1)/x) q
     Bessel K:     p' = (nu/x) p - b q,        q' = -b p - ((nu+1)/x) q
     Sin, Cos:     s' = b c,                   c' = -b s
     Airy Ai, Bi:  p = Ai(z), q = Ai'(z):      p' = b q,   q' = b^2 x p
   Orders nu0 + k (k a non-negative integer) of one family are rewritten in the
   basis {Z_nu0, Z_(nu0+1)} by the three-term recurrence, so J_nu J_(nu+1) is a
   symmetric power rather than two unrelated kernels.  Products of kernels then
   live in the space of monomials of a fixed multidegree in the family symbols,
   which d/dx maps into itself.

   ANSATZ.  For an integrand  x^alpha Sum_j Sum_M f_(M,j) x^j M  (alpha the
   common non-integer part of the exponents, M the monomials) seek
     F = x^alpha Sum_M Sum_(k = lo..hi) c_(M,k) x^k M,
   and equate coefficients of x^j M in dF/dx - f: a linear system in the
   c_(M,k) with coefficients rational in the parameters.  Any solution is an
   exact antiderivative (the derivation rules are identities), and it is
   re-checked symbolically before it is returned.

   WINDOW.  At x = 0 the module matrix has a simple pole whose residue is
   diagonal in the monomial basis (e_M = Sum of the nu, -(nu+1) entries), so a
   homogeneous solution starts at x^(-alpha - e_M): lo is the least of
   jmin + 1 and every such exponent that is an integer -- exact.  Infinity is
   an irregular singular point, where a complete degree bound needs Barkatou's
   super-reduction (not implemented); hi is a heuristic, deg f + 1 (or, for an
   integer alpha, D - alpha with D the total kernel degree, because kernel
   combinations such as J_nu^2 + J_(nu+1)^2 ~ 2/(Pi x) decay faster than their
   coefficients grow) plus the dimension.  The window only limits what is
   found: every returned answer is exact.  An integrand whose antiderivative is
   not in the module (x J_1, which needs Struve H) has no solution and
   declines. *)

Holonomic`LinQ[e_, x_] := PolynomialQ[e, x] && Exponent[e, x] == 1 && Coefficient[e, x, 0] === 0;

(* x-dependent kernels of f: {kind, head, order, b} *)
Holonomic`Kernels[f_, x_] := DeleteDuplicates[Join[
  Cases[f, (h : BesselJ | BesselY | BesselI | BesselK)[nu_, z_] /;
      FreeQ[nu, x] && Holonomic`LinQ[z, x] :> {"bessel", h, nu, Coefficient[z, x, 1]}, {0, Infinity}],
  Cases[f, (h : Sin | Cos)[z_] /; Holonomic`LinQ[z, x] :> {"trig", Sin, 0, Coefficient[z, x, 1]}, {0, Infinity}],
  Cases[f, (h : AiryAi | AiryAiPrime)[z_] /; Holonomic`LinQ[z, x] :> {"airy", AiryAi, 0, Coefficient[z, x, 1]}, {0, Infinity}],
  Cases[f, (h : AiryBi | AiryBiPrime)[z_] /; Holonomic`LinQ[z, x] :> {"airy", AiryBi, 0, Coefficient[z, x, 1]}, {0, Infinity}]]];

(* Group kernels into families: same kind, head and scale, orders differing by
   integers.  Returns {kind, head, nu0, b} per family. *)
Holonomic`Families[ks_] := Module[{fams = {}, k, pos, d},
  Do[
    pos = 0;
    Do[If[fams[[i, 1]] === k[[1]] && fams[[i, 2]] === k[[2]] &&
          Simplify[fams[[i, 4]] - k[[4]]] === 0 &&
          IntegerQ[d = Simplify[k[[3]] - fams[[i, 3]]]],
         pos = i;
         If[d < 0, fams[[i, 3]] = k[[3]]]; Break[]], {i, Length[fams]}];
    If[pos === 0, AppendTo[fams, k]],
    {k, ks}];
  fams];

(* Z_(nu0+k)(b x) in the basis {p, q} = {Z_nu0, Z_(nu0+1)}, k >= 0. *)
Holonomic`Shift[h_, nu0_, b_, k_, p_, q_, x_] := Module[{z0 = p, z1 = q, z2, n},
  If[k == 0, Return[p]];
  Do[
    n = nu0 + j;  (* z0 = Z_n, z1 = Z_(n+1) -> z2 = Z_(n+2) *)
    z2 = Which[
      h === BesselJ || h === BesselY, (2 (n + 1)/(b x)) z1 - z0,
      h === BesselI, z0 - (2 (n + 1)/(b x)) z1,
      h === BesselK, z0 + (2 (n + 1)/(b x)) z1];
    z0 = z1; z1 = Expand[z2],
    {j, 0, k - 2}];
  z1];

(* Derivatives of the two family symbols: {D p, D q}. *)
Holonomic`Deriv[{kind_, h_, nu_, b_}, p_, q_, x_] := Which[
  kind === "bessel" && (h === BesselJ || h === BesselY), {(nu/x) p - b q, b p - ((nu + 1)/x) q},
  kind === "bessel" && h === BesselI, {(nu/x) p + b q, b p - ((nu + 1)/x) q},
  kind === "bessel" && h === BesselK, {(nu/x) p - b q, -b p - ((nu + 1)/x) q},
  kind === "trig", {b q, -b p},
  kind === "airy", {b q, b^2 x p}];

(* The kernel-to-symbol rules.  The pattern names are context-qualified: a
   plain nu_ would capture the user's own nu inside the family order nu0. *)
Holonomic`Rules[fam : {kind_, h_, nu0_, b_}, p_, q_, x_] := Which[
  kind === "bessel",
    {h[Holonomic`nn_, Holonomic`zz_] /; Simplify[Holonomic`zz - b x] === 0 && IntegerQ[Simplify[Holonomic`nn - nu0]] :>
       Holonomic`Shift[h, nu0, b, Simplify[Holonomic`nn - nu0], p, q, x]},
  kind === "trig",
    {Sin[Holonomic`zz_] /; Simplify[Holonomic`zz - b x] === 0 :> p, Cos[Holonomic`zz_] /; Simplify[Holonomic`zz - b x] === 0 :> q},
  kind === "airy" && h === AiryAi,
    {AiryAi[Holonomic`zz_] /; Simplify[Holonomic`zz - b x] === 0 :> p, AiryAiPrime[Holonomic`zz_] /; Simplify[Holonomic`zz - b x] === 0 :> q},
  kind === "airy" && h === AiryBi,
    {AiryBi[Holonomic`zz_] /; Simplify[Holonomic`zz - b x] === 0 :> p, AiryBiPrime[Holonomic`zz_] /; Simplify[Holonomic`zz - b x] === 0 :> q}];

Holonomic`Back[{kind_, h_, nu0_, b_}, x_] := Which[
  kind === "bessel", {h[nu0, b x], h[nu0 + 1, b x]},
  kind === "trig", {Sin[b x], Cos[b x]},
  kind === "airy" && h === AiryAi, {AiryAi[b x], AiryAiPrime[b x]},
  kind === "airy" && h === AiryBi, {AiryBi[b x], AiryBiPrime[b x]}];

(* d/dx of a polynomial in the symbols with coefficients in x. *)
Holonomic`Der[e_, syms_, dsyms_, x_] := D[e, x] + Total[MapThread[D[e, #1] #2 &, {syms, dsyms}]];

(* The diagonal residue at x = 0 of D acting on a monomial with exponent
   vector ev: Sum over symbols of exponent * (coefficient of sym/x in D sym). *)
Holonomic`Residue0[ev_, syms_, dsyms_, x_] :=
  Total[MapThread[#1 Coefficient[Expand[x #3] /. x -> 0, #2] &, {ev, syms, dsyms}]];

(* Split c x^e (or a sum of such) into {e, c} pairs. *)
Holonomic`PowerTerms[r_, x_] := Module[{ts = If[Head[r] === Plus, List @@ r, {r}]},
  Map[Function[t, Module[{e = Exponent[t, x], c},
      c = Simplify[t/x^e];
      If[!FreeQ[c, x], Throw[$Failed, Holonomic`fail]];
      {e, c}]], ts]];

(* Integrate one exponent class: g = x^alpha Sum (rules) with integer offsets. *)
Holonomic`Solve1[terms_, alpha_, syms_, dsyms_, degs_, x_] := Module[
  {mons, js, jmin, jmax, lo, hi, dim, ks, F, DF, eqs, unk, sol, rhs, cr, ev},
  (* the monomials of the multidegree *)
  mons = Holonomic`Monomials[syms, degs];
  dim = Length[mons];
  js = terms[[All, 2]];
  jmin = Min[js]; jmax = Max[js];
  ks = Select[Function[mon, Simplify[-alpha -
           Holonomic`Residue0[(Exponent[mon, #] &) /@ syms, syms, dsyms, x]]] /@ mons, IntegerQ];
  lo = Min[Join[{jmin + 1}, ks]];
  (* At Infinity a combination of kernel monomials can decay faster than its
     coefficients grow (J_nu^2 + J_(nu+1)^2 ~ 2/(Pi x)), so with an integer
     alpha the antiderivative of a decaying integrand can need coefficients up
     to x^(D - alpha), D the total kernel degree. *)
  hi = If[IntegerQ[alpha], Max[jmax + 1, Total[degs] - alpha], jmax + 1] + dim;
  unk = Flatten[Table[Holonomic`c[m, k], {m, dim}, {k, lo, hi}]];
  (* F / x^alpha and d/dx F / x^alpha *)
  F = Sum[Holonomic`c[m, k] x^k mons[[m]], {m, dim}, {k, lo, hi}];
  DF = Sum[Holonomic`c[m, k] ((alpha + k) x^(k - 1) mons[[m]] +
        x^k Holonomic`Der[mons[[m]], syms, dsyms, x]), {m, dim}, {k, lo, hi}];
  rhs = Sum[t[[3]] x^t[[2]] t[[1]], {t, terms}];
  cr = CoefficientRules[Expand[x^(1 - lo + 2) (DF - rhs)], Append[syms, x]];
  eqs = (#[[2]] == 0 &) /@ cr;
  If[eqs === {}, Throw[$Failed, Holonomic`fail]];
  (* a particular solution; unknowns left free are set to 0 *)
  sol = Quiet[Solve[eqs, unk]];
  If[!MatchQ[sol, {{___Rule}, ___}], Throw[$Failed, Holonomic`fail]];
  F = (F /. First[sol]) /. Holonomic`c[_, _] -> 0;
  (* exact re-check of the derivative *)
  If[Simplify[Together[Expand[(alpha F + x Holonomic`Der[F, syms, dsyms, x])/x - rhs]]] =!= 0,
    Throw[$Failed, Holonomic`fail]];
  F];

Holonomic`Monomials[syms_, degs_] := Module[{groups},
  (* syms come in pairs per family; degs is the total degree per family *)
  groups = Table[Table[syms[[2 i - 1]]^a syms[[2 i]]^(degs[[i]] - a), {a, 0, degs[[i]]}],
                 {i, Length[degs]}];
  Flatten[Outer[Times, Sequence @@ groups]]];

(* WRONSKIAN ROUTE.  For two solutions y1, y2 of one second-order equation,
   t = y2/y1 has dt/dx = W/y1^2 (W = y1 y2' - y1' y2), so an integrand
   W(x) g(y1, y2) with g homogeneous of degree -2 integrates as
     INT W g(y1, y2) dx = INT g(1, t) dt,   t = y2/y1.
   Pairs: AiryAi/AiryBi (W = b/Pi), BesselJ/BesselY (W = 2/(Pi x)) and
   BesselI/BesselK (W = -1/x) of one order and argument b x.  This closes
   INT dx/AiryAi[x]^2 = Pi AiryBi/AiryAi and INT dx/(x BesselJ[nu, x]^2), which
   the module ansatz cannot (negative powers). *)
Holonomic`WPairs[f_, x_] := DeleteDuplicates[Join[
  Cases[f, (AiryAi | AiryBi)[z_] /; Holonomic`LinQ[z, x] :>
     {AiryAi[z], AiryBi[z], Coefficient[z, x, 1]/Pi}, {0, Infinity}],
  Cases[f, (BesselJ | BesselY)[nu_, z_] /; FreeQ[nu, x] && Holonomic`LinQ[z, x] :>
     {BesselJ[nu, z], BesselY[nu, z], 2/(Pi x)}, {0, Infinity}],
  Cases[f, (BesselI | BesselK)[nu_, z_] /; FreeQ[nu, x] && Holonomic`LinQ[z, x] :>
     {BesselI[nu, z], BesselK[nu, z], -1/x}, {0, Infinity}]]];

Holonomic`Wronskian[f_, x_Symbol] := Module[{prs, y1, y2, w, h, g, F},
  prs = Holonomic`WPairs[f, x];
  If[Length[prs] != 1, Throw[$Failed, Holonomic`fail]];
  {y1, y2, w} = First[prs];
  h = (f /. {y1 -> Holonomic`y1, y2 -> Holonomic`y2})/w;
  If[!FreeQ[h, x], Throw[$Failed, Holonomic`fail]];
  (* homogeneous of degree -2 iff y1^2 h(y1, y1 t) is free of y1; PowerExpand
     is legitimate in this formal identity (symbolic exponents) *)
  g = Simplify[PowerExpand[Holonomic`y1^2 (h /. Holonomic`y2 -> Holonomic`y1 Holonomic`t)]];
  If[!FreeQ[g, Holonomic`y1], Throw[$Failed, Holonomic`fail]];
  F = Integrate[g, Holonomic`t];
  If[!FreeQ[F, Integrate], Throw[$Failed, Holonomic`fail]];
  (* t^e -> y2^e y1^-e rather than (y2/y1)^e: d(u^c) = c u^(c-1) u' holds on
     every interval whatever the branch, so the monomial form differentiates
     back to the integrand for a non-integer e too *)
  F /. {Holonomic`t^e_ :> y2^e y1^(-e)} /. Holonomic`t -> y2/y1];

Holonomic`Integrate[f_, x_Symbol] := Module[{r},
  r = Catch[Holonomic`Ansatz[f, x], Holonomic`fail];
  If[r === $Failed, r = Catch[Holonomic`Wronskian[f, x], Holonomic`fail]];
  r];

Holonomic`Ansatz[f_, x_Symbol] := Catch[Module[
  {ks, fams, n, syms, dsyms, rules, back, g, cr, groups, keyed, res = 0, alpha, terms, degs, F},
  ks = Holonomic`Kernels[f, x];
  (* the gate: at least one Bessel or Airy kernel *)
  If[!MemberQ[ks[[All, 1]], "bessel" | "airy"], Throw[$Failed, Holonomic`fail]];
  fams = Holonomic`Families[ks];
  n = Length[fams];
  syms = Flatten[Table[{Holonomic`p[i], Holonomic`q[i]}, {i, n}]];
  dsyms = Flatten[Table[Holonomic`Deriv[fams[[i]], Holonomic`p[i], Holonomic`q[i], x], {i, n}]];
  rules = Flatten[Table[Holonomic`Rules[fams[[i]], Holonomic`p[i], Holonomic`q[i], x], {i, n}]];
  g = Expand[f /. rules];
  (* every x-dependence must now be through the symbols and powers of x *)
  If[!FreeQ[g, BesselJ | BesselY | BesselI | BesselK | AiryAi | AiryBi | AiryAiPrime | AiryBiPrime | Sin | Cos],
    Throw[$Failed, Holonomic`fail]];
  cr = CoefficientRules[g, syms];
  (* a polynomial in the kernels: negative powers (1/AiryAi^2) are the
     Wronskian route's, and the module ansatz has no monomials for them *)
  If[cr === {} || !AllTrue[Flatten[cr[[All, 1]]], NonNegative], Throw[$Failed, Holonomic`fail]];
  (* key each (monomial, power) term by {family multidegree, exponent class} *)
  keyed = Flatten[Map[Function[rule, Module[{ev = rule[[1]], mon, dg},
      mon = Times @@ (syms^ev);
      dg = Table[ev[[2 i - 1]] + ev[[2 i]], {i, n}];
      If[Total[dg] == 0, Throw[$Failed, Holonomic`fail]];
      ({dg, #[[1]], #[[2]], mon} &) /@ Holonomic`PowerTerms[rule[[2]], x]]], cr], 1];
  (* group by multidegree, then by exponent class (integer differences) *)
  groups = GatherBy[keyed, First];
  Do[
    Module[{cls = {}, placed},
      Do[placed = False;
        Do[If[IntegerQ[Simplify[t[[2]] - cls[[i, 1, 2]]]],
             cls[[i]] = Append[cls[[i]], t]; placed = True; Break[]], {i, Length[cls]}];
        If[!placed, AppendTo[cls, {t}]], {t, grp}];
      Do[
        alpha = cl[[1, 2]];
        terms = ({#[[4]], Simplify[#[[2]] - alpha], #[[3]]} &) /@ cl;
        degs = cl[[1, 1]];
        F = Holonomic`Solve1[terms, alpha, syms, dsyms, degs, x];
        res = res + x^alpha F,
        {cl, cls}]],
    {grp, groups}];
  If[res === 0, Throw[$Failed, Holonomic`fail]];
  back = Flatten[Table[Holonomic`Back[fams[[i]], x], {i, n}]];
  res = res /. Thread[syms -> back];
  Collect[Expand[res], back, Together]],
  Holonomic`fail];
