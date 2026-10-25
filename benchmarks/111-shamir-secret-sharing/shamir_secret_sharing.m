(* ==========================================================================
   Experiment 111 -- Shamir secret sharing (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  (t, n) threshold sharing over the prime field F_p, p a
   256-bit prime: a secret s is the constant term of a random degree t-1
   polynomial f, share i is (i, f(i) mod p), and any t shares recover
   s = f(0) by Lagrange interpolation at 0.

   "app:" cases are the protocol written in the language: Horner evaluation
   mod p (split), Lagrange-at-0 with one PowerMod[d, -1, p] inverse per term
   (reconstruct), and the batch form a share server uses -- the Lagrange
   weights depend only on the share abscissae, so they are computed once and
   every further secret is a dot product mod p.

   "prim:" is the WL route a user reaches for first: InterpolatingPolynomial
   over the t shares, evaluated at 0, then reduced mod p.  WL has no
   finite-field interpolation, so this runs over Q (rational coefficients,
   huge denominators) and is expected to be slow -- that is the finding the
   case exists to measure.  Python's counterpart is the same Q-route
   (sympy.polys.polyfuncs.interpolate, sympy has no GF(p) interpolation);
   PARI's is polinterpolate on Mod(y, p) entries, which stays in F_p.

   Inputs come from crypto_data.m's shared LCG; every check is a sum of
   field elements folded mod $CheckMod.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod", "InterpolatingPolynomial", "Fold"}];

$P = primeOf[1111, 256];
$T = 50;
$N = 200;

(* Coefficients {s, a1, ..., a_{t-1}} of the k-th secret's polynomial. *)
polyOf[k_] := Table[Mod[bigOf[100000 k + j, 256], $P], {j, 0, $T - 1}];

(* f(x) mod p by Horner, highest coefficient first. *)
hornerMod[cs_, x_] := Fold[Mod[#1 x + #2, $P] &, 0, Reverse[cs]];

split[cs_] := Table[{i, hornerMod[cs, i]}, {i, 1, $N}];

(* f(0) from t shares {x_i, y_i}: sum y_i prod_{j != i} x_j / (x_j - x_i). *)
lagrange0[sh_] := Module[{xs = sh[[All, 1]], num, den},
  Mod[Sum[
      num = 1; den = 1;
      Do[If[j != i,
            num = Mod[num xs[[j]], $P];
            den = Mod[den (xs[[j]] - xs[[i]]), $P]], {j, Length[xs]}];
      Mod[sh[[i, 2]] num PowerMod[den, -1, $P], $P],
      {i, Length[xs]}], $P]];

(* The share subset every reconstruction uses: shares 3, 7, 11, ... (t of them). *)
pick[sh_] := sh[[Range[3, 3 + 4 ($T - 1), 4]]];

polys = Table[polyOf[k], {k, 1, 20}];

(* ---- app: split ------------------------------------------------------- *)

bench["app: split 256-bit n=200 t=50 x20", shares = split /@ polys;];
check["app: split 256-bit n=200 t=50 x20", Mod[Total[Flatten[shares[[All, All, 2]]]], $CheckMod]];

(* ---- app: reconstruct ------------------------------------------------- *)

bench["app: reconstruct Lagrange t=50 x20", recs = lagrange0[pick[#]] & /@ shares;];
check["app: reconstruct Lagrange t=50 x20", Mod[Total[recs], $CheckMod]];

(* ---- app: batch reconstruct with precomputed weights ------------------- *)

(* Weights w_i = prod_{j != i} x_j / (x_j - x_i) mod p for the fixed abscissae;
   each secret is then sum w_i y_i mod p.  200 secrets, all split from
   polys[[k]] shifted by k so they differ. *)
batchShares = Table[{#[[1]], Mod[#[[2]] + k, $P]} & /@ pick[shares[[Mod[k, 20] + 1]]], {k, 1, 2000}];

bench["app: batch reconstruct t=50 x2000",
  xs0 = batchShares[[1, All, 1]];
  w = Table[Mod[Product[If[j != i, xs0[[j]], 1], {j, $T}]
            PowerMod[Product[If[j != i, xs0[[j]] - xs0[[i]], 1], {j, $T}], -1, $P], $P], {i, $T}];
  brec = Mod[#[[All, 2]] . w, $P] & /@ batchShares;];
check["app: batch reconstruct t=50 x2000", Mod[Total[brec], $CheckMod]];

(* ---- prim: InterpolatingPolynomial over Q, then mod p ------------------ *)

(* The rational value r = f_Q(0) agrees with s mod p because every
   denominator is a product of share-abscissa differences, coprime to p. *)
qmod[r_] := Mod[Numerator[r] PowerMod[Denominator[r], -1, $P], $P];

bench["prim: InterpolatingPolynomial t=50 mod p x20",
  irecs = (qmod[InterpolatingPolynomial[pick[#], z] /. z -> 0]) & /@ shares;];
check["prim: InterpolatingPolynomial t=50 mod p x20", Mod[Total[irecs], $CheckMod]];

(* All n = 200 shares of secrets 1..4 (over-determined: the degree-49 polynomial
   is recovered exactly), to show how the Q-route scales with the point count. *)
benchOnce["prim: InterpolatingPolynomial n=200 mod p x4",
  Do[irec200 = qmod[InterpolatingPolynomial[shares[[r]], z] /. z -> 0], {r, 4}];];
check["prim: InterpolatingPolynomial n=200 mod p x4", Mod[irec200, $CheckMod]];
