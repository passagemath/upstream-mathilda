(* ==========================================================================
   Experiment 113 -- lattice attacks (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  LLL lattice reduction as the cryptanalyst uses it:

   (a) breaking a Merkle-Hellman knapsack public key with the CJLOSS lattice
       (Coster-Joux-LaMacchia-Odlyzko-Schnorr-Stern, the improved
       Lagarias-Odlyzko attack): a low-density knapsack is built from a
       superincreasing sequence, a message is encrypted as a subset sum, and
       the ciphertext is broken by reducing the (n+1)-dimensional basis
           b_i = (2 e_i, N a_i),   b_{n+1} = (1, ..., 1, N s)
       whose shortest vector is (2 x_i - 1 ..., 0), all entries +-1;
   (b) a hidden-number / approximate-GCD attack: from t+1 noisy multiples
       x_i = p q_i + r_i of a secret 128-bit p (16-bit noise) the lattice
           (2^(rho+1), x_1, ..., x_t),  -x_0 e_i  (i = 1..t)
       has the short vector (q_0 2^(rho+1), q_0 x_i - q_i x_0), giving q_0 and
       then p = Quotient[x_0, q_0];
   (c) plain LLL timing on a square knapsack-type (q-ary) basis
           (q, 0, ..., 0),  (a_i, e_i)        n = 30, 50.

   The legitimate Merkle-Hellman decryption (modular inverse + greedy) is the
   "app:" control: the key holder's cost against the attacker's.

   LLL output is not unique across implementations (Mathilda: exact LLL with
   delta = 3/4; PARI qflll: delta = 0.99; sympy: delta = 3/4), so the checks
   are invariants only: the recovered message / secret integer for (a), (b),
   and |det| of the reduced basis (which must equal q: the transform is
   unimodular) for (c).

   All inputs come from crypto_data.m's shared LCG.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"LatticeReduce", "PowerMod", "FromDigits", "Det"}];

(* ---- input vectors: NOT raw LCG output -------------------------------- *)

(* Raw LCG integers are NOT usable as lattice inputs, because the LCG is
   linear and LLL is precisely the tool that sees linearity:
     - word k of bigOf[s, .] is (A^k s + c_k) mod 2^64, so bigOf[seed + i, .]
       for consecutive i has words in exact arithmetic progression; the
       q-ary basis below then had a shortest vector of squared norm 4;
     - chunks of ONE lcgStates stream are still affine in the start state
       mod 2^64; the 30-dim q-ary basis then reduced to squared norm 120,
       against ~3.8e6 for truly random entries (the Gaussian heuristic).
   mixBigs therefore pushes each LCG word through a modular exponentiation
   g^e mod P (P a (bits+64)-bit prime from primeOf) before truncating to a
   bits-bit integer with its top bit set.  Reduction of the result is then
   indistinguishable from the random case (squared norms ~1e6). *)
mixBigs[seed_, cnt_, bits_] := Module[{p = primeOf[seed, bits + 64]},
  Mod[PowerMod[3, #, p] & /@ lcgStates[seed, cnt], 2^(bits - 1)] + 2^(bits - 1)];

(* ---- (a) Merkle-Hellman knapsack ---------------------------------------- *)

(* {w, q, m, a}: superincreasing w_i = (sum_{j<i} w_j) + r_i with k-bit
   LCG gaps r_i, a prime modulus q > sum w, multiplier m, public key
   a_i = m w_i mod q.  Density n / (k + n): 0.30 at (n, k) = (24, 56),
   0.44 at (32, 40).  At (32, 24), density 0.57, the delta = 3/4 reducers
   (Mathilda, sympy) miss 1 message in 20 that delta = 0.99 (PARI)
   recovers, so the sizes are chosen where every reducer succeeds -- see
   README. *)
mhKey[seed_, n_, k_] := Module[{r, w = {}, acc = 0, wi, q, m},
  r = mixBigs[seed, n, k];
  Do[wi = acc + r[[i]]; acc += wi; AppendTo[w, wi], {i, 1, n}];
  q = primeOf[seed, k + n + 1];
  m = Mod[bigOf[seed + 1000, k + n - 1], q];
  {w, q, m, Mod[m w, q]}];

mhMessage[seed_, n_] := smallList[seed, n, 2];          (* bit vector *)
mhEncrypt[a_, x_] := x . a;

(* The key holder: undo the multiplier, then the greedy superincreasing
   solve, most significant item last. *)
mhDecrypt[{w_, q_, m_, _}, c_] := Module[{t = Mod[c PowerMod[m, -1, q], q],
                                          n = Length[w], x},
  x = ConstantArray[0, n];
  Do[If[t >= w[[i]], x[[i]] = 1; t -= w[[i]]], {i, n, 1, -1}];
  FromDigits[x, 2]];

(* The attacker: CJLOSS basis, LatticeReduce, then scan the reduced rows for
   a +-1 vector with last coordinate 0 whose bits re-encrypt to s.
   Returns the message integer, or -1 on failure. *)
cjlossBasis[a_, s_] := Module[{n = Length[a]},
  Append[Table[Append[2 UnitVector[n, i], n a[[i]]], {i, 1, n}],
         Append[ConstantArray[1, n], n s]]];

cjlossRecover[a_, s_, red_] := Module[{found = -1, v, x},
  Do[v = red[[j]];
     If[found < 0 && Last[v] == 0 && Union[Abs[Most[v]]] === {1},
       x = (1 + Most[v])/2;
       If[x . a == s, found = FromDigits[x, 2],
         x = (1 - Most[v])/2;
         If[x . a == s, found = FromDigits[x, 2]]]],
     {j, 1, Length[red]}];
  found];

mhBreak[a_, s_] := cjlossRecover[a, s, LatticeReduce[cjlossBasis[a, s]]];

key24 = mhKey[11, 24, 56];
key32 = mhKey[13, 32, 40];
msgs24 = Table[mhMessage[2000 + j, 24], {j, 1, 20}];
msgs32 = Table[mhMessage[3000 + j, 32], {j, 1, 12}];
cts24 = mhEncrypt[key24[[4]], #] & /@ msgs24;
cts32 = mhEncrypt[key32[[4]], #] & /@ msgs32;

benchOnce["prim: LatticeReduce CJLOSS knapsack n=24 x20",
      brk24 = mhBreak[key24[[4]], #] & /@ cts24;];
check["prim: LatticeReduce CJLOSS knapsack n=24 x20", Mod[Total[brk24], $CheckMod]];

benchOnce["prim: LatticeReduce CJLOSS knapsack n=32 x12",
      brk32 = mhBreak[key32[[4]], #] & /@ cts32;];
check["prim: LatticeReduce CJLOSS knapsack n=32 x12", Mod[Total[brk32], $CheckMod]];

(* the key holder's decryption of 2000 ciphertexts under key32 *)
msgsDec = Table[mhMessage[5000 + j, 32], {j, 1, 2000}];
ctsDec = mhEncrypt[key32[[4]], #] & /@ msgsDec;
bench["app: Merkle-Hellman decrypt n=32 x2000",
      dec32 = mhDecrypt[key32, #] & /@ ctsDec;];
check["app: Merkle-Hellman decrypt n=32 x2000", Mod[Total[dec32], $CheckMod]];

(* ---- (b) approximate-GCD / hidden-number secret recovery --------------- *)

(* t+1 samples x_i = p q_i + r_i: p a 128-bit prime, q_i 128-bit, r_i in
   [2^15, 2^16).  rho = 16. *)
agcdSamples[seed_, t_] := Module[{p = primeOf[seed, 128]},
  {p, p mixBigs[seed + 10, t + 1, 128] + mixBigs[seed + 500, t + 1, 16]}];

agcdRecover[xs_] := Module[{t = Length[xs] - 1, x0 = First[xs], basis, red, q0},
  basis = Prepend[Table[-x0 UnitVector[t + 1, i + 1], {i, 1, t}],
                  Prepend[Rest[xs], 2^17]];
  red = LatticeReduce[basis];
  q0 = Abs[red[[1, 1]]]/2^17;
  If[IntegerQ[q0] && q0 > 0, Quotient[x0, q0], -1]];

agcd = Table[agcdSamples[7000 + 1000 j, 8], {j, 1, 40}];
benchOnce["app: approximate-GCD secret recovery dim=9 x40",
      ps = agcdRecover[#[[2]]] & /@ agcd;];
check["app: approximate-GCD secret recovery dim=9 x40",
      Mod[Total[ps], $CheckMod] + If[ps === (First /@ agcd), 0, 1]];

(* ---- (c) plain LLL on a square knapsack-type basis --------------------- *)

qaryBasis[seed_, n_, bits_] := Module[{q = primeOf[seed, bits]},
  Prepend[MapThread[Prepend, {IdentityMatrix[n - 1], Mod[mixBigs[seed + 1, n - 1, bits], q]}],
          Prepend[ConstantArray[0, n - 1], q]]];

qa30 = qaryBasis[200000, 30, 300];
qa50 = qaryBasis[300000, 50, 200];

benchOnce["prim: LatticeReduce q-ary basis n=30 300-bit", red30 = LatticeReduce[qa30];];
check["prim: LatticeReduce q-ary basis n=30 300-bit", Mod[Abs[Det[red30]], $CheckMod]];

benchOnce["prim: LatticeReduce q-ary basis n=50 200-bit", red50 = LatticeReduce[qa50];];
check["prim: LatticeReduce q-ary basis n=50 200-bit", Mod[Abs[Det[red50]], $CheckMod]];
