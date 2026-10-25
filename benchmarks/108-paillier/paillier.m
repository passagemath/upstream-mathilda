(* ==========================================================================
   Experiment 108 -- Paillier additively homomorphic encryption
   (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  Paillier with 1024-bit primes (n ~ 2048 bits, n^2 ~
   4096 bits), written in the language ("app:" cases), g = n + 1:

       keygen:   lambda = lcm(p-1, q-1),  mu = L(g^lambda mod n^2)^-1 mod n,
                 L(u) = (u - 1) / n
       encrypt:  c = g^m r^n mod n^2         (r^n mod n^2 dominates)
                 fast form: g^m = 1 + m n (mod n^2) since g = n + 1
       decrypt:  m = L(c^lambda mod n^2) mu mod n
       add:      E(m1) E(m2) = E(m1 + m2)    (product mod n^2)
       scalar:   E(m)^k = E(k m)

   Almost all of it is exponentiation modulo a 4096-bit n^2, so the "prim:"
   case is that kernel batched: PowerMod / gmpy2.powmod / Mod(b, n^2)^e.

   Primes are primeOf[seed, 1024] (setup, untimed); messages, randomisers and
   scalars come from bigOf.  The 1000 ciphertexts summed homomorphically are
   built with the fast form from 20 precomputed r^n values (setup, untimed);
   the sum + decrypt is repeated 10 times to clear the PARI timer's floor.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod", "LCM", "Quotient"}];

(* ---- keys ------------------------------------------------------------- *)

pqs = Table[{primeOf[2 s + 3001, 1024], primeOf[2 s + 3002, 1024]}, {s, 1, 8}];

plL[u_, n_] := Quotient[u - 1, n];

plKeygen[{p_, q_}] := Module[{n = p q, n2, lam, mu},
  n2 = n^2;
  lam = LCM[p - 1, q - 1];
  mu = PowerMod[plL[PowerMod[n + 1, lam, n2], n], -1, n];
  {n, lam, mu}];

bench["app: Paillier keygen x8 (given p, q)", keys = plKeygen /@ pqs;];
check["app: Paillier keygen x8 (given p, q)", Mod[Total[Flatten[keys]], $CheckMod]];

{nn, lam, mu} = First[keys];
nn2 = nn^2;

(* ---- prim: the core cost ---------------------------------------------- *)

pmBases = Table[Mod[bigOf[s, 4096], nn2], {s, 3101, 3150}];
bench["prim: PowerMod 4096-bit modulus x50",
      pms = PowerMod[#, nn, nn2] & /@ pmBases;];
check["prim: PowerMod 4096-bit modulus x50", Mod[Total[pms], $CheckMod]];

(* ---- app: encrypt / decrypt ------------------------------------------- *)

plEncrypt[m_, r_] := Mod[PowerMod[nn + 1, m, nn2] PowerMod[r, nn, nn2], nn2];
plEncryptFast[m_, r_] := Mod[(1 + m nn) PowerMod[r, nn, nn2], nn2];
plDecrypt[c_] := Mod[plL[PowerMod[c, lam, nn2], nn] mu, nn];

msgs = Table[Mod[bigOf[s, 2048], nn], {s, 3201, 3220}];
rs = Table[Mod[bigOf[s, 2048], nn], {s, 3301, 3320}];

bench["app: Paillier encrypt x20", cts = MapThread[plEncrypt, {msgs, rs}];];
check["app: Paillier encrypt x20", Mod[Total[cts], $CheckMod]];

bench["app: Paillier encrypt x20 (g=n+1: 1+mn)",
      ctsF = MapThread[plEncryptFast, {msgs, rs}];];
check["app: Paillier encrypt x20 (g=n+1: 1+mn)", Mod[Total[ctsF], $CheckMod]];

bench["app: Paillier decrypt x20", dec = plDecrypt /@ cts;];
check["app: Paillier decrypt x20", If[dec === msgs, Mod[Total[dec], $CheckMod], -1]];

(* ---- app: homomorphic sum of 1000 ciphertexts ------------------------- *)

rns = PowerMod[#, nn, nn2] & /@ rs;
sumMsgs = Table[bigOf[s, 2000], {s, 4001, 5000}];
sumCts = Table[Mod[(1 + sumMsgs[[i]] nn) rns[[Mod[i, 20] + 1]], nn2], {i, 1, 1000}];

bench["app: homomorphic sum of 1000 ciphertexts + decrypt x10",
      Do[cSum = Fold[Mod[#1 #2, nn2] &, 1, sumCts]; mSum = plDecrypt[cSum], {10}];];
check["app: homomorphic sum of 1000 ciphertexts + decrypt x10",
      If[mSum == Mod[Total[sumMsgs], nn], Mod[mSum, $CheckMod], -1]];

(* ---- app: scalar multiplication E(m)^k = E(k m) ----------------------- *)

scalars = Table[Mod[bigOf[s, 2048], nn], {s, 3401, 3420}];
bench["app: Paillier scalar mult c^k x20 (k < n)",
      scs = MapThread[PowerMod[#1, #2, nn2] &, {cts, scalars}];];
check["app: Paillier scalar mult c^k x20 (k < n)",
      If[plDecrypt[First[scs]] == Mod[First[scalars] First[msgs], nn],
         Mod[Total[scs], $CheckMod], -1]];
