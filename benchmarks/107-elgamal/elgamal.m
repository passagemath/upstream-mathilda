(* ==========================================================================
   Experiment 107 -- ElGamal over the RFC 3526 MODP-2048 group
   (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  Textbook ElGamal in Z_p^* (p = $MODP2048, g = 2),
   written in the language ("app:" cases):

       keygen:   x in [1, p-2],  h = g^x
       encrypt:  y in [1, p-2],  (c1, c2) = (g^y, m h^y)
       decrypt:  s = c1^x,  m = c2 s^-1   -- the inverse two ways:
                 PowerMod[s, -1, p] (extended Euclid) and s^(p-2) (Fermat)
       homomorphism:     prod (c1_i, c2_i) decrypts to prod m_i
       re-randomisation: (c1 g^r, c2 h^r) decrypts to the same m

   All of it is 2048-bit modular exponentiation, so the "prim:" case is that
   one kernel batched: PowerMod / gmpy2.powmod / Mod(b, p)^e, x200.

   Exponents, messages, and randomisers all come from crypto_data's bigOf;
   checks are sums (or products) folded mod $CheckMod.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod"}];

p = $MODP2048;
g = 2;

(* ---- prim: the core cost ---------------------------------------------- *)

pmBases = Table[Mod[bigOf[s, 2048], p], {s, 2001, 2200}];
pmExps = Table[bigOf[s, 2048], {s, 2201, 2400}];
bench["prim: PowerMod 2048-bit x200",
      pms = MapThread[PowerMod[#1, #2, p] &, {pmBases, pmExps}];];
check["prim: PowerMod 2048-bit x200", Mod[Total[pms], $CheckMod]];

(* ---- app: ElGamal ----------------------------------------------------- *)

egKeygen[seed_] := Module[{x = Mod[bigOf[seed, 2048], p - 2] + 1},
  {x, PowerMod[g, x, p]}];

egEncrypt[m_, h_, y_] := {PowerMod[g, y, p], Mod[m PowerMod[h, y, p], p]};

egDecryptInv[{c1_, c2_}, x_] := Mod[c2 PowerMod[PowerMod[c1, x, p], -1, p], p];

egDecryptFermat[{c1_, c2_}, x_] :=
  Mod[c2 PowerMod[PowerMod[c1, x, p], p - 2, p], p];

bench["app: ElGamal keygen x50", keys = Table[egKeygen[s], {s, 901, 950}];];
check["app: ElGamal keygen x50", Mod[Total[Flatten[keys]], $CheckMod]];

{privX, pubH} = First[keys];
msgs = Table[Mod[bigOf[s, 2048], p - 1] + 1, {s, 1001, 1050}];
ephs = Table[Mod[bigOf[s, 2048], p - 2] + 1, {s, 1101, 1150}];

bench["app: ElGamal encrypt x50",
      cts = MapThread[egEncrypt[#1, pubH, #2] &, {msgs, ephs}];];
check["app: ElGamal encrypt x50", Mod[Total[Flatten[cts]], $CheckMod]];

bench["app: ElGamal decrypt x50 (PowerMod[s,-1,p])",
      dec1 = egDecryptInv[#, privX] & /@ cts;];
check["app: ElGamal decrypt x50 (PowerMod[s,-1,p])",
      If[dec1 === msgs, Mod[Total[dec1], $CheckMod], -1]];

bench["app: ElGamal decrypt x50 (s^(p-2))",
      dec2 = egDecryptFermat[#, privX] & /@ cts;];
check["app: ElGamal decrypt x50 (s^(p-2))",
      If[dec2 === msgs, Mod[Total[dec2], $CheckMod], -1]];

(* product of the 50 ciphertexts, decrypted; repeated 20 times so the case
   clears the PARI timer's resolution *)
bench["app: homomorphic product of 50 ciphertexts + decrypt x20",
      Do[prodC = {Mod[Times @@ cts[[All, 1]], p], Mod[Times @@ cts[[All, 2]], p]};
         prodM = egDecryptInv[prodC, privX], {20}];];
check["app: homomorphic product of 50 ciphertexts + decrypt x20",
      If[prodM == Mod[Times @@ msgs, p], Mod[prodM, $CheckMod], -1]];

rands = Table[Mod[bigOf[s, 2048], p - 2] + 1, {s, 1201, 1250}];
bench["app: ElGamal re-randomise x50",
      rcts = MapThread[{Mod[#1[[1]] PowerMod[g, #2, p], p],
                        Mod[#1[[2]] PowerMod[pubH, #2, p], p]} &, {cts, rands}];];
check["app: ElGamal re-randomise x50",
      If[(egDecryptInv[#, privX] & /@ rcts) === msgs,
         Mod[Total[Flatten[rcts]], $CheckMod], -1]];
