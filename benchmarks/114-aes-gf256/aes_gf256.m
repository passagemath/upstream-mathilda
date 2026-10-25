(* ==========================================================================
   Experiment 114 -- AES-128 and GF(2^8) (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  AES-128 (FIPS-197) written in the language from its
   field arithmetic up: GF(2^8) multiply by repeated xtime, the inverse as
   a^254 (seven squarings), the S-box as inverse + affine map, the key
   expansion, and full block encryption (SubBytes, ShiftRows, MixColumns,
   AddRoundKey).  Bytes are integers 0..255; a block / state is a list of
   16 bytes in FIPS column order (state[r + 4c] = byte r of column c).

   Two spellings of the SAME algorithm:
     * idiomatic   -- BitXor / BitAnd / BitOr / BitShiftLeft / BitShiftRight,
                      as any WL programmer writes it.  Mathilda has no
                      bitwise heads, so these cases are benchIf'd (ABSENT).
     * (no-bitwise) -- a portable variant using only arithmetic and table
                      lookup: byte XOR is a lookup in a 256x256 table built
                      ONCE (untimed) from IntegerDigits[., 2, 8], xtime is
                      If[b >= 128, xor(2 b - 256, 27), 2 b], rotations are
                      Mod/Quotient.  Python and PARI run the same variant, so
                      Mathilda gets a timed, checked row.

   The GF(2)[x] route ("prim:") inverts all 255 nonzero field elements with
   PolynomialExtendedGCD[., x^8+x^4+x^3+x+1, x, Modulus -> 2].  The library
   route ("prim: ... library") is WL's Encrypt with an AES-128 ECB
   SymmetricKey -- absent in Mathilda (benchIf'd).

   Checks: byte lists are folded as FromDigits[bytes, 256] mod $CheckMod.
   The FIPS-197 Appendix C.1 vector (key 000102..0f, plaintext 00112233..ff)
   must give 69c4e0d86a7b0430d8cdb78070b4c55a, i.e. CHECK 413762112.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"BitXor", "BitAnd", "BitOr", "BitShiftLeft", "BitShiftRight",
         "PolynomialExtendedGCD", "Encrypt", "SymmetricKey", "ByteArray",
         "IntegerDigits", "FromDigits"}];

hexBytes[s_] := FromDigits[#, 16] & /@ StringPartition[s, 2];
fold[bytes_] := Mod[FromDigits[bytes, 256], $CheckMod];

fipsKey = Range[0, 15];
fipsPt  = hexBytes["00112233445566778899aabbccddeeff"];

(* ShiftRows as a permutation of state indices: new[r + 4c] = old[r + 4((c + r) mod 4)] *)
shiftIdx = Table[r + 4 Mod[c + r, 4] + 1, {c, 0, 3}, {r, 0, 3}] // Flatten;

(* ECB data: one LCG key, 128 LCG blocks *)
ecbKey = smallList[1401, 16, 256];
ecbBlocks = Partition[smallList[1402, 128*16, 256], 16];

(* =================== idiomatic (BitXor & friends) ======================= *)

xtime[b_] := BitXor[BitAnd[BitShiftLeft[b, 1], 255], If[b >= 128, 27, 0]];

gmul[a_, b_] := Module[{p = 0, x = a, y = b},
  While[y > 0, If[OddQ[y], p = BitXor[p, x]]; x = xtime[x]; y = Quotient[y, 2]];
  p];

(* a^254 = a^-1: a^(2 + 4 + ... + 128) *)
ginv[a_] := Module[{r = 1, s = a}, Do[s = gmul[s, s]; r = gmul[r, s], {7}]; r];

rotl8[b_, k_] := BitAnd[BitOr[BitShiftLeft[b, k], BitShiftRight[b, 8 - k]], 255];

makeSbox[] := Table[With[{b = ginv[a]},
    BitXor[b, rotl8[b, 1], rotl8[b, 2], rotl8[b, 3], rotl8[b, 4], 99]], {a, 0, 255}];

(* 11 round keys of 16 bytes each, from the 44-word schedule *)
expandKey[key_, sb_] := Module[{w = Partition[key, 4], t, rc = 1},
  Do[t = w[[i]];
     If[Mod[i, 4] == 0,
       t = sb[[RotateLeft[t] + 1]];
       t[[1]] = BitXor[t[[1]], rc]; rc = xtime[rc]];
     AppendTo[w, BitXor[w[[i - 3]], t]], {i, 4, 43}];
  Flatten /@ Partition[w, 4]];

mixColumns[s_] := Flatten[Function[col, Module[{t = BitXor @@ col},
    Table[BitXor[col[[r]], t, xtime[BitXor[col[[r]], col[[Mod[r, 4] + 1]]]]], {r, 1, 4}]]] /@
  Partition[s, 4]];

encryptBlock[pt_, rk_, sb_] := Module[{s = BitXor[pt, rk[[1]]]},
  Do[s = BitXor[mixColumns[sb[[s + 1]][[shiftIdx]]], rk[[r]]], {r, 2, 10}];
  BitXor[sb[[s + 1]][[shiftIdx]], rk[[11]]]];

(* =================== (no-bitwise): arithmetic + table lookup ============ *)

xorTab = Table[FromDigits[Mod[IntegerDigits[a, 2, 8] + IntegerDigits[b, 2, 8], 2], 2],
               {a, 0, 255}, {b, 0, 255}];
x8[a_, b_] := xorTab[[a + 1, b + 1]];
x8v[u_, v_] := MapThread[x8, {u, v}];                   (* byte-wise on lists *)

xtimeN[b_] := If[b >= 128, x8[2 b - 256, 27], 2 b];

gmulN[a_, b_] := Module[{p = 0, x = a, y = b},
  While[y > 0, If[OddQ[y], p = x8[p, x]]; x = xtimeN[x]; y = Quotient[y, 2]];
  p];

ginvN[a_] := Module[{r = 1, s = a}, Do[s = gmulN[s, s]; r = gmulN[r, s], {7}]; r];

rotl8N[b_, k_] := Mod[b 2^k, 256] + Quotient[b, 2^(8 - k)];

makeSboxN[] := Table[With[{b = ginvN[a]},
    Fold[x8, b, {rotl8N[b, 1], rotl8N[b, 2], rotl8N[b, 3], rotl8N[b, 4], 99}]], {a, 0, 255}];

expandKeyN[key_, sb_] := Module[{w = Partition[key, 4], t, rc = 1},
  Do[t = w[[i]];
     If[Mod[i, 4] == 0,
       t = sb[[RotateLeft[t] + 1]];
       t[[1]] = x8[t[[1]], rc]; rc = xtimeN[rc]];
     AppendTo[w, x8v[w[[i - 3]], t]], {i, 4, 43}];
  Flatten /@ Partition[w, 4]];

mixColumnsN[s_] := Flatten[Function[col, Module[{t = Fold[x8, col]},
    Table[x8[x8[col[[r]], t], xtimeN[x8[col[[r]], col[[Mod[r, 4] + 1]]]]], {r, 1, 4}]]] /@
  Partition[s, 4]];

encryptBlockN[pt_, rk_, sb_] := Module[{s = x8v[pt, rk[[1]]]},
  Do[s = x8v[mixColumnsN[sb[[s + 1]][[shiftIdx]]], rk[[r]]], {r, 2, 10}];
  x8v[sb[[s + 1]][[shiftIdx]], rk[[11]]]];

(* =================== cases ============================================== *)

(* ---- S-box generation from the field ---- *)

benchIf["app: AES S-box generation x4", "BitXor", Do[sbox = makeSbox[], {4}];];
checkIf["app: AES S-box generation x4", "BitXor", fold[sbox]];

bench["app: AES S-box generation x4 (no-bitwise)", Do[sboxN = makeSboxN[], {4}];];
check["app: AES S-box generation x4 (no-bitwise)", fold[sboxN]];

(* ---- FIPS-197 C.1: key expansion + one block, repeated ---- *)

benchIf["app: AES-128 FIPS-197 C.1 x100", "BitXor",
  Do[fipsCt = encryptBlock[fipsPt, expandKey[fipsKey, sbox], sbox], {100}];];
checkIf["app: AES-128 FIPS-197 C.1 x100", "BitXor", fold[fipsCt]];

bench["app: AES-128 FIPS-197 C.1 x100 (no-bitwise)",
  Do[fipsCtN = encryptBlockN[fipsPt, expandKeyN[fipsKey, sboxN], sboxN], {100}];];
check["app: AES-128 FIPS-197 C.1 x100 (no-bitwise)", fold[fipsCtN]];

(* ---- ECB: 128 LCG blocks under one LCG key ---- *)

benchIf["app: AES-128 ECB encrypt 128 blocks", "BitXor",
  ecbCt = With[{rk = expandKey[ecbKey, sbox]}, encryptBlock[#, rk, sbox] & /@ ecbBlocks];];
checkIf["app: AES-128 ECB encrypt 128 blocks", "BitXor", Mod[Total[fold /@ ecbCt], $CheckMod]];

bench["app: AES-128 ECB encrypt 128 blocks (no-bitwise)",
  ecbCtN = With[{rk = expandKeyN[ecbKey, sboxN]}, encryptBlockN[#, rk, sboxN] & /@ ecbBlocks];];
check["app: AES-128 ECB encrypt 128 blocks (no-bitwise)", Mod[Total[fold /@ ecbCtN], $CheckMod]];

(* ---- GF(2)[x] route: all 255 inverses by the extended Euclidean algorithm ---- *)

aesPoly = x^8 + x^4 + x^3 + x + 1;
toPoly[a_] := IntegerDigits[a, 2, 8] . x^Range[7, 0, -1];
polyInv[a_] := PolynomialExtendedGCD[toPoly[a], aesPoly, x, Modulus -> 2][[2, 1]] /. x -> 2;

bench["prim: GF(2^8) inverse via PolynomialExtendedGCD Modulus->2 x7650",
      Do[invs = polyInv /@ Range[255], {30}];];
check["prim: GF(2^8) inverse via PolynomialExtendedGCD Modulus->2 x7650",
      Mod[Range[255] . invs, $CheckMod]];

(* ---- the library: Encrypt with an AES-128 ECB key, 65536 blocks ---- *)

(* The cipher pads (PKCS#7) to one extra block; only the first 16 N bytes
   are the ECB image of the data. *)
libData = Flatten[ConstantArray[ecbBlocks, 512]];
benchIf["prim: AES-128 ECB library 65536 blocks", "Encrypt",
  libCt = Normal[Encrypt[SymmetricKey[<|"Cipher" -> "AES128", "BlockMode" -> "ECB",
                                        "Key" -> ByteArray[ecbKey]|>],
                         ByteArray[libData]]["Data"]][[;; 16*65536]];];
checkIf["prim: AES-128 ECB library 65536 blocks", "Encrypt",
  Mod[Total[fold /@ Partition[libCt, 16]], $CheckMod]];
