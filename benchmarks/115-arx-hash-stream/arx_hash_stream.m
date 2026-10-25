(* ==========================================================================
   Experiment 115 -- ARX hashing and stream ciphers (group F)
   ==========================================================================
   WHAT IT MEASURES.  The add-rotate-xor primitives on 32-bit words, written
   in the language: SHA-256 (FIPS 180-4: padding, message schedule, 64-round
   compression), HMAC-SHA256 (RFC 2104), and the ChaCha20 block function
   (RFC 7539: quarter rounds, 10 double rounds, little-endian serialisation).
   Words are integers in [0, 2^32); addition is Mod[., 2^32].

   Two spellings of the SAME algorithm:
     * idiomatic   -- BitXor / BitAnd / BitOr / BitShiftLeft / BitShiftRight /
                      BitNot, as a WL programmer writes it.  Mathilda has no
                      bitwise heads: these rows are benchIf'd (ABSENT).
     * (no-bitwise) -- 32-bit XOR / AND as four 8-bit lookups in 256x256
                      tables built once (untimed) from IntegerDigits[., 2, 8];
                      bytes split off with Quotient/Mod; rotations and shifts
                      by Quotient/Mod with powers of two; NOT x = 2^32 - 1 - x.
                      Python and PARI run the identical variant.

   The native primitive ("prim:") is Hash[string, "SHA256"]; Mathilda's Hash
   returns unevaluated for that form, so the case is probed (SKIP).

   Correctness anchors, verified by the checks:
     SHA-256("abc")         = ba7816bf...f20015ad      -> CHECK 997922219
     HMAC-SHA256 RFC 4231 #2 ("Jefe", "what do ya want for nothing?")
                            = 5bdcc146...64ec3843
     ChaCha20 block, RFC 7539 2.3.2 (key 00..1f, nonce 000000090000004a00000000,
                              counter 1) = 10f1e7e4 d13b5915 500fdd1f ...
   Digests are folded as FromDigits[., 2^32 or 256] mod $CheckMod.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"BitXor", "BitAnd", "BitOr", "BitShiftLeft", "BitShiftRight", "BitNot",
         "Hash", "IntegerDigits", "FromDigits", "ToCharacterCode"}];

hexWords[s_] := FromDigits[#, 16] & /@ StringSplit[s];
hexBytes[s_] := FromDigits[#, 16] & /@ StringPartition[s, 2];

$M32 = 2^32;
shaK = hexWords["428a2f98 71374491 b5c0fbcf e9b5dba5 3956c25b 59f111f1 923f82a4 ab1c5ed5 d807aa98 12835b01 243185be 550c7dc3 72be5d74 80deb1fe 9bdc06a7 c19bf174 e49b69c1 efbe4786 0fc19dc6 240ca1cc 2de92c6f 4a7484aa 5cb0a9dc 76f988da 983e5152 a831c66d b00327c8 bf597fc7 c6e00bf3 d5a79147 06ca6351 14292967 27b70a85 2e1b2138 4d2c6dfc 53380d13 650a7354 766a0abb 81c2c92e 92722c85 a2bfe8a1 a81a664b c24b8b70 c76c51a3 d192e819 d6990624 f40e3585 106aa070 19a4c116 1e376c08 2748774c 34b0bcb5 391c0cb3 4ed8aa4a 5b9cca4f 682e6ff3 748f82ee 78a5636f 84c87814 8cc70208 90befffa a4506ceb bef9a3f7 c67178f2"];
shaH0 = hexWords["6a09e667 bb67ae85 3c6ef372 a54ff53a 510e527f 9b05688c 1f83d9ab 5be0cd19"];

(* FIPS 180-4 padding: 0x80, zeros, 64-bit big-endian bit length *)
shaPad[msg_] := Join[msg, {128}, ConstantArray[0, Mod[55 - Length[msg], 64]],
                     IntegerDigits[8 Length[msg], 256, 8]];

(* printable-ASCII messages from the LCG, so the bytes are also a WL string *)
asciiMsg[seed_, len_] := smallList[seed, len, 94] + 32;
msgs = Table[asciiMsg[1500 + j, 256], {j, 1, 8}];       (* 8 x 256 bytes *)
abc = ToCharacterCode["abc"];
hmacKey = ToCharacterCode["Jefe"];
hmacMsg = ToCharacterCode["what do ya want for nothing?"];
chachaKey = Range[0, 31];
chachaNonce = hexBytes["000000090000004a00000000"];

(* =================== idiomatic (BitXor & friends) ======================= *)

rotr[x_, n_] := BitOr[BitShiftRight[x, n], BitAnd[BitShiftLeft[x, 32 - n], $M32 - 1]];
rotl[x_, n_] := rotr[x, 32 - n];

shaCompress[h_, blk_] := Module[{w = Join[FromDigits[#, 256] & /@ Partition[blk, 4],
                                         ConstantArray[0, 48]],
                                 a, b, c, d, e, f, g, hh, s0, s1, t1, t2},
  Do[s0 = BitXor[rotr[w[[t - 15]], 7], rotr[w[[t - 15]], 18], BitShiftRight[w[[t - 15]], 3]];
     s1 = BitXor[rotr[w[[t - 2]], 17], rotr[w[[t - 2]], 19], BitShiftRight[w[[t - 2]], 10]];
     w[[t]] = Mod[w[[t - 16]] + s0 + w[[t - 7]] + s1, $M32], {t, 17, 64}];
  {a, b, c, d, e, f, g, hh} = h;
  Do[t1 = Mod[hh + BitXor[rotr[e, 6], rotr[e, 11], rotr[e, 25]] +
              BitXor[BitAnd[e, f], BitAnd[BitXor[e, $M32 - 1], g]] + shaK[[t]] + w[[t]], $M32];
     t2 = Mod[BitXor[rotr[a, 2], rotr[a, 13], rotr[a, 22]] +
              BitXor[BitAnd[a, b], BitAnd[a, c], BitAnd[b, c]], $M32];
     hh = g; g = f; f = e; e = Mod[d + t1, $M32];
     d = c; c = b; b = a; a = Mod[t1 + t2, $M32], {t, 1, 64}];
  Mod[h + {a, b, c, d, e, f, g, hh}, $M32]];

sha256[msg_] := FromDigits[Fold[shaCompress, shaH0, Partition[shaPad[msg], 64]], $M32];

hmacSha256[key_, msg_] := Module[{k = PadRight[key, 64]},
  sha256[Join[BitXor[k, 92], IntegerDigits[sha256[Join[BitXor[k, 54], msg]], 256, 32]]]];

quarter[s_, {i_, j_, k_, l_}] := Module[{a = s[[i]], b = s[[j]], c = s[[k]], d = s[[l]], r = s},
  a = Mod[a + b, $M32]; d = rotl[BitXor[d, a], 16];
  c = Mod[c + d, $M32]; b = rotl[BitXor[b, c], 12];
  a = Mod[a + b, $M32]; d = rotl[BitXor[d, a], 8];
  c = Mod[c + d, $M32]; b = rotl[BitXor[b, c], 7];
  r[[{i, j, k, l}]] = {a, b, c, d}; r];

(* ---- shared by both variants: state layout and the round schedule ---- *)
le32[bytes_] := FromDigits[Reverse[#], 256] & /@ Partition[bytes, 4];
chachaInit[key_, ctr_, nonce_] :=
  Join[{1634760805, 857760878, 2036477234, 1797285236}, le32[key], {ctr}, le32[nonce]];
chachaRounds = {{1, 5, 9, 13}, {2, 6, 10, 14}, {3, 7, 11, 15}, {4, 8, 12, 16},
                {1, 6, 11, 16}, {2, 7, 12, 13}, {3, 8, 9, 14}, {4, 5, 10, 15}};

chachaBlock[key_, ctr_, nonce_] := Module[{s0 = chachaInit[key, ctr, nonce], s},
  s = s0;
  Do[s = Fold[quarter, s, chachaRounds], {10}];
  Flatten[Reverse[IntegerDigits[#, 256, 4]] & /@ Mod[s + s0, $M32]]];

(* =================== (no-bitwise): byte tables + arithmetic ============= *)

xorTab = Table[FromDigits[Mod[IntegerDigits[a, 2, 8] + IntegerDigits[b, 2, 8], 2], 2],
               {a, 0, 255}, {b, 0, 255}];
andTab = Table[FromDigits[IntegerDigits[a, 2, 8] IntegerDigits[b, 2, 8], 2],
               {a, 0, 255}, {b, 0, 255}];
bytePow = {1, 256, 65536, 16777216};

x32[a_, b_] := bytePow . Extract[xorTab,
  Transpose[{Mod[Quotient[a, bytePow], 256] + 1, Mod[Quotient[b, bytePow], 256] + 1}]];
a32[a_, b_] := bytePow . Extract[andTab,
  Transpose[{Mod[Quotient[a, bytePow], 256] + 1, Mod[Quotient[b, bytePow], 256] + 1}]];
n32[a_] := $M32 - 1 - a;
rotrN[x_, n_] := Quotient[x, 2^n] + Mod[x, 2^n] 2^(32 - n);
rotlN[x_, n_] := rotrN[x, 32 - n];

shaCompressN[h_, blk_] := Module[{w = Join[FromDigits[#, 256] & /@ Partition[blk, 4],
                                          ConstantArray[0, 48]],
                                  a, b, c, d, e, f, g, hh, s0, s1, t1, t2},
  Do[s0 = x32[x32[rotrN[w[[t - 15]], 7], rotrN[w[[t - 15]], 18]], Quotient[w[[t - 15]], 8]];
     s1 = x32[x32[rotrN[w[[t - 2]], 17], rotrN[w[[t - 2]], 19]], Quotient[w[[t - 2]], 1024]];
     w[[t]] = Mod[w[[t - 16]] + s0 + w[[t - 7]] + s1, $M32], {t, 17, 64}];
  {a, b, c, d, e, f, g, hh} = h;
  Do[t1 = Mod[hh + x32[x32[rotrN[e, 6], rotrN[e, 11]], rotrN[e, 25]] +
              x32[a32[e, f], a32[n32[e], g]] + shaK[[t]] + w[[t]], $M32];
     t2 = Mod[x32[x32[rotrN[a, 2], rotrN[a, 13]], rotrN[a, 22]] +
              x32[x32[a32[a, b], a32[a, c]], a32[b, c]], $M32];
     hh = g; g = f; f = e; e = Mod[d + t1, $M32];
     d = c; c = b; b = a; a = Mod[t1 + t2, $M32], {t, 1, 64}];
  Mod[h + {a, b, c, d, e, f, g, hh}, $M32]];

sha256N[msg_] := FromDigits[Fold[shaCompressN, shaH0, Partition[shaPad[msg], 64]], $M32];

hmacSha256N[key_, msg_] := Module[{k = PadRight[key, 64]},
  sha256N[Join[xorTab[[k + 1, 93]],
               IntegerDigits[sha256N[Join[xorTab[[k + 1, 55]], msg]], 256, 32]]]];

quarterN[s_, {i_, j_, k_, l_}] := Module[{a = s[[i]], b = s[[j]], c = s[[k]], d = s[[l]], r = s},
  a = Mod[a + b, $M32]; d = rotlN[x32[d, a], 16];
  c = Mod[c + d, $M32]; b = rotlN[x32[b, c], 12];
  a = Mod[a + b, $M32]; d = rotlN[x32[d, a], 8];
  c = Mod[c + d, $M32]; b = rotlN[x32[b, c], 7];
  r[[{i, j, k, l}]] = {a, b, c, d}; r];

chachaBlockN[key_, ctr_, nonce_] := Module[{s0 = chachaInit[key, ctr, nonce], s},
  s = s0;
  Do[s = Fold[quarterN, s, chachaRounds], {10}];
  Flatten[Reverse[IntegerDigits[#, 256, 4]] & /@ Mod[s + s0, $M32]]];

(* =================== cases ============================================== *)

(* ---- SHA-256("abc"), the FIPS 180-4 one-block vector ---- *)
benchIf["app: SHA-256 abc x40", "BitXor", Do[dAbc = sha256[abc], {40}];];
checkIf["app: SHA-256 abc x40", "BitXor", Mod[dAbc, $CheckMod]];

bench["app: SHA-256 abc x40 (no-bitwise)", Do[dAbcN = sha256N[abc], {40}];];
check["app: SHA-256 abc x40 (no-bitwise)", Mod[dAbcN, $CheckMod]];

(* ---- SHA-256 of 8 LCG messages of 256 bytes (5 blocks each) ---- *)
benchIf["app: SHA-256 8x256-byte messages", "BitXor", dMsgs = sha256 /@ msgs;];
checkIf["app: SHA-256 8x256-byte messages", "BitXor", Mod[Total[dMsgs], $CheckMod]];

bench["app: SHA-256 8x256-byte messages (no-bitwise)", dMsgsN = sha256N /@ msgs;];
check["app: SHA-256 8x256-byte messages (no-bitwise)", Mod[Total[dMsgsN], $CheckMod]];

(* ---- HMAC-SHA256, RFC 4231 test case 2 ---- *)
benchIf["app: HMAC-SHA256 RFC4231 #2 x12", "BitXor", Do[tag = hmacSha256[hmacKey, hmacMsg], {12}];];
checkIf["app: HMAC-SHA256 RFC4231 #2 x12", "BitXor", Mod[tag, $CheckMod]];

bench["app: HMAC-SHA256 RFC4231 #2 x12 (no-bitwise)", Do[tagN = hmacSha256N[hmacKey, hmacMsg], {12}];];
check["app: HMAC-SHA256 RFC4231 #2 x12 (no-bitwise)", Mod[tagN, $CheckMod]];

(* ---- ChaCha20 keystream: RFC 7539 key/nonce, counters 1..48 ---- *)
benchIf["app: ChaCha20 keystream 48 blocks", "BitXor",
  ks = Table[chachaBlock[chachaKey, ctr, chachaNonce], {ctr, 1, 48}];];
checkIf["app: ChaCha20 keystream 48 blocks", "BitXor",
  Mod[Total[FromDigits[#, 256] & /@ ks], $CheckMod]];

bench["app: ChaCha20 keystream 48 blocks (no-bitwise)",
  ksN = Table[chachaBlockN[chachaKey, ctr, chachaNonce], {ctr, 1, 48}];];
check["app: ChaCha20 keystream 48 blocks (no-bitwise)",
  Mod[Total[FromDigits[#, 256] & /@ ksN], $CheckMod]];

(* ---- the native primitive: Hash[string, "SHA256"] on the same messages, x100 ---- *)
hashOkQ = IntegerQ[Hash["abc", "SHA256"]];
msgStrs = FromCharacterCode /@ msgs;
If[hashOkQ,
  bench["prim: Hash SHA256 8x256-byte messages x100", Do[dPrim = Hash[#, "SHA256"] & /@ msgStrs, {100}];];
  check["prim: Hash SHA256 8x256-byte messages x100", Mod[Total[dPrim], $CheckMod]],
  Print["SKIP\tprim: Hash SHA256 8x256-byte messages x100\tHash[SHA256]"]];
