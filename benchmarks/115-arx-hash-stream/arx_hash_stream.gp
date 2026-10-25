\\ Experiment 115 -- ARX hashing and stream ciphers (group F).  See arx_hash_stream.m.
\\
\\ PARI/GP baseline.  SHA-256, HMAC-SHA256 and the ChaCha20 block function in
\\ GP, twice: idiomatically with bitxor / bitand / bitor / shift, and as the
\\ "(no-bitwise)" variant (32-bit XOR/AND as four 8-bit lookups in 256x256
\\ matrices built untimed from base-2 digits; rotations/shifts by \ and %) --
\\ the identical algorithm to the .m / .py halves.  PARI has no hash
\\ primitive: the native case is skip()'d.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["bitxor", "bitand", "bitor", "shift", "binary", "fromdigits"]);

hexWords(s) = apply(t -> eval(Str("0x", t)), strsplit(s, " "));
hexBytes(s) = my(c = Vec(s)); vector(#c / 2, i, eval(Str("0x", c[2 * i - 1], c[2 * i])));
toDigits(v, b, n) = vector(n, i, (v \ b^(n - i)) % b);

M32 = 2^32;
SHAK = hexWords("428a2f98 71374491 b5c0fbcf e9b5dba5 3956c25b 59f111f1 923f82a4 ab1c5ed5 d807aa98 12835b01 243185be 550c7dc3 72be5d74 80deb1fe 9bdc06a7 c19bf174 e49b69c1 efbe4786 0fc19dc6 240ca1cc 2de92c6f 4a7484aa 5cb0a9dc 76f988da 983e5152 a831c66d b00327c8 bf597fc7 c6e00bf3 d5a79147 06ca6351 14292967 27b70a85 2e1b2138 4d2c6dfc 53380d13 650a7354 766a0abb 81c2c92e 92722c85 a2bfe8a1 a81a664b c24b8b70 c76c51a3 d192e819 d6990624 f40e3585 106aa070 19a4c116 1e376c08 2748774c 34b0bcb5 391c0cb3 4ed8aa4a 5b9cca4f 682e6ff3 748f82ee 78a5636f 84c87814 8cc70208 90befffa a4506ceb bef9a3f7 c67178f2");
SHAH0 = hexWords("6a09e667 bb67ae85 3c6ef372 a54ff53a 510e527f 9b05688c 1f83d9ab 5be0cd19");

shaPad(msg) = concat([msg, [128], vector((55 - #msg) % 64), toDigits(8 * #msg, 256, 8)]);

msgs = vector(8, j, apply(c -> c + 32, smallList(1500 + j, 256, 94)));
abc = Vec(Vecsmall("abc"));
hmacKey = Vec(Vecsmall("Jefe"));
hmacMsg = Vec(Vecsmall("what do ya want for nothing?"));
chachaKey = vector(32, i, i - 1);
chachaNonce = hexBytes("000000090000004a00000000");

\\ =================== idiomatic ============================================

rotr(x, n) = bitor(shift(x, -n), bitand(shift(x, 32 - n), M32 - 1));
rotl(x, n) = rotr(x, 32 - n);

shaCompress(h, blk) =
{
  my(w = concat(vector(16, i, fromdigits(blk[4 * i - 3 .. 4 * i], 256)), vector(48)), a, b, c, d, e, f, g, hh, s0, s1, t1, t2);
  for (t = 17, 64,
    s0 = bitxor(bitxor(rotr(w[t - 15], 7), rotr(w[t - 15], 18)), shift(w[t - 15], -3));
    s1 = bitxor(bitxor(rotr(w[t - 2], 17), rotr(w[t - 2], 19)), shift(w[t - 2], -10));
    w[t] = (w[t - 16] + s0 + w[t - 7] + s1) % M32);
  [a, b, c, d, e, f, g, hh] = h;
  for (t = 1, 64,
    t1 = (hh + bitxor(bitxor(rotr(e, 6), rotr(e, 11)), rotr(e, 25)) + bitxor(bitand(e, f), bitand(bitxor(e, M32 - 1), g)) + SHAK[t] + w[t]) % M32;
    t2 = (bitxor(bitxor(rotr(a, 2), rotr(a, 13)), rotr(a, 22)) + bitxor(bitxor(bitand(a, b), bitand(a, c)), bitand(b, c))) % M32;
    hh = g; g = f; f = e; e = (d + t1) % M32;
    d = c; c = b; b = a; a = (t1 + t2) % M32);
  (h + [a, b, c, d, e, f, g, hh]) % M32;
}

sha256(msg) =
{
  my(p = shaPad(msg), h = SHAH0);
  for (i = 0, #p / 64 - 1, h = shaCompress(h, p[64 * i + 1 .. 64 * i + 64]));
  fromdigits(h, M32);
}

hmacSha256(key, msg) =
{
  my(k = concat(key, vector(64 - #key)));
  sha256(concat(apply(b -> bitxor(b, 92), k), toDigits(sha256(concat(apply(b -> bitxor(b, 54), k), msg)), 256, 32)));
}

le32(bs) = vector(#bs / 4, i, fromdigits(Vecrev(bs[4 * i - 3 .. 4 * i]), 256));
chachaInit(key, ctr, nonce) = concat([[1634760805, 857760878, 2036477234, 1797285236], le32(key), [ctr], le32(nonce)]);
CHACHAROUNDS = [[1, 5, 9, 13], [2, 6, 10, 14], [3, 7, 11, 15], [4, 8, 12, 16], [1, 6, 11, 16], [2, 7, 12, 13], [3, 8, 9, 14], [4, 5, 10, 15]];

quarter(s, q) =
{
  my(a = s[q[1]], b = s[q[2]], c = s[q[3]], d = s[q[4]]);
  a = (a + b) % M32; d = rotl(bitxor(d, a), 16);
  c = (c + d) % M32; b = rotl(bitxor(b, c), 12);
  a = (a + b) % M32; d = rotl(bitxor(d, a), 8);
  c = (c + d) % M32; b = rotl(bitxor(b, c), 7);
  s[q[1]] = a; s[q[2]] = b; s[q[3]] = c; s[q[4]] = d;
  s;
}

chachaBlock(key, ctr, nonce) =
{
  my(s0 = chachaInit(key, ctr, nonce), s);
  s = s0;
  for (r = 1, 10, for (q = 1, 8, s = quarter(s, CHACHAROUNDS[q])));
  concat(vector(16, i, Vecrev(toDigits((s[i] + s0[i]) % M32, 256, 4))));
}

\\ =================== (no-bitwise) =========================================

bits8(a) = binary(a + 256)[2 .. 9];
XT = matrix(256, 256, a, b, fromdigits((bits8(a - 1) + bits8(b - 1)) % 2, 2));
AT = matrix(256, 256, a, b, fromdigits(vector(8, i, bits8(a - 1)[i] * bits8(b - 1)[i]), 2));
BYTEPOW = [1, 256, 65536, 16777216];

x32(a, b) = sum(i = 1, 4, my(p = BYTEPOW[i]); p * XT[(a \ p) % 256 + 1, (b \ p) % 256 + 1]);
a32(a, b) = sum(i = 1, 4, my(p = BYTEPOW[i]); p * AT[(a \ p) % 256 + 1, (b \ p) % 256 + 1]);
n32(a) = M32 - 1 - a;
rotrN(x, n) = x \ 2^n + (x % 2^n) * 2^(32 - n);
rotlN(x, n) = rotrN(x, 32 - n);

shaCompressN(h, blk) =
{
  my(w = concat(vector(16, i, fromdigits(blk[4 * i - 3 .. 4 * i], 256)), vector(48)), a, b, c, d, e, f, g, hh, s0, s1, t1, t2);
  for (t = 17, 64,
    s0 = x32(x32(rotrN(w[t - 15], 7), rotrN(w[t - 15], 18)), w[t - 15] \ 8);
    s1 = x32(x32(rotrN(w[t - 2], 17), rotrN(w[t - 2], 19)), w[t - 2] \ 1024);
    w[t] = (w[t - 16] + s0 + w[t - 7] + s1) % M32);
  [a, b, c, d, e, f, g, hh] = h;
  for (t = 1, 64,
    t1 = (hh + x32(x32(rotrN(e, 6), rotrN(e, 11)), rotrN(e, 25)) + x32(a32(e, f), a32(n32(e), g)) + SHAK[t] + w[t]) % M32;
    t2 = (x32(x32(rotrN(a, 2), rotrN(a, 13)), rotrN(a, 22)) + x32(x32(a32(a, b), a32(a, c)), a32(b, c))) % M32;
    hh = g; g = f; f = e; e = (d + t1) % M32;
    d = c; c = b; b = a; a = (t1 + t2) % M32);
  (h + [a, b, c, d, e, f, g, hh]) % M32;
}

sha256N(msg) =
{
  my(p = shaPad(msg), h = SHAH0);
  for (i = 0, #p / 64 - 1, h = shaCompressN(h, p[64 * i + 1 .. 64 * i + 64]));
  fromdigits(h, M32);
}

hmacSha256N(key, msg) =
{
  my(k = concat(key, vector(64 - #key)));
  sha256N(concat(apply(b -> XT[b + 1, 93], k), toDigits(sha256N(concat(apply(b -> XT[b + 1, 55], k), msg)), 256, 32)));
}

quarterN(s, q) =
{
  my(a = s[q[1]], b = s[q[2]], c = s[q[3]], d = s[q[4]]);
  a = (a + b) % M32; d = rotlN(x32(d, a), 16);
  c = (c + d) % M32; b = rotlN(x32(b, c), 12);
  a = (a + b) % M32; d = rotlN(x32(d, a), 8);
  c = (c + d) % M32; b = rotlN(x32(b, c), 7);
  s[q[1]] = a; s[q[2]] = b; s[q[3]] = c; s[q[4]] = d;
  s;
}

chachaBlockN(key, ctr, nonce) =
{
  my(s0 = chachaInit(key, ctr, nonce), s);
  s = s0;
  for (r = 1, 10, for (q = 1, 8, s = quarterN(s, CHACHAROUNDS[q])));
  concat(vector(16, i, Vecrev(toDigits((s[i] + s0[i]) % M32, 256, 4))));
}

\\ =================== cases ================================================

bench("app: SHA-256 abc x40", () -> for (k = 1, 40, dAbc = sha256(abc)));
check("app: SHA-256 abc x40", dAbc % CHECKMOD);
bench("app: SHA-256 abc x40 (no-bitwise)", () -> for (k = 1, 40, dAbcN = sha256N(abc)));
check("app: SHA-256 abc x40 (no-bitwise)", dAbcN % CHECKMOD);

bench("app: SHA-256 8x256-byte messages", () -> dMsgs = apply(sha256, msgs));
check("app: SHA-256 8x256-byte messages", vecsum(dMsgs) % CHECKMOD);
bench("app: SHA-256 8x256-byte messages (no-bitwise)", () -> dMsgsN = apply(sha256N, msgs));
check("app: SHA-256 8x256-byte messages (no-bitwise)", vecsum(dMsgsN) % CHECKMOD);

bench("app: HMAC-SHA256 RFC4231 #2 x12", () -> for (k = 1, 12, tag = hmacSha256(hmacKey, hmacMsg)));
check("app: HMAC-SHA256 RFC4231 #2 x12", tag % CHECKMOD);
bench("app: HMAC-SHA256 RFC4231 #2 x12 (no-bitwise)", () -> for (k = 1, 12, tagN = hmacSha256N(hmacKey, hmacMsg)));
check("app: HMAC-SHA256 RFC4231 #2 x12 (no-bitwise)", tagN % CHECKMOD);

bench("app: ChaCha20 keystream 48 blocks", () -> ks = vector(48, c, chachaBlock(chachaKey, c, chachaNonce)));
check("app: ChaCha20 keystream 48 blocks", vecsum(apply(b -> fromdigits(b, 256), ks)) % CHECKMOD);
bench("app: ChaCha20 keystream 48 blocks (no-bitwise)", () -> ksN = vector(48, c, chachaBlockN(chachaKey, c, chachaNonce)));
check("app: ChaCha20 keystream 48 blocks (no-bitwise)", vecsum(apply(b -> fromdigits(b, 256), ksN)) % CHECKMOD);

skip("prim: Hash SHA256 8x256-byte messages x100", "SHA256");

quit();
