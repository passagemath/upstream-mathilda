\\ Experiment 114 -- AES-128 and GF(2^8) (group F).  See aes_gf256.m.
\\
\\ PARI/GP baseline.  The same AES-128 twice in GP: idiomatically with
\\ bitxor / bitand / bitor / shift, and as the "(no-bitwise)" variant that
\\ XORs bytes by lookup in a 256x256 matrix built (untimed) from base-2
\\ digits.  The GF(2)[x] route inverts Mod(Mod(1,2)*poly, aes poly) -- PARI's
\\ native finite-field arithmetic.  PARI has no AES primitive: the library
\\ case is skip()'d.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["bitxor", "bitand", "bitor", "shift", "Mod(Pol)^-1", "binary", "fromdigits"]);

hexBytes(s) = my(c = Vec(s)); vector(#c / 2, i, eval(Str("0x", c[2 * i - 1], c[2 * i])));
foldB(bs) = fromdigits(bs, 256) % CHECKMOD;

fipsKey = vector(16, i, i - 1);
fipsPt = hexBytes("00112233445566778899aabbccddeeff");
shiftIdx = concat(vector(4, c, vector(4, r, (r - 1) + 4 * ((c - 1 + r - 1) % 4) + 1)));
ecbKey = smallList(1401, 16, 256);
ecbFlat = smallList(1402, 128 * 16, 256);
ecbBlocks = vector(128, i, ecbFlat[16 * (i - 1) + 1 .. 16 * i]);

\\ =================== idiomatic (bitxor & friends) ========================

xtime(b) = bitxor(bitand(shift(b, 1), 255), if (b >= 128, 27, 0));

gmul(a, b) =
{
  my(p = 0);
  while (b > 0, if (bitand(b, 1), p = bitxor(p, a)); a = xtime(a); b = shift(b, -1));
  p;
}

ginv(a) = my(r = 1, s = a); for (i = 1, 7, s = gmul(s, s); r = gmul(r, s)); r;

rotl8(b, k) = bitand(bitor(shift(b, k), shift(b, k - 8)), 255);

makeSbox() = vector(256, i, my(b = ginv(i - 1)); bitxor(bitxor(bitxor(bitxor(bitxor(b, rotl8(b, 1)), rotl8(b, 2)), rotl8(b, 3)), rotl8(b, 4)), 99));

expandKey(key, sb) =
{
  my(w = vector(44), t, rc = 1);
  for (i = 1, 4, w[i] = key[4 * i - 3 .. 4 * i]);
  for (i = 4, 43,
    t = w[i];
    if (i % 4 == 0,
      t = vector(4, j, sb[t[j % 4 + 1] + 1]);
      t[1] = bitxor(t[1], rc); rc = xtime(rc));
    w[i + 1] = vector(4, j, bitxor(w[i - 3][j], t[j])));
  vector(11, r, concat(w[4 * r - 3 .. 4 * r]));
}

mixColumns(s) =
{
  my(out = vector(16), col, t);
  for (c = 0, 3,
    col = s[4 * c + 1 .. 4 * c + 4];
    t = bitxor(bitxor(bitxor(col[1], col[2]), col[3]), col[4]);
    for (r = 1, 4, out[4 * c + r] = bitxor(bitxor(col[r], t), xtime(bitxor(col[r], col[r % 4 + 1])))));
  out;
}

encryptBlock(pt, rk, sb) =
{
  my(s = vector(16, i, bitxor(pt[i], rk[1][i])), u);
  for (r = 2, 10,
    u = vector(16, i, sb[s[shiftIdx[i]] + 1]);
    u = mixColumns(u);
    s = vector(16, i, bitxor(u[i], rk[r][i])));
  vector(16, i, bitxor(sb[s[shiftIdx[i]] + 1], rk[11][i]));
}

\\ =================== (no-bitwise): arithmetic + table lookup ==============

XT = matrix(256, 256, a, b, fromdigits((binary(a - 1 + 256)[2..9] + binary(b - 1 + 256)[2..9]) % 2, 2));
x8(a, b) = XT[a + 1, b + 1];
x8v(u, v) = vector(#u, i, XT[u[i] + 1, v[i] + 1]);

xtimeN(b) = if (b >= 128, x8(2 * b - 256, 27), 2 * b);

gmulN(a, b) =
{
  my(p = 0);
  while (b > 0, if (b % 2 == 1, p = x8(p, a)); a = xtimeN(a); b = b \ 2);
  p;
}

ginvN(a) = my(r = 1, s = a); for (i = 1, 7, s = gmulN(s, s); r = gmulN(r, s)); r;

rotl8N(b, k) = (b * 2^k) % 256 + b \ 2^(8 - k);

makeSboxN() = vector(256, i, my(b = ginvN(i - 1)); x8(x8(x8(x8(x8(b, rotl8N(b, 1)), rotl8N(b, 2)), rotl8N(b, 3)), rotl8N(b, 4)), 99));

expandKeyN(key, sb) =
{
  my(w = vector(44), t, rc = 1);
  for (i = 1, 4, w[i] = key[4 * i - 3 .. 4 * i]);
  for (i = 4, 43,
    t = w[i];
    if (i % 4 == 0,
      t = vector(4, j, sb[t[j % 4 + 1] + 1]);
      t[1] = x8(t[1], rc); rc = xtimeN(rc));
    w[i + 1] = x8v(w[i - 3], t));
  vector(11, r, concat(w[4 * r - 3 .. 4 * r]));
}

mixColumnsN(s) =
{
  my(out = vector(16), col, t);
  for (c = 0, 3,
    col = s[4 * c + 1 .. 4 * c + 4];
    t = x8(x8(x8(col[1], col[2]), col[3]), col[4]);
    for (r = 1, 4, out[4 * c + r] = x8(x8(col[r], t), xtimeN(x8(col[r], col[r % 4 + 1])))));
  out;
}

encryptBlockN(pt, rk, sb) =
{
  my(s = x8v(pt, rk[1]), u);
  for (r = 2, 10,
    u = vector(16, i, sb[s[shiftIdx[i]] + 1]);
    s = x8v(mixColumnsN(u), rk[r]));
  x8v(vector(16, i, sb[s[shiftIdx[i]] + 1]), rk[11]);
}

\\ =================== cases ================================================

bench("app: AES S-box generation x4", () -> for (k = 1, 4, sbox = makeSbox()));
check("app: AES S-box generation x4", foldB(sbox));

bench("app: AES S-box generation x4 (no-bitwise)", () -> for (k = 1, 4, sboxN = makeSboxN()));
check("app: AES S-box generation x4 (no-bitwise)", foldB(sboxN));

bench("app: AES-128 FIPS-197 C.1 x100", () -> for (k = 1, 100, fipsCt = encryptBlock(fipsPt, expandKey(fipsKey, sbox), sbox)));
check("app: AES-128 FIPS-197 C.1 x100", foldB(fipsCt));

bench("app: AES-128 FIPS-197 C.1 x100 (no-bitwise)", () -> for (k = 1, 100, fipsCtN = encryptBlockN(fipsPt, expandKeyN(fipsKey, sboxN), sboxN)));
check("app: AES-128 FIPS-197 C.1 x100 (no-bitwise)", foldB(fipsCtN));

bench("app: AES-128 ECB encrypt 128 blocks", () -> my(rk = expandKey(ecbKey, sbox)); ecbCt = apply(b -> encryptBlock(b, rk, sbox), ecbBlocks));
check("app: AES-128 ECB encrypt 128 blocks", vecsum(apply(foldB, ecbCt)) % CHECKMOD);

bench("app: AES-128 ECB encrypt 128 blocks (no-bitwise)", () -> my(rk = expandKeyN(ecbKey, sboxN)); ecbCtN = apply(b -> encryptBlockN(b, rk, sboxN), ecbBlocks));
check("app: AES-128 ECB encrypt 128 blocks (no-bitwise)", vecsum(apply(foldB, ecbCtN)) % CHECKMOD);

\\ ---- GF(2)[x] route: PARI's native inverse in F_2[x]/(aes poly) ----------

aesPoly = Mod(1, 2) * (x^8 + x^4 + x^3 + x + 1);
toPoly(a) = Mod(1, 2) * Pol(binary(a), x);
polyInv(a) = subst(lift(lift(Mod(toPoly(a), aesPoly)^-1)), x, 2);

bench("prim: GF(2^8) inverse via PolynomialExtendedGCD Modulus->2 x7650", () -> for (k = 1, 30, invs = vector(255, a, polyInv(a))));
check("prim: GF(2^8) inverse via PolynomialExtendedGCD Modulus->2 x7650", sum(a = 1, 255, a * invs[a]) % CHECKMOD);

skip("prim: AES-128 ECB library 65536 blocks", "AES");

quit();
