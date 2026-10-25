\\ Experiment 106 -- ECDSA sign / verify over secp256k1 (group F).  See ecdsa.m.
\\
\\ PARI/GP baseline.  "app:" cases are the same Jacobian arithmetic,
\\ Shamir's-trick verifier and deterministic (NOT RFC 6979) LCG nonces as the
\\ .m half, on plain GP integers with %, inverse via Mod(x, m)^-1.
\\ "prim:" ECDSA cases use PARI's native curve group (ellmul / elladd on
\\ ellinit([0,0,0,0,7], p)).  PARI has no SHA-256, so that case is skip()'d.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Mod^-1", "binary", "ellinit", "ellmul", "elladd"]);

ecP = SECP256K1_P;
ecN = SECP256K1_N;
ecG = [SECP256K1_GX, SECP256K1_GY];
INFJ = [1, 1, 0];

\\ ---- Jacobian point arithmetic (as in 105-ec-arithmetic) -------------------

toJ(P) = if (#P == 0, INFJ, [P[1], P[2], 1]);

fromJ(P) =
{
  my(zi);
  if (P[3] == 0, return([]));
  zi = lift(Mod(P[3], ecP)^-1);
  [(P[1] * zi^2) % ecP, (P[2] * zi^3) % ecP];
}

ecDblJ(P) =
{
  my(x, y, z, a, b, c, d, e, x3);
  x = P[1]; y = P[2]; z = P[3];
  if (z == 0 || y == 0, return(INFJ));
  a = x^2 % ecP; b = y^2 % ecP;
  c = b^2 % ecP;
  d = (2 * ((x + b)^2 - a - c)) % ecP;
  e = 3 * a;
  x3 = (e^2 - 2 * d) % ecP;
  [x3, (e * (d - x3) - 8 * c) % ecP, (2 * y * z) % ecP];
}

ecAddJ(P1, P2) =
{
  my(z1z1, z2z2, u1, u2, s1, s2, h, r, h2, h3, u1h2, x3);
  if (P1[3] == 0, return(P2));
  if (P2[3] == 0, return(P1));
  z1z1 = P1[3]^2 % ecP; z2z2 = P2[3]^2 % ecP;
  u1 = (P1[1] * z2z2) % ecP; u2 = (P2[1] * z1z1) % ecP;
  s1 = (P1[2] * P2[3] * z2z2) % ecP; s2 = (P2[2] * P1[3] * z1z1) % ecP;
  h = (u2 - u1) % ecP; r = (s2 - s1) % ecP;
  if (h == 0, return(if (r == 0, ecDblJ(P1), INFJ)));
  h2 = h^2 % ecP; h3 = (h * h2) % ecP; u1h2 = (u1 * h2) % ecP;
  x3 = (r^2 - h3 - 2 * u1h2) % ecP;
  [x3, (r * (u1h2 - x3) - s1 * h3) % ecP, (P1[3] * P2[3] * h) % ecP];
}

ecMulJ(k, pt) =
{
  my(r = INFJ, pj = toJ(pt), bits = binary(k));
  for (i = 1, #bits, r = ecDblJ(r); if (bits[i], r = ecAddJ(r, pj)));
  fromJ(r);
}

\\ Shamir's trick: u1 G + u2 Q in one pass, G + Q precomputed.
shamirJ(u1, u2, g, q) =
{
  my(gj = toJ(g), qj = toJ(q), gq, len, b1, b2, r = INFJ);
  gq = ecAddJ(gj, qj);
  len = max(#binary(u1), #binary(u2));
  b1 = binary(u1); b1 = concat(vector(len - #b1), b1);
  b2 = binary(u2); b2 = concat(vector(len - #b2), b2);
  for (i = 1, len, r = ecDblJ(r); if (b1[i] && b2[i], r = ecAddJ(r, gq), if (b1[i], r = ecAddJ(r, gj), if (b2[i], r = ecAddJ(r, qj)))));
  fromJ(r);
}

\\ ---- ECDSA ----------------------------------------------------------------

ecdsaSign(z, d, k) =
{
  my(r, s);
  r = ecMulJ(k, ecG)[1] % ecN;
  s = (lift(Mod(k, ecN)^-1) * (z + r * d)) % ecN;
  [r, s];
}

ecdsaVerify(z, sig, q) =
{
  my(r = sig[1], s = sig[2], w, x);
  if (!(1 <= r && r < ecN && 1 <= s && s < ecN), return(0));
  w = lift(Mod(s, ecN)^-1);
  x = shamirJ((z * w) % ecN, (r * w) % ecN, ecG, q);
  if (#x == 0, return(0));
  if (x[1] % ecN == r, 1, 0);
}

privD = bigOf(701, 256) % ecN;
pubQ = ecMulJ(privD, ecG);

zs = vector(48, i, bigOf(800 + i, 256));
nonces = vector(48, i, bigOf(1800 + i, 256) % ecN);

bench("app: ECDSA sign x48", () -> sigs = vector(48, i, ecdsaSign(zs[i], privD, nonces[i])));
check("app: ECDSA sign x48", sum(i = 1, 48, sigs[i][1] + sigs[i][2]) % CHECKMOD);

bench("app: ECDSA verify x48 (Shamir)", () -> nValid = sum(i = 1, 48, ecdsaVerify(zs[i], sigs[i], pubQ)));
check("app: ECDSA verify x48 (Shamir)", nValid);

bench("app: ECDSA verify tampered x48 (Shamir)", () -> nBad = sum(i = 1, 48, ecdsaVerify(zs[i] + 1, sigs[i], pubQ)));
check("app: ECDSA verify tampered x48 (Shamir)", nBad);

\\ ---- prim -----------------------------------------------------------------

skip("prim: SHA-256 of 1000 messages", "sha256");

E = ellinit([0, 0, 0, 0, 7], ecP);

primSign(z, d, k) =
{
  my(r, s);
  r = lift(ellmul(E, ecG, k)[1]) % ecN;
  s = (lift(Mod(k, ecN)^-1) * (z + r * d)) % ecN;
  [r, s];
}

primVerify(z, sig, q) =
{
  my(r = sig[1], s = sig[2], w, x);
  if (!(1 <= r && r < ecN && 1 <= s && s < ecN), return(0));
  w = lift(Mod(s, ecN)^-1);
  x = elladd(E, ellmul(E, ecG, (z * w) % ecN), ellmul(E, q, (r * w) % ecN));
  if (x == [0], return(0));
  if (lift(x[1]) % ecN == r, 1, 0);
}

bench("prim: ECDSA sign x48", () -> psigs = vector(48, i, primSign(zs[i], privD, nonces[i])));
check("prim: ECDSA sign x48", sum(i = 1, 48, psigs[i][1] + psigs[i][2]) % CHECKMOD);

bench("prim: ECDSA verify x48", () -> pValid = sum(i = 1, 48, primVerify(zs[i], sigs[i], pubQ)));
check("prim: ECDSA verify x48", pValid);

quit();
