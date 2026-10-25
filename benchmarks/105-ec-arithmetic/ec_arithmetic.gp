\\ Experiment 105 -- elliptic-curve point arithmetic on secp256k1 (group F).
\\ See ec_arithmetic.m.
\\
\\ PARI/GP baseline.  The "app:" cases are the same affine / Jacobian / ladder
\\ code as the .m half, on plain GP integers reduced with %, inverse via
\\ Mod(x, p)^-1.  The "prim:" cases are PARI's native curve group: ellmul /
\\ elladd on ellinit([0,0,0,0,7], p) -- the cross-check for the app cases.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Mod^-1", "binary", "ellinit", "ellmul", "elladd"]);

ecP = SECP256K1_P;
ecN = SECP256K1_N;
ecG = [SECP256K1_GX, SECP256K1_GY];
INFJ = [1, 1, 0];

\\ ---- affine coordinates ([] is the point at infinity) ---------------------

ecDblA(P) =
{
  my(x, y, l, x3);
  if (#P == 0 || P[2] == 0, return([]));
  x = P[1]; y = P[2];
  l = (3 * x^2 * lift(Mod(2 * y, ecP)^-1)) % ecP;
  x3 = (l^2 - 2 * x) % ecP;
  [x3, (l * (x - x3) - y) % ecP];
}

ecAddA(P, Q) =
{
  my(l, x3);
  if (#P == 0, return(Q));
  if (#Q == 0, return(P));
  if (P[1] == Q[1], return(if ((P[2] + Q[2]) % ecP == 0, [], ecDblA(P))));
  l = ((Q[2] - P[2]) * lift(Mod(Q[1] - P[1], ecP)^-1)) % ecP;
  x3 = (l^2 - P[1] - Q[1]) % ecP;
  [x3, (l * (P[1] - x3) - P[2]) % ecP];
}

ecMulA(k, pt) =
{
  my(r = [], bits = binary(k));
  for (i = 1, #bits, r = ecDblA(r); if (bits[i], r = ecAddA(r, pt)));
  r;
}

\\ ---- Jacobian coordinates -------------------------------------------------

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

ecLadderJ(k, pt) =
{
  my(r0 = INFJ, r1 = toJ(pt), bits = binary(k));
  for (i = 1, #bits, if (bits[i] == 0, r1 = ecAddJ(r0, r1); r0 = ecDblJ(r0), r0 = ecAddJ(r0, r1); r1 = ecDblJ(r1)));
  fromJ(r0);
}

xDigest(pts) = sum(i = 1, #pts, lift(pts[i][1])) % CHECKMOD;

ks = vector(64, i, bigOf(500 + i, 256) % ecN);

\\ ---- app: single operations, batched -------------------------------------

g2 = ecDblA(ecG);

bench("app: affine point add x10000", () -> qa = g2; for (i = 1, 10000, qa = ecAddA(qa, ecG)));
check("app: affine point add x10000", qa[1] % CHECKMOD);

gJ = toJ(ecG);
bench("app: Jacobian point add x10000", () -> qj = toJ(g2); for (i = 1, 10000, qj = ecAddJ(qj, gJ)));
check("app: Jacobian point add x10000", fromJ(qj)[1] % CHECKMOD);

bench("app: affine point double x10000", () -> da = ecG; for (i = 1, 10000, da = ecDblA(da)));
check("app: affine point double x10000", da[1] % CHECKMOD);

bench("app: Jacobian point double x10000", () -> dj = gJ; for (i = 1, 10000, dj = ecDblJ(dj)));
check("app: Jacobian point double x10000", fromJ(dj)[1] % CHECKMOD);

\\ ---- app: scalar multiplication -------------------------------------------

bench("app: scalar mult affine double-and-add x64", () -> ma = apply(k -> ecMulA(k, ecG), ks));
check("app: scalar mult affine double-and-add x64", xDigest(ma));

bench("app: scalar mult Jacobian double-and-add x64", () -> mj = apply(k -> ecMulJ(k, ecG), ks));
check("app: scalar mult Jacobian double-and-add x64", xDigest(mj));

bench("app: scalar mult Jacobian Montgomery ladder x64", () -> ml = apply(k -> ecLadderJ(k, ecG), ks));
check("app: scalar mult Jacobian Montgomery ladder x64", xDigest(ml));

\\ ---- prim: PARI's native curve group --------------------------------------

E = ellinit([0, 0, 0, 0, 7], ecP);

bench("prim: point add x10000", () -> qe = g2; for (i = 1, 10000, qe = elladd(E, qe, ecG)));
check("prim: point add x10000", lift(qe[1]) % CHECKMOD);

bench("prim: scalar mult k*G x64", () -> me = apply(k -> ellmul(E, ecG, k), ks));
check("prim: scalar mult k*G x64", xDigest(me));

quit();
