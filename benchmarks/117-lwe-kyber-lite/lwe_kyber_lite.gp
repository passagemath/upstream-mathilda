\\ Experiment 117 -- Module-LWE / Kyber-like encryption (group F).  See lwe_kyber_lite.m.
\\
\\ PARI/GP baseline.  A ring element of R_q = Z_q[x]/(x^256 + 1) is a native
\\ POLMOD over INTMOD coefficients, Mod(Pol(Mod(c, q)), x^256 + 1): that is how
\\ a GP programmer writes this ring, and PARI multiplies it with its Fp[x]
\\ kernels (Kronecker substitution).  Plain LWE uses integer matrices then
\\ % q.  The modular inverse / solve are the Mod-matrix ^-1 and matsolve
\\ (matsolvemod is ~20x slower here, so the fast native spelling is timed).
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Mod(Pol, x^256+1)", "matrix * col % q", "Mod-matrix ^-1", "matsolve"]);

q = 3329;
NN = 256;
HALF = 1665;
RMOD = x^256 + 1;

cbd(seed, len) = my(b = smallList(seed, 4 * len, 2)); vector(len, i, b[4*i-3] + b[4*i-2] - b[4*i-1] - b[4*i]);
uniform(seed, len) = smallList(seed, len, q);

\\ coefficient vector <-> ring element
toR(v) = Mod(Polrev(v) * Mod(1, q), RMOD);
coeffs(p) = Vecrev(lift(lift(p)), NN);

toBits(c) = apply(d -> ((4 * d + q) \ (2 * q)) % 2, c);

\\ position-weighted fold of a flat vector
foldHash(f) = sum(i = 1, #f, i * f[i]) % CHECKMOD;

\\ ---- prim: the ring product ------------------------------------------------

pa = toR(uniform(11, NN));
pb = toR(cbd(12, NN));
bench("prim: negacyclic polymul n=256 x200", () -> for (i = 1, 200, pc = pa * pb));
check("prim: negacyclic polymul n=256 x200", foldHash(coeffs(pc)));

\\ ---- prim: plain LWE on a 512 x 512 matrix ----------------------------------

lweL = uniform(21, 512 * 512);
lweA = matrix(512, 512, i, j, lweL[512 * (i - 1) + j]);
lweS = cbd(22, 512)~;
lweE = cbd(23, 512)~;
bench("prim: LWE Mod[A.s+e,q] 512x512 x200", () -> for (i = 1, 200, lweB = (lweA * lweS + lweE) % q));
check("prim: LWE Mod[A.s+e,q] 512x512 x200", foldHash(lweB));

lweRL = smallList(24, 256 * 512, 2);
lweR = matrix(256, 512, i, j, lweRL[512 * (i - 1) + j]);
bench("prim: LWE batch Mod[R.A,q] 256x512x512", () -> lweU = (lweR * lweA) % q);
check("prim: LWE batch Mod[R.A,q] 256x512x512", foldHash(concat(Vec(lweU~))));

\\ ---- prim: modular matrix inverse / solve ----------------------------------

invL = uniform(31, 128 * 128);
invM = matrix(128, 128, i, j, invL[128 * (i - 1) + j]);
invB = uniform(32, 128)~;
bench("prim: Inverse[M,Modulus->q] 128x128 x4", () -> for (i = 1, 4, invMi = lift((invM * Mod(1, q))^-1)));
check("prim: Inverse[M,Modulus->q] 128x128 x4", foldHash(concat(Vec(invMi~))));
bench("prim: LinearSolve[M,b,Modulus->q] 128x128 x16", () -> for (i = 1, 16, invX = lift(matsolve(invM * Mod(1, q), invB))));
check("prim: LinearSolve[M,b,Modulus->q] 128x128 x16", foldHash(invX));

\\ ---- app: Kyber-lite --------------------------------------------------------

\\ Inputs for rank k and N messages, as ring elements; seeds offset by 1000 k.
kyberData(k, nMsg) =
{
  my(b = 1000 * k);
  [matrix(k, k, i, j, toR(uniform(b + 10 * i + j, NN))),
   vector(k, i, toR(cbd(b + 100 + i, NN))),
   vector(k, i, toR(cbd(b + 150 + i, NN))),
   vector(nMsg, t, vector(k, i, toR(cbd(b + 200 + 10 * t + i, NN)))),
   vector(nMsg, t, vector(k, i, toR(cbd(b + 500 + 10 * t + i, NN)))),
   vector(nMsg, t, toR(cbd(b + 800 + t, NN))),
   vector(nMsg, t, smallList(b + 900 + t, NN, 2))];
}

kyberKeygen(a, s, e) = my(k = #s); vector(k, i, sum(j = 1, k, a[i, j] * s[j]) + e[i]);

kyberEncrypt(a, t, r, e1, e2, m) =
{
  my(k = #r, u, v);
  u = vector(k, i, sum(j = 1, k, a[j, i] * r[j]) + e1[i]);
  v = sum(j = 1, k, t[j] * r[j]) + e2 + HALF * toR(m);
  concat(u, [v]);
}

kyberDecrypt(s, ct) = my(k = #s); toBits(coeffs(ct[k + 1] - sum(j = 1, k, s[j] * ct[j])));

kyberRun(d) =
{
  my(t, cts, msg, correct = 0, flat = []);
  t = kyberKeygen(d[1], d[2], d[3]);
  cts = vector(#d[7], j, kyberEncrypt(d[1], t, d[4][j], d[5][j], d[6][j], d[7][j]));
  msg = vector(#cts, j, kyberDecrypt(d[2], cts[j]));
  for (j = 1, #cts, for (i = 1, NN, if (msg[j][i] == d[7][j][i], correct++)));
  for (j = 1, #cts, for (i = 1, #cts[j], flat = concat(flat, coeffs(cts[j][i]))));
  correct * 10^10 + foldHash(flat);
}

ky2 = kyberData(2, 8);
bench("app: Kyber-lite k=2 keygen+enc+dec x8", () -> kres2 = kyberRun(ky2));
check("app: Kyber-lite k=2 keygen+enc+dec x8", kres2);

ky3 = kyberData(3, 8);
bench("app: Kyber-lite k=3 keygen+enc+dec x8", () -> kres3 = kyberRun(ky3));
check("app: Kyber-lite k=3 keygen+enc+dec x8", kres3);

\\ GP has one spelling of the ring-vector sum; same program, under the label
\\ of the .m half's MapThread variant.
bench("app: Kyber-lite k=2 MapThread spelling x8", () -> kres2m = kyberRun(ky2));
check("app: Kyber-lite k=2 MapThread spelling x8", kres2m);

\\ ---- app: Regev LWE, one bit per 512-vector ---------------------------------

lweMsg = smallList(25, 256, 2);

lweRun() =
{
  my(b, cts, dec, correct = 0, flat = List());
  b = (lweA * lweS + lweE) % q;
  cts = vector(256, t, my(r = lweR[t, ]); [(r * lweA) % q, (r * b + HALF * lweMsg[t]) % q]);
  dec = vector(256, t, toBits([(cts[t][2] - cts[t][1] * lweS) % q])[1]);
  for (t = 1, 256, if (dec[t] == lweMsg[t], correct++));
  for (t = 1, 256, listput(flat, cts[t][1]); listput(flat, [cts[t][2]]));
  correct * 10^10 + foldHash(concat(Vec(flat)));
}

bench("app: Regev LWE 512x512 enc+dec 256 bits", () -> lres = lweRun());
check("app: Regev LWE 512x512 enc+dec 256 bits", lres);

quit();
