\\ Experiment 108 -- Paillier additively homomorphic encryption (group F).
\\ See paillier.m.
\\
\\ PARI/GP baseline.  "app:" cases: the same Paillier (g = n + 1) on plain GP
\\ integers, exponentiation as lift(Mod(b, n^2)^e), inverse as Mod(x, n)^-1,
\\ L(u) = (u - 1) \ n.  "prim:" is Mod(b, n^2)^n batched x50.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Mod^", "Mod^-1", "lcm"]);

\\ ---- keys ---------------------------------------------------------------

pqs = vector(8, s, [primeOf(2 * s + 3001, 1024), primeOf(2 * s + 3002, 1024)]);

plL(u, n) = (u - 1) \ n;

plKeygen(pq) =
{
  my(p = pq[1], q = pq[2], n, n2, lam, mu);
  n = p * q; n2 = n^2;
  lam = lcm(p - 1, q - 1);
  mu = lift(Mod(plL(lift(Mod(n + 1, n2)^lam), n), n)^-1);
  [n, lam, mu];
}

bench("app: Paillier keygen x8 (given p, q)", () -> keys = apply(plKeygen, pqs));
check("app: Paillier keygen x8 (given p, q)", sum(i = 1, 8, vecsum(keys[i])) % CHECKMOD);

nn = keys[1][1];
lam = keys[1][2];
mu = keys[1][3];
nn2 = nn^2;

\\ ---- prim ---------------------------------------------------------------

pmBases = vector(50, i, bigOf(3100 + i, 4096) % nn2);
bench("prim: PowerMod 4096-bit modulus x50", () -> pms = apply(b -> lift(Mod(b, nn2)^nn), pmBases));
check("prim: PowerMod 4096-bit modulus x50", vecsum(pms) % CHECKMOD);

\\ ---- app: encrypt / decrypt -------------------------------------------

pw(b, e) = lift(Mod(b, nn2)^e);

plEncrypt(m, r) = (pw(nn + 1, m) * pw(r, nn)) % nn2;

plEncryptFast(m, r) = ((1 + m * nn) * pw(r, nn)) % nn2;

plDecrypt(c) = (plL(pw(c, lam), nn) * mu) % nn;

msgs = vector(20, i, bigOf(3200 + i, 2048) % nn);
rs = vector(20, i, bigOf(3300 + i, 2048) % nn);

bench("app: Paillier encrypt x20", () -> cts = vector(20, i, plEncrypt(msgs[i], rs[i])));
check("app: Paillier encrypt x20", vecsum(cts) % CHECKMOD);

bench("app: Paillier encrypt x20 (g=n+1: 1+mn)", () -> ctsF = vector(20, i, plEncryptFast(msgs[i], rs[i])));
check("app: Paillier encrypt x20 (g=n+1: 1+mn)", vecsum(ctsF) % CHECKMOD);

bench("app: Paillier decrypt x20", () -> dec = apply(plDecrypt, cts));
check("app: Paillier decrypt x20", if (dec == msgs, vecsum(dec) % CHECKMOD, -1));

\\ ---- app: homomorphic sum of 1000 ciphertexts --------------------------

rns = apply(r -> pw(r, nn), rs);
sumMsgs = vector(1000, i, bigOf(4000 + i, 2000));
sumCts = vector(1000, i, ((1 + sumMsgs[i] * nn) * rns[i % 20 + 1]) % nn2);

bench("app: homomorphic sum of 1000 ciphertexts + decrypt x10", () -> for (r = 1, 10, cSum = 1; for (i = 1, 1000, cSum = (cSum * sumCts[i]) % nn2); mSum = plDecrypt(cSum)));
check("app: homomorphic sum of 1000 ciphertexts + decrypt x10", if (mSum == vecsum(sumMsgs) % nn, mSum % CHECKMOD, -1));

\\ ---- app: scalar multiplication ----------------------------------------

scalars = vector(20, i, bigOf(3400 + i, 2048) % nn);
bench("app: Paillier scalar mult c^k x20 (k < n)", () -> scs = vector(20, i, pw(cts[i], scalars[i])));
check("app: Paillier scalar mult c^k x20 (k < n)", if (plDecrypt(scs[1]) == (scalars[1] * msgs[1]) % nn, vecsum(scs) % CHECKMOD, -1));

quit();
