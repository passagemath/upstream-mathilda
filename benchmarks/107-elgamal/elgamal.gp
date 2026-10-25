\\ Experiment 107 -- ElGamal over the RFC 3526 MODP-2048 group (group F).
\\ See elgamal.m.
\\
\\ PARI/GP baseline.  "app:" cases: the same textbook ElGamal on plain GP
\\ integers, exponentiation as lift(Mod(b, p)^e), inverse as Mod(s, p)^-1.
\\ Products of all ciphertexts are formed whole (vecprod) and reduced once,
\\ as the .m half's Times @@ does.  "prim:" is Mod(b, p)^e batched x200.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Mod^", "Mod^-1"]);

p = MODP2048;
g = 2;

\\ ---- prim ---------------------------------------------------------------

pmBases = vector(200, i, bigOf(2000 + i, 2048) % p);
pmExps = vector(200, i, bigOf(2200 + i, 2048));
bench("prim: PowerMod 2048-bit x200", () -> pms = vector(200, i, lift(Mod(pmBases[i], p)^pmExps[i])));
check("prim: PowerMod 2048-bit x200", vecsum(pms) % CHECKMOD);

\\ ---- app ----------------------------------------------------------------

pw(b, e) = lift(Mod(b, p)^e);

egKeygen(seed) = my(x = bigOf(seed, 2048) % (p - 2) + 1); [x, pw(g, x)];

egEncrypt(m, h, y) = [pw(g, y), (m * pw(h, y)) % p];

egDecryptInv(c, x) = (c[2] * lift(Mod(pw(c[1], x), p)^-1)) % p;

egDecryptFermat(c, x) = (c[2] * pw(pw(c[1], x), p - 2)) % p;

bench("app: ElGamal keygen x50", () -> keys = vector(50, i, egKeygen(900 + i)));
check("app: ElGamal keygen x50", sum(i = 1, 50, keys[i][1] + keys[i][2]) % CHECKMOD);

privX = keys[1][1];
pubH = keys[1][2];
msgs = vector(50, i, bigOf(1000 + i, 2048) % (p - 1) + 1);
ephs = vector(50, i, bigOf(1100 + i, 2048) % (p - 2) + 1);

bench("app: ElGamal encrypt x50", () -> cts = vector(50, i, egEncrypt(msgs[i], pubH, ephs[i])));
check("app: ElGamal encrypt x50", sum(i = 1, 50, cts[i][1] + cts[i][2]) % CHECKMOD);

bench("app: ElGamal decrypt x50 (PowerMod[s,-1,p])", () -> dec1 = vector(50, i, egDecryptInv(cts[i], privX)));
check("app: ElGamal decrypt x50 (PowerMod[s,-1,p])", if (dec1 == msgs, vecsum(dec1) % CHECKMOD, -1));

bench("app: ElGamal decrypt x50 (s^(p-2))", () -> dec2 = vector(50, i, egDecryptFermat(cts[i], privX)));
check("app: ElGamal decrypt x50 (s^(p-2))", if (dec2 == msgs, vecsum(dec2) % CHECKMOD, -1));

bench("app: homomorphic product of 50 ciphertexts + decrypt x20", () -> for (r = 1, 20, prodC = [vecprod(vector(50, i, cts[i][1])) % p, vecprod(vector(50, i, cts[i][2])) % p]; prodM = egDecryptInv(prodC, privX)));
check("app: homomorphic product of 50 ciphertexts + decrypt x20", if (prodM == vecprod(msgs) % p, prodM % CHECKMOD, -1));

rands = vector(50, i, bigOf(1200 + i, 2048) % (p - 2) + 1);
bench("app: ElGamal re-randomise x50", () -> rcts = vector(50, i, [(cts[i][1] * pw(g, rands[i])) % p, (cts[i][2] * pw(pubH, rands[i])) % p]));
check("app: ElGamal re-randomise x50", if (vector(50, i, egDecryptInv(rcts[i], privX)) == msgs, sum(i = 1, 50, rcts[i][1] + rcts[i][2]) % CHECKMOD, -1));

quit();
