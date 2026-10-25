\\ Experiment 101 -- RSA-2048 encrypt / decrypt, plain vs CRT (group F).  See rsa_crt_decrypt.m.
\\
\\ PARI/GP baseline.  "prim:" = Mod(c, n)^e (PARI's modular power) per element
\\ and over the list via apply, chinese() for the recombination.  "app:" CRT
\\ decryption / Garner recombination written with plain GP integers and Mod.
\\ The key is generated untimed, exactly as in 100-rsa-keygen.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Mod^", "chinese", "ispseudoprime", "lcm", "gcd"]);

E = 65537;

\\ ---- the key (untimed) ---------------------------------------------------

genPrime(start) =
{
  my(c = start + 1 - start % 2);
  while (!(ispseudoprime(c) && gcd(c - 1, E) == 1), c += 2);
  c;
}

rsaP = genPrime(bigOf(1101, 1024));
rsaQ = genPrime(bigOf(1102, 1024));
rsaN = rsaP * rsaQ;
rsaD = lift(Mod(E, lcm(rsaP - 1, rsaQ - 1))^-1);
rsaDp = rsaD % (rsaP - 1);
rsaDq = rsaD % (rsaQ - 1);
rsaQinv = lift(Mod(rsaQ, rsaP)^-1);

msgs = vector(4000, i, bigOf(2000 + i, 2000));
msgs200 = msgs[1..200];
cts = apply(m -> lift(Mod(m, rsaN)^E), msgs);
cts50 = cts[1..50];
cts200 = cts[1..200];

\\ ---- prim: encryption ----------------------------------------------------

bench("prim: PowerMod encrypt e=65537 2048-bit x4000", () -> enc = vector(#msgs, i, lift(Mod(msgs[i], rsaN)^E)));
check("prim: PowerMod encrypt e=65537 2048-bit x4000", vecsum(enc) % CHECKMOD);

bench("prim: PowerMod[list] encrypt e=65537 2048-bit x4000", () -> encL = apply(m -> lift(Mod(m, rsaN)^E), msgs));
check("prim: PowerMod[list] encrypt e=65537 2048-bit x4000", vecsum(encL) % CHECKMOD);

\\ ---- prim: plain decryption ---------------------------------------------

bench("prim: PowerMod decrypt c^d mod n 2048-bit x50", () -> dec = apply(c -> lift(Mod(c, rsaN)^rsaD), cts50));
check("prim: PowerMod decrypt c^d mod n 2048-bit x50", vecsum(dec) % CHECKMOD);

\\ ---- app: CRT (Garner) decryption --------------------------------------

crtDecrypt(c) =
{
  my(m1, m2, h);
  m1 = lift(Mod(c, rsaP)^rsaDp);
  m2 = lift(Mod(c, rsaQ)^rsaDq);
  h = (rsaQinv * (m1 - m2)) % rsaP;
  m2 + h * rsaQ;
}

bench("app: RSA-2048 CRT decrypt x200", () -> decC = apply(crtDecrypt, cts200));
check("app: RSA-2048 CRT decrypt x200", vecsum(decC) % CHECKMOD);

\\ ---- recombination only -------------------------------------------------

m1s = apply(c -> lift(Mod(c, rsaP)^rsaDp), cts);
m2s = apply(c -> lift(Mod(c, rsaQ)^rsaDq), cts);

bench("prim: ChineseRemainder recombine 2048-bit x20000", () -> for (r = 1, 5, recC = vector(#m1s, i, lift(chinese(Mod(m1s[i], rsaP), Mod(m2s[i], rsaQ))))));
check("prim: ChineseRemainder recombine 2048-bit x20000", vecsum(recC) % CHECKMOD);

bench("app: Garner recombine 2048-bit x20000", () -> for (r = 1, 5, recG = vector(#m1s, i, m2s[i] + ((rsaQinv * (m1s[i] - m2s[i])) % rsaP) * rsaQ)));
check("app: Garner recombine 2048-bit x20000", vecsum(recG) % CHECKMOD);

\\ ---- app: full round trip -----------------------------------------------

bench("app: RSA-2048 encrypt+CRT decrypt round trip x200", () -> rt = #select(m -> crtDecrypt(lift(Mod(m, rsaN)^E)) == m, msgs200));
check("app: RSA-2048 encrypt+CRT decrypt round trip x200", rt);

quit();
