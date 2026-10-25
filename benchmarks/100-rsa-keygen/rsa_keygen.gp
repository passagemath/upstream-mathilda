\\ Experiment 100 -- RSA key generation (group F).  See rsa_keygen.m.
\\
\\ PARI/GP baseline.  nextprime / ispseudoprime (BPSW) are the primitives a GP
\\ keygen uses; Mod(e, l)^-1 is the inverse.  The "app:" keygen is the same
\\ incremental search as the .m half, driven by ispseudoprime.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["nextprime", "ispseudoprime", "Mod^-1", "lcm", "gcd"]);

E = 65537;

\\ ---- prim ---------------------------------------------------------------

starts1024 = vector(4, i, bigOf(100 + i, 1024));
bench("prim: NextPrime 1024-bit x4", () -> nps = apply(x -> nextprime(x + 1), starts1024));
check("prim: NextPrime 1024-bit x4", vecsum(nps) % CHECKMOD);

primes2048 = vector(8, i, primeOf(200 + i, 2048));
bench("prim: PrimeQ 2048-bit prime x8", () -> cnt = #select(ispseudoprime, primes2048));
check("prim: PrimeQ 2048-bit prime x8", cnt);

odds2048 = vector(400, i, my(b = bigOf(300 + i, 2048)); b + 1 - b % 2);
bench("prim: PrimeQ 2048-bit odd x400", () -> cnt = #select(ispseudoprime, odds2048));
check("prim: PrimeQ 2048-bit odd x400", cnt);

lams = select(x -> gcd(x, E) == 1, vector(2100, i, bigOf(1000 + i, 2048)))[1..2000];
bench("prim: PowerMod[e,-1,lambda] 2048-bit x20000", () -> for (r = 1, 10, ds = apply(l -> lift(Mod(E, l)^-1), lams)));
check("prim: PowerMod[e,-1,lambda] 2048-bit x20000", vecsum(ds) % CHECKMOD);

bench("prim: ModularInverse 2048-bit x20000", () -> for (r = 1, 10, ds2 = apply(l -> lift(Mod(E, l)^-1), lams)));
check("prim: ModularInverse 2048-bit x20000", vecsum(ds2) % CHECKMOD);

\\ ---- app ----------------------------------------------------------------

genPrime(start) =
{
  my(c = start + 1 - start % 2);
  while (!(ispseudoprime(c) && gcd(c - 1, E) == 1), c += 2);
  c;
}

rsaKey(bits, s) =
{
  my(p, q, lam, d);
  p = genPrime(bigOf(s, bits / 2));
  q = genPrime(bigOf(s + 1, bits / 2));
  lam = lcm(p - 1, q - 1);
  d = lift(Mod(E, lam)^-1);
  [p * q, d, d % (p - 1), d % (q - 1), lift(Mod(q, p)^-1)];
}

keyDigest(keys) = sum(i = 1, #keys, vecsum(keys[i])) % CHECKMOD;

bench("app: RSA-1024 keygen x8", () -> k1 = vector(8, s, rsaKey(1024, 2 * s)));
check("app: RSA-1024 keygen x8", keyDigest(k1));

benchOnce("app: RSA-2048 keygen x2", () -> k2 = vector(2, s, rsaKey(2048, 2 * s + 50)));
check("app: RSA-2048 keygen x2", keyDigest(k2));

quit();
