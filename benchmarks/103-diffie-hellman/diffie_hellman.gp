\\ Experiment 103 -- finite-field Diffie-Hellman (group F).  See diffie_hellman.m.
\\
\\ PARI/GP baseline.  "prim:" = znprimroot, Mod(2, p)^a, kronecker.  "app:" =
\\ the same safe-prime search (driven by ispseudoprime), textbook primitive
\\ root, DH exchange and Euler-criterion validation as the .m half, in plain
\\ GP integers and Mod().
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["ispseudoprime", "znprimroot", "Mod^", "kronecker"]);

P = MODP2048;
Q = (P - 1) / 2;

\\ ---- app: safe-prime generation -----------------------------------------

safePrime(start) =
{
  my(q = start + 1 - start % 2);
  while (!(ispseudoprime(q) && ispseudoprime(2*q + 1)), q += 2);
  2*q + 1;
}

bench("app: safe-prime search 256-bit x4", () -> sp256 = vector(4, s, safePrime(bigOf(s, 255))));
check("app: safe-prime search 256-bit x4", vecsum(sp256) % CHECKMOD);

benchOnce("app: safe-prime search 512-bit x2", () -> sp512 = vector(2, i, safePrime(bigOf(8 + i, 511))));
check("app: safe-prime search 512-bit x2", vecsum(sp512) % CHECKMOD);

\\ ---- primitive roots of safe primes --------------------------------------

safes = vector(8, s, safePrime(bigOf(s, 255)));

bench("prim: PrimitiveRoot 256-bit safe prime x800", () -> for (r = 1, 100, roots = apply(p -> lift(znprimroot(p)), safes)));
check("prim: PrimitiveRoot 256-bit safe prime x800", vecsum(roots));

safeRoot(p) =
{
  my(q = (p - 1) / 2, g = 2);
  while (g^2 % p == 1 || Mod(g, p)^q == 1, g++);
  g;
}

bench("app: primitive root of safe prime 256-bit x800", () -> for (r = 1, 100, roots2 = apply(safeRoot, safes)));
check("app: primitive root of safe prime 256-bit x800", vecsum(roots2));

\\ ---- MODP-2048 key exchange ----------------------------------------------

privA = vector(50, i, bigOf(6000 + i, 256));
privB = vector(50, i, bigOf(7000 + i, 256));

bench("prim: PowerMod 2^a mod MODP-2048 (256-bit a) x100", () -> pubs = apply(a -> lift(Mod(2, P)^a), concat(privA, privB)));
check("prim: PowerMod 2^a mod MODP-2048 (256-bit a) x100", vecsum(pubs) % CHECKMOD);

dhExchange(a, b) =
{
  my(pa, pb, sa, sb);
  pa = lift(Mod(2, P)^a);
  pb = lift(Mod(2, P)^b);
  sa = lift(Mod(pb, P)^a);
  sb = lift(Mod(pa, P)^b);
  if (sa == sb, sa, 0);
}

bench("app: MODP-2048 DH exchange x50", () -> secrets = vector(50, i, dhExchange(privA[i], privB[i])));
check("app: MODP-2048 DH exchange x50", vecsum(secrets) % CHECKMOD);

\\ ---- peer-key validation -------------------------------------------------

cands = concat(apply(a -> lift(Mod(2, P)^a), privA), vector(50, i, bigOf(8000 + i, 2048) % P));

validQ(y) = 1 < y && y < P - 1 && Mod(y, P)^Q == 1;

bench("app: subgroup validation y^q mod p MODP-2048 x100", () -> nv = #select(validQ, cands));
check("app: subgroup validation y^q mod p MODP-2048 x100", nv);

bench("prim: JacobiSymbol validation MODP-2048 x10000", () -> for (r = 1, 100, nj = #select(y -> 1 < y && y < P - 1 && kronecker(y, P) == 1, cands)));
check("prim: JacobiSymbol validation MODP-2048 x10000", nj);

quit();
