\\ Experiment 102 -- primality testing (group F).  See primality_testing.m.
\\
\\ PARI/GP baseline.  "prim:" = ispseudoprime (BPSW) for PrimeQ and isprime
\\ (APRCL / ECPP proof) for the proven row.  "app:" = the same Miller-Rabin
\\ and Baillie-PSW code as the .m half, plain GP integers, Mod() powers and
\\ kronecker() for the Jacobi symbol.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["ispseudoprime", "isprime", "kronecker", "valuation"]);

\\ ---- inputs (untimed) -----------------------------------------------------

prs = concat([vector(10, i, primeOf(400 + i, 512)), vector(10, i, primeOf(410 + i, 1024)), vector(10, i, primeOf(420 + i, 2048))]);
odds1024 = vector(2000, i, my(b = bigOf(5000 + i, 1024)); b + 1 - b % 2);
odds500 = odds1024[1..500];

chernick(start) =
{
  my(k = start);
  while (!(ispseudoprime(6*k + 1) && ispseudoprime(12*k + 1) && ispseudoprime(18*k + 1)), k++);
  (6*k + 1) * (12*k + 1) * (18*k + 1);
}

carmichaels = vector(8, i, chernick(bigOf(3000 + i, 64)));
spsps = [3215031751, 2152302898747, 3474749660383, 341550071728321, 3825123056546413051, 318665857834031151167461, 3317044064679887385961981];
mersenneC = apply(p -> 2^p - 1, select(isprime, [1000..1199]));
fermatC = vector(4, i, 2^(2^(7 + i)) + 1);
pseudos = concat([carmichaels, spsps, mersenneC, fermatC]);

\\ ---- app: the tests written in the language ------------------------------

MRBASES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

mrRound(n, a, d, s) =
{
  my(x = lift(Mod(a, n)^d));
  if (x == 1 || x == n - 1, return(1));
  for (r = 1, s - 1, x = x^2 % n; if (x == n - 1, return(1)));
  0;
}

millerRabin(n) =
{
  my(s, d);
  if (n < 2, return(0));
  if (n % 2 == 0, return(n == 2));
  s = valuation(n - 1, 2); d = (n - 1) >> s;
  for (i = 1, #MRBASES, if (MRBASES[i] % n != 0 && !mrRound(n, MRBASES[i], d, s), return(0)));
  1;
}

half(x, n) = if (x % 2, (x + n) / 2, x / 2);

strongLucas(n) =
{
  my(dd = 5, j, p = 1, q, d, s, bits, u, v, qk, t);
  while (1,
    j = kronecker(dd, n);
    if (j == -1, break);
    if (j == 0 && abs(dd) != n, return(0));
    dd = if (dd > 0, -dd - 2, -dd + 2));
  q = (1 - dd) / 4;
  s = valuation(n + 1, 2); d = (n + 1) >> s;
  bits = binary(d);
  u = 1; v = p; qk = q % n;
  for (i = 2, #bits,
    t = (u * v) % n; v = (v^2 - 2 * qk) % n; u = t;
    qk = qk^2 % n;
    if (bits[i] == 1,
      t = half((p * u + v) % n, n); v = half((dd * u + p * v) % n, n); u = t;
      qk = (qk * q) % n));
  if (u == 0 || v == 0, return(1));
  for (r = 1, s - 1,
    v = (v^2 - 2 * qk) % n;
    if (v == 0, return(1));
    qk = qk^2 % n);
  0;
}

bpsw(n) =
{
  my(s, d);
  if (n < 2, return(0));
  if (n % 2 == 0, return(n == 2));
  s = valuation(n - 1, 2); d = (n - 1) >> s;
  if (!mrRound(n, 2, d, s), return(0));
  if (issquare(n), return(0));
  strongLucas(n);
}

\\ The .m "no-Return" variants: early exits as loop conditions.  Same verdicts.

mrRound2(n, a, d, s) =
{
  my(x = lift(Mod(a, n)^d), r = 1);
  if (x == 1 || x == n - 1, return(1));
  while (r < s && x != n - 1, x = x^2 % n; r++);
  x == n - 1;
}

millerRabin2(n) =
{
  my(s, d);
  if (n < 2, return(0));
  if (n % 2 == 0, return(n == 2));
  s = valuation(n - 1, 2); d = (n - 1) >> s;
  for (i = 1, #MRBASES, if (MRBASES[i] % n != 0 && !mrRound2(n, MRBASES[i], d, s), return(0)));
  1;
}

strongLucas2(n) =
{
  my(dd = 5, p = 1, q, d, s, bits, u, v, qk, t, r = 1);
  while (kronecker(dd, n) == 1, dd = if (dd > 0, -dd - 2, -dd + 2));
  if (kronecker(dd, n) == 0 && abs(dd) != n, return(0));
  q = (1 - dd) / 4;
  s = valuation(n + 1, 2); d = (n + 1) >> s;
  bits = binary(d);
  u = 1; v = p; qk = q % n;
  for (i = 2, #bits,
    t = (u * v) % n; v = (v^2 - 2 * qk) % n; u = t;
    qk = qk^2 % n;
    if (bits[i] == 1,
      t = half((p * u + v) % n, n); v = half((dd * u + p * v) % n, n); u = t;
      qk = (qk * q) % n));
  if (u == 0 || v == 0, return(1));
  while (r < s && v != 0, v = (v^2 - 2 * qk) % n; qk = qk^2 % n; r++);
  v == 0;
}

bpsw2(n) =
{
  my(s, d);
  if (n < 2, return(0));
  if (n % 2 == 0, return(n == 2));
  s = valuation(n - 1, 2); d = (n - 1) >> s;
  if (!mrRound2(n, 2, d, s), return(0));
  if (issquare(n), return(0));
  strongLucas2(n);
}

cnt(f, v) = #select(f, v);

\\ ---- prim -------------------------------------------------------------------

bench("prim: PrimeQ 512/1024/2048-bit primes x30", () -> c1 = cnt(ispseudoprime, prs));
check("prim: PrimeQ 512/1024/2048-bit primes x30", c1);

bench("prim: PrimeQ 1024-bit odd x2000", () -> c2 = cnt(ispseudoprime, odds1024));
check("prim: PrimeQ 1024-bit odd x2000", c2);

bench("prim: PrimeQ pseudoprimes (Carmichael/spsp/Mersenne/Fermat) x47", () -> c3 = cnt(ispseudoprime, pseudos));
check("prim: PrimeQ pseudoprimes (Carmichael/spsp/Mersenne/Fermat) x47", c3);

bench("prim: PrimeQ Mersenne 2^4423-1", () -> c4 = ispseudoprime(2^4423 - 1));
check("prim: PrimeQ Mersenne 2^4423-1", c4);

proofPrimes = vector(4, i, primeOf(430 + i, 256));
benchOnce("prim: proven primality 256-bit prime x4", () -> c5 = cnt(isprime, proofPrimes));
check("prim: proven primality 256-bit prime x4", c5);

\\ ---- app ----------------------------------------------------------------------

bench("app: Miller-Rabin 12 bases primes x30", () -> a1 = cnt(millerRabin, prs));
check("app: Miller-Rabin 12 bases primes x30", a1);

bench("app: Miller-Rabin 12 bases pseudoprimes x47", () -> a2 = cnt(millerRabin, pseudos));
check("app: Miller-Rabin 12 bases pseudoprimes x47", a2);

bench("app: Baillie-PSW primes x30", () -> a3 = cnt(bpsw, prs));
check("app: Baillie-PSW primes x30", a3);

bench("app: Baillie-PSW pseudoprimes x47", () -> a4 = cnt(bpsw, pseudos));
check("app: Baillie-PSW pseudoprimes x47", a4);

bench("app: Baillie-PSW 1024-bit odd x500", () -> a5 = cnt(bpsw, odds500));
check("app: Baillie-PSW 1024-bit odd x500", a5);

bench("app: Miller-Rabin 12 bases (no-Return) primes x30", () -> b1 = cnt(millerRabin2, prs));
check("app: Miller-Rabin 12 bases (no-Return) primes x30", b1);

bench("app: Miller-Rabin 12 bases (no-Return) pseudoprimes x47", () -> b2 = cnt(millerRabin2, pseudos));
check("app: Miller-Rabin 12 bases (no-Return) pseudoprimes x47", b2);

bench("app: Baillie-PSW (no-Return) primes x30", () -> b3 = cnt(bpsw2, prs));
check("app: Baillie-PSW (no-Return) primes x30", b3);

bench("app: Baillie-PSW (no-Return) 1024-bit odd x500", () -> b4 = cnt(bpsw2, odds500));
check("app: Baillie-PSW (no-Return) 1024-bit odd x500", b4);

quit();
