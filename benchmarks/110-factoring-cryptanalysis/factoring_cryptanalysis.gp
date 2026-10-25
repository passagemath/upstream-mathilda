\\ Experiment 110 -- Factoring cryptanalysis (group F).  See factoring_cryptanalysis.m.
\\
\\ PARI/GP baseline.  The "app:" algorithms are written in plain GP (Pollard
\\ rho/Brent, Pollard p-1 stage 1, Fermat), identical to the .m/.py halves; the
\\ "prim:" case is factor(), PARI's native factoring.  One statement per
\\ top-level line; multi-line bodies in { } braces.

read("../harness.gp");
read("../crypto_data.gp");

require(["factor", "Mod^k", "gcd", "nextprime", "isprime"]);

smaller(n, f) = min(f, n / f);

\\ ---- app: Pollard rho (Brent) -------------------------------------------

rhoBrent(n) =
{
  my(y = 2, c = 1, m = 128, g = 1, r = 1, q = 1, x, ys, k, lim);
  if (n % 2 == 0, return(2));
  while (g == 1,
    x = y; for (j = 1, r, y = (y^2 + c) % n);
    k = 0;
    while (k < r && g == 1,
      ys = y; lim = min(m, r - k);
      for (j = 1, lim, y = (y^2 + c) % n; q = (q * abs(x - y)) % n);
      g = gcd(q, n); k += m);
    r *= 2);
  if (g == n, g = 1; while (g == 1, ys = (ys^2 + c) % n; g = gcd(abs(x - ys), n)));
  g;
}

weakRho(s, small, bits) = nextprime(bigOf(s, small) + 1) * nextprime(bigOf(s + 1, bits - small) + 1);

rhoNs = [weakRho(2501, 28, 64), weakRho(2511, 28, 80), weakRho(2521, 28, 96), weakRho(2531, 28, 112)];

bench("app: Pollard rho (Brent) 28-bit factor x4", () -> rhoF = apply(rhoBrent, rhoNs));
check("app: Pollard rho (Brent) 28-bit factor x4", sum(i = 1, 4, smaller(rhoNs[i], rhoF[i])) % CHECKMOD + CHECKMOD * (sum(i = 1, 4, rhoNs[i] % rhoF[i] == 0)));

\\ ---- app: Pollard p-1 (stage 1) -----------------------------------------

PRIMES = primes(primepi(20000));

pm1(n, B1) =
{
  my(a = 2, e, pr, g);
  for (i = 1, #PRIMES, pr = PRIMES[i]; if (pr > B1, break);
    e = floor(log(B1) / log(pr));
    a = lift(Mod(a, n)^(pr^e)));
  gcd(a - 1, n);
}

smoothPrime(b, t0) =
{
  my(base = 2 * prod(i = 1, primepi(b), prime(i)), t = t0, p);
  p = base * t + 1;
  while (!isprime(p), t++; p = base * t + 1);
  p;
}

pm1Keys = [[smoothPrime(40, 1) * nextprime(bigOf(2601, 45) + 1), 60], [smoothPrime(50, 3) * nextprime(bigOf(2602, 55) + 1), 80], [smoothPrime(60, 2) * nextprime(bigOf(2603, 60) + 1), 100], [smoothPrime(70, 5) * nextprime(bigOf(2604, 64) + 1), 120]];

bench("app: Pollard p-1 stage 1 x120", () -> for (r = 1, 30, pm1F = apply(key -> pm1(key[1], key[2]), pm1Keys)));
check("app: Pollard p-1 stage 1 x120", sum(i = 1, 4, smaller(pm1Keys[i][1], pm1F[i])) % CHECKMOD + CHECKMOD * (sum(i = 1, 4, pm1F[i] > 1 && pm1F[i] < pm1Keys[i][1] && pm1Keys[i][1] % pm1F[i] == 0)));

\\ ---- app: Fermat (close primes) -----------------------------------------

fermat(n) =
{
  my(a = ceil(sqrt(n)), b2);
  b2 = a^2 - n;
  while (!issquare(b2), a++; b2 = a^2 - n);
  a - sqrtint(b2);
}

closeSemiprime(s, bits, gap) = my(p = nextprime(bigOf(s, bits) + 1)); p * nextprime(p + gap + 1);

fermatNs = [closeSemiprime(2701, 40, 2^25), closeSemiprime(2702, 40, 2^25), closeSemiprime(2703, 40, 2^26), closeSemiprime(2704, 40, 2^24)];

bench("app: Fermat (close primes) x4", () -> ferF = apply(fermat, fermatNs));
check("app: Fermat (close primes) x4", sum(i = 1, 4, smaller(fermatNs[i], ferF[i])) % CHECKMOD + CHECKMOD * (sum(i = 1, 4, ferF[i] > 1 && ferF[i] < fermatNs[i] && fermatNs[i] % ferF[i] == 0)));

\\ ---- prim: factor on balanced semiprimes --------------------------------

bal(s, bits) = nextprime(bigOf(s, bits \ 2) + 1) * nextprime(bigOf(s + 1, bits \ 2) + 1);
smallestFactor(n) = factor(n)[1, 1];

bal80 = vector(25, k, bal(2800 + 2 * (k - 1), 80));
benchOnce("prim: FactorInteger balanced 80-bit x25", () -> sf80 = apply(smallestFactor, bal80));
check("prim: FactorInteger balanced 80-bit x25", vecsum(sf80) % CHECKMOD);

bal100 = vector(10, k, bal(2850 + 2 * (k - 1), 100));
benchOnce("prim: FactorInteger balanced 100-bit x10", () -> sf100 = apply(smallestFactor, bal100));
check("prim: FactorInteger balanced 100-bit x10", vecsum(sf100) % CHECKMOD);

bal120 = vector(3, k, bal(2880 + 2 * (k - 1), 120));
benchOnce("prim: FactorInteger balanced 120-bit x3", () -> sf120 = apply(smallestFactor, bal120));
check("prim: FactorInteger balanced 120-bit x3", vecsum(sf120) % CHECKMOD);

n160 = bal(2900, 160);
benchOnce("prim: FactorInteger balanced 160-bit x1", () -> sf160 = smallestFactor(n160));
check("prim: FactorInteger balanced 160-bit x1", sf160 % CHECKMOD);

quit();
