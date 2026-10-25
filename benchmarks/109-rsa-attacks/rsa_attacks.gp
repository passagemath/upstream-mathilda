\\ Experiment 109 -- RSA attacks (group F).  See rsa_attacks.m.
\\
\\ PARI/GP baseline.  The "app:" attacks use plain GP integers, Mod() and
\\ bezout()/lift as a GP programmer writes them; the "prim:" continued-fraction
\\ cases use contfrac (ContinuedFraction) and contfracpnqn (Convergents, the
\\ head Mathilda lacks); the Hastad cube root uses sqrtnint, PARI's exact
\\ integer n-th root.  One statement per top-level line.

read("../harness.gp");
read("../crypto_data.gp");

require(["contfrac", "contfracpnqn", "bezout", "chinese", "sqrtnint"]);

E1 = 65537;
E2 = 65539;

\\ ---- (a) Wiener ---------------------------------------------------------

wienerKey(s) =
{
  my(p = nextprime(bigOf(s, 512)), q = nextprime(bigOf(s + 1, 512)), n, phi, d, e);
  n = p * q; phi = (p - 1) * (q - 1);
  d = bitor(bigOf(s + 2, 120), 1);
  while (gcd(d, phi) != 1, d += 2);
  e = lift(Mod(d, phi)^-1);
  [e, n, d];
}

wiener(e, n) =
{
  my(cf = contfrac(e / n), hm1 = 1, hm2 = 0, km1 = 0, km2 = 1, h, k, phi, b, disc, found = 0);
  for (i = 1, #cf,
    h = cf[i] * hm1 + hm2; k = cf[i] * km1 + km2;
    if (h > 0 && k > 0 && (e * k - 1) % h == 0,
      phi = (e * k - 1) / h; b = n - phi + 1; disc = b^2 - 4 * n;
      if (disc >= 0 && issquare(disc), found = k));
    hm2 = hm1; hm1 = h; km2 = km1; km1 = k);
  found;
}

wkeys = vector(20, s, wienerKey(1700 + 3 * (s - 1)));

bench("prim: ContinuedFraction e/N 1024-bit x600", () -> for (r = 1, 30, cfs = apply(w -> contfrac(w[1] / w[2]), wkeys)));
check("prim: ContinuedFraction e/N 1024-bit x600", (sum(i = 1, 20, vecsum(cfs[i])) + sum(i = 1, 20, #cfs[i])) % CHECKMOD);

wienerFromConv(cv, e, n) =
{
  my(found = 0, h, k, phi, b, disc, nc = matsize(cv)[2]);
  for (j = 1, nc,
    h = cv[1, j]; k = cv[2, j];
    if (h > 0 && k > 0 && (e * k - 1) % h == 0,
      phi = (e * k - 1) / h; b = n - phi + 1; disc = b^2 - 4 * n;
      if (disc >= 0 && issquare(disc), found = k)));
  found;
}

bench("prim: Convergents e/N 1024-bit x40", () -> for (r = 1, 2, cvrec = apply(w -> my(cf = contfrac(w[1] / w[2])); wienerFromConv(contfracpnqn(cf, #cf), w[1], w[2]), wkeys)));
check("prim: Convergents e/N 1024-bit x40", vecsum(cvrec) % CHECKMOD);

bench("app: Wiener attack 1024-bit x20", () -> wrec = apply(w -> wiener(w[1], w[2]), wkeys));
check("app: Wiener attack 1024-bit x20", vecsum(wrec) % CHECKMOD + CHECKMOD * (sum(i = 1, 20, wrec[i] == wkeys[i][3])));

\\ ---- (b) common modulus -------------------------------------------------

cmN = nextprime(bigOf(1801, 512)) * nextprime(bigOf(1802, 512));

commonModulus(m) =
{
  my(c1 = lift(Mod(m, cmN)^E1), c2 = lift(Mod(m, cmN)^E2), g = bezout(E1, E2));
  lift(Mod(c1, cmN)^g[1] * Mod(c2, cmN)^g[2]);
}

cmMsgs = vector(100, k, bigOf(1809 + k, 900) % cmN);
bench("app: common-modulus attack 1024-bit x600", () -> for (r = 1, 6, cmRec = apply(commonModulus, cmMsgs)));
check("app: common-modulus attack 1024-bit x600", vecsum(cmRec) % CHECKMOD + CHECKMOD * (sum(k = 1, 100, cmRec[k] == cmMsgs[k])));

\\ ---- (c) Hastad broadcast (e = 3) ---------------------------------------

hMods = [nextprime(bigOf(1900, 512)) * nextprime(bigOf(1901, 512)), nextprime(bigOf(1902, 512)) * nextprime(bigOf(1903, 512)), nextprime(bigOf(1904, 512)) * nextprime(bigOf(1905, 512))];

hastad(m) =
{
  my(cs = apply(n -> Mod(m, n)^3, hMods));
  sqrtnint(lift(chinese(cs)), 3);
}

hMsgs = vector(50, k, bigOf(1949 + k, 300));
bench("app: Hastad broadcast e=3 x1500", () -> for (r = 1, 30, hRec = apply(hastad, hMsgs)));
check("app: Hastad broadcast e=3 x1500", vecsum(hRec) % CHECKMOD + CHECKMOD * (sum(k = 1, 50, hRec[k] == hMsgs[k])));
quit();
