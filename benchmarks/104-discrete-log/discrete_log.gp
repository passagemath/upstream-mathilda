\\ Experiment 104 -- discrete logarithms (group F).  See discrete_log.m.
\\
\\ PARI/GP baseline.  "prim:" = znlog(h, Mod(g, p), r) (PARI's generic
\\ discrete log: Pohlig-Hellman, then BSGS / rho / index calculus per prime
\\ factor; the subgroup order r is passed when known, as a GP user would).
\\ "app:" = the same BSGS (Map(), mapput / mapisdefined, and a bulk Map(matrix)),
\\ Pollard rho and Pohlig-Hellman (chinese) as the .m half, in plain GP integers.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["znlog", "chinese", "Map", "factor"]);

\\ ---- instances (untimed) --------------------------------------------------

subgroup(bits, s) =
{
  my(r = primeOf(s, bits), k = 1, p, b = 2);
  while (!ispseudoprime(2*k*r + 1), k++);
  p = 2*k*r + 1;
  while (lift(Mod(b, p)^((p - 1) / r)) == 1, b++);
  [p, r, lift(Mod(b, p)^((p - 1) / r))];
}

instance(bits, s, t) =
{
  my(v = subgroup(bits, s), x);
  x = bigOf(t, 64) % v[2];
  [v[1], v[2], v[3], lift(Mod(v[3], v[1])^x), x];
}

inst24 = vector(4, i, instance(24, 8999 + i, 9099 + i));
inst32 = vector(2, i, instance(32, 9199 + i, 9299 + i));

qs = vector(8, i, nextprime(40000 + 5000 * (i - 1) + 1));
smoothL = vecprod(qs);
smoothK = 1;
while (!ispseudoprime(2 * smoothK * smoothL + 1), smoothK++);
pS = 2 * smoothK * smoothL + 1;
fS = factor(pS - 1);
gS = 2;
while (sum(i = 1, #fS~, lift(Mod(gS, pS)^((pS - 1) / fS[i, 1])) == 1) > 0, gS++);
xsSmooth = vector(10, i, bigOf(9400 + i, 128) % (pS - 1));
hsSmooth = apply(x -> lift(Mod(gS, pS)^x), xsSmooth);

\\ ---- app: baby-step giant-step --------------------------------------------

bsgsInc(g, h, p, r) =
{
  my(m = sqrtint(r) + 1, tab = Map(), y = 1, c, j = 0, v);
  for (i = 0, m - 1, mapput(tab, y, i); y = y * g % p);
  c = lift(Mod(g, p)^-m);
  y = h;
  while (!mapisdefined(tab, y, &v), y = y * c % p; j++);
  j * m + v;
}

bsgs(g, h, p, r) =
{
  my(m = sqrtint(r) + 1, pw, tab, c, y = h, j = 0, v);
  pw = vector(m); pw[1] = 1;
  for (i = 2, m, pw[i] = pw[i - 1] * g % p);
  tab = Map(matconcat([pw~, [0 .. m - 1]~]));
  c = lift(Mod(g, p)^-m);
  while (!mapisdefined(tab, y, &v), y = y * c % p; j++);
  j * m + v;
}

bench("app: BSGS 24-bit subgroup (incremental table) x4", () -> b1 = apply(v -> bsgsInc(v[3], v[4], v[1], v[2]), inst24));
check("app: BSGS 24-bit subgroup (incremental table) x4", vecsum(b1) % CHECKMOD);

bench("app: BSGS 24-bit subgroup (bulk table) x4", () -> b2 = apply(v -> bsgs(v[3], v[4], v[1], v[2]), inst24));
check("app: BSGS 24-bit subgroup (bulk table) x4", vecsum(b2) % CHECKMOD);

bench("app: BSGS 32-bit subgroup (bulk table) x2", () -> b3 = apply(v -> bsgs(v[3], v[4], v[1], v[2]), inst32));
check("app: BSGS 32-bit subgroup (bulk table) x2", vecsum(b3) % CHECKMOD);

\\ ---- app: Pollard rho for logarithms ---------------------------------------

rhoStep(st, g, h, p, r) =
{
  my(x = st[1], a = st[2], b = st[3], s = x % 3);
  if (s == 0, [x^2 % p, 2*a % r, 2*b % r],
      s == 1, [x*g % p, (a + 1) % r, b],
              [x*h % p, a, (b + 1) % r]);
}

pollardRho(g, h, p, r) =
{
  my(t = 0, ans = -1, tort, hare, den);
  while (ans < 0,
    tort = [lift(Mod(g, p)^t), t, 0];
    hare = tort;
    tort = rhoStep(tort, g, h, p, r);
    hare = rhoStep(rhoStep(hare, g, h, p, r), g, h, p, r);
    while (tort[1] != hare[1],
      tort = rhoStep(tort, g, h, p, r);
      hare = rhoStep(rhoStep(hare, g, h, p, r), g, h, p, r));
    den = (tort[3] - hare[3]) % r;
    if (den != 0, ans = ((hare[2] - tort[2]) * lift(Mod(den, r)^-1)) % r, t++));
  ans;
}

benchOnce("app: Pollard rho 32-bit subgroup x2", () -> r1 = apply(v -> pollardRho(v[3], v[4], v[1], v[2]), inst32));
check("app: Pollard rho 32-bit subgroup x2", vecsum(r1) % CHECKMOD);

\\ ---- app: Pohlig-Hellman over a smooth p ------------------------------------

pohligHellman(g, h, p, fs) =
{
  my(n = p - 1, res = Mod(0, 1), q, e, gq, xk, hk);
  for (k = 1, #fs~,
    q = fs[k, 1]; e = fs[k, 2];
    gq = lift(Mod(g, p)^(n / q));
    xk = 0;
    for (j = 0, e - 1,
      hk = lift(Mod(lift(Mod(g, p)^-xk) * h, p)^(n / q^(j + 1)));
      xk += bsgs(gq, hk, p, q) * q^j);
    res = chinese(res, Mod(xk, q^e)));
  lift(res);
}

bench("app: Pohlig-Hellman 131-bit smooth p x10", () -> ph = apply(h -> pohligHellman(gS, h, pS, fS), hsSmooth));
check("app: Pohlig-Hellman 131-bit smooth p x10", vecsum(ph) % CHECKMOD);

\\ ---- prim: znlog ------------------------------------------------------------

bench("prim: discrete log 24-bit subgroup x40", () -> for (rr = 1, 10, d1 = apply(v -> znlog(v[4], Mod(v[3], v[1]), v[2]), inst24)));
check("prim: discrete log 24-bit subgroup x40", vecsum(d1) % CHECKMOD);

bench("prim: discrete log 32-bit subgroup x10", () -> for (rr = 1, 5, d2 = apply(v -> znlog(v[4], Mod(v[3], v[1]), v[2]), inst32)));
check("prim: discrete log 32-bit subgroup x10", vecsum(d2) % CHECKMOD);

bench("prim: discrete log 131-bit smooth p x40", () -> for (rr = 1, 4, d3 = apply(h -> znlog(h, Mod(gS, pS)), hsSmooth)));
check("prim: discrete log 131-bit smooth p x40", vecsum(d3) % CHECKMOD);

quit();
