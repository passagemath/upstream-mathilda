\\ Experiment 112 -- Rabin cryptosystem & modular square roots (group F).
\\ See rabin_sqrt_mod.m.
\\
\\ PARI/GP baseline.  Tonelli-Shanks and Cipolla as plain GP ("app:"); the
\\ native primitive ("prim:") is sqrt(Mod(a, p)), PARI's modular square root.
\\ The "prim: PowerMod[a,1/2,p]" and "prim: SqrtMod" labels both call
\\ sqrt(Mod(a, p)) here -- PowerMod[a,1/2,p] is the WL spelling of this very
\\ primitive, SqrtMod a name WL lacks.  Every root is canonicalised to
\\ min(r, p-r).  One statement per top-level line.

read("../harness.gp");
read("../crypto_data.gp");

require(["sqrt(Mod)", "ispseudoprime", "nextprime"]);

primeMod(seed, bits, r, m) =
{
  my(p = primeOf(seed, bits));
  while (p % m != r, p = nextprime(p + 1));
  p;
}

canon(r, p) = min(r, p - r);

tonelli(n, p) =
{
  my(Q = p - 1, S = 0, z, M, c, t, R, i, b, tt);
  while (Q % 2 == 0, Q /= 2; S++);
  if (S == 1, return(lift(Mod(n, p)^((p + 1) / 4))));
  z = 2; while (lift(Mod(z, p)^((p - 1) / 2)) != p - 1, z++);
  M = S; c = lift(Mod(z, p)^Q); t = lift(Mod(n, p)^Q); R = lift(Mod(n, p)^((Q + 1) / 2));
  while (t != 1,
    i = 0; tt = t; while (tt != 1, tt = (tt * tt) % p; i++);
    b = lift(Mod(c, p)^(2^(M - i - 1))); M = i; c = (b * b) % p;
    t = (t * c) % p; R = (R * b) % p);
  R;
}

qr(x, p) = x == 0 || lift(Mod(x, p)^((p - 1) / 2)) == 1;

fmul(a, b, w, p) = [(a[1] * b[1] + a[2] * b[2] * w) % p, (a[1] * b[2] + a[2] * b[1]) % p];

fpow(base, e, w, p) =
{
  my(result = [1, 0], b = base);
  while (e > 0,
    if (e % 2 == 1, result = fmul(result, b, w, p));
    b = fmul(b, b, w, p); e = e \ 2);
  result;
}

cipolla(n, p) =
{
  my(a0 = 1, w);
  while (qr((a0 * a0 - n) % p, p), a0++);
  w = (a0 * a0 - n) % p;
  fpow([a0, 1], (p + 1) / 2, w, p)[1];
}

P256 = primeMod(1201, 256, 1, 8);
as256 = vector(100, k, lift(Mod(bigOf(k + 9, 256), P256)^2));
P1024 = primeMod(1251, 1024, 1, 8);
as1024 = vector(10, k, lift(Mod(bigOf(k + 9, 1024), P1024)^2));

bench("app: Tonelli-Shanks 256-bit p=1mod8 x100", () -> ts256 = apply(a -> canon(tonelli(a, P256), P256), as256));
check("app: Tonelli-Shanks 256-bit p=1mod8 x100", vecsum(ts256) % CHECKMOD);

benchOnce("app: Tonelli-Shanks 1024-bit p=1mod8 x10", () -> ts1024 = apply(a -> canon(tonelli(a, P1024), P1024), as1024));
check("app: Tonelli-Shanks 1024-bit p=1mod8 x10", vecsum(ts1024) % CHECKMOD);

bench("app: Cipolla 256-bit p=1mod8 x100", () -> cp256 = apply(a -> canon(cipolla(a, P256), P256), as256));
check("app: Cipolla 256-bit p=1mod8 x100", vecsum(cp256) % CHECKMOD);

benchOnce("app: Cipolla 1024-bit p=1mod8 x10", () -> cp1024 = apply(a -> canon(cipolla(a, P1024), P1024), as1024));
check("app: Cipolla 1024-bit p=1mod8 x10", vecsum(cp1024) % CHECKMOD);

bench("prim: PowerMod[a,1/2,p] 256-bit x500", () -> for (r = 1, 5, pm256 = apply(a -> canon(lift(sqrt(Mod(a, P256))), P256), as256)));
check("prim: PowerMod[a,1/2,p] 256-bit x500", vecsum(pm256) % CHECKMOD);

bench("prim: PowerMod[a,1/2,p] 1024-bit x40", () -> for (r = 1, 4, pm1024 = apply(a -> canon(lift(sqrt(Mod(a, P1024))), P1024), as1024)));
check("prim: PowerMod[a,1/2,p] 1024-bit x40", vecsum(pm1024) % CHECKMOD);

bench("prim: SqrtMod 256-bit x500", () -> for (r = 1, 5, sm256 = apply(a -> canon(lift(sqrt(Mod(a, P256))), P256), as256)));
check("prim: SqrtMod 256-bit x500", vecsum(sm256) % CHECKMOD);

PR = primeMod(1301, 256, 3, 4);
QR = primeMod(1302, 256, 3, 4);
NR = PR * QR;
YPR = lift(Mod(PR, QR)^-1);
YQR = lift(Mod(QR, PR)^-1);

rabinRoots(c) =
{
  my(mp = lift(Mod(c, PR)^((PR + 1) / 4)), mq = lift(Mod(c, QR)^((QR + 1) / 4)), u, v);
  u = (mp * QR * YQR + mq * PR * YPR) % NR;
  v = (mp * QR * YQR - mq * PR * YPR) % NR;
  (canon(u, NR) + canon(v, NR)) % NR;
}

msgs = vector(500, k, bigOf(k + 499, 256) % NR);
cts = apply(m -> (m * m) % NR, msgs);
bench("app: Rabin decrypt 512-bit n x500", () -> rr = apply(rabinRoots, cts));
check("app: Rabin decrypt 512-bit n x500", vecsum(rr) % CHECKMOD);

quit();
