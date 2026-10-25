\\ Experiment 113 -- lattice attacks (group F).  See lattice_attacks.m.
\\
\\ PARI/GP baseline.  qflll is PARI's LLL (delta = 0.99, fplll-style floating
\\ point with exact fallback).  qflll reduces the COLUMNS of its argument, so
\\ every row basis B is passed as B~ and the reduced basis read back as
\\ qflll(B~, 3)~ (flag 3: integral input, returns the reduced basis itself).
\\ Everything around the reduction is plain GP integers / Mod.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["qflll", "Mod^-1", "fromdigits", "matdet"]);

lll(B) = qflll(B~, 3)~;

\\ ---- input vectors: NOT raw LCG output ------------------------------------
\\ (see lattice_attacks.m: LCG output is linear and LLL sees it, so each LCG
\\ word goes through g^e mod P before truncation to a bits-bit integer)

mixBigs(seed, cnt, bits) =
{
  my(p = primeOf(seed, bits + 64));
  apply(e -> lift(Mod(3, p)^e) % 2^(bits - 1) + 2^(bits - 1), lcgStates(seed, cnt));
}

\\ ---- (a) Merkle-Hellman knapsack ------------------------------------------

mhKey(seed, n, k) =
{
  my(r, w = vector(n), acc = 0, q, m);
  r = mixBigs(seed, n, k);
  for (i = 1, n, w[i] = acc + r[i]; acc += w[i]);
  q = primeOf(seed, k + n + 1);
  m = bigOf(seed + 1000, k + n - 1) % q;
  [w, q, m, apply(x -> m * x % q, w)];
}

mhMessage(seed, n) = smallList(seed, n, 2);
mhEncrypt(a, x) = sum(i = 1, #a, x[i] * a[i]);

mhDecrypt(key, c) =
{
  my(w = key[1], q = key[2], m = key[3], t, n = #key[1], x);
  t = lift(Mod(c, q) * Mod(m, q)^-1);
  x = vector(n);
  forstep (i = n, 1, -1, if (t >= w[i], x[i] = 1; t -= w[i]));
  fromdigits(x, 2);
}

cjlossBasis(a, s) =
{
  my(n = #a, B = matrix(n + 1, n + 1));
  for (i = 1, n, B[i, i] = 2; B[i, n + 1] = n * a[i]);
  for (j = 1, n, B[n + 1, j] = 1);
  B[n + 1, n + 1] = n * s;
  B;
}

cjlossRecover(a, s, red) =
{
  my(n = #a, v, x, found = -1);
  for (j = 1, #red[, 1],
    v = red[j, ];
    if (found < 0 && v[n + 1] == 0 && Set(apply(abs, v[1..n])) == [1],
      x = vector(n, i, (1 + v[i]) / 2);
      if (mhEncrypt(a, x) == s, found = fromdigits(x, 2),
        x = vector(n, i, (1 - v[i]) / 2);
        if (mhEncrypt(a, x) == s, found = fromdigits(x, 2)))));
  found;
}

mhBreak(a, s) = cjlossRecover(a, s, lll(cjlossBasis(a, s)));

key24 = mhKey(11, 24, 56);
key32 = mhKey(13, 32, 40);
cts24 = vector(20, j, mhEncrypt(key24[4], mhMessage(2000 + j, 24)));
cts32 = vector(12, j, mhEncrypt(key32[4], mhMessage(3000 + j, 32)));

benchOnce("prim: LatticeReduce CJLOSS knapsack n=24 x20", () -> brk24 = apply(s -> mhBreak(key24[4], s), cts24));
check("prim: LatticeReduce CJLOSS knapsack n=24 x20", vecsum(brk24) % CHECKMOD);

benchOnce("prim: LatticeReduce CJLOSS knapsack n=32 x12", () -> brk32 = apply(s -> mhBreak(key32[4], s), cts32));
check("prim: LatticeReduce CJLOSS knapsack n=32 x12", vecsum(brk32) % CHECKMOD);

ctsDec = vector(2000, j, mhEncrypt(key32[4], mhMessage(5000 + j, 32)));
bench("app: Merkle-Hellman decrypt n=32 x2000", () -> dec32 = apply(c -> mhDecrypt(key32, c), ctsDec));
check("app: Merkle-Hellman decrypt n=32 x2000", vecsum(dec32) % CHECKMOD);

\\ ---- (b) approximate-GCD / hidden-number secret recovery ------------------

agcdSamples(seed, t) =
{
  my(p = primeOf(seed, 128));
  [p, p * mixBigs(seed + 10, t + 1, 128) + mixBigs(seed + 500, t + 1, 16)];
}

agcdRecover(xs) =
{
  my(t = #xs - 1, x0 = xs[1], B = matrix(t + 1, t + 1), red, v);
  B[1, 1] = 2^17;
  for (j = 2, t + 1, B[1, j] = xs[j]);
  for (i = 2, t + 1, B[i, i] = -x0);
  red = lll(B);
  v = abs(red[1, 1]);
  if (v % 2^17 == 0 && v > 0, x0 \ (v / 2^17), -1);
}

agcd = vector(40, j, agcdSamples(7000 + 1000 * j, 8));
benchOnce("app: approximate-GCD secret recovery dim=9 x40", () -> ps = apply(s -> agcdRecover(s[2]), agcd));
check("app: approximate-GCD secret recovery dim=9 x40", vecsum(ps) % CHECKMOD + if (ps == apply(s -> s[1], agcd), 0, 1));

\\ ---- (c) plain LLL on a square knapsack-type basis ------------------------

qaryBasis(seed, n, bits) =
{
  my(q = primeOf(seed, bits), B = matrix(n, n), a = mixBigs(seed + 1, n - 1, bits));
  B[1, 1] = q;
  for (i = 1, n - 1, B[i + 1, 1] = a[i] % q; B[i + 1, i + 1] = 1);
  B;
}

qa30 = qaryBasis(200000, 30, 300);
qa50 = qaryBasis(300000, 50, 200);

benchOnce("prim: LatticeReduce q-ary basis n=30 300-bit", () -> red30 = lll(qa30));
check("prim: LatticeReduce q-ary basis n=30 300-bit", abs(matdet(red30)) % CHECKMOD);

benchOnce("prim: LatticeReduce q-ary basis n=50 200-bit", () -> red50 = lll(qa50));
check("prim: LatticeReduce q-ary basis n=50 200-bit", abs(matdet(red50)) % CHECKMOD);

quit();
