\\ Experiment 111 -- Shamir secret sharing (group F).  See shamir_secret_sharing.m.
\\
\\ PARI/GP baseline.  The "app:" protocol uses plain GP integers with % p and
\\ Mod(d, p)^-1, as a GP programmer writes it.  The "prim:" interpolation is
\\ polinterpolate(xs, ys, 0) on Mod(y, p) entries, so -- unlike the WL and
\\ sympy routes, which go through Q -- it stays in F_p throughout.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Mod^-1", "polinterpolate"]);

P = primeOf(1111, 256);
T = 50;
N = 200;

polyOf(k) = vector(T, j, bigOf(100000 * k + j - 1, 256) % P);

hornerMod(cs, x) =
{
  my(acc = 0);
  forstep (i = #cs, 1, -1, acc = (acc * x + cs[i]) % P);
  acc;
}

split(cs) = vector(N, i, [i, hornerMod(cs, i)]);

lagrange0(sh) =
{
  my(xs = vector(#sh, i, sh[i][1]), total = 0, num, den);
  for (i = 1, #xs,
    num = 1; den = 1;
    for (j = 1, #xs,
      if (j != i,
        num = (num * xs[j]) % P;
        den = (den * (xs[j] - xs[i])) % P));
    total += (sh[i][2] * num * lift(Mod(den, P)^-1)) % P);
  total % P;
}

pick(sh) = vector(T, k, sh[3 + 4 * (k - 1)]);

polys = vector(20, k, polyOf(k));

\\ ---- app ----------------------------------------------------------------

bench("app: split 256-bit n=200 t=50 x20", () -> shares = apply(split, polys));
check("app: split 256-bit n=200 t=50 x20", sum(a = 1, 20, sum(b = 1, N, shares[a][b][2])) % CHECKMOD);

bench("app: reconstruct Lagrange t=50 x20", () -> recs = vector(20, k, lagrange0(pick(shares[k]))));
check("app: reconstruct Lagrange t=50 x20", vecsum(recs) % CHECKMOD);

batchShares = vector(2000, k, apply(s -> [s[1], (s[2] + k) % P], pick(shares[k % 20 + 1])));

batchRec() =
{
  my(xs0 = vector(T, i, batchShares[1][i][1]), w);
  w = vector(T, i, (prod(j = 1, T, if (j != i, xs0[j], 1)) * lift(Mod(prod(j = 1, T, if (j != i, xs0[j] - xs0[i], 1)), P)^-1)) % P);
  brec = vector(2000, k, sum(i = 1, T, batchShares[k][i][2] * w[i]) % P);
}

bench("app: batch reconstruct t=50 x2000", batchRec);
check("app: batch reconstruct t=50 x2000", vecsum(brec) % CHECKMOD);

\\ ---- prim ---------------------------------------------------------------

interp0(sh) = lift(polinterpolate(vector(#sh, i, sh[i][1]), vector(#sh, i, Mod(sh[i][2], P)), 0));

bench("prim: InterpolatingPolynomial t=50 mod p x20", () -> irecs = vector(20, k, interp0(pick(shares[k]))));
check("prim: InterpolatingPolynomial t=50 mod p x20", vecsum(irecs) % CHECKMOD);

benchOnce("prim: InterpolatingPolynomial n=200 mod p x4", () -> for (r = 1, 4, i200 = interp0(shares[r])));
check("prim: InterpolatingPolynomial n=200 mod p x4", i200 % CHECKMOD);

quit();
