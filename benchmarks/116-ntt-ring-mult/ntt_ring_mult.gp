\\ Experiment 116 -- NTT ring multiplication (group F).  See ntt_ring_mult.m.
\\
\\ PARI/GP baseline.  "app:" the same iterative radix-2 Cooley-Tukey NTT as
\\ the .m / .py halves, written with GP vectors and integer % q.  "prim:"
\\ PARI's polynomial arithmetic: the ListConvolve case multiplies integer
\\ polynomials (Polrev) and folds by x^n = -1; the PolynomialRemainder case
\\ multiplies Mod(., q)-coefficient polynomials and reduces % (x^n + 1).
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["Polrev", "znprimroot", "Mod^-1", "binary"]);

nttSetup(q, n) =
{
  my(g = lift(znprimroot(q)), psi, psiInv, w, wInv, k = logint(n, 2), nInv);
  psi = lift(Mod(g, q)^((q - 1) / (2 * n))); psiInv = lift(Mod(psi, q)^-1);
  w = psi^2 % q; wInv = lift(Mod(w, q)^-1); nInv = lift(Mod(n, q)^-1);
  [q, n,
   vector(n, i, fromdigits(Vecrev(binary(i - 1 + 2^k))[1 .. k], 2) + 1),
   vector(n / 2, j, lift(Mod(w, q)^(j - 1))),
   vector(n / 2, j, lift(Mod(wInv, q)^(j - 1))),
   vector(n, j, lift(Mod(psi, q)^(j - 1))),
   vector(n, j, nInv * lift(Mod(psiInv, q)^(j - 1)) % q)];
}

\\ ctx = [q, n, perm, tw, twInv, psiPow, postPow]
ntt(a, tw, ctx) =
{
  my(q = ctx[1], n = ctx[2], perm = ctx[3], x, h = 1, step, u, v);
  x = vector(n, i, a[perm[i]]);
  while (h < n,
    step = n / (2 * h);
    forstep (s = 0, n - 1, 2 * h,
      for (j = 0, h - 1,
        u = x[s + j + 1];
        v = x[s + j + h + 1] * tw[j * step + 1] % q;
        x[s + j + 1] = (u + v) % q;
        x[s + j + h + 1] = (u - v) % q));
    h *= 2);
  x;
}

nttMul(a, b, ctx) =
{
  my(q = ctx[1], n = ctx[2], fa, fb, c);
  fa = ntt(vector(n, i, a[i] * ctx[6][i] % q), ctx[4], ctx);
  fb = ntt(vector(n, i, b[i] * ctx[6][i] % q), ctx[4], ctx);
  c = ntt(vector(n, i, fa[i] * fb[i] % q), ctx[5], ctx);
  vector(n, i, c[i] * ctx[7][i] % q);
}

convMul(a, b, q, n) =
{
  my(full = Vecrev(Polrev(a) * Polrev(b)), m = #full);
  vector(n, i, ((if (i <= m, full[i], 0)) - (if (i + n <= m, full[i + n], 0))) % q);
}

polyMul(a, b, q, n) =
{
  my(r = (Mod(1, q) * Polrev(a) * Polrev(b)) % (x^n + 1), v = lift(Vecrev(r)));
  concat(v, vector(n - #v));
}

digest(cs) = sum(k = 1, #cs, sum(i = 1, #cs[k], i * cs[k][i])) % CHECKMOD;

ctxA = nttSetup(12289, 1024);
ctxB = nttSetup(7681, 256);
pairsA = vector(20, k, [smallList(1600 + 2 * k, 1024, 12289), smallList(1601 + 2 * k, 1024, 12289)]);
pairsB = vector(100, k, [smallList(1700 + 2 * k, 256, 7681), smallList(1701 + 2 * k, 256, 7681)]);

bench("app: NTT negacyclic mult q=12289 n=1024 x20", () -> cA1 = apply(p -> nttMul(p[1], p[2], ctxA), pairsA));
check("app: NTT negacyclic mult q=12289 n=1024 x20", digest(cA1));

bench("prim: ListConvolve negacyclic q=12289 n=1024 x4", () -> cA2 = apply(p -> convMul(p[1], p[2], 12289, 1024), pairsA[1 .. 4]));
check("prim: ListConvolve negacyclic q=12289 n=1024 x4", digest(cA2));

bench("app: NTT negacyclic mult q=7681 n=256 x100", () -> cB1 = apply(p -> nttMul(p[1], p[2], ctxB), pairsB));
check("app: NTT negacyclic mult q=7681 n=256 x100", digest(cB1));

bench("prim: ListConvolve negacyclic q=7681 n=256 x100", () -> cB2 = apply(p -> convMul(p[1], p[2], 7681, 256), pairsB));
check("prim: ListConvolve negacyclic q=7681 n=256 x100", digest(cB2));

pairsC = vector(2, k, [smallList(1900 + 2 * k, 128, 7681), smallList(1901 + 2 * k, 128, 7681)]);
benchOnce("prim: PolynomialRemainder negacyclic q=7681 n=128 x2", () -> cC3 = apply(p -> polyMul(p[1], p[2], 7681, 128), pairsC));
check("prim: PolynomialRemainder negacyclic q=7681 n=128 x2", digest(cC3));

quit();
