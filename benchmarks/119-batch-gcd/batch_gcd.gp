\\ Experiment 119 -- batch GCD over RSA moduli (group F).  See batch_gcd.m.
\\
\\ PARI/GP baseline.  The "app:" scan is written with plain GP integers
\\ (pairwise products, %, \, gcd), exactly the .m algorithm.  The "prim:"
\\ product is vecprod (PARI's own balanced product), and the GCD vector is
\\ gcd() applied pairwise.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["vecprod", "gcd", "nextprime"]);

MM = 1024;
bgP = vector(MM, i, primeOf(5000 + 2 * i, 512));
bgQ = vector(MM, i, if (i % 37 == 0, bgP[i - 17], primeOf(5001 + 2 * i, 512)));
moduli = vector(MM, i, bgP[i] * bgQ[i]);
moduli256 = moduli[1..256];

\\ ---- the algorithm, written in the language ---------------------------------

productTree(ns) =
{
  my(tree = List([ns]));
  while (#ns > 1,
    ns = vector((#ns + 1) \ 2, i, if (2 * i <= #ns, ns[2*i-1] * ns[2*i], ns[2*i-1]));
    listput(tree, ns));
  Vec(tree);
}

remainderTree(tree) =
{
  my(r = tree[#tree], level);
  forstep (k = #tree - 1, 1, -1,
    level = tree[k];
    r = vector(#level, j, r[(j + 1) \ 2] % level[j]^2));
  r;
}

batchGCD(ns) = my(z = remainderTree(productTree(ns))); vector(#ns, i, gcd(ns[i], z[i] \ ns[i]));

weakDigest(g) = my(w = select(x -> x > 1, g)); #w * 10^10 + vecsum(w) % CHECKMOD;

\\ ---- prim: the big-integer kernels ------------------------------------------

bench("prim: Times@@ 1024 moduli (1M-bit product) x4", () -> for (r = 1, 4, bigProd = vecprod(moduli)));
check("prim: Times@@ 1024 moduli (1M-bit product) x4", bigProd % CHECKMOD);

tree1024 = productTree(moduli);
z1024 = remainderTree(tree1024);
zq1024 = vector(MM, i, z1024[i] \ moduli[i]);
bench("prim: GCD 1024-bit pairs x1024x8", () -> for (r = 1, 8, gFinal = vector(MM, i, gcd(moduli[i], zq1024[i]))));
check("prim: GCD 1024-bit pairs x1024x8", weakDigest(gFinal));

\\ ---- app: the two tree phases, then the whole scan ---------------------------

bench("app: product tree M=1024 x4", () -> for (r = 1, 4, tr = productTree(moduli)));
check("app: product tree M=1024 x4", tr[#tr][1] % CHECKMOD + 10^10 * #tr);

bench("app: remainder tree mod N^2 M=1024", () -> zz = remainderTree(tree1024));
check("app: remainder tree mod N^2 M=1024", vecsum(zz) % CHECKMOD);

bench("app: batch GCD M=256 x4", () -> for (r = 1, 4, g256 = batchGCD(moduli256)));
check("app: batch GCD M=256 x4", weakDigest(g256));

bench("app: batch GCD M=1024", () -> g1024 = batchGCD(moduli));
check("app: batch GCD M=1024", weakDigest(g1024));

naiveGCD(ns) =
{
  my(m = #ns, g = vector(#ns, i, 1), d);
  for (i = 1, m - 1, for (j = i + 1, m,
    d = gcd(ns[i], ns[j]);
    if (d > 1, g[i] *= d; g[j] *= d)));
  g;
}

benchOnce("app: naive pairwise GCD M=256 (32640 GCDs)", () -> gNaive = naiveGCD(moduli256));
check("app: naive pairwise GCD M=256 (32640 GCDs)", weakDigest(gNaive));

quit();
