\\ Experiment 118 -- quadratic-residuosity cryptography (group F).  See qr_prngs.m.
\\
\\ PARI/GP baseline.  The "app:" cases use plain GP integers (x^2 % n) and
\\ Mod(c, p)^e for Euler's criterion; the symbol "prim:" cases use
\\ kronecker(), PARI's one symbol builtin (Jacobi = Kronecker on odd moduli).
\\ The Blum primes come from the same deterministic search driven by
\\ ispseudoprime (BPSW).
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");
read("../crypto_data.gp");

require(["kronecker", "ispseudoprime", "Mod^e", "fromdigits"]);

blumPrime(start) =
{
  my(c = start - start % 4 + 3);
  while (!ispseudoprime(c), c += 4);
  c;
}

bbsP = blumPrime(bigOf(1181, 512));
bbsQ = blumPrime(bigOf(1182, 512));
bbsN = bbsP * bbsQ;
bbsX0 = bigOf(1183, 1000)^2 % bbsN;

\\ ---- app: Blum-Blum-Shub ----------------------------------------------------

bbsParity(x, n, steps) =
{
  my(h = 0);
  for (i = 1, steps, x = x^2 % n; h = (2 * h + x % 2) % CHECKMOD);
  h;
}

bbsLow10(x, n, steps) =
{
  my(h = 0);
  for (i = 1, steps, x = x^2 % n; h = (1024 * h + x % 1024) % CHECKMOD);
  h;
}

bench("app: BBS-1024 parity x20000 bits", () -> bbs1 = bbsParity(bbsX0, bbsN, 20000));
check("app: BBS-1024 parity x20000 bits", bbs1);

bench("app: BBS-1024 10 bits/step x20000 steps", () -> bbs2 = bbsLow10(bbsX0, bbsN, 20000));
check("app: BBS-1024 10 bits/step x20000 steps", bbs2);

\\ ---- app / prim: Goldwasser-Micali -----------------------------------------

gmY = bbsN - 1;
gmMsg = smallList(1184, 2000, 2);
gmR = vector(2000, i, bigOf(2000 + i, 1000));

bench("app: GM-1024 encrypt x2000 bits x10", () -> for (r = 1, 10, gmC = vector(2000, i, gmR[i]^2 * gmY^gmMsg[i] % bbsN)));
check("app: GM-1024 encrypt x2000 bits x10", vecsum(gmC) % CHECKMOD);

bench("prim: JacobiSymbol GM-1024 decrypt x2000 bits x4", () -> for (r = 1, 4, gmD1 = apply(c -> (1 - kronecker(c, bbsP)) / 2, gmC)));
check("prim: JacobiSymbol GM-1024 decrypt x2000 bits x4", fromdigits(gmD1, 2) % CHECKMOD + if (gmD1 == gmMsg, 10^10, 0));

gmE = (bbsP - 1) / 2;
bench("app: GM-1024 decrypt Euler criterion x2000 bits", () -> gmD2 = apply(c -> if (Mod(c, bbsP)^gmE == 1, 0, 1), gmC));
check("app: GM-1024 decrypt Euler criterion x2000 bits", fromdigits(gmD2, 2) % CHECKMOD + if (gmD2 == gmMsg, 10^10, 0));

\\ ---- prim: the symbol kernels on 1024-bit operands -------------------------

jacA = vector(5000, i, bigOf(3000 + i, 1024));
jacM = vector(5000, i, my(b = bigOf(9000 + i, 1024)); b + 1 - b % 2);
bench("prim: JacobiSymbol 1024-bit x5000x4", () -> for (r = 1, 4, jac = vector(5000, i, kronecker(jacA[i], jacM[i]))));
check("prim: JacobiSymbol 1024-bit x5000x4", vecsum(jac) % CHECKMOD);

kroA = vector(5000, i, (-1)^i * jacA[i]);
kroM = vector(5000, i, jacM[i] - i % 2);
bench("prim: KroneckerSymbol signed/even 1024-bit x5000x4", () -> for (r = 1, 4, kro = vector(5000, i, kronecker(kroA[i], kroM[i]))));
check("prim: KroneckerSymbol signed/even 1024-bit x5000x4", vecsum(kro) % CHECKMOD);

quit();
