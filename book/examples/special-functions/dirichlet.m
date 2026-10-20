# The non-principal character mod 4
Table[DirichletCharacter[4, 2, n], {n, 1, 8}]
# A complex character mod 5: values are roots of unity
Table[DirichletCharacter[5, 2, n], {n, 1, 5}]
# The real character mod 5 is the Legendre symbol
Table[DirichletCharacter[5, 3, n], {n, 1, 5}] == Table[JacobiSymbol[n, 5], {n, 1, 5}]
# Exact values at s = 1 (non-principal) and s <= 0
DirichletL[4, 2, 1]
DirichletL[3, 2, 1]
DirichletL[4, 2, -2]
# The principal character mod 1 is the Riemann zeta function
DirichletL[1, 1, 2]
# The principal character mod 4 is zeta with the factor 2 removed
DirichletL[4, 1, -1]
# L(2, chi_4) is Catalan's constant
N[DirichletL[4, 2, 2], 30]
N[DirichletL[4, 2, 2], 30] == N[Catalan, 30]
# The real character mod 5 at s = 1: 2 Log[GoldenRatio]/Sqrt[5]
DirichletL[5, 3, 1]
N[DirichletL[5, 3, 1], 30] == N[2 Log[GoldenRatio]/Sqrt[5], 30]
# Machine and arbitrary precision, on and off the critical line
DirichletL[5, 3, 0.5]
N[DirichletL[5, 3, 1/2], 40]
DirichletL[7, 2, 1/2 + 3. I]
# A rigorous enclosure of L(s, chi_4) for s in [2, 3]
DirichletL[4, 2, Interval[{2, 3}]]
