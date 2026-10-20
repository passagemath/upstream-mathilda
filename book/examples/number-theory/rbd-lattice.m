# Section 5.6.4: a planted miss, beyond the gcd window, found by the lattice detector.
p = NextPrime[Floor[3 (22/7)^59] + 2^28];
q = NextPrime[2 10^29];
n = p q
Floor[(7/22)^59 n] - 3 q
Quiet[FactorInteger`BlakeRationalBaseDescentClassic[n, "Base" -> 22/7]]
FactorInteger`BlakeRationalBaseDescentLattice[n, "Base" -> 22/7]
FactorInteger[n, Method -> {"BlakeRationalBaseDescent", "Variant" -> "Coppersmith", "Base" -> 22/7, "Coefficient" -> 3, "Exponent" -> 59}]
