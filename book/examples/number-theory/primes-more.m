# Section 5.6.2: random primes, prime powers, and fast prime counting.
SeedRandom[2026]; RandomPrime[1000, 6]
SeedRandom[2026]; RandomPrime[1000, 6]
PrimeQ[RandomPrime[{10^99, 10^100}]]
PrimePowerQ[{8, 9, 12, 1, -8, 1/49}]
PrimePowerQ[(2 + I)^3]
PrimePowerQ[5, GaussianIntegers -> True]
{t, v} = AbsoluteTiming[PrimePi[10^12]]; v
t < 1
Prime[10^9]
