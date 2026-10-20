# Section 5.6.3: the factoring methods, each callable on its own.
Names["FactorInteger`*"]
FactorInteger`ECM[2^128 + 1]
FactorInteger[2^128 + 1, Method -> {"ECM", "B1" -> 5*10^4}]
FactorInteger`Fermat[NextPrime[10^30] NextPrime[10^30 + 10^6]]
Quiet[FactorInteger`Fermat[(2^61 - 1) NextPrime[10^40]]]
FactorInteger`PollardPMinus1[(2^61 - 1) NextPrime[10^40]]
FactorInteger`QuadraticSieve[NextPrime[10^29] NextPrime[3 10^29]]
FactorInteger[2^256 + 1]
