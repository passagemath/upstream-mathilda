# Section 5.6.4: Blake's rational base descent.
P = NextPrime[Floor[3 (22/7)^150]];
m = NextPrime[7 10^29];
P - Floor[3 (22/7)^150]
IntegerLength[P m]
FactorInteger`BlakeRationalBaseDescentClassic[P m, "Base" -> 22/7]
FactorInteger`BlakeRationalBaseDescent[P m]
TimeConstrained[FactorInteger[P m], 2]
