# A machine argument runs the stable three-term recurrence
HermiteH[100, 5.]
# Substituting into the expanded polynomial cancels catastrophically
HermiteH[100, x] /. x -> 5.
# The exact value, to twenty digits
N[HermiteH[100, 5], 20]
LaguerreL[100, 5.]
LaguerreL[100, x] /. x -> 5.
# Non-integer and complex order
ChebyshevT[1/2, 0.3]
N[ChebyshevT[1/7, 1/3], 30]
ChebyshevT[2 + I, 0.5]
# Derivatives close within the family
D[ChebyshevT[n, x], x]
D[HermiteH[n, x], x]
D[LaguerreL[n, x], x]
# Associated Legendre functions, exact and over a list of machine reals
LegendreP[3, 1, x]
LegendreP[30, 4, {0.1, 0.3, 0.5}]
