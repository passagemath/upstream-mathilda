# Classical orthogonal polynomials: exact integer-coefficient polynomials
ChebyshevT[5, x]
ChebyshevU[4, x]
HermiteH[5, x]
LaguerreL[4, x]
# The generalized Laguerre polynomial carries its parameter symbolically
LaguerreL[2, a, x]
# Chebyshev's defining identity T_n(cos t) = cos(n t)
Simplify[ChebyshevT[3, Cos[t]]]
# Orthogonality against the Chebyshev and Hermite weights
Integrate[ChebyshevT[3, x] ChebyshevT[5, x]/Sqrt[1 - x^2], {x, -1, 1}]
Integrate[HermiteH[3, x]^2 Exp[-x^2], {x, -Infinity, Infinity}]
