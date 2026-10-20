# Tricomi's U: exact reductions
HypergeometricU[a, a + 1, z]
HypergeometricU[2, 3, z]
D[HypergeometricU[a, b, z], z]
N[HypergeometricU[1/3, 2, 5], 30]
# The Whittaker functions, real and complex
WhittakerM[2, 3, 4.5]
WhittakerW[2, 3, 4.5]
N[WhittakerW[1/5, 2 - I, 2]]
# The value at the origin follows the leading exponent
WhittakerW[k, 1/2, 0]
D[WhittakerW[k, m, z], z]
# M is a scaled Kummer 1F1, checked at a machine point
Chop[WhittakerM[k, m, z] - Exp[-z/2] z^(m + 1/2) Hypergeometric1F1[m - k + 1/2, 2 m + 1, z] /. {k -> 0.3, m -> 0.7, z -> 2.1}]
