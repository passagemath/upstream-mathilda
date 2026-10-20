# One inexact argument numericalizes the exact ones
Zeta[2., Sqrt[2]]
Zeta[2.`30, Sqrt[2]]
BesselJ[Sqrt[2], 2.]
# q-Pochhammer: exact, then a non-integer n
QPochhammer[1/2, 1/2, 3]
QPochhammer[0.5, 0.5, 2.5]
