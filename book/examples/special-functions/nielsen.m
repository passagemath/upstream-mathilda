# Nielsen's generalized polylogarithm S_{n,p}(z)
PolyLog[1, 2, 0.5]
N[PolyLog[1, 2, 1/2], 30]
N[Zeta[3]/8 - Log[2]^3/6, 30]
N[PolyLog[2, 2, 1], 30]
N[Zeta[4]/4, 30]
# Beyond the unit disk for integer n and p
PolyLog[1, 2, 3.]
# A non-integer order off the disk stays unevaluated
PolyLog[1.5, 2, 3.]
PolyLog[2, 2, {0.1, 0.5, -0.9}]
# LerchPhi summed over all integers k
LerchPhi[1, 2, 1/2, DoublyInfinite -> True]
