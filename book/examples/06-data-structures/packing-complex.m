# 11.5 Complex buffers and n-ary special functions over packed arrays
z = RandomComplex[{-1 - I, 1 + I}, 1000];
{NDArrayQ[z], DataType[z]}
c = Table[Complex[k, 1.], {k, 1., 10.}];
{NDArrayQ[c], DataType[c]}
NDArrayQ[Exp[I c]]
NDArrayQ[Table[If[k == 3., k, k + I], {k, 1., 10.}]]
NDArrayQ[Table[k + I, {k, 10}]]
x = Range[0., 0.9, 0.1];
NDArrayQ[Hypergeometric2F1[1, 2, 3, x]]
NDArrayQ[LerchPhi[x, 2, 1]]
NDArrayQ[Hypergeometric2F1[x, 2, 3, 0.5]]
Hypergeometric2F1[x, 2, 3, 0.5] == Table[Hypergeometric2F1[a, 2, 3, 0.5], {a, x}]
