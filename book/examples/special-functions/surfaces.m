ks = Range[0.05, 1., 0.05];
v = WhittakerM[ks, 0.7, 2.1];
PackedArrayQ[v]
# The same bits from a scalar call...
WhittakerM[ks[[7]], 0.7, 2.1] === v[[7]]
# ...compiled, with the array in the first operand...
f = Compile[{{k, _Real, 1}}, WhittakerM[k, 0.7, 2.1]];
f[ks] === v
# ...and over a visible NDArray
Normal[WhittakerM[NDArray[ks], 0.7, 2.1]] === v
CompileDiagnostics[{{k, _Real, 1}}, WhittakerM[k, 0.7, 2.1]]
