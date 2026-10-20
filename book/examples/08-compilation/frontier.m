# Section 13.4: the frontier moves. Heads that recently gained a machine lowering.
CompileDiagnostics[{{n, _Integer}}, MoebiusMu[n] + EulerPhi[n] + DivisorSigma[1, n]]
CompileDiagnostics[{{n, _Integer}}, If[PrimeQ[n], PowerMod[2, n - 1, n], JacobiSymbol[n, 7]]]
CompileDiagnostics[{{v, _Real, 1}}, Total[Clip[Chop[v], {-1., 1.}]]]
CompileDiagnostics[{{a, _Real, 1}}, Hypergeometric2F1[a, 1.5, 2.5, 0.3]]
lerch = Compile[{{v, _Real, 1}}, LerchPhi[v, 2., 0.5]];
lerch[{0.1, 0.2, 0.3}] === LerchPhi[{0.1, 0.2, 0.3}, 2., 0.5]
