# 5.3 Definite integration -- piecewise integrands and symbolic endpoints
Integrate[Abs[x - a], {x, 0, 2 a}, Assumptions -> a > 0]
Integrate[Min[x, a], {x, 0, 2 a}, Assumptions -> a > 0]
Integrate[Max[Sin[x], Cos[x]], {x, 0, Pi}]
Integrate[Piecewise[{{x^2, x < 1}}, 2 - x], {x, 0, 2}]
Integrate[Sqrt[1 - Cos[2 x]], {x, 0, 2 Pi}]
Integrate[1/Sqrt[a^2 - x^2], {x, 0, a}]
Integrate[1/Sqrt[a^2 - x^2], {x, 0, a}, Assumptions -> a > 0]
Integrate[x/Sqrt[a^2 - x^2], {x, 0, a}, Assumptions -> a > 0]
