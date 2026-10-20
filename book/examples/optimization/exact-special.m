# Piecewise (Abs) objectives, separable objectives, and the time budget.
Minimize[Sum[Abs[x - i^2], {i, 1, 10}], x]
Minimize[x^2 + Abs[x - 2], x]
Minimize[2 x + Abs[x], x]
st = Sum[v^4 - 16 v^2 + 5 v, {v, {a, b, c, d, e, f}}]
m = Minimize[st, {a, b, c, d, e, f}];
N[m]
NMinimize[st, {a, b, c, d, e, f}]
Minimize[{x^(2/3) + y^(2/3), x^4 + y^4 <= 1}, {x, y}, TimeConstraint -> 1]
