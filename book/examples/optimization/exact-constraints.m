# Constrained exact optimization: Lagrange/KKT candidates, certified globally.
Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]
r = Minimize[{x^2 + y^2, (x - 2)^2 + (y - 3)^2 <= 1}, {x, y}]
Expand[(Sqrt[13] - 1)^2]
N[First[r]]
Minimize[{x^2 + y^2, x y == 1}, {x, y}]
Minimize[{Sqrt[x] + Sqrt[y], x + y == 2}, {x, y}]
Minimize[{(x^2 + 1)/(x + 1), x >= 0}, x]
Minimize[{(x^2 + 1)/(x - 1), x > 1}, x]
