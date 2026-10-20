# Exact optimization: Minimize and Maximize return proven global optima in closed form.
Minimize[2 x^2 - 3 x + 5, x]
Maximize[x (10 - x), x]
Minimize[Expand[(x - 1) (x - 2) (x - 4) (x - 7) (x - 9) (x - 10)], x]
q = Minimize[x^4 - 3 x^2 + x, x]
RootReduce[First[q]]
N[q]
Minimize[(x y - 3)^2 + 1, {x, y}]
