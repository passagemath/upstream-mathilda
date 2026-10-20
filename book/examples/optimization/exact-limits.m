# What a proof can say that a search cannot: unbounded, infeasible, and declined problems.
Minimize[x^3 - x, x]
NMinimize[x^3 - x, x]
Minimize[{x + y, x <= y^2}, {x, y}]
NMinimize[{x + y, x <= y^2}, {x, y}]
Minimize[{x + y, x^2 < -1}, {x, y}]
Check[Minimize[x^3 - x, x], "no minimum"]
Minimize[x^2 + Sin[x], x]
Minimize[2.5 x^2 - 3 x, x]
