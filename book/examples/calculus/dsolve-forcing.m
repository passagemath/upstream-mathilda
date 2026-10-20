# 5.3 DSolve -- step and impulse forcing, special-function and inert answers
DSolve[{y'[t] + y[t] == UnitStep[t - 1], y[0] == 0}, y[t], t]
DSolve[{y'[t] + 2 y[t] == DiracDelta[t - 3], y[0] == 1}, y[t], t]
DSolve[{y''[t] + y[t] == DiracDelta[t - Pi], y[0] == 0, y'[0] == 0}, y[t], t]
DSolve[x y''[x] + (1 - x) y'[x] + 2 y[x] == 0, y[x], x]
DSolve[y'[x] == f[x] + g[x] y[x], y[x], x]
DSolve[y'[x] == y[x]^2 - 2/x^2, y[x], x]
