# Reduce without an explicit variable list
Reduce[x^2 == 1]
Reduce[x^2 < 1, Reals]
Reduce[x + y + z == 1 && x >= 0 && y >= 0 && z >= 0 && x^2 + y^2 + z^2 > 4, Reals]
