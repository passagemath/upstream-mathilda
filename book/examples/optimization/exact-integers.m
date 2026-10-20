# Exact linear and integer programs: the relaxation and the integer optimum, both proven.
Maximize[{x + y, -2 x + 2 y >= 1 && -8 x + 10 y <= 13 && x >= 0 && y >= 0}, {x, y}]
Maximize[{x + y, -2 x + 2 y >= 1 && -8 x + 10 y <= 13 && x >= 0 && y >= 0}, {x, y}, Integers]
Minimize[{x + y, x^2 + y^2 == 25}, {x, y}, Integers]
Minimize[{x + y + z, x^2 + y^2 + z^2 == 9}, {x, y, z}, Integers]
Minimize[{x + y, x^2 + y^2 <= 4 && Element[{x, y}, Integers]}, {x, y}]
Maximize[{60 a + 100 b + 120 c, 10 a + 20 b + 30 c <= 50 && 0 <= a <= 1 && 0 <= b <= 1 && 0 <= c <= 1}, {a, b, c}, Integers]
