# 5.3 Series under Assumptions -- branch selection
Series[Abs[x - a], {x, 0, 2}]
Series[Abs[x - a], {x, 0, 2}, Assumptions -> a > 0]
Series[Sqrt[1 - Cos[x]], {x, 0, 3}, Assumptions -> x < 0]
Series[Sqrt[a - x], {x, 0, 1}, Assumptions -> a < 0]
Series[ArcTan[a x], {x, Infinity, 1}, Assumptions -> a < 0]
