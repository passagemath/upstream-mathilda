# The simplest G-functions are elementary
MeijerG[{{}, {}}, {{0}, {}}, z]
MeijerG[{{1/2}, {}}, {{0}, {}}, z]
MeijerG[{{1, 1}, {}}, {{1}, {0}}, z]
# Bessel functions are G-functions of x^2/4
MeijerG[{{}, {}}, {{v}, {-v}}, z]
MeijerG[{{}, {}}, {{1/3, 0}, {}}, z]
# Three b's: no single classical form, so G stays
MeijerG[{{}, {}}, {{1/3, 0, 1/2}, {}}, z]
# FunctionExpand writes out its Slater expansion
FunctionExpand[MeijerG[{{}, {}}, {{1/3, 0, 1/2}, {}}, z]]
# Numerics: machine, arbitrary precision, a list
MeijerG[{{}, {}}, {{1/3, 0, 1/2}, {}}, {0.1, 0.7, 2.}]
N[MeijerG[{{}, {}}, {{1/3, 0, 1/2}, {}}, 7/10], 25]
# The log case (b's an integer apart) has no Slater expansion
MeijerG[{{}, {2}}, {{1/2, 3/2}, {}}, z]
N[MeijerG[{{}, {2}}, {{1/2, 3/2}, {}}, 3], 30]
