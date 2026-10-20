# Ch.12 -- Flat heads: a rule may match part of a longer sum, product, or ** chain
1 + x + y /. x + y -> z
MatchQ[1 + x + y, x + y]
a b c /. a c -> s
x ** y ** z
x ** y ** z /. y ** z -> w
x ** y ** z /. x ** z -> w
Attributes[NonCommutativeMultiply]
FullForm[Hold[a . b ** c^2]]
op /: op[i_] ** op[j_] := op[i + j]
op[1] ** op[2] ** op[3]
a ** op[1] ** op[2] ** b
op[1] ** a ** op[2]
