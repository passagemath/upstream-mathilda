# InputForm prints machine reals with the shortest digits that read back exactly.
0.1 + 0.2
InputForm[0.1 + 0.2]
InputForm[{0.1, 1/3., 1.5*^-5}]
ToExpression[ToString[0.1 + 0.2, InputForm]] === 0.1 + 0.2
