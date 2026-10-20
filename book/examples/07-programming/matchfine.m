# Ch.12 -- matcher fine print: sym:patt binds one atom; unbound alternatives; contiguous OrderlessPatternSequence
FullForm[Hold[n_. a:_g]]
3 g[x] /. n_. a:_g :> {n, a}
g[x] /. n_. a:_g :> {n, a}
f[1] /. f[x_] | h[y_] :> {x, y}
MatchQ[f[a, b, c], f[___, OrderlessPatternSequence[c, b], ___]]
MatchQ[f[a, b, c], f[___, OrderlessPatternSequence[a, c], ___]]
