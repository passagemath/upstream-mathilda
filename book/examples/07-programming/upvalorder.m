# Ch.12 -- evaluation order: UpValues before DownValues; HoldAllComplete switches them off
f[x_] := "general"
g /: f[g[x_]] := "from g"
{f[1], f[g[1]]}
SetAttributes[hc, HoldAllComplete]
g /: hc[g[x_]] := "from g"
hc[g[1]]
SetAttributes[ho, HoldAll]
g /: ho[g[x_]] := "from g"
ho[g[1]]
Definition[g]
DownValues[f] = {HoldPattern[f[x_]] :> x^2}
{f[3], f[g[3]]}
