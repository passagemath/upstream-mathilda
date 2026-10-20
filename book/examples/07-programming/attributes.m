# Ch.7 -- attributes drive evaluation: Orderless, Flat, Listable
Attributes[Plus]
SetAttributes[circ, Orderless]
circ[c, a, b]
SetAttributes[gp, Flat]
gp[a, gp[b, c], d]
SetAttributes[hl, Listable]
hl[{1, 2, 3}]
hl[{1, 2}, {3, 4}]
ClearAttributes[hl, Listable]
hl[{1, 2, 3}]
# Changing a builtin's attributes (v0.480, issue #103); state restored at the end
Attributes[Sin]
ClearAttributes[Sin, Listable]
Sin[{0, Pi/2}]
SetAttributes[Sin, Listable]
Sin[{0, Pi/2}]
Unprotect[Sin]
Sin[t0] = 3/5
Sin[t0]^2 + Cos[t0]^2
Sin[t0] =.
Protect[Sin]
Attributes[Sin]
t1 /: Sin[t1] = 3/5
{Sin[t1], Sin[t2]}
