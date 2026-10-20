# Ch.12 -- UpValues on a Protected head: integers modulo p as a new number type
mod /: mod[a_, p_] + mod[b_, p_] := mod[Mod[a + b, p], p]
mod /: mod[a_, p_] mod[b_, p_] := mod[Mod[a b, p], p]
mod[1, 5] + mod[2, 5] + mod[3, 5]
mod[1, 5] + mod[2, 5] + x
mod[2, 5] + mod[2, 5]
mod[3, 7] mod[5, 7]
mod[1, 5] + mod[2, 7]
Attributes[Plus]
DownValues[Plus]
