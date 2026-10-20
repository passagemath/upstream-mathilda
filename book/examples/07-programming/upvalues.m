# Ch.12 -- UpValues: attaching a rule to the argument's symbol (UpSet, TagSetDelayed, TagUnset)
area[sq[s_]] ^= s^2
area[sq[3]]
area[sq[a + b]]
UpValues[sq]
DownValues[area]
sq /: perimeter[sq[s_]] := 4 s
perimeter[sq[5]]
area[r_] := Pi r^2
{area[2], area[sq[2]]}
sq /: perimeter[sq[s_]] =.
perimeter[sq[5]]
UpValues[sq]
Definition[sq]
