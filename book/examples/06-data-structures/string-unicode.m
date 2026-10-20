# 11.4 Strings are sequences of characters, not bytes (outputs kept ASCII for the listings)
g = "\[Alpha]\[Beta]\[Gamma]";
StringLength[g]
ToCharacterCode[g]
ToCharacterCode[StringTake[g, {2}]]
ToCharacterCode[ToUpperCase[g]]
"\:03b1" === "\[Alpha]" === FromCharacterCode[945]
StringLength["caf\[EAcute]"]
StringMatchQ[g, LetterCharacter ..]
