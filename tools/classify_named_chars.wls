(* Classify every named character in src/strings/named_chars.inc by how the
   reference parser treats it in source code (not inside a string).

   Usage: wolframscript -file tools/classify_named_chars.wls src/strings/named_chars.inc > out.tsv

   One TSV line per name: name, code point, class, head, associativity, precedence.
     LETTER    forms identifiers (a<c>b is one symbol); the head column holds
               the symbol a one-character name aliases (\[Pi] alone is Pi)
     SPACE     whitespace; NEWLINE: a line break
     INFIX / PREFIX / POSTFIX  an operator with the given head
     OTHER     anything else (rejected, or brackets/structural)
   tools/gen_named_chars.py merges the result into the .inc table. *)
file = Last[$ScriptCommandLine];
lines = Select[ReadList[file, String], StringStartsQ[#, "    {\""] &];
entries = StringCases[#, "{\"" ~~ n : WordCharacter .. ~~ "\", 0x" ~~ h : HexadecimalCharacter .. ~~ "u}" :> {n, FromDigits[h, 16]}][[1]] & /@ lines;
hp[s_] := Quiet@Check[ToExpression[s, InputForm, Hold], $Failed];
symName[Hold[x_Symbol]] := SymbolName[Unevaluated[x]];
symName[_] := None;
heads[Hold[h_[___]]] := ToString[Unevaluated[h]];
classify[{name_, cp_}] := Module[{c = FromCharacterCode[cp], alone, inf, pre, post, tri, h, assoc = "", prec = ""},
  alone = hp[c]; inf = hp["a" <> c <> "b"]; pre = hp[c <> "b"]; post = hp["a" <> c];
  Which[
    symName[inf] === "a" <> c <> "b",
      {"LETTER", If[symName[alone] =!= c && symName[alone] =!= None, symName[alone], ""], "", ""},
    inf === Hold[a, b], {"NEWLINE", "", "", ""},
    MatchQ[alone, Hold[] | Hold[Null]] && inf === Hold[a b], {"SPACE", "", "", ""},
    MatchQ[alone, Hold[_Symbol]] && symName[alone] =!= c, {"CONSTANT", symName[alone], "", ""},
    MatchQ[inf, Hold[_[a, b]]],
      h = heads[inf]; tri = hp["a" <> c <> "b" <> c <> "d"];
      assoc = Which[tri === ToExpression["Hold[" <> h <> "[a,b,d]]"], "FLAT",
                    tri === ToExpression["Hold[" <> h <> "[" <> h <> "[a,b],d]]"], "LEFT",
                    tri === ToExpression["Hold[" <> h <> "[a," <> h <> "[b,d]]]"], "RIGHT",
                    True, "OTHER"];
      prec = Quiet@Check[Precedence[ToExpression[h]], ""];
      {"INFIX", h, assoc, ToString[prec]},
    MatchQ[pre, Hold[_[b]]], h = heads[pre]; {"PREFIX", h, "", ToString[Quiet@Check[Precedence[ToExpression[h]], ""]]},
    MatchQ[post, Hold[_[a]]], h = heads[post]; {"POSTFIX", h, "", ToString[Quiet@Check[Precedence[ToExpression[h]], ""]]},
    True, {"OTHER", "", "", ""}]];
Do[WriteString["stdout", StringRiffle[Join[{e[[1]], IntegerString[e[[2]], 16, 4]}, classify[e]], "\t"], "\n"], {e, entries}];
