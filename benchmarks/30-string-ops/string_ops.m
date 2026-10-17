(* Experiment 30 -- String operations and regular expressions.
   MEASURES src/strings/ (17 files) and src/strings/regex/ (11 files, PCRE2-backed
   when USE_REGEX is compiled in).

   NOT COVERED BY THE EXISTING CORPUS AT ALL.  Every one of the twenty existing
   experiments is numeric; the string subsystem has never been timed against
   anything.  Python's `re` is a mature C engine and its str methods are heavily
   optimised, so this is a real execution baseline rather than a pure-Python one.

   Checks are LENGTHS and COUNTS on small deterministic inputs. *)

Get["../harness.m"];
Get["../data.m"];

require[{"StringJoin", "StringSplit", "StringReplace", "StringLength",
         "StringTake", "StringPosition", "StringCases", "StringCount",
         "StringReverse", "StringRiffle", "RegularExpression", "StringMatchQ",
         "ToUpperCase", "Characters"}];

(* A deterministic ~200000-character subject built by repetition. *)
unit = "the quick brown fox jumps over the lazy dog ";
big = StringJoin[Table[unit, {5000}]];

bench["StringLength of 200k chars", StringLength[big];];
check["StringLength of 200k chars", StringLength[unit]];

bench["StringReplace literal, 200k chars",
  StringReplace[big, "quick" -> "slow"];];
check["StringReplace literal, 200k chars",
  StringLength[StringReplace[unit, "quick" -> "slow"]]];

bench["StringSplit on space, 200k chars", StringSplit[big];];
check["StringSplit on space, 200k chars", Length[StringSplit[unit]]];

bench["StringCount substring, 200k chars", StringCount[big, "the"];];
check["StringCount substring, 200k chars", StringCount[unit, "the"]];

bench["StringReverse 200k chars", StringReverse[big];];
check["StringReverse 200k chars", StringLength[StringReverse[unit]]];

benchIf["StringCases regex, 200k chars", "RegularExpression",
  StringCases[big, RegularExpression["[aeiou]+"]], 1];
checkIf["StringCases regex, 200k chars", "RegularExpression",
  Length[StringCases[unit, RegularExpression["[aeiou]+"]]]];

benchIf["StringReplace regex, 200k chars", "RegularExpression",
  StringReplace[big, RegularExpression["[aeiou]"] -> "*"], 1];
checkIf["StringReplace regex, 200k chars", "RegularExpression",
  StringLength[StringReplace[unit, RegularExpression["[aeiou]"] -> "*"]]];

bench["Characters of 200k chars", Characters[big], 1];
check["Characters of 200k chars", Length[Characters[unit]]];

(* ---- String-subsystem review (v0.455+): the cases the review targeted. ----
   words: the 40000 words of `big`; bigU: 200000 characters of 2-byte UTF-8. *)
words = StringSplit[big];
bigU = StringRepeat["h\:00e9llo w\:00f6rld ", 16667];

(* One pattern over many short strings: pays a pattern compile per call unless
   compiled patterns are cached. *)
bench["StringMatchQ over 40k words",
  Map[StringMatchQ[#, __ ~~ "o" ~~ __] &, words];, 1];
check["StringMatchQ over 40k words",
  Count[Map[StringMatchQ[#, __ ~~ "o" ~~ __] &, StringSplit[unit]], True]];

(* Several rules, one of which never matches: re-searching every rule at every
   step made this quadratic. *)
bench["StringReplace 3 rules, 200k chars",
  StringReplace[big, {"quick" -> "slow", "fox" -> "cat", "zzz" -> "y"}];, 1];
check["StringReplace 3 rules, 200k chars",
  StringLength[StringReplace[unit, {"quick" -> "slow", "fox" -> "cat", "zzz" -> "y"}]]];

(* A named capture bound into an evaluated right-hand side (issue #98). *)
bench["StringCases x:p :> f[x], 200k chars",
  StringCases[big, x:LetterCharacter.. ~~ " " :> StringLength[x]];, 1];
check["StringCases x:p :> f[x], 200k chars",
  Total[StringCases[unit, x:LetterCharacter.. ~~ " " :> StringLength[x]]]];

bench["StringPosition overlapping, 200k chars", StringPosition[big, "o"];, 1];
check["StringPosition overlapping, 200k chars", Length[StringPosition[unit, "o"]]];

bench["StringLength UTF-8, 200k chars", StringLength[bigU];];
check["StringLength UTF-8, 200k chars", StringLength["h\:00e9llo w\:00f6rld "]];

bench["Characters UTF-8, 200k chars", Characters[bigU];, 1];
check["Characters UTF-8, 200k chars", Length[Characters["h\:00e9llo w\:00f6rld "]]];

bench["StringTake over 40k words", StringTake[words, UpTo[3]];, 1];
check["StringTake over 40k words", StringLength[StringJoin[StringTake[StringSplit[unit], UpTo[3]]]]];

bench["ToUpperCase 200k chars", ToUpperCase[big];];
check["ToUpperCase 200k chars", StringLength[ToUpperCase[unit]]];

bench["StringJoin of 40k words", StringJoin[words];];
check["StringJoin of 40k words", StringLength[StringJoin[StringSplit[unit]]]];

bench["StringRiffle of 40k words", StringRiffle[words, ","];, 1];
check["StringRiffle of 40k words", StringLength[StringRiffle[StringSplit[unit], ","]]];
