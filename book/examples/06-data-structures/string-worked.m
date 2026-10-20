# 11.4 Worked example: from raw text to structured data
config = "host = example.org\nport=8080\n# a comment\ndebug = true";
lines = Select[StringSplit[config, "\n"], ! StringStartsQ[#, "#"] &]
kv = Association[StringCases[lines, k : WordCharacter .. ~~ Whitespace ... ~~ "=" ~~ Whitespace ... ~~ v__ :> k -> v]]
kv["port"]
ToExpression[kv["port"]] + 1
text = "The cat sat on the mat. The mat was flat, and the cat was fat!";
words = ToLowerCase[StringCases[text, WordCharacter ..]]
Counts[words]
Select[Counts[words], # > 1 &]
