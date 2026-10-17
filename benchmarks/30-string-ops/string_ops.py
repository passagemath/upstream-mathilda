#!/usr/bin/env python3
"""Experiment 30 -- String operations and regular expressions (Python column).

Same eight kernels as ``string_ops.m``, same order and sizes.

NOT COVERED BY THE EXISTING CORPUS AT ALL -- every one of the twenty existing
experiments is numeric, and Mathilda's string subsystem (28 files) has never been
timed against anything.

Python's ``re`` is a mature C engine and ``str`` methods are heavily optimised, so
unlike the sympy and networkx columns this IS a real compiled baseline.  Mathilda
routes through PCRE2 when USE_REGEX is compiled in, so the regex rows compare two
C engines.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import re

from harness import bench, check, require

require(["re:sub", "re:findall"])

unit = "the quick brown fox jumps over the lazy dog "
big = unit * 5000

VOWELS = re.compile("[aeiou]+")
VOWEL = re.compile("[aeiou]")

bench("StringLength of 200k chars", lambda: len(big))
check("StringLength of 200k chars", len(unit))

bench("StringReplace literal, 200k chars",
      lambda: big.replace("quick", "slow"))
check("StringReplace literal, 200k chars", len(unit.replace("quick", "slow")))

bench("StringSplit on space, 200k chars", lambda: big.split())
check("StringSplit on space, 200k chars", len(unit.split()))

bench("StringCount substring, 200k chars", lambda: big.count("the"))
check("StringCount substring, 200k chars", unit.count("the"))

bench("StringReverse 200k chars", lambda: big[::-1])
check("StringReverse 200k chars", len(unit[::-1]))

bench("StringCases regex, 200k chars", lambda: VOWELS.findall(big), reps=1)
check("StringCases regex, 200k chars", len(VOWELS.findall(unit)))

bench("StringReplace regex, 200k chars",
      lambda: VOWEL.sub("*", big), reps=1)
check("StringReplace regex, 200k chars", len(VOWEL.sub("*", unit)))

bench("Characters of 200k chars", lambda: list(big), reps=1)
check("Characters of 200k chars", len(list(unit)))

# ---- String-subsystem review (v0.455+): the cases the review targeted. ----
words = big.split()
bigU = "héllo wörld " * 16667
O_PAT = re.compile(r".+o.+", re.S)

bench("StringMatchQ over 40k words",
      lambda: [O_PAT.fullmatch(w) is not None for w in words], reps=1)
check("StringMatchQ over 40k words",
      sum(O_PAT.fullmatch(w) is not None for w in unit.split()))

THREE = re.compile("quick|fox|zzz")
REPL = {"quick": "slow", "fox": "cat", "zzz": "y"}
bench("StringReplace 3 rules, 200k chars",
      lambda: THREE.sub(lambda m: REPL[m.group(0)], big), reps=1)
check("StringReplace 3 rules, 200k chars",
      len(THREE.sub(lambda m: REPL[m.group(0)], unit)))

WORD_SP = re.compile(r"([^\W\d_]+) ")
bench("StringCases x:p :> f[x], 200k chars",
      lambda: [len(m.group(1)) for m in WORD_SP.finditer(big)], reps=1)
check("StringCases x:p :> f[x], 200k chars",
      sum(len(m.group(1)) for m in WORD_SP.finditer(unit)))

bench("StringPosition overlapping, 200k chars",
      lambda: [(m.start() + 1, m.end()) for m in re.finditer("o", big)], reps=1)
check("StringPosition overlapping, 200k chars", unit.count("o"))

bench("StringLength UTF-8, 200k chars", lambda: len(bigU))
check("StringLength UTF-8, 200k chars", len("héllo wörld "))

bench("Characters UTF-8, 200k chars", lambda: list(bigU), reps=1)
check("Characters UTF-8, 200k chars", len(list("héllo wörld ")))

bench("StringTake over 40k words", lambda: [w[:3] for w in words], reps=1)
check("StringTake over 40k words", len("".join(w[:3] for w in unit.split())))

bench("ToUpperCase 200k chars", lambda: big.upper())
check("ToUpperCase 200k chars", len(unit.upper()))

bench("StringJoin of 40k words", lambda: "".join(words))
check("StringJoin of 40k words", len("".join(unit.split())))

bench("StringRiffle of 40k words", lambda: ",".join(words), reps=1)
check("StringRiffle of 40k words", len(",".join(unit.split())))
