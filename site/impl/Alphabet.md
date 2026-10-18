---
source: src/strings/charcode.c
---
**Algorithm.** `builtin_alphabet` accepts only the zero-argument call and returns `char_list('a', 'z')` — the same helper as `CharacterRange`, which encodes each code point `0x61`–`0x7A` as a one-character `String`. Any argument (a language name, a type spec) makes it return `NULL`, so the call stays unevaluated.

**Data structures.** A fresh 26-element `List` of one-character strings on every call.

**Complexity / limits.** Constant time. Mathematica's language-specific forms (`Alphabet["Greek"]`, `Alphabet[lang, type]`) and the language-data they rely on are not implemented.
