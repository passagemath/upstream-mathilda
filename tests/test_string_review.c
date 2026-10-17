/*
 * test_string_review.c - the String* review (v0.456): every expected value
 * below is Mathematica 14's answer to the same input (generated from the
 * differential corpus with wolframscript; never written by hand).
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include <stdio.h>

static void test_unicode_code_points(void) {
    /* Unicode code-point semantics: lengths, characters, reversal, padding,
     * repetition, positions, Unicode classes and case folding in patterns. */
    assert_eval_eq("StringLength[\"h\\:00e9llo\"]",
                   "5", 0);
    assert_eval_eq("StringLength[\"\\:65e5\\:672c\\:8a9e\"]",
                   "3", 0);
    assert_eval_eq("StringLength[\"a\\|01f600b\"]",
                   "3", 0);
    assert_eval_eq("StringLength[\"\"]",
                   "0", 0);
    assert_eval_eq("Characters[\"h\\:00e9llo\"]",
                   "{\"h\", \"é\", \"l\", \"l\", \"o\"}", 0);
    assert_eval_eq("Characters[\"\\:65e5\\:672c\"]",
                   "{\"日\", \"本\"}", 0);
    assert_eval_eq("Characters[\"\"]",
                   "{}", 0);
    assert_eval_eq("StringReverse[\"h\\:00e9llo\"]",
                   "\"olléh\"", 0);
    assert_eval_eq("StringReverse[\"a\\|01f600b\"]",
                   "\"b😀a\"", 0);
    assert_eval_eq("StringPadLeft[\"\\:00e9\", 3]",
                   "\"  é\"", 0);
    assert_eval_eq("StringPadRight[\"\\:00e9\", 3, \"\\:00e8\"]",
                   "\"éèè\"", 0);
    assert_eval_eq("StringPadLeft[\"abc\", 2]",
                   "\"bc\"", 0);
    assert_eval_eq("StringPadLeft[{\"a\", \"bcd\"}]",
                   "{\"  a\", \"bcd\"}", 0);
    assert_eval_eq("StringPadRight[\"ab\", 5, \"xy\"]",
                   "\"abxyx\"", 0);
    assert_eval_eq("StringPadLeft[\"ab\", 5, \"xy\"]",
                   "\"yxyab\"", 0);
    assert_eval_eq("StringRepeat[\"\\:00e9a\", 3]",
                   "\"éaéaéa\"", 0);
    assert_eval_eq("StringRepeat[\"\\:00e9a\", 3, 5]",
                   "\"éaéaé\"", 0);
    assert_eval_eq("StringRepeat[\"ab\", 0]",
                   "\"\"", 0);
    assert_eval_eq("StringPosition[\"h\\:00e9llo\", \"l\"]",
                   "{{3, 3}, {4, 4}}", 0);
    assert_eval_eq("StringPosition[\"\\:00e9\\:00e9\", \"\\:00e9\"]",
                   "{{1, 1}, {2, 2}}", 0);
    assert_eval_eq("StringCases[\"\\:00fc1\\:03b1\", LetterCharacter]",
                   "{\"ü\", \"α\"}", 0);
    assert_eval_eq("StringCases[\"\\:00dcber\", \"\\:00fc\", IgnoreCase -> True]",
                   "{\"Ü\"}", 0);
    assert_eval_eq("StringMatchQ[\"\\:00e9\", _]",
                   "True", 0);
    assert_eval_eq("StringSplit[\"h\\:00e9\", \"\"]",
                   "{\"h\", \"é\"}", 0);
    assert_eval_eq("StringCount[\"\\:00e9\\:00e9\\:00e9\", _]",
                   "3", 0);
    assert_eval_eq("StringReplace[\"h\\:00e9llo\", \"\\:00e9\" -> \"e\"]",
                   "\"hello\"", 0);
    assert_eval_eq("StringJoin[\"\\:00e9\", \"a\"]",
                   "\"éa\"", 0);
    assert_eval_eq("StringRiffle[{\"\\:00e9\", \"b\"}, \"-\"]",
                   "\"é-b\"", 0);
    assert_eval_eq("StringContainsQ[\"\\:00c9t\\:00e9\", \"\\:00e9t\\:00e9\", IgnoreCase -> True]",
                   "True", 0);
    assert_eval_eq("StringCases[\"a b\\:00a0c\", WhitespaceCharacter]",
                   "{\" \", \" \"}", 0);
    assert_eval_eq("StringCases[\"x\\:0662y\", DigitCharacter]",
                   "{\"٢\"}", 0);
    assert_eval_eq("StringTrim[\"\\:00a0a \"]",
                   "\"a\"", 0);
}

static void test_pad_repeat_join_riffle(void) {
    /* StringPadLeft/Right cyclic alignment, StringRepeat, StringJoin,
     * StringRiffle, StringExtract, and the Listable heads over lists. */
    assert_eval_eq("StringPadLeft[\"ab\", 6, \"xyz\"]",
                   "\"xyzxab\"", 0);
    assert_eval_eq("StringPadRight[\"ab\", 6, \"xyz\"]",
                   "\"abzxyz\"", 0);
    assert_eval_eq("StringPadLeft[\"ab\", 7, \"xyz\"]",
                   "\"zxyzxab\"", 0);
    assert_eval_eq("StringPadRight[\"ab\", 7, \"xyz\"]",
                   "\"abzxyzx\"", 0);
    assert_eval_eq("StringPadLeft[\"ab\", 5, \"xy\"]",
                   "\"yxyab\"", 0);
    assert_eval_eq("StringPadRight[\"ab\", 5, \"xy\"]",
                   "\"abxyx\"", 0);
    assert_eval_eq("StringPadLeft[{\"a\", \"bcd\"}, 4, \"xy\"]",
                   "{\"xyxa\", \"xbcd\"}", 0);
    assert_eval_eq("StringPadRight[{\"a\", \"bcd\"}, 4, \"xy\"]",
                   "{\"ayxy\", \"bcdy\"}", 0);
    assert_eval_eq("StringPadLeft[\"ab\", 5, \"\"]",
                   "StringPadLeft[\"ab\", 5, \"\"]", 0);
    assert_eval_eq("StringPadLeft[\"abcdef\", 3]",
                   "\"def\"", 0);
    assert_eval_eq("StringPadRight[\"abcdef\", 3]",
                   "\"abc\"", 0);
    assert_eval_eq("StringPadLeft[\"ab\", 0]",
                   "\"\"", 0);
    assert_eval_eq("StringPadLeft[\"ab\", -1]",
                   "StringPadLeft[\"ab\", -1]", 0);
    assert_eval_eq("StringPadLeft[{\"a\", \"bc\"}, 3, {\"<\", \">\"}]",
                   "StringPadLeft[{\"a\", \"bc\"}, 3, {\"<\", \">\"}]", 0);
    assert_eval_eq("StringPadLeft[{}, 3]",
                   "\"\"", 0);
    assert_eval_eq("StringPadLeft[{\"a\", 1}, 3]",
                   "StringPadLeft[{\"a\", 1}, 3]", 0);
    assert_eval_eq("StringPadLeft[\"abc\"]",
                   "StringPadLeft[\"abc\"]", 0);
    assert_eval_eq("StringPadRight[{\"a\", \"bcd\", \"\"}]",
                   "{\"a  \", \"bcd\", \"   \"}", 0);
    assert_eval_eq("StringPadLeft[\"\\:00e9\", 4, \"\\:00e8x\"]",
                   "\"èxèé\"", 0);
    assert_eval_eq("StringRepeat[\"ab\", 2.5]",
                   "StringRepeat[\"ab\", 2.5]", 0);
    assert_eval_eq("StringRepeat[\"ab\", -1]",
                   "StringRepeat[\"ab\", -1]", 0);
    assert_eval_eq("StringRepeat[\"\", 3]",
                   "\"\"", 0);
    assert_eval_eq("StringRepeat[\"abc\", 2, 4]",
                   "\"abca\"", 0);
    assert_eval_eq("StringRepeat[\"abc\", 10, 2]",
                   "\"ab\"", 0);
    assert_eval_eq("StringRepeat[\"a\\:00e9\", 2, 3]",
                   "\"aéa\"", 0);
    assert_eval_eq("StringRepeat[{\"a\"}, 2]",
                   "StringRepeat[{\"a\"}, 2]", 0);
    assert_eval_eq("StringJoin[\"a\", {\"b\", {\"c\"}}, \"d\"]",
                   "\"abcd\"", 0);
    assert_eval_eq("StringJoin[{}]",
                   "\"\"", 0);
    assert_eval_eq("StringJoin[]",
                   "\"\"", 0);
    assert_eval_eq("StringJoin[\"a\", 1]",
                   "StringJoin[\"a\", 1]", 0);
    assert_eval_eq("StringJoin[{\"a\", \"b\"}, {\"c\"}]",
                   "\"abc\"", 0);
    assert_eval_eq("StringRiffle[{\"a\", \"b\", \"c\"}]",
                   "\"a b c\"", 0);
    assert_eval_eq("StringRiffle[{\"a\", \"b\", \"c\"}, \", \"]",
                   "\"a, b, c\"", 0);
    assert_eval_eq("StringRiffle[{\"a\", \"b\", \"c\"}, {\"(\", \", \", \")\"}]",
                   "\"(a, b, c)\"", 0);
    assert_eval_eq("StringRiffle[{{\"a\", \"b\"}, {\"c\", \"d\"}}]",
                   "\"a b\\nc d\"", 0);
    assert_eval_eq("StringRiffle[{{\"a\", \"b\"}, {\"c\", \"d\"}}, \";\", \",\"]",
                   "\"a,b;c,d\"", 0);
    assert_eval_eq("StringRiffle[{1, 2, 3}, \"-\"]",
                   "\"1-2-3\"", 0);
    assert_eval_eq("StringRiffle[{x, \"b\"}]",
                   "\"x b\"", 0);
    assert_eval_eq("StringRiffle[{}, \"-\"]",
                   "\"\"", 0);
    assert_eval_eq("StringReverse[{\"ab\", \"cd\"}]",
                   "{\"ba\", \"dc\"}", 0);
    assert_eval_eq("StringLength[{\"ab\", \"\\:00e9\"}]",
                   "{2, 1}", 0);
    assert_eval_eq("Characters[{\"ab\", \"c\"}]",
                   "{{\"a\", \"b\"}, {\"c\"}}", 0);
    assert_eval_eq("StringExtract[\"a b c\", 2]",
                   "\"b\"", 0);
    assert_eval_eq("StringExtract[\"a b c\", -1]",
                   "\"c\"", 0);
    assert_eval_eq("StringExtract[\"a:b c:d\", \" \" -> 2]",
                   "\"c:d\"", 0);
    assert_eval_eq("StringExtract[\"a:b c:d\", {\" \", \":\"} -> {2, 1}]",
                   "{\"b\", \"a\"}", 0);
    assert_eval_eq("StringExtract[\"a\\:00e9 b\", 1]",
                   "\"aé\"", 0);
    assert_eval_eq("StringExtract[\"a b\", 5]",
                   "Missing[\"PartAbsent\", 5]", 0);
}

static void test_new_builtins_and_forms(void) {
    /* StringDelete, StringTakeDrop, StringRotateLeft/Right, StringReplaceList,
     * case mapping, character codes, the *Q predicates, CharacterRange,
     * Alphabet, operator forms and occurrence limits. */
    assert_eval_eq("StringDelete[\"abcabc\", \"b\"]",
                   "\"acac\"", 0);
    assert_eval_eq("StringDelete[\"a1b22c\", DigitCharacter..]",
                   "\"abc\"", 0);
    assert_eval_eq("StringDelete[{\"ab\", \"cb\"}, \"b\"]",
                   "{\"a\", \"c\"}", 0);
    assert_eval_eq("StringDelete[\"ABab\", \"a\", IgnoreCase -> True]",
                   "\"Bb\"", 0);
    assert_eval_eq("StringDelete[\"b\"][\"abc\"]",
                   "\"ac\"", 0);
    assert_eval_eq("StringTakeDrop[\"abcdef\", 2]",
                   "{\"ab\", \"cdef\"}", 0);
    assert_eval_eq("StringTakeDrop[\"abcdef\", -2]",
                   "{\"ef\", \"abcd\"}", 0);
    assert_eval_eq("StringTakeDrop[\"abcdef\", {2, 4}]",
                   "{\"bcd\", \"aef\"}", 0);
    assert_eval_eq("StringTakeDrop[{\"abc\", \"de\"}, 1]",
                   "{{\"a\", \"bc\"}, {\"d\", \"e\"}}", 0);
    assert_eval_eq("StringRotateLeft[\"abcdef\", 2]",
                   "\"cdefab\"", 0);
    assert_eval_eq("StringRotateLeft[\"abcdef\"]",
                   "\"bcdefa\"", 0);
    assert_eval_eq("StringRotateLeft[\"abcdef\", -1]",
                   "\"fabcde\"", 0);
    assert_eval_eq("StringRotateLeft[\"abc\", 7]",
                   "\"bca\"", 0);
    assert_eval_eq("StringRotateRight[\"abcdef\", 2]",
                   "\"efabcd\"", 0);
    assert_eval_eq("StringRotateRight[\"abcdef\"]",
                   "\"fabcde\"", 0);
    assert_eval_eq("StringRotateLeft[\"\", 3]",
                   "\"\"", 0);
    assert_eval_eq("StringRotateLeft[\"h\\:00e9llo\", 1]",
                   "\"élloh\"", 0);
    assert_eval_eq("StringReplaceList[\"aaa\", \"a\" -> \"b\"]",
                   "{\"baa\", \"aba\", \"aab\"}", 0);
    assert_eval_eq("StringReplaceList[\"abc\", {\"a\" -> \"x\", \"b\" -> \"y\"}]",
                   "{\"xbc\", \"ayc\"}", 0);
    assert_eval_eq("StringReplaceList[\"aaa\", \"aa\" -> \"b\"]",
                   "{\"ba\", \"ab\"}", 0);
    assert_eval_eq("StringReplaceList[\"abc\", \"a\" -> \"x\", 0]",
                   "{}", 0);
    assert_eval_eq("StringReplaceList[\"aaaa\", \"a\" -> \"x\", 2]",
                   "{\"xaaa\", \"axaa\"}", 0);
    assert_eval_eq("StringReplaceList[{\"ab\", \"a\"}, \"a\" -> \"x\"]",
                   "{{\"xb\"}, {\"x\"}}", 0);
    assert_eval_eq("ToUpperCase[\"Hello W\\:00f6rld\"]",
                   "\"HELLO WÖRLD\"", 0);
    assert_eval_eq("ToLowerCase[\"Hello W\\:00d6RLD\"]",
                   "\"hello wörld\"", 0);
    assert_eval_eq("ToUpperCase[{\"a\", \"b\"}]",
                   "{\"A\", \"B\"}", 0);
    assert_eval_eq("ToUpperCase[\"\\:00df\"]",
                   "\"SS\"", 0);
    assert_eval_eq("ToUpperCase[x]",
                   "ToUpperCase[x]", 0);
    assert_eval_eq("ToLowerCase[\"\\[CapitalSigma]\"]",
                   "\"σ\"", 0);
    assert_eval_eq("Capitalize[\"hello world\"]",
                   "\"Hello world\"", 0);
    assert_eval_eq("Capitalize[\"\\:00e9t\\:00e9\"]",
                   "\"Été\"", 0);
    assert_eval_eq("Capitalize[\"\"]",
                   "\"\"", 0);
    assert_eval_eq("Capitalize[{\"ab\", \"cd\"}]",
                   "{\"Ab\", \"Cd\"}", 0);
    assert_eval_eq("Decapitalize[\"Hello World\"]",
                   "\"hello World\"", 0);
    assert_eval_eq("Decapitalize[\"HELLO\"]",
                   "\"hELLO\"", 0);
    assert_eval_eq("ToCharacterCode[\"abc\"]",
                   "{97, 98, 99}", 0);
    assert_eval_eq("ToCharacterCode[\"h\\:00e9\"]",
                   "{104, 233}", 0);
    assert_eval_eq("ToCharacterCode[{\"ab\", \"c\"}]",
                   "{{97, 98}, {99}}", 0);
    assert_eval_eq("ToCharacterCode[\"\"]",
                   "{}", 0);
    assert_eval_eq("ToCharacterCode[\"a\\|01f600\"]",
                   "{97, 128512}", 0);
    assert_eval_eq("FromCharacterCode[97]",
                   "\"a\"", 0);
    assert_eval_eq("FromCharacterCode[{97, 98}]",
                   "\"ab\"", 0);
    assert_eval_eq("FromCharacterCode[{{97, 98}, {99}}]",
                   "{\"ab\", \"c\"}", 0);
    assert_eval_eq("FromCharacterCode[233]",
                   "\"é\"", 0);
    assert_eval_eq("FromCharacterCode[128512]",
                   "\"😀\"", 0);
    assert_eval_eq("FromCharacterCode[{}]",
                   "\"\"", 0);
    assert_eval_eq("FromCharacterCode[-1]",
                   "FromCharacterCode[-1]", 0);
    assert_eval_eq("FromCharacterCode[1114112]",
                   "FromCharacterCode[1114112]", 0);
    assert_eval_eq("LetterQ[\"abc\"]",
                   "True", 0);
    assert_eval_eq("LetterQ[\"ab1\"]",
                   "False", 0);
    assert_eval_eq("LetterQ[\"\"]",
                   "True", 0);
    assert_eval_eq("LetterQ[\"\\:00e9\"]",
                   "True", 0);
    assert_eval_eq("LetterQ[1]",
                   "LetterQ[1]", 0);
    assert_eval_eq("DigitQ[\"123\"]",
                   "True", 0);
    assert_eval_eq("DigitQ[\"12a\"]",
                   "False", 0);
    assert_eval_eq("DigitQ[\"\"]",
                   "True", 0);
    assert_eval_eq("DigitQ[\"\\:0662\"]",
                   "False", 0);
    assert_eval_eq("UpperCaseQ[\"ABC\"]",
                   "True", 0);
    assert_eval_eq("UpperCaseQ[\"AbC\"]",
                   "False", 0);
    assert_eval_eq("UpperCaseQ[\"AB1\"]",
                   "False", 0);
    assert_eval_eq("UpperCaseQ[\"\"]",
                   "True", 0);
    assert_eval_eq("LowerCaseQ[\"abc\"]",
                   "True", 0);
    assert_eval_eq("LowerCaseQ[\"\\:00e9\"]",
                   "True", 0);
    assert_eval_eq("LowerCaseQ[\"ab c\"]",
                   "False", 0);
    assert_eval_eq("CharacterRange[\"a\", \"e\"]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\"}", 0);
    assert_eval_eq("CharacterRange[97, 100]",
                   "{\"a\", \"b\", \"c\", \"d\"}", 0);
    assert_eval_eq("CharacterRange[\"e\", \"a\"]",
                   "{}", 0);
    assert_eval_eq("CharacterRange[\"\\[Alpha]\", \"\\[Delta]\"]",
                   "{\"α\", \"β\", \"γ\", \"δ\"}", 0);
    assert_eval_eq("Alphabet[]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"g\", \"h\", \"i\", \"j\", \"k\", \"l\", \"m\", \"n\", \"o\", \"p\", \"q\", \"r\", \"s\", \"t\", \"u\", \"v\", \"w\", \"x\", \"y\", \"z\"}", 0);
    assert_eval_eq("StringCases[\"abc\", \"a\" | \"b\", 1]",
                   "{\"a\"}", 0);
    assert_eval_eq("StringCases[\"abcabc\", \"a\" ~~ _, 1]",
                   "{\"ab\"}", 0);
    assert_eval_eq("StringCases[\"a\"][\"cab\"]",
                   "{\"a\"}", 0);
    assert_eval_eq("StringReplace[\"aaaa\", \"a\" -> \"b\", 2]",
                   "\"bbaa\"", 0);
    assert_eval_eq("StringReplace[\"a\" -> \"b\"][\"aaa\"]",
                   "\"bbb\"", 0);
    assert_eval_eq("StringReplace[\"ABab\", \"a\" -> \"x\", IgnoreCase -> True]",
                   "\"xBxb\"", 0);
    assert_eval_eq("StringMatchQ[\"a\" ~~ ___][\"abc\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"A\" ~~ ___, IgnoreCase -> True]",
                   "True", 0);
    assert_eval_eq("StringSplit[\",\"][\"a,b\"]",
                   "{\",\"}[\"a,b\"]", 0);
    assert_eval_eq("StringPosition[\"a\"][\"banana\"]",
                   "{{2, 2}, {4, 4}, {6, 6}}", 0);
    assert_eval_eq("StringCount[\"a\"][\"banana\"]",
                   "StringCount[\"a\"][\"banana\"]", 0);
    assert_eval_eq("StringPosition[\"banana\", \"a\", 2]",
                   "{{2, 2}, {4, 4}}", 0);
    assert_eval_eq("StringCount[\"banana\", \"a\", IgnoreCase -> True]",
                   "3", 0);
    assert_eval_eq("StringFreeQ[\"abc\", \"B\", IgnoreCase -> True]",
                   "False", 0);
    assert_eval_eq("StringQ[\"a\"]",
                   "True", 0);
    assert_eval_eq("StringQ[1]",
                   "False", 0);
}

/*
 * test_strings_utf8_and_specs: StringTake/StringDrop/StringPart/StringInsert/
 * StringReplacePart/StringPartition index characters (code points), and take
 * every Mathematica 14 spec form. Expected values come from Mathematica 14.0
 * (wolframscript); non-ASCII expectations are spelled as UTF-8 byte escapes.
 */
static void test_strings_utf8_and_specs(void) {
    assert_eval_eq("StringTake[\"h\\:00e9llo\", 2]",
                   "\"h\xc3\xa9\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", -3]",
                   "\"llo\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", {2}]",
                   "\"\xc3\xa9\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", {2, 4}]",
                   "\"\xc3\xa9""ll\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", {1, -1, 2}]",
                   "\"hlo\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", {-1, 1, -2}]",
                   "\"olh\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", UpTo[10]]",
                   "\"h\xc3\xa9""llo\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", All]",
                   "\"h\xc3\xa9""llo\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", None]",
                   "\"\"", 0);
    assert_eval_eq("StringTake[\"h\\:00e9llo\", {{1, 2}, {-2, -1}, All, None, UpTo[2], {}}]",
                   "{\"h\xc3\xa9\", \"lo\", \"h\xc3\xa9""llo\", \"\", \"h\xc3\xa9\", \"\"}", 0);
    assert_eval_eq("StringTake[\"abc\", {1, 2, 3, 4}]",
                   "{\"a\", \"ab\", \"abc\", StringTake[\"abc\", 4]}", 0);
    assert_eval_eq("StringTake[\"abc\", {1, 3, 0}]",
                   "{\"a\", \"abc\", \"\"}", 0);
    assert_eval_eq("StringTake[\"abc\", {{1, 2}, {2, 5}}]",
                   "{\"ab\", StringTake[\"abc\", {2, 5}]}", 0);
    assert_eval_eq("StringTake[\"abc\", {3, 2}]",
                   "\"\"", 0);
    assert_eval_eq("StringTake[\"abc\", {4, 3}]",
                   "\"\"", 0);
    assert_eval_eq("StringTake[\"abc\", {0, 3}]",
                   "\"\"", 0);
    assert_eval_eq("StringTake[\"abc\", {5, 3}]",
                   "StringTake[\"abc\", {5, 3}]", 0);
    assert_eval_eq("StringTake[\"abc\", 5]",
                   "StringTake[\"abc\", 5]", 0);
    assert_eval_eq("StringTake[\"abc\", {1, 3, -1}]",
                   "StringTake[\"abc\", {1, 3, -1}]", 0);
    assert_eval_eq("StringTake[\"abc\", 10^30]",
                   "StringTake[\"abc\", 1000000000000000000000000000000]", 0);
    assert_eval_eq("StringTake[\"abc\", 1.]",
                   "\"a\"", 0);
    assert_eval_eq("StringTake[{\"h\\:00e9llo\", \"\\[Alpha]\"}, 2]",
                   "{\"h\xc3\xa9\", StringTake[\"\xce\xb1\", 2]}", 0);
    assert_eval_eq("StringTake[{\"abc\", x}, 1]",
                   "StringTake[{\"abc\", x}, 1]", 0);
    assert_eval_eq("StringTake[{\"h\\:00e9llo\", \"\\[Alpha]\\[Beta]\"}, {{1}, {-1}}]",
                   "{{\"h\", \"o\"}, {\"\xce\xb1\", \"\xce\xb2\"}}", 0);
    assert_eval_eq("StringTake[\"\\:00e9\", {1, 1, 9223372036854775807}]",
                   "\"\xc3\xa9\"", 0);
    assert_eval_eq("StringDrop[\"h\\:00e9llo\", 2]",
                   "\"llo\"", 0);
    assert_eval_eq("StringDrop[\"h\\:00e9llo\", -2]",
                   "\"h\xc3\xa9""l\"", 0);
    assert_eval_eq("StringDrop[\"h\\:00e9llo\", {2}]",
                   "\"hllo\"", 0);
    assert_eval_eq("StringDrop[\"h\\:00e9llo\", {1, -1, 2}]",
                   "\"\xc3\xa9""l\"", 0);
    assert_eval_eq("StringDrop[\"h\\:00e9llo w\\:00f6rld\", {2, -2, 3}]",
                   "\"hll wrld\"", 0);
    assert_eval_eq("StringDrop[\"abc\", All]",
                   "\"\"", 0);
    assert_eval_eq("StringDrop[\"abc\", None]",
                   "\"abc\"", 0);
    assert_eval_eq("StringDrop[\"abc\", {}]",
                   "\"abc\"", 0);
    assert_eval_eq("StringDrop[\"abc\", UpTo[5]]",
                   "\"\"", 0);
    assert_eval_eq("StringDrop[\"abc\", {5, 4}]",
                   "StringDrop[\"abc\", {5, 4}]", 0);
    assert_eval_eq("StringDrop[\"abcdef\", {4, 2}]",
                   "StringDrop[\"abcdef\", {4, 2}]", 0);
    assert_eval_eq("StringDrop[\"abc\", {1, 3, 0}]",
                   "StringDrop[\"abc\", {1, 3, 0}]", 0);
    assert_eval_eq("StringDrop[{\"h\\:00e9llo\", \"\\[Alpha]\"}, 2]",
                   "{\"llo\", StringDrop[\"\xce\xb1\", 2]}", 0);
    assert_eval_eq("StringPart[\"h\\:00e9llo\", 2]",
                   "\"\xc3\xa9\"", 0);
    assert_eval_eq("StringPart[\"h\\:00e9llo\", {1, 2, -1}]",
                   "{\"h\", \"\xc3\xa9\", \"o\"}", 0);
    assert_eval_eq("StringPart[\"h\\:00e9llo\", All]",
                   "{\"h\", \"\xc3\xa9\", \"l\", \"l\", \"o\"}", 0);
    assert_eval_eq("StringPart[\"h\\:00e9llo\", 2 ;; 3]",
                   "{\"\xc3\xa9\", \"l\"}", 0);
    assert_eval_eq("StringPart[\"h\\:00e9llo\", -1 ;; 1 ;; -2]",
                   "{\"o\", \"l\", \"h\"}", 0);
    assert_eval_eq("StringPart[\"abc\", 0 ;; 3]",
                   "{}", 0);
    assert_eval_eq("StringPart[\"abc\", 4]",
                   "StringPart[\"abc\", 4]", 0);
    assert_eval_eq("StringPart[\"abc\", {1, 4}]",
                   "StringPart[\"abc\", {1, 4}]", 0);
    assert_eval_eq("StringPart[\"\", All]",
                   "{}", 0);
    assert_eval_eq("StringPart[{\"h\\:00e9llo\", \"\\[Alpha]\"}, 2]",
                   "{\"\xc3\xa9\", StringPart[\"\xce\xb1\", 2]}", 0);
    assert_eval_eq("StringInsert[\"h\\:00e9llo\", \"X\", 3]",
                   "\"h\xc3\xa9""Xllo\"", 0);
    assert_eval_eq("StringInsert[\"h\\:00e9llo\", \"X\", -1]",
                   "\"h\xc3\xa9""lloX\"", 0);
    assert_eval_eq("StringInsert[\"h\\:00e9llo\", \"\\[Alpha]\", {2, -2}]",
                   "\"h\xce\xb1""\xc3\xa9""ll\xce\xb1""o\"", 0);
    assert_eval_eq("StringInsert[\"abc\", \"X\", {4, 1}]",
                   "\"XabcX\"", 0);
    assert_eval_eq("StringInsert[\"abc\", \"X\", -4]",
                   "\"Xabc\"", 0);
    assert_eval_eq("StringInsert[\"abc\", \"X\", 5]",
                   "StringInsert[\"abc\", \"X\", 5]", 0);
    assert_eval_eq("StringInsert[\"abc\", \"X\", x]",
                   "StringInsert[\"abc\", \"X\", x]", 0);
    assert_eval_eq("StringInsert[{\"abc\", \"d\"}, \"X\", 3]",
                   "{\"abXc\", StringInsert[\"d\", \"X\", 3]}", 0);
    assert_eval_eq("StringReplacePart[\"h\\:00e9llo\", \"X\", {2, 3}]",
                   "\"hXlo\"", 0);
    assert_eval_eq("StringReplacePart[\"h\\:00e9llo\", {\"X\", \"Y\"}, {{1, 1}, {5, 5}}]",
                   "\"X\xc3\xa9""llY\"", 0);
    assert_eval_eq("StringReplacePart[\"abc\", \"X\", 2]",
                   "\"Xc\"", 0);
    assert_eval_eq("StringReplacePart[\"abc\", \"X\", 0]",
                   "\"Xabc\"", 0);
    assert_eval_eq("StringReplacePart[\"abc\", \"X\", {3, 2}]",
                   "\"abXc\"", 0);
    assert_eval_eq("StringReplacePart[\"abc\", \"X\", {{1, 2}, {2, 3}}]",
                   "\"Xc\"", 0);
    assert_eval_eq("StringReplacePart[\"abcde\", {\"X\", \"Y\", \"Z\"}, {{1, 5}, {2, 1}, {3, 2}}]",
                   "\"Xde\"", 0);
    assert_eval_eq("StringReplacePart[\"abc\", {\"X\", \"Y\"}, {{2, 1}, {2, 1}}]",
                   "\"aXYbc\"", 0);
    assert_eval_eq("StringReplacePart[\"abc\", \"X\", {1, 4}]",
                   "StringReplacePart[\"abc\", \"X\", {1, 4}]", 0);
    assert_eval_eq("StringReplacePart[{\"abc\", \"d\"}, \"X\", {2, 2}]",
                   "StringReplacePart[{\"abc\", \"d\"}, \"X\", {2, 2}]", 0);
    assert_eval_eq("StringPartition[\"h\\:00e9llo\", 2]",
                   "{\"h\xc3\xa9\", \"ll\"}", 0);
    assert_eval_eq("StringPartition[\"h\\:00e9llo\", 2, 1]",
                   "{\"h\xc3\xa9\", \"\xc3\xa9""l\", \"ll\", \"lo\"}", 0);
    assert_eval_eq("StringPartition[\"h\\:00e9llo\", UpTo[2]]",
                   "{\"h\xc3\xa9\", \"ll\", \"o\"}", 0);
    assert_eval_eq("StringPartition[\"abcdefg\", UpTo[3], 2]",
                   "{\"abc\", \"cde\", \"efg\", \"g\"}", 0);
    assert_eval_eq("StringPartition[\"abcd\", 0, 2]",
                   "{\"\", \"\", \"\"}", 0);
    assert_eval_eq("StringPartition[\"\", 0, 1]",
                   "{\"\"}", 0);
    assert_eval_eq("StringPartition[{\"abcd\", \"e\"}, 2]",
                   "StringPartition[{\"abcd\", \"e\"}, 2]", 0);
    assert_eval_eq("StringPartition[\"abcd\", 0]",
                   "StringPartition[\"abcd\", 0]", 0);
}

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_REGEX
    TEST(test_unicode_code_points);
    TEST(test_pad_repeat_join_riffle);
    TEST(test_new_builtins_and_forms);
    TEST(test_strings_utf8_and_specs);
#else
    printf("USE_REGEX not defined; skipping string review tests\n");
#endif

    printf("All string review tests passed!\n");
    return 0;
}
