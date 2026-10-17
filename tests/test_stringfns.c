/*
 * test_stringfns.c - regression tests for the regex-aware string functions
 * StringMatchQ / StringCases / StringReplace / StringSplit, driven directly
 * from the worked examples in the RegularExpression specification.
 *
 * NOTE on escaping: Mathilda's string lexer consumes one backslash from every
 * "\x" escape, so a regex backslash must be written as "\\" in Mathilda source
 * and therefore as four backslashes ("\\\\") in this C source.  Real newlines
 * are embedded as literal '\n' bytes (Mathilda has no \n string escape).
 *
 * The whole suite is a no-op when built without PCRE2 (USE_REGEX undefined).
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"

#ifdef USE_REGEX

/* ============================= StringCases ============================= */

static void test_cases_charclass(void) {
    assert_eval_eq("StringCases[\"adefgh12c34\", RegularExpression[\"[a-e]+\"]]",
                   "{\"ade\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"a13b12c1da32efg\", RegularExpression[\"[ab]\"]]",
                   "{\"a\", \"b\", \"a\"}", 0);
    assert_eval_eq("StringCases[\"adefgh12c34\", RegularExpression[\"[a-e]\"]]",
                   "{\"a\", \"d\", \"e\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"a13b12c17a32\", RegularExpression[\"[^a1]\"]]",
                   "{\"3\", \"b\", \"2\", \"c\", \"7\", \"3\", \"2\"}", 0);
    assert_eval_eq("StringCases[\"aaaa bbbb 1234\", RegularExpression[\"[a-z]+\"]]",
                   "{\"aaaa\", \"bbbb\"}", 0);
}

static void test_cases_quantifiers(void) {
    assert_eval_eq("StringCases[\"a23b4222c63333d80\", RegularExpression[\"\\\\d+\"]]",
                   "{\"23\", \"4222\", \"63333\", \"80\"}", 0);
    assert_eval_eq("StringCases[\"aabc1aaaagh2ade\", RegularExpression[\"a{2,3}\"]]",
                   "{\"aa\", \"aaa\"}", 0);
    /* non-greedy */
    assert_eval_eq("StringCases[\"abc1agh2cde\", RegularExpression[\"a.+?\\\\d\"]]",
                   "{\"abc1\", \"agh2\"}", 0);
}

static void test_cases_shorthand_classes(void) {
    assert_eval_eq("StringCases[\"a2322c63333d80\", RegularExpression[\"\\\\d\"]]",
                   "{\"2\", \"3\", \"2\", \"2\", \"6\", \"3\", \"3\", \"3\", \"3\", \"8\", \"0\"}", 0);
    assert_eval_eq("StringCases[\"a2322c63333d80\", RegularExpression[\"\\\\D\"]]",
                   "{\"a\", \"c\", \"d\"}", 0);
    assert_eval_eq("StringCases[\"a23b42c63,d80\", RegularExpression[\"\\\\w\"]]",
                   "{\"a\", \"2\", \"3\", \"b\", \"4\", \"2\", \"c\", \"6\", \"3\", \"d\", \"8\", \"0\"}", 0);
    assert_eval_eq("StringCases[\"a23b:42c63;d80\", RegularExpression[\"\\\\W\"]]",
                   "{\":\", \";\"}", 0);
}

static void test_cases_named_class(void) {
    assert_eval_eq("StringCases[\"AaBBccDDeefG\", RegularExpression[\"[[:upper:]]+\"]]",
                   "{\"A\", \"BB\", \"DD\", \"G\"}", 0);
}

static void test_cases_literal_and_dot(void) {
    /* A bare string is a LITERAL pattern (dot matches only a dot). */
    assert_eval_eq("StringCases[\"a.b.c\", \".\"]", "{\".\", \".\"}", 0);
    /* RegularExpression["."] matches any character. */
    assert_eval_eq("StringCases[\"a1b2\", RegularExpression[\".\"]]",
                   "{\"a\", \"1\", \"b\", \"2\"}", 0);
}

static void test_cases_inline_option(void) {
    assert_eval_eq("StringCases[\"AaBbCc\", RegularExpression[\"(?i)[a-c]\"]]",
                   "{\"A\", \"a\", \"B\", \"b\", \"C\", \"c\"}", 0);
}

static void test_cases_empty(void) {
    assert_eval_eq("StringCases[\"\", RegularExpression[\"\\\\d\"]]", "{}", 0);
}

static void test_cases_thread(void) {
    assert_eval_eq("StringCases[{\"a1\", \"b22\"}, RegularExpression[\"\\\\d+\"]]",
                   "{{\"1\"}, {\"22\"}}", 0);
}

/* Overlaps / IgnoreCase, shared with StringCount and StringPosition via
 * regex_scan().  Exhaustive coverage lives in test_stringcount.c; these guard
 * the StringCases-specific bits, in particular that $n replacement templates
 * still resolve when the scan carries capture offsets. */
static void test_cases_overlaps(void) {
    /* The default is Overlaps -> False: non-overlapping, left to right. */
    assert_eval_eq("StringCases[\"AAAA\", \"AA\"]", "{\"AA\", \"AA\"}", 0);
    assert_eval_eq("StringCases[\"AAAA\", \"AA\", Overlaps -> True]",
                   "{\"AA\", \"AA\", \"AA\"}", 0);
    assert_eval_eq("StringCases[\"AAA\", x__, Overlaps -> All]",
                   "{\"AAA\", \"AA\", \"A\", \"AA\", \"A\", \"A\"}", 0);
}

static void test_cases_ignorecase(void) {
    assert_eval_eq("StringCases[\"aAbB\", \"a\", IgnoreCase -> True]",
                   "{\"a\", \"A\"}", 0);
    assert_eval_eq("StringCases[\"aAbB\", \"a\", IgnoreCase -> False]", "{\"a\"}", 0);
    assert_eval_eq("StringCases[\"ABC\", \"b\" -> \"X\", IgnoreCase -> True]",
                   "{\"X\"}", 0);
}

static void test_cases_rule_with_options(void) {
    /* $0 (whole match) under an overlapping scan. $n is a template only for a
     * RegularExpression pattern; for any other pattern "$0" is literal text. */
    assert_eval_eq("StringCases[\"AAAA\", RegularExpression[\"AA\"] -> \"$0!\", Overlaps -> True]",
                   "{\"AA!\", \"AA!\", \"AA!\"}", 0);
    assert_eval_eq("StringCases[\"AAAA\", \"AA\" -> \"$0!\", Overlaps -> True]",
                   "{\"$0!\", \"$0!\", \"$0!\"}", 0);
    /* Numbered capture groups survive the shared capture pool. */
    assert_eval_eq("StringCases[\"a1b22\", "
                   "RegularExpression[\"([a-z])(\\\\d+)\"] -> \"$2$1\"]",
                   "{\"1a\", \"22b\"}", 0);
    /* ... including in the All policy, where matches are found against
     * substrings and the capture offsets must be rebased to the subject. */
    assert_eval_eq("StringCases[\"ab\", "
                   "RegularExpression[\"([a-z])\"] -> \"<$1>\", Overlaps -> All]",
                   "{\"<a>\", \"<b>\"}", 0);
}

static void test_cases_arity(void) {
    /* Below the required arity: StringCases::argrx, call left unevaluated. */
    assert_eval_eq("StringCases[]", "StringCases[]", 0);
    assert_eval_eq("StringCases[\"abc\"]", "StringCases[\"abc\"]", 0);
    /* The occurrence-limit form keeps the first n matches. */
    assert_eval_eq("StringCases[\"aaa\", \"a\", 2]", "{\"a\", \"a\"}", 0);
    /* Non-string subject leaves the call unevaluated. */
    assert_eval_eq("StringCases[123, \"a\"]", "StringCases[123, \"a\"]", 0);
}

/* ============================ StringMatchQ ============================= */

static void test_matchq_basic(void) {
    assert_eval_eq("StringMatchQ[\"12345\", RegularExpression[\"\\\\d+\"]]", "True", 0);
    assert_eval_eq("StringMatchQ[\"12a45\", RegularExpression[\"\\\\d+\"]]", "False", 0);
    /* whole-string anchoring: a partial match is not enough */
    assert_eval_eq("StringMatchQ[\"abc\", \"ab\"]", "False", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"abc\"]", "True", 0);
}

static void test_matchq_newline(void) {
    /* Real newline bytes; \\w|\\s covers letters, digits and whitespace. */
    assert_eval_eq("StringMatchQ[\"abcd\nefgh\n1234\", RegularExpression[\"(\\\\w|\\\\s)*\"]]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"abcd\nefgh\n1234\", RegularExpression[\"(.*|\\\\s*)*\"]]",
                   "True", 0);
}

static void test_matchq_thread(void) {
    assert_eval_eq("StringMatchQ[{\"12\", \"x\"}, RegularExpression[\"\\\\d+\"]]",
                   "{True, False}", 0);
}

static void test_matchq_unevaluated(void) {
    /* Non-string subject leaves the call unevaluated. */
    assert_eval_eq("StringMatchQ[xyz, RegularExpression[\"a\"]]",
                   "StringMatchQ[xyz, RegularExpression[\"a\"]]", 0);
}

/* ============================ StringReplace =========================== */

static void test_replace_backref(void) {
    assert_eval_eq("StringReplace[\"a13b12c1da32efg\", RegularExpression[\"(\\\\d+)\"] -> \"[$1]\"]",
                   "\"a[13]b[12]c[1]da[32]efg\"", 0);
}

static void test_replace_zero_width(void) {
    /* \\b is a zero-width word boundary; :> keeps the RHS a constant. */
    assert_eval_eq("StringReplace[\"123 45 6 789\", RegularExpression[\"\\\\b\"] :> \"X\"]",
                   "\"X123X X45X X6X X789X\"", 0);
}

static void test_replace_literal(void) {
    assert_eval_eq("StringReplace[\"hello\", \"l\" -> \"L\"]", "\"heLLo\"", 0);
    assert_eval_eq("StringReplace[\"a+b\", \"+\" -> \"plus\"]", "\"aplusb\"", 0);
}

static void test_replace_multi_rule(void) {
    assert_eval_eq("StringReplace[\"abcde\", {\"a\" -> \"1\", \"c\" -> \"3\"}]",
                   "\"1b3de\"", 0);
}

/* ============================= StringSplit ============================ */

static void test_split_delimiter(void) {
    assert_eval_eq("StringSplit[\"1.23, 4.56  7.89\", RegularExpression[\"(\\\\s|,)+\"]]",
                   "{\"1.23\", \"4.56\", \"7.89\"}", 0);
    assert_eval_eq("StringSplit[\"a,b,c\", \",\"]", "{\"a\", \"b\", \"c\"}", 0);
    /* Literal delimiter: the dot is not a metacharacter here. */
    assert_eval_eq("StringSplit[\"a.b.c\", \".\"]", "{\"a\", \"b\", \"c\"}", 0);
}

static void test_split_zero_width(void) {
    /* Zero-width delimiters at line starts / ends; result printed via InputForm.
     * Mathilda prints newlines raw, so the expected strings embed real '\n'. */
    assert_eval_eq("InputForm[StringSplit[\"line1\nline2\nline3\", RegularExpression[\"(?m)^\"]]]",
                   "{\"line1\\n\", \"line2\\n\", \"line3\"}", 0);
    assert_eval_eq("InputForm[StringSplit[\"line1\nline2\nline3\", RegularExpression[\"(?m)$\"]]]",
                   "{\"line1\", \"\\nline2\", \"\\nline3\"}", 0);
}

static void test_split_no_match(void) {
    /* Delimiter absent -> single piece; empty subject -> empty list. */
    assert_eval_eq("StringSplit[\"abc\", RegularExpression[\"x\"]]", "{\"abc\"}", 0);
    assert_eval_eq("StringSplit[\"\", RegularExpression[\"x\"]]", "{}", 0);
}

static void test_split_whitespace_default(void) {
    /* One-argument form splits at runs of whitespace. */
    assert_eval_eq("StringSplit[\"a bbb  cccc aa   d\"]",
                   "{\"a\", \"bbb\", \"cccc\", \"aa\", \"d\"}", 0);
    assert_eval_eq("StringSplit[\"This is a sentence--which goes on.\"]",
                   "{\"This\", \"is\", \"a\", \"sentence--which\", \"goes\", \"on.\"}", 0);
}

static void test_split_literal_alternatives(void) {
    /* Literal delimiter (dot is not a metacharacter). */
    assert_eval_eq("StringSplit[\"192.168.0.1\", \".\"]",
                   "{\"192\", \"168\", \"0\", \"1\"}", 0);
    /* A list of delimiters and the equivalent Alternatives (|) form. */
    assert_eval_eq("StringSplit[\"a-b:c-d:e-f-g\", {\":\", \"-\"}]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"g\"}", 0);
    assert_eval_eq("StringSplit[\"a-b:c-d:e-f-g\", \":\" | \"-\"]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"g\"}", 0);
}

static void test_split_charclass(void) {
    assert_eval_eq("StringSplit[\"123  2.3  4  6\", WhitespaceCharacter ..]",
                   "{\"123\", \"2.3\", \"4\", \"6\"}", 0);
    assert_eval_eq("StringSplit[\"11a22b3\", _?LetterQ]",
                   "{\"11\", \"22\", \"3\"}", 0);
    assert_eval_eq("StringSplit[\"This is a sentence, which goes on.\", Except[WordCharacter] ..]",
                   "{\"This\", \"is\", \"a\", \"sentence\", \"which\", \"goes\", \"on\"}", 0);
}

static void test_split_stringexpression(void) {
    /* ~~ concatenates string patterns: whitespace, a digit, whitespace. */
    assert_eval_eq("StringSplit[\"primes: 2 two 3 three 5 five ...\", "
                   "Whitespace ~~ RegularExpression[\"\\\\d\"] ~~ Whitespace]",
                   "{\"primes:\", \"two\", \"three\", \"five ...\"}", 0);
}

static void test_split_rules(void) {
    /* Insert a value at each delimiter (patt -> val). */
    assert_eval_eq("StringSplit[\"a b::c d::e f g\", \"::\" -> \"--\"]",
                   "{\"a b\", \"--\", \"c d\", \"--\", \"e f g\"}", 0);
    /* :> with a named pattern keeps the delimiter itself. */
    assert_eval_eq("StringSplit[\"a--b c--d e\", x : \"--\" :> x]",
                   "{\"a\", \"--\", \"b c\", \"--\", \"d e\"}", 0);
    /* Equivalent to StringReplace when the pieces are rejoined. */
    assert_eval_eq("StringJoin[StringSplit[\"ab::c::d::ef\", \"::\" -> \"X\"]]",
                   "\"abXcXdXef\"", 0);
}

static void test_split_n_and_all(void) {
    /* At most n substrings; the remainder is the last piece. */
    assert_eval_eq("StringSplit[\"a b c d\", \" \", 2]", "{\"a\", \"b c d\"}", 0);
    /* All keeps the leading/trailing empty substrings. */
    assert_eval_eq("InputForm[StringSplit[\":a:b:c:\", \":\", All]]",
                   "{\"\", \"a\", \"b\", \"c\", \"\"}", 0);
    /* Interior empties are kept by default; boundary empties dropped. */
    assert_eval_eq("InputForm[StringSplit[\"a,,b\", \",\"]]",
                   "{\"a\", \"\", \"b\"}", 0);
    assert_eval_eq("InputForm[StringSplit[\":a:b:c:\", \":\"]]",
                   "{\"a\", \"b\", \"c\"}", 0);
}

static void test_split_null_and_ignorecase(void) {
    /* A null delimiter splits at every character. */
    assert_eval_eq("StringSplit[\"abcdefg\", \"\"]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"g\"}", 0);
    /* IgnoreCase treats c and C as equivalent. */
    assert_eval_eq("StringSplit[\"cat Cat hat CAT\", \"c\", IgnoreCase -> True]",
                   "{\"at \", \"at hat \", \"AT\"}", 0);
}

static void test_split_thread_and_regex(void) {
    assert_eval_eq("StringSplit[{\"a:b:c:d\", \"listable:element\"}, \":\"]",
                   "{{\"a\", \"b\", \"c\", \"d\"}, {\"listable\", \"element\"}}", 0);
    assert_eval_eq("StringSplit[\"A tree, an apple, four pears. And more: two sacks\", "
                   "RegularExpression[\"\\\\W+\"]]",
                   "{\"A\", \"tree\", \"an\", \"apple\", \"four\", \"pears\", "
                   "\"And\", \"more\", \"two\", \"sacks\"}", 0);
}

static void test_split_unevaluated(void) {
    /* A non-string subject leaves the call unevaluated. */
    assert_eval_eq("StringSplit[xyz, \":\"]", "StringSplit[xyz, \":\"]", 0);
}

/* The WL string-pattern vocabulary is shared: StringCases gains it too. */
static void test_cases_wl_pattern(void) {
    assert_eval_eq("StringCases[\"a bbb  cccc aa   d\", Except[WhitespaceCharacter] ..]",
                   "{\"a\", \"bbb\", \"cccc\", \"aa\", \"d\"}", 0);
    assert_eval_eq("StringReplace[\"a1b2c3\", DigitCharacter -> \"#\"]",
                   "\"a#b#c#\"", 0);
}

/* StringTrim - trims matching substrings from both ends of a string. */
static void test_trim_whitespace_default(void) {
    /* Default pattern trims leading/trailing whitespace runs. */
    assert_eval_eq("StringTrim[\"   aaa bbb ccc   \"]",
                   "\"aaa bbb ccc\"", 0);
    /* No trimmable ends: string is returned unchanged. */
    assert_eval_eq("StringTrim[\"hello\"]", "\"hello\"", 0);
    /* Tabs and real newlines count as whitespace (embedded as literal bytes). */
    assert_eval_eq("StringTrim[\"\t\n  hi there  \n\t\"]", "\"hi there\"", 0);
    /* An all-whitespace string trims to empty; empty input stays empty. */
    assert_eval_eq("StringTrim[\"   \"]", "\"\"", 0);
    assert_eval_eq("StringTrim[\"\"]", "\"\"", 0);
}

static void test_trim_pattern(void) {
    /* Multi-character alternative repeated at both ends: ("+" | "-") ... */
    assert_eval_eq("StringTrim[\"++++aaa bbb ccc----\", (\"+\" | \"-\") ...]",
                   "\"aaa bbb ccc\"", 0);
    /* One match is removed at each end (Mathematica), not a fixed point. */
    assert_eval_eq("StringTrim[\"xxabcxx\", \"x\"]", "\"xabcx\"", 0);
    /* Character-class Repeated strips digit runs from both ends. */
    assert_eval_eq("StringTrim[\"007bond007\", DigitCharacter ..]", "\"bond\"", 0);
    /* Interior matches are left untouched; only the ends are trimmed. */
    assert_eval_eq("StringTrim[\"xxaxbxx\", \"x\"]", "\"xaxbx\"", 0);
}

static void test_trim_regex_front_only(void) {
    /* RegularExpression["^ *"] anchors to the start, so only the front trims. */
    assert_eval_eq("StringTrim[\"   aaa bbb ccc   \", RegularExpression[\"^ *\"]]",
                   "\"aaa bbb ccc   \"", 0);
}

static void test_trim_thread(void) {
    /* Threads over a list of strings; a non-string element is an error
     * (StringTrim::strse) and leaves the call unevaluated, as in Mathematica. */
    assert_eval_eq("StringTrim[{\"  a  \", \"  b  \"}]", "{\"a\", \"b\"}", 0);
    assert_eval_eq("StringTrim[{\"  a  \", 5, x}]", "StringTrim[{\"  a  \", 5, x}]", 0);
}

static void test_trim_unevaluated(void) {
    /* A non-string subject leaves the call unevaluated. */
    assert_eval_eq("StringTrim[x]", "StringTrim[x]", 0);
    /* Zero arguments: unevaluated (a StringTrim::argt message is written). */
    assert_eval_eq("StringTrim[]", "StringTrim[]", 0);
}

/* ============================ StringExtract =========================== */

/* Single-level whitespace position specs: n, -n, {..}, span, All. */
static void test_extract_positions(void) {
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", 2]", "\"bbb\"", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", -1]", "\"d\"", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", 2 ;; 4]",
                   "{\"bbb\", \"cccc\", \"aa\"}", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", {1, 3}]",
                   "{\"a\", \"cccc\"}", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", {1, 4, 5}]",
                   "{\"a\", \"aa\", \"d\"}", 0);
}

/* {pos..} is equivalent to Part[StringSplit[..], {pos..}]. */
static void test_extract_part_equivalence(void) {
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", {1, 4, 5}] === "
                   "Part[StringSplit[\"a bbb  cccc aa   d\"], {1, 4, 5}]", "True", 0);
}

/* All extracts every whitespace-delimited token (== StringSplit). */
static void test_extract_all(void) {
    assert_eval_eq("StringExtract[\"A tree, an apple, four pears. And more: two sacks\", All]",
                   "{\"A\", \"tree,\", \"an\", \"apple,\", \"four\", \"pears.\", "
                   "\"And\", \"more:\", \"two\", \"sacks\"}", 0);
    assert_eval_eq("StringExtract[\"A tree, an apple, four pears. And more: two sacks\", All] === "
                   "StringSplit[\"A tree, an apple, four pears. And more: two sacks\"]", "True", 0);
}

/* sep -> pos with a literal separator and a pattern separator. */
static void test_extract_separator(void) {
    assert_eval_eq("StringExtract[\"a--bbb--ccc--dddd\", \"--\" -> 3]", "\"ccc\"", 0);
    /* sep -> All is exactly StringSplit[string, sep]. */
    assert_eval_eq("StringExtract[\"a--bbb---ccc--dddd\", \"--\" -> All] === "
                   "StringSplit[\"a--bbb---ccc--dddd\", \"--\"]", "True", 0);
    /* A string pattern as separator: runs of non-word characters. */
    assert_eval_eq("StringExtract[\"A tree, an apple, four pears. And more: two sacks\", "
                   "Except[WordCharacter] .. -> All]",
                   "{\"A\", \"tree\", \"an\", \"apple\", \"four\", \"pears\", "
                   "\"And\", \"more\", \"two\", \"sacks\"}", 0);
}

/* Multi-level extraction with whitespace/newline depth defaults. */
static void test_extract_multilevel(void) {
    /* table = "a 1\nb 2\nc 3 x" with real newline bytes. */
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", All, 1]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", 3, 1]", "\"c\"", 0);
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", 3, All]",
                   "{\"c\", \"3\", \"x\"}", 0);
    /* Three-level: "\n\n" -> 2, "\n" -> 2, Whitespace -> 3. */
    assert_eval_eq("StringExtract[\"a b c\nd e f\ng h i\n\nj k l\nm n o\np q r s"
                   "\n\nt u v\nw x y z\", 2, 2, 3]", "\"o\"", 0);
}

/* Out-of-range blocks yield Missing["PartAbsent", pos]. */
static void test_extract_missing(void) {
    assert_eval_eq("StringExtract[\"a b c\", 5]", "Missing[\"PartAbsent\", 5]", 0);
    assert_eval_eq("StringExtract[\"a b c\", -9]", "Missing[\"PartAbsent\", -9]", 0);
    /* Mapped across lines: absent third word becomes Missing per line. */
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", All, 3]",
                   "{Missing[\"PartAbsent\", 3], Missing[\"PartAbsent\", 3], \"x\"}", 0);
}

/* A list of strings threads. */
static void test_extract_thread(void) {
    assert_eval_eq("StringExtract[{\"a b c\", \"x y z\"}, 2]", "{\"b\", \"y\"}", 0);
    assert_eval_eq("StringExtract[{\"a b c\", \"x y z\"}, {1, 3}]",
                   "{{\"a\", \"c\"}, {\"x\", \"z\"}}", 0);
}

static void test_extract_unevaluated(void) {
    /* A non-string, non-list subject leaves the call unevaluated. */
    assert_eval_eq("StringExtract[xyz, 2]", "StringExtract[xyz, 2]", 0);
    /* Zero / one argument: unevaluated (a StringExtract::argm message is written). */
    assert_eval_eq("StringExtract[]", "StringExtract[]", 0);
    assert_eval_eq("StringExtract[\"abc\"]", "StringExtract[\"abc\"]", 0);
}

#endif /* USE_REGEX */

#ifdef USE_REGEX
static void test_issue98_named_captures(void) {
    /* Issue #98 and its family: pattern names bound into a rule RHS
     * (StringCases / StringReplace / StringSplit); a non-string StringReplace
     * value gives a StringExpression; $n templates only for RegularExpression.
     * Expected values from Mathematica 14. */
    assert_eval_eq("StringCases[\"abc123\", (x__ ~~ \"123\") :> x]",
                   "{\"abc\"}", 0);
    assert_eval_eq("StringCases[\"abc123\", a__ ~~ \"3\" :> {a}]",
                   "{{\"abc12\"}}", 0);
    assert_eval_eq("StringCases[\"xabc123\", \"x\" ~~ a__ ~~ \"123\" :> {a}]",
                   "{{\"abc\"}}", 0);
    assert_eval_eq("StringCases[\"abc123\", x_ ~~ \"123\" :> x]",
                   "{\"c\"}", 0);
    assert_eval_eq("StringCases[\"abc123\", a:LetterCharacter.. :> f[a]]",
                   "{f[\"abc\"]}", 0);
    assert_eval_eq("StringCases[\"abc123\", a:LetterCharacter.. :> 7]",
                   "{7}", 0);
    assert_eval_eq("StringCases[\"abc123\", a:LetterCharacter.. -> {a}]",
                   "{{\"abc\"}}", 0);
    assert_eval_eq("StringCases[\"a1b22c333\", x:DigitCharacter.. :> ToExpression[x]]",
                   "{1, 22, 333}", 0);
    assert_eval_eq("StringCases[\"a1b22c333\", x:DigitCharacter.. :> StringLength[x]]",
                   "{1, 2, 3}", 0);
    assert_eval_eq("StringReplace[\"abc123\", a__ ~~ \"123\" :> a <> \"!\"]",
                   "\"abc!\"", 0);
    assert_eval_eq("StringReplace[\"abc123\", x:DigitCharacter :> f[x]]",
                   "\"abc\"~~f[\"1\"]~~f[\"2\"]~~f[\"3\"]", 0);
    assert_eval_eq("StringReplace[\"abc123\", x:DigitCharacter :> ToExpression[x]]",
                   "\"abc\"~~1~~2~~3", 0);
    assert_eval_eq("StringReplace[\"abc\", \"b\" -> \"$0\"]",
                   "\"a$0c\"", 0);
    assert_eval_eq("StringReplace[\"abc\", RegularExpression[\"(b)\"] -> \"<$1>\"]",
                   "\"a<b>c\"", 0);
    assert_eval_eq("StringReplace[\"abc\", RegularExpression[\"b\"] -> \"$0$0\"]",
                   "\"abbc\"", 0);
    assert_eval_eq("StringCases[\"abc\", x_ ~~ \"c\" -> \"$1\"]",
                   "{\"$1\"}", 0);
    assert_eval_eq("StringCases[\"abc\", RegularExpression[\"(.)c\"] -> \"$1\"]",
                   "{\"b\"}", 0);
    assert_eval_eq("StringCases[\"abc\", RegularExpression[\"(.)c\"] :> \"$1\"]",
                   "{\"b\"}", 0);
    assert_eval_eq("StringCases[\"ab\", (x:\"a\") | \"b\" :> {x}]",
                   "{{\"a\"}, {\"\"}}", 0);
    assert_eval_eq("StringCases[\"ab\", (x:\"a\") | (y:\"b\") :> {x, y}]",
                   "{{\"a\", \"\"}, {\"\", \"b\"}}", 0);
    assert_eval_eq("StringSplit[\"a-1b-2c\", \"-\" ~~ x:DigitCharacter :> x]",
                   "{\"a\", \"1\", \"b\", \"2\", \"c\"}", 0);
    assert_eval_eq("StringSplit[\"a-1b-2c\", \"-\" ~~ x:DigitCharacter :> f[x]]",
                   "{\"a\", f[\"1\"], \"b\", f[\"2\"], \"c\"}", 0);
    assert_eval_eq("StringSplit[\"a1b2c\", x:DigitCharacter :> ToExpression[x]]",
                   "{\"a\", 1, \"b\", 2, \"c\"}", 0);
    assert_eval_eq("StringCases[\"abc123def45\", x:LetterCharacter.. ~~ y:DigitCharacter.. :> {x, y}]",
                   "{{\"abc\", \"123\"}, {\"def\", \"45\"}}", 0);
    assert_eval_eq("StringReplace[\"aXbXc\", \"X\" :> f[1]]",
                   "\"a\"~~f[1]~~\"b\"~~f[1]~~\"c\"", 0);
    assert_eval_eq("StringCases[\"ab\", x_ :> {x}]",
                   "{{\"a\"}, {\"b\"}}", 0);
    assert_eval_eq("StringReplace[\"1\", \"1\" :> 5]",
                   "StringExpression[5]", 0);
    assert_eval_eq("StringReplace[\"a1\", \"1\" :> 5]",
                   "\"a\"~~5", 0);
    assert_eval_eq("StringReplace[\"1\", \"1\" :> \"x\"]",
                   "\"x\"", 0);
    assert_eval_eq("StringReplace[\"12\", x:DigitCharacter :> ToExpression[x]]",
                   "1~~2", 0);
    assert_eval_eq("StringReplace[\"ab\", \"a\" :> {1}]",
                   "{1}~~\"b\"", 0);
    assert_eval_eq("StringReplace[\"ab\", \"a\" :> \"\"]",
                   "\"b\"", 0);
    assert_eval_eq("StringReplace[\"ab\", \"a\" :> Sequence[]]",
                   "\"b\"", 0);
    assert_eval_eq("StringReplace[{\"a1\", \"b\"}, \"1\" :> 5]",
                   "{\"a\"~~5, \"b\"}", 0);
    assert_eval_eq("StringSplit[\"a1b\", \"1\" :> 5]",
                   "{\"a\", 5, \"b\"}", 0);
    assert_eval_eq("StringSplit[\"a1b\", x:\"1\" :> {x}]",
                   "{\"a\", {\"1\"}, \"b\"}", 0);
    assert_eval_eq("StringCases[\"a1b\", \"1\" :> Sequence[]]",
                   "{}", 0);
    assert_eval_eq("StringCases[\"a1b2\", x:DigitCharacter :> Sequence[x, x]]",
                   "{\"1\", \"1\", \"2\", \"2\"}", 0);
    assert_eval_eq("StringCases[\"ab\", x_ :> If[x == \"a\", 1, Nothing]]",
                   "{1}", 0);
}
#endif

#ifdef USE_REGEX
static void test_string_pattern_language(void) {
    /* The string-pattern language against Mathematica 14: repeated names
     * (backreferences, also across | and inside ..), Shortest/Longest,
     * Repeated[p, spec], Condition and PatternTest via callouts, blanks match
     * newlines (RegularExpression keeps PCRE's .), line anchors, character
     * classes, StringMatchQ wildcards * and @, StringTrim one match per end,
     * IgnoreCase. */
    assert_eval_eq("StringCases[\"ab\", (x:\"a\") | (x:\"b\") :> {x}]",
                   "{{\"a\"}}", 0);
    assert_eval_eq("StringMatchQ[\"aa\", (x:\"a\"|\"b\") ~~ x_]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"ab\", (x:\"a\"|\"b\") ~~ x_]",
                   "False", 0);
    assert_eval_eq("StringCases[\"aaabbc\", (x:LetterCharacter)..]",
                   "{\"aaa\", \"bb\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"aaabbc\", (x:LetterCharacter).. :> x]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"aaabbc\", x:(LetterCharacter..)]",
                   "{\"aaabbc\"}", 0);
    assert_eval_eq("StringMatchQ[\"abab\", (x:\"ab\")..]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"abab\", x__ ~~ x__]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"abcabc\", x__ ~~ x__]",
                   "True", 0);
    assert_eval_eq("StringCases[\"abcabcx\", x__ ~~ x__ :> x]",
                   "{\"abc\"}", 0);
    assert_eval_eq("StringCases[\"aXbXc\", Shortest[\"a\" ~~ __ ~~ \"X\"]]",
                   "{\"aXbX\"}", 0);
    assert_eval_eq("StringCases[\"aXbXc\", \"a\" ~~ __ ~~ \"X\"]",
                   "{\"aXbX\"}", 0);
    assert_eval_eq("StringCases[\"aXbXc\", Longest[\"a\" ~~ __ ~~ \"X\"]]",
                   "{\"aXbX\"}", 0);
    assert_eval_eq("StringCases[\"<a><b>\", Shortest[\"<\" ~~ x__ ~~ \">\"] :> x]",
                   "{\"a\", \"b\"}", 0);
    assert_eval_eq("StringCases[\"aaaa\", Repeated[\"a\", 2]]",
                   "{\"aa\", \"aa\"}", 0);
    assert_eval_eq("StringCases[\"aaaaa\", Repeated[\"a\", {2, 3}]]",
                   "{\"aaa\", \"aa\"}", 0);
    assert_eval_eq("StringCases[\"aaaaa\", Repeated[\"a\", {2}]]",
                   "{\"aa\", \"aa\"}", 0);
    assert_eval_eq("StringCases[\"aaaaa\", RepeatedNull[\"a\", 2]]",
                   "{\"aa\", \"aa\", \"a\", \"\"}", 0);
    assert_eval_eq("StringCases[\"a1b22c333\", x:DigitCharacter.. /; StringLength[x] > 1]",
                   "{\"22\", \"333\"}", 0);
    assert_eval_eq("StringCases[\"a1b22c333\", (x:DigitCharacter.. /; StringLength[x] > 1) :> ToExpression[x]]",
                   "{22, 333}", 0);
    assert_eval_eq("StringCases[\"a1b22c333\", x:DigitCharacter.. :> ToExpression[x] /; StringLength[x] > 1]",
                   "{1 /; StringLength[\"1\"] > 1, 22 /; StringLength[\"22\"] > 1, 333 /; StringLength[\"333\"] > 1}", 0);
    assert_eval_eq("StringCases[\"abc123\", _?DigitQ]",
                   "{\"1\", \"2\", \"3\"}", 0);
    assert_eval_eq("StringCases[\"abc123\", __?DigitQ]",
                   "{\"123\"}", 0);
    assert_eval_eq("StringCases[\"abc123\", x_?(# == \"b\" &)]",
                   "{\"b\"}", 0);
    assert_eval_eq("StringCases[\"abc123\", x_?(# == \"b\" &) :> x <> x]",
                   "{\"bb\"}", 0);
    assert_eval_eq("StringCases[\"abc123\", (x__?LetterQ) ~~ y__ :> {x, y}]",
                   "{{\"abc\", \"123\"}}", 0);
    assert_eval_eq("StringCases[\"aB3\", _?UpperCaseQ]",
                   "{\"B\"}", 0);
    assert_eval_eq("StringCases[\"aB3\", _?LowerCaseQ]",
                   "{\"a\"}", 0);
    assert_eval_eq("StringMatchQ[\"a\\nb\", ___]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"a\\nb\", \"a\" ~~ _ ~~ \"b\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"a\\nb\", RegularExpression[\"a.b\"]]",
                   "False", 0);
    assert_eval_eq("StringCases[\"ab\\ncd\", StartOfLine ~~ x_ :> x]",
                   "{\"a\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"ab\\ncd\", x_ ~~ EndOfLine :> x]",
                   "{\"b\", \"d\"}", 0);
    assert_eval_eq("StringCases[\"one two\", WordBoundary ~~ x_ :> x]",
                   "{\"o\", \" \", \"t\"}", 0);
    assert_eval_eq("StringCases[\"abcxyz\", CharacterRange[\"a\", \"c\"]..]",
                   "{\"abc\"}", 0);
    assert_eval_eq("StringCases[\"0x1fZ\", HexadecimalCharacter..]",
                   "{\"0\", \"1f\"}", 0);
    assert_eval_eq("StringCases[\"a,b.c!\", PunctuationCharacter]",
                   "{\",\", \".\", \"!\"}", 0);
    assert_eval_eq("StringCases[\"x 1.5 -2 3e5 .7\", NumberString]",
                   "{\"1.5\", \"-2\", \"3\", \"5\", \".7\"}", 0);
    assert_eval_eq("StringCases[\"AbC\", \"b\", IgnoreCase -> True]",
                   "{\"b\"}", 0);
    assert_eval_eq("StringCases[\"a b\", Except[WhitespaceCharacter]..]",
                   "{\"a\", \"b\"}", 0);
    assert_eval_eq("StringCases[\"abc\", Except[\"b\"]]",
                   "{\"a\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"abc\", Except[\"b\", LetterCharacter]]",
                   "{\"a\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"a1b2\", Except[DigitCharacter]..]",
                   "{\"a\", \"b\"}", 0);
    assert_eval_eq("StringCases[\"hello world\", WordCharacter..]",
                   "{\"hello\", \"world\"}", 0);
    assert_eval_eq("StringCases[\"a  b\\tc\", Whitespace]",
                   "{\"  \", \"\\t\"}", 0);
    assert_eval_eq("StringCases[\"abc\", _ ~~ _]",
                   "{\"ab\"}", 0);
    assert_eval_eq("StringCases[\"abcd\", {\"ab\", \"bc\", \"cd\"}]",
                   "{\"ab\", \"cd\"}", 0);
    assert_eval_eq("StringCases[\"abcd\", \"ab\" | \"bc\" | \"cd\", Overlaps -> True]",
                   "{\"ab\", \"bc\", \"cd\"}", 0);
    assert_eval_eq("StringCases[\"abc\", __, Overlaps -> All]",
                   "{\"abc\", \"ab\", \"a\", \"bc\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringPosition[\"abcabc\", \"bc\"]",
                   "{{2, 3}, {5, 6}}", 0);
    assert_eval_eq("StringCount[\"aaaa\", \"aa\"]",
                   "2", 0);
    assert_eval_eq("StringCount[\"aaaa\", \"aa\", Overlaps -> True]",
                   "3", 0);
    assert_eval_eq("StringCases[\"The cat\", StartOfString ~~ x_ :> x]",
                   "{\"T\"}", 0);
    assert_eval_eq("StringCases[\"The cat\", x_ ~~ EndOfString :> x]",
                   "{\"t\"}", 0);
    assert_eval_eq("StringReplace[\"a.b\", \".\" -> \"!\"]",
                   "\"a!b\"", 0);
    assert_eval_eq("StringReplace[\"a+b\", \"+\" -> \"-\"]",
                   "\"a-b\"", 0);
    assert_eval_eq("StringCases[\"a(b)c\", \"(\" ~~ x_ ~~ \")\" :> x]",
                   "{\"b\"}", 0);
    assert_eval_eq("StringCases[\"x^2\", \"^\"]",
                   "{\"^\"}", 0);
    assert_eval_eq("StringCases[\"a\\\\b\", \"\\\\\"]",
                   "{\"\\\\\"}", 0);
    assert_eval_eq("StringCases[\"abc\", Alternatives[]]",
                   "{}", 0);
    assert_eval_eq("StringCases[\"abc\", RegularExpression[\"(?i)B\"]]",
                   "{\"b\"}", 0);
    assert_eval_eq("StringCases[\"ABC\", RegularExpression[\"b\"], IgnoreCase -> True]",
                   "{\"B\"}", 0);
    assert_eval_eq("StringCases[\"aaa\", x:\"a\" ~~ y:\"a\" :> {x, y}]",
                   "{{\"a\", \"a\"}}", 0);
    assert_eval_eq("StringCases[\"abcd\", (x:\"a\"|\"c\") ~~ y_ :> x <> y]",
                   "{\"ab\", \"cd\"}", 0);
    assert_eval_eq("StringReplace[\"hello\", Shortest[x__] ~~ \"l\" :> x]",
                   "\"helo\"", 0);
    assert_eval_eq("StringCases[\"aa\", x_ ~~ x_ :> x]",
                   "{\"a\"}", 0);
    assert_eval_eq("StringCases[\"ab\", x_ ~~ x_ :> x]",
                   "{}", 0);
    assert_eval_eq("StringCases[\"abc\", (x_ ~~ y_) /; x != y :> x <> y]",
                   "{\"ab\"}", 0);
    assert_eval_eq("StringCases[\"aabc\", (x_ ~~ y_) /; x != y :> x <> y]",
                   "{\"ab\"}", 0);
    assert_eval_eq("StringCases[\"a1\", LetterCharacter ~~ DigitCharacter]",
                   "{\"a1\"}", 0);
    assert_eval_eq("StringCases[\"!\\\"#$%&'()*+,-./:;<=>?@[\\\\]^_`{|}~\", PunctuationCharacter]",
                   "{\"!\", \"\\\"\", \"#\", \"$\", \"%\", \"&\", \"'\", \"(\", \")\", \"*\", \"+\", \",\", \"-\", \".\", \"/\", \":\", \";\", \"<\", \"=\", \">\", \"?\", \"@\", \"[\", \"\\\\\", \"]\", \"^\", \"_\", \"`\", \"{\", \"|\", \"}\", \"~\"}", 0);
    assert_eval_eq("StringCases[\"a_b-c\", WordCharacter..]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"a b\\tc\\nd\", WhitespaceCharacter]",
                   "{\" \", \"\\t\", \"\\n\"}", 0);
    assert_eval_eq("StringCases[\"0123456789abcdefABCDEFgG\", HexadecimalCharacter..]",
                   "{\"0123456789abcdefABCDEF\"}", 0);
    assert_eval_eq("StringCases[\"+5 -.5 1. 007 1e3\", NumberString]",
                   "{\"+5\", \"-.5\", \"1.\", \"007\", \"1\", \"3\"}", 0);
    assert_eval_eq("StringCases[\"aaa\", \"\"]",
                   "{\"\", \"\", \"\", \"\"}", 0);
    assert_eval_eq("StringCases[\"abc\", \"\" ~~ x_ :> x]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"ab\", x___ :> {x}]",
                   "{{\"ab\"}, {\"\"}}", 0);
    assert_eval_eq("StringCases[\"aab\", Shortest[x__] :> x]",
                   "{\"a\", \"a\", \"b\"}", 0);
    assert_eval_eq("StringCases[\"aabb\", Shortest[a__ ~~ b__] :> {a, b}]",
                   "{{\"a\", \"a\"}, {\"b\", \"b\"}}", 0);
    assert_eval_eq("StringCases[\"abc\", x:Repeated[_, {2}] :> x]",
                   "{\"ab\"}", 0);
    assert_eval_eq("StringCases[\"aaaa\", (x:\"a\") ~~ (x..) :> x]",
                   "StringCases[\"aaaa\", x:\"a\"~~x.. :> x]", 0);
    assert_eval_eq("StringCases[\"abcabc\", (x:(\"a\"|\"b\")) ~~ __ ~~ x :> x]",
                   "StringCases[\"abcabc\", x:\"a\" | \"b\"~~__~~x :> x]", 0);
    assert_eval_eq("StringCases[\"x1y22\", (l:LetterCharacter) ~~ (d:DigitCharacter..) /; StringLength[d] == 2 :> l]",
                   "{\"y\"}", 0);
    assert_eval_eq("StringCases[\"a1b2\", LetterCharacter ~~ d:DigitCharacter /; ToExpression[d] > 1]",
                   "{\"b2\"}", 0);
    assert_eval_eq("StringCases[\"abc\", Except[\"b\"]..]",
                   "{\"a\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"one two  three\", Except[\" \"]..]",
                   "{\"one\", \"two\", \"three\"}", 0);
    assert_eval_eq("StringReplace[\"a\\nb\", _ -> \"x\"]",
                   "\"xxx\"", 0);
    assert_eval_eq("StringSplit[\"a\\nb\\nc\", \"\\n\"]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringSplit[\"a1b2c3\", DigitCharacter]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringSplit[\"a b  c\"]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringSplit[\" a b \", \" \", All]",
                   "{\"\", \"a\", \"b\", \"\"}", 0);
    assert_eval_eq("StringSplit[\"a,b;c\", {\",\", \";\"}]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringSplit[\"aXbxc\", \"x\", IgnoreCase -> True]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringCount[\"abcabc\", \"a\" | \"b\"]",
                   "4", 0);
    assert_eval_eq("StringPosition[\"abcabc\", \"a\" | \"bc\"]",
                   "{{1, 1}, {2, 3}, {4, 4}, {5, 6}}", 0);
    assert_eval_eq("StringPosition[\"aaa\", \"aa\"]",
                   "{{1, 2}, {2, 3}}", 0);
    assert_eval_eq("StringPosition[\"aaa\", \"aa\", Overlaps -> False]",
                   "{{1, 2}}", 0);
    assert_eval_eq("StringPosition[\"abcabc\", \"bc\", 1]",
                   "{{2, 3}}", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"a*\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"a@\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"ABC\", \"abc\", IgnoreCase -> True]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[{\"a1\", \"b\"}, LetterCharacter ~~ DigitCharacter]",
                   "{True, False}", 0);
    assert_eval_eq("StringContainsQ[\"abc\", \"B\", IgnoreCase -> True]",
                   "True", 0);
    assert_eval_eq("StringStartsQ[\"abc\", \"ab\"]",
                   "True", 0);
    assert_eval_eq("StringEndsQ[\"abc\", \"bc\"]",
                   "True", 0);
    assert_eval_eq("StringFreeQ[\"abc\", \"d\"]",
                   "True", 0);
    assert_eval_eq("StringTrim[\"  ab  \"]",
                   "\"ab\"", 0);
    assert_eval_eq("StringTrim[\"xxabxx\", \"x\"]",
                   "\"xabx\"", 0);
    assert_eval_eq("StringTrim[\"xxabxx\", \"x\"..]",
                   "\"ab\"", 0);
    assert_eval_eq("StringCases[\"The Quick\", CharacterRange[\"A\", \"Z\"]]",
                   "{\"T\", \"Q\"}", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"a\" ~~ \"*\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"a*\", \"a\\\\*\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"ab\", \"a\\\\*\"]",
                   "False", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"a@\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"aBc\", \"a@\"]",
                   "False", 0);
    assert_eval_eq("StringMatchQ[\"a\", \"a@\"]",
                   "False", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"*c\" | \"x\"]",
                   "True", 0);
    assert_eval_eq("StringContainsQ[\"abc\", \"a*\"]",
                   "False", 0);
    assert_eval_eq("StringCases[\"a*b\", \"a*\"]",
                   "{\"a*\"}", 0);
    assert_eval_eq("StringPosition[\"abc\", \"a*\"]",
                   "{}", 0);
    assert_eval_eq("StringFreeQ[\"abc\", \"a*\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"\", \"*\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"ABC\", \"a*\", IgnoreCase -> True]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"a\\nb\", \"a*\"]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"x\", RegularExpression[\"x*\"]]",
                   "True", 0);
    assert_eval_eq("StringTrim[\"xxabcxx\", \"x\"]",
                   "\"xabcx\"", 0);
    assert_eval_eq("StringTrim[\"x\", \"x\"]",
                   "\"\"", 0);
    assert_eval_eq("StringTrim[\"xx\", \"x\"]",
                   "\"\"", 0);
    assert_eval_eq("StringTrim[\"xxx\", \"xx\"]",
                   "\"x\"", 0);
    assert_eval_eq("StringTrim[\"aba\", \"a\" ~~ ___]",
                   "\"\"", 0);
    assert_eval_eq("StringTrim[\" a\\n\"]",
                   "\"a\"", 0);
    assert_eval_eq("StringTrim[{\" a \", \"b  \"}]",
                   "{\"a\", \"b\"}", 0);
    assert_eval_eq("StringTrim[\"abcab\", \"ab\"]",
                   "\"c\"", 0);
    assert_eval_eq("StringTrim[\"abc\", \"\"]",
                   "\"abc\"", 0);
    assert_eval_eq("StringTrim[\"123abc456\", DigitCharacter..]",
                   "\"abc\"", 0);
    assert_eval_eq("StringTrim[\"abc\", x_ /; x == \"a\"]",
                   "\"bc\"", 0);
    assert_eval_eq("StringTrim[\"ABcab\", \"ab\", IgnoreCase -> True]",
                   "\"c\"", 0);
}
#endif

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_REGEX
    TEST(test_cases_charclass);
    TEST(test_cases_quantifiers);
    TEST(test_cases_shorthand_classes);
    TEST(test_cases_named_class);
    TEST(test_cases_literal_and_dot);
    TEST(test_cases_inline_option);
    TEST(test_cases_empty);
    TEST(test_cases_thread);
    TEST(test_cases_overlaps);
    TEST(test_cases_ignorecase);
    TEST(test_cases_rule_with_options);
    TEST(test_cases_arity);

    TEST(test_matchq_basic);
    TEST(test_matchq_newline);
    TEST(test_matchq_thread);
    TEST(test_matchq_unevaluated);

    TEST(test_replace_backref);
    TEST(test_replace_zero_width);
    TEST(test_replace_literal);
    TEST(test_replace_multi_rule);

    TEST(test_split_delimiter);
    TEST(test_split_zero_width);
    TEST(test_split_no_match);
    TEST(test_split_whitespace_default);
    TEST(test_split_literal_alternatives);
    TEST(test_split_charclass);
    TEST(test_split_stringexpression);
    TEST(test_split_rules);
    TEST(test_split_n_and_all);
    TEST(test_split_null_and_ignorecase);
    TEST(test_split_thread_and_regex);
    TEST(test_split_unevaluated);
    TEST(test_cases_wl_pattern);
    TEST(test_issue98_named_captures);
    TEST(test_string_pattern_language);

    TEST(test_trim_whitespace_default);
    TEST(test_trim_pattern);
    TEST(test_trim_regex_front_only);
    TEST(test_trim_thread);
    TEST(test_trim_unevaluated);

    TEST(test_extract_positions);
    TEST(test_extract_part_equivalence);
    TEST(test_extract_all);
    TEST(test_extract_separator);
    TEST(test_extract_multilevel);
    TEST(test_extract_missing);
    TEST(test_extract_thread);
    TEST(test_extract_unevaluated);
#else
    printf("USE_REGEX not defined; skipping regex string-function tests\n");
#endif

    printf("All string-function tests passed!\n");
    return 0;
}
