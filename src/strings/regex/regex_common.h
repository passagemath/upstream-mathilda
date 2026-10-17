#ifndef MATHILDA_REGEX_COMMON_H
#define MATHILDA_REGEX_COMMON_H

/*
 * regex_common.h - Expr-aware helpers shared by the regex string builtins
 * (StringMatchQ / StringCases / StringCount / StringReplace / StringSplit /
 * StringPosition).
 *
 * This layer sits above the pure PCRE2 wrapper (regex_engine.h): it turns a
 * Mathilda pattern argument (RegularExpression[...], a literal string, a
 * Rule/RuleDelayed, or a List of those) into an array of compiled rules,
 * enumerates the matches of that rule set under an Overlaps policy, and
 * provides $n replacement-template expansion, trailing-option decoding, and a
 * small growable string buffer.  It compiles unconditionally; when PCRE2 is
 * absent every build fails cleanly via regex_available().
 */

#include <stddef.h>
#include "expr.h"
#include "regex_engine.h"

/*
 * One compiled pattern, optionally paired with a replacement RHS taken from a
 * Rule/RuleDelayed.  `rhs` is BORROWED from the caller's expression (never
 * owned/freed here); it is NULL for a bare pattern.
 */
/* A point in a translated pattern where Mathilda code must run mid-match
 * (emitted as a PCRE2 callout): Condition[p, c] tests `c` with the names bound
 * so far; PatternTest[p, f] tests f[text p matched], `group` capturing it. */
enum { SP_CALLOUT_CONDITION = 1, SP_CALLOUT_TEST = 2 };
typedef struct {
    int   kind;     /* SP_CALLOUT_* */
    Expr* expr;     /* borrowed from the pattern: the condition, or the test f */
    int   group;    /* SP_CALLOUT_TEST: capture group holding the tested text */
} SpCallout;

/*
 * What the string-pattern translator learned about one pattern besides its
 * regex source: which capture group binds which pattern name, so a rule RHS
 * can be instantiated with the named pieces of each match (x__ ~~ "1" :> x),
 * and the callouts (Condition / PatternTest) the matcher must evaluate.
 */
typedef struct {
    int          ngroups;     /* capture groups emitted by the translator */
    const char** group_name;  /* [0..ngroups]; interned name bound by group g,
                                 NULL for an anonymous group (owned array,
                                 borrowed names); NULL when ngroups == 0 */
    int          nnamed;      /* how many groups carry a name */
    int          has_regex;   /* pattern contains RegularExpression[...]: a
                                 string RHS is a $n template */
    int          is_null;     /* pattern is the empty string "" */
    int          raw_regex;   /* the pattern IS RegularExpression["re"]: source
                                 is the user's regex verbatim (compile without
                                 DOTALL, so `.` keeps PCRE's meaning) */
    SpCallout*   callouts;    /* [ncallouts], owned; callout k is (?C<k+1>) */
    int          ncallouts;
} SpInfo;

typedef struct {
    RegexProgram* prog;   /* owned; free via regex_rules_free */
    Expr*         rhs;    /* borrowed; NULL when the element was a bare pattern */
    Expr*         lhs;    /* borrowed; the delimiter pattern (rule LHS or bare pattern). */
    SpInfo        info;   /* owned; capture-name table for binding the RHS */
} RegexRule;

/*
 * Build an array of RegexRule from a pattern argument.  Accepts
 * RegularExpression["re"], a literal string (matched literally), Rule/
 * RuleDelayed (rhs recorded), or a List of any of those (alternatives / rule
 * set).  `anchored` is a bit set: REGEX_RULE_ANCHORED wraps each pattern
 * \A(?:...)\z for a whole-string match (StringMatchQ), REGEX_RULE_WILDCARDS
 * makes `*` and `@` in literal strings metacharacters (StringMatchQ only:
 * `*` is any run of characters, `@` one or more non-uppercase characters,
 * `\\*` / `\\@` literal). Historically `anchored` was a boolean, so 1 still
 * means REGEX_RULE_ANCHORED.  When anchored, each pattern is wrapped \A(?:...)\z for a
 * whole-string match (StringMatchQ).
 *
 * Returns the rule count (>= 1) with a malloc'd array stored in *out, or -1 on
 * an unsupported pattern, a compile error, or when regex support is absent
 * (a diagnostic keyed on `head` is printed to stderr in the latter two cases).
 * A successful result must be released with regex_rules_free().
 */
enum { REGEX_RULE_ANCHORED = 1, REGEX_RULE_WILDCARDS = 2 };
int  regex_rules_build(Expr* patt, int anchored, RegexRule** out, const char* head);
/*
 * Extended form: when `caseless` is nonzero each pattern matches case-insensitively
 * (a `(?i)` inline modifier is prepended to the compiled source). `regex_rules_build`
 * is the thin `caseless == 0` wrapper. Used by StringSplit's IgnoreCase option.
 */
int  regex_rules_build_ex(Expr* patt, int anchored, int caseless,
                          RegexRule** out, const char* head);
/*
 * General form: `rule_flags` are REGEX_RULE_* bits (as `anchored` above) and
 * `regex_flags` are engine REGEX_* compile flags passed straight through
 * (REGEX_CASELESS, REGEX_ANCHORED, REGEX_ENDANCHORED, ...).
 */
int  regex_rules_build_opts(Expr* patt, unsigned rule_flags, unsigned regex_flags,
                            RegexRule** out, const char* head);
void regex_rules_free(RegexRule* rules, int n);

/*
 * Translate a Wolfram string pattern into malloc'd PCRE source (caller frees), or
 * return NULL if the pattern is unsupported. Handles literal strings,
 * RegularExpression["re"], the character-class heads (Whitespace, WhitespaceCharacter,
 * LetterCharacter, DigitCharacter, WordCharacter, NumberString), StringExpression,
 * Alternatives, Repeated/RepeatedNull, Except, Blank, Pattern (capture group), and
 * PatternTest with a known predicate (LetterQ/DigitQ/UpperCaseQ/LowerCaseQ).
 *
 * When `is_null` is non-NULL it is set to 1 iff `patt` is the empty string "" (the
 * "split at every character" null delimiter); callers that care handle it specially.
 */
char* wl_pattern_to_regex(Expr* patt, int* is_null);

/* As wl_pattern_to_regex, additionally filling *info (when non-NULL) with the
 * capture-group name table. Release *info with sp_info_free (also safe after a
 * NULL return, when *info is zeroed). */
char* wl_pattern_to_regex_ex(Expr* patt, SpInfo* info);
/* As wl_pattern_to_regex_ex with REGEX_RULE_WILDCARDS-style literal wildcards
 * when `wildcards` is nonzero. */
char* wl_pattern_to_regex_wild(Expr* patt, SpInfo* info, int wildcards);
void  sp_info_free(SpInfo* info);

/*
 * Expand a replacement template (`$0`..`$N`, and `$$` -> literal `$`) using the
 * capture offsets `ov` (npairs pairs of byte offsets) against subject `subj`.
 * Unknown/unset groups expand to the empty string.  Returns a malloc'd
 * NUL-terminated string (caller frees), or NULL on OOM.
 */
char* regex_expand_template(const char* tpl, const char* subj,
                            const size_t* ov, size_t npairs);

/*
 * Instantiate rule `r`'s right-hand side for one match and return the new,
 * evaluated Expr (caller owns), or NULL for a bare pattern (no RHS).
 *
 *   - A string RHS of a rule whose pattern contains RegularExpression[...] is a
 *     $n template, expanded against the captures (Mathematica semantics); any
 *     other string RHS is returned literally.
 *   - Otherwise every pattern name of the LHS is bound to the substring its
 *     capture group matched ("" when that group did not participate -- e.g. the
 *     untaken branch of an Alternatives), substituted into the RHS with the
 *     matcher's own scoping-aware replace_bindings, and the result evaluated.
 *
 * `ov` holds `npairs` (start, end) byte-offset pairs (pair 0 = whole match).
 */
Expr* regex_rule_eval_rhs(const RegexRule* r, const char* subj,
                          const size_t* ov, size_t npairs);

/* ------------------------------------------------------------------ */
/* Shared match scanner                                               */
/* ------------------------------------------------------------------ */

/*
 * How overlapping matches are treated, mirroring the Overlaps option of
 * StringCases / StringCount / StringPosition.
 *
 *   REGEX_OV_TRUE   Overlapping substrings count separately, but only the first
 *                   (natural, leftmost-longest per the engine) match starting at
 *                   each position is reported.
 *   REGEX_OV_FALSE  No overlaps: greedy left-to-right, each match resumes the
 *                   scan at the end of the previous one.
 *   REGEX_OV_ALL    Every matching substring at every start, all lengths. The
 *                   caller MUST have built `rules` anchored (\A(?:...)\z), since
 *                   this mode tests whole candidate substrings.
 */
typedef enum { REGEX_OV_TRUE = 0, REGEX_OV_FALSE = 1, REGEX_OV_ALL = 2 }
        RegexOverlapMode;

/*
 * One recorded match. [ms, me) are half-open BYTE offsets into the subject;
 * `rule` is the index of the RegexRule that produced it (so a caller can find
 * the matching rule's rhs). When the scan was asked for captures,
 * caps[caps_off .. caps_off + 2*npairs) holds that match's capture offsets in
 * the same (start, end) layout regex_match writes, always subject-relative and
 * therefore directly usable with regex_expand_template; npairs is 0 otherwise.
 */
typedef struct { size_t ms, me; int rule; size_t caps_off, npairs; } RegexSpan;

/* Result of one scan. Zero-initialise (RegexScan s = {0};) and release with
 * regex_scan_free(). `caps` is NULL unless captures were requested. */
typedef struct { RegexSpan* spans; size_t count; size_t* caps; } RegexScan;

/*
 * Enumerate the matches of `rules` (nr of them) in subj[0, len) under `mode`,
 * writing them to *out ordered by start offset. Ties at the same start keep
 * rule order. This is the single shared implementation behind StringCases,
 * StringCount and StringPosition, so those three agree by construction.
 *
 * `want_captures` nonzero allocates the capture pool (only callers that expand
 * $n replacement templates need it; leaving it 0 makes a pure-counting scan
 * allocate nothing per match beyond the span array).
 *
 * Returns the match count (>= 0), or -1 on allocation failure. A successful
 * result must be released with regex_scan_free().
 */
long regex_scan(const char* subj, size_t len, RegexRule* rules, int nr,
                RegexOverlapMode mode, int want_captures, RegexScan* out);
void regex_scan_free(RegexScan* s);

/*
 * Recognise a trailing option Rule/RuleDelayed[opt_sym, value] and decode its
 * value into *value: with `overlaps` nonzero, the RegexOverlapMode (All ->
 * REGEX_OV_ALL, False -> REGEX_OV_FALSE, anything else -> REGEX_OV_TRUE);
 * otherwise a boolean, 1 iff the value is the symbol True. Returns 1 when `e`
 * is an option for `opt_sym`, 0 when it is not (and *value is untouched).
 */
int regex_match_opt(const Expr* e, const char* opt_sym, int* value, int overlaps);

/* Small growable byte buffer used to assemble result strings. Zero-initialise
 * (RegexBuf b = {0};) and release with free(b.p). */
typedef struct { char* p; size_t len, cap; } RegexBuf;

/* Append n bytes; keeps the buffer NUL-terminated. Returns 0, or -1 on OOM. */
int regexbuf_add(RegexBuf* b, const char* s, size_t n);

#endif /* MATHILDA_REGEX_COMMON_H */
