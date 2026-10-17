#ifndef MATHILDA_REGEX_ENGINE_H
#define MATHILDA_REGEX_ENGINE_H

/*
 * regex_engine.h - thin, Expr-agnostic wrapper around PCRE2.
 *
 * This is the ONLY translation unit in the project that includes <pcre2.h>.
 * The wrapper is compiled unconditionally; when Mathilda is built without
 * PCRE2 (USE_REGEX undefined) every entry point degrades to a stub so the
 * rest of the string subsystem links and runs, with regex_available()
 * returning 0 and regex_compile() returning NULL.  This mirrors the
 * USE_MPFR / USE_LAPACK / USE_GRAPHICS graceful-degrade policy.
 *
 * Byte semantics: all offsets are byte offsets into the (raw char*) subject,
 * also under REGEX_UTF (callers convert to character positions).
 *
 * Re-entrancy: a match may run Mathilda code through a callout, and that code
 * may match with the same program; such a nested match uses its own match
 * data, so the outer match's state is never clobbered.
 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque compiled program; owns the underlying pcre2_code + match data. */
typedef struct RegexProgram RegexProgram;

/* Returns 1 when built with PCRE2 support, 0 otherwise. */
int regex_available(void);

/*
 * Compile `pattern` (a NUL-terminated PCRE regular expression).
 *
 * On success returns a non-NULL owned program; free it with regex_free().
 * On failure (syntax error, OOM, or built without USE_REGEX) returns NULL and,
 * when errbuf != NULL and errcap > 0, writes a NUL-terminated diagnostic into
 * errbuf.
 */
RegexProgram* regex_compile(const char* pattern, char* errbuf, size_t errcap);

/* Compile flags for regex_compile_ex. */
enum {
    REGEX_CASELESS = 1u << 0,   /* case-insensitive (PCRE2_CASELESS) */
    REGEX_DOTALL   = 1u << 1,   /* `.` also matches newline (PCRE2_DOTALL) */
    REGEX_UTF      = 1u << 2,   /* UTF-8 subject/pattern, Unicode classes
                                   (PCRE2_UTF | PCRE2_UCP | MATCH_INVALID_UTF) */
    REGEX_ANCHORED = 1u << 3,   /* match only at the start offset */
    REGEX_ENDANCHORED = 1u << 4 /* match must end at the subject's end */
};

/* Largest number of (start, end) offset pairs exposed per match: the whole
 * match plus up to 63 capture groups. */
#define REGEX_MAX_PAIRS 64

/*
 * Callout hook, run when the engine reaches a (?C<n>) point. `ov` holds
 * `npairs` (start, end) byte offsets of the captures set so far (pair 0 is the
 * match in progress: its start and the current position). Return 0 to
 * continue, nonzero to fail at this point (the engine backtracks).
 */
typedef int (*RegexCalloutFn)(int number, const char* subj, size_t len,
                              const size_t* ov, size_t npairs, void* data);

/*
 * As regex_compile, with REGEX_* `flags` and an optional callout hook (NULL:
 * none). With a hook, auto-possessification is disabled so a failing callout
 * can backtrack into the preceding quantifier, as a failed Condition must.
 */
RegexProgram* regex_compile_ex(const char* pattern, unsigned flags,
                               RegexCalloutFn callout, void* callout_data,
                               char* errbuf, size_t errcap);

/* Free a program returned by regex_compile(). NULL is ignored. */
void regex_free(RegexProgram* prog);

/* Number of capture groups (excluding whole-match group 0). */
int regex_group_count(const RegexProgram* prog);

/*
 * Attempt one match at or after byte offset `start` in subj[0..len).
 *
 * `ov` receives (start, end) byte offsets for group 0..ngroups, i.e. it must
 * have room for at least ov_pairs*2 size_t entries; unset/absent groups get
 * (REGEX_UNSET, REGEX_UNSET).  Returns 1 on match, 0 on no match, -1 on error.
 */
int regex_match(RegexProgram* prog, const char* subj, size_t len, size_t start,
                size_t* ov, size_t ov_pairs);

/* Sentinel written into `ov` for an unset capture group. */
#define REGEX_UNSET ((size_t)-1)

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_REGEX_ENGINE_H */
