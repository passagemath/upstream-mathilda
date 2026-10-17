/*
 * regex_common.c - implementation of the Expr-aware regex helpers.
 */

#include "regex_common.h"
#include "sym_names.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "match.h"     /* MatchEnv / replace_bindings: RHS instantiation */
#include "eval.h"

/* ------------------------------------------------------------------ */
/* Growable byte buffer                                               */
/* ------------------------------------------------------------------ */

int regexbuf_add(RegexBuf* b, const char* s, size_t n) {
    if (b->len + n + 1 > b->cap) {
        size_t nc = b->cap ? b->cap : 32;
        while (nc < b->len + n + 1) nc *= 2;
        char* np = realloc(b->p, nc);
        if (!np) return -1;
        b->p = np;
        b->cap = nc;
    }
    memcpy(b->p + b->len, s, n);
    b->len += n;
    b->p[b->len] = '\0';
    return 0;
}

/* ------------------------------------------------------------------ */
/* Small local helpers                                                */
/* ------------------------------------------------------------------ */

/* C99-safe strdup replacement. */
static char* rc_strdup(const char* s) {
    size_t n = strlen(s) + 1;
    char* d = malloc(n);
    if (d) memcpy(d, s, n);
    return d;
}

/* Length of the leading PCRE start-of-pattern verbs -- (*UTF), (*UCP), (*CR),
 * ... -- which must stay at the very front of the compiled source. */
static size_t leading_verbs_len(const char* src) {
    size_t i = 0;
    while (src[i] == '(' && src[i + 1] == '*') {
        const char* close = strchr(src + i, ')');
        if (!close) break;
        i = (size_t)(close - src) + 1;
    }
    return i;
}

/* Wrap regex source as \A(?:src)\z for a whole-string (anchored) match,
 * keeping any leading (*VERB)s in front. */
static char* wrap_anchored(const char* src) {
    size_t v = leading_verbs_len(src);
    size_t n = strlen(src);
    char* out = malloc(n + 9);              /* "\A(?:" (5) + src + ")\z" (3) + NUL */
    if (!out) return NULL;
    memcpy(out, src, v);
    memcpy(out + v, "\\A(?:", 5);
    memcpy(out + v + 5, src + v, n - v);
    memcpy(out + 5 + n, ")\\z", 3);
    out[8 + n] = '\0';
    return out;
}

/* New EXPR_STRING holding subj[s, e). */
static Expr* rc_substr(const char* subj, size_t s, size_t e) {
    size_t n = e - s;
    char* buf = malloc(n + 1);
    if (!buf) return expr_new_string("");
    memcpy(buf, subj + s, n);
    buf[n] = '\0';
    Expr* r = expr_new_string(buf);
    free(buf);
    return r;
}

/*
 * A MatchEnv binding each pattern name of rule `r` to the text of the first
 * participating group that carries it (a name repeated across Alternatives
 * branches owns a group in the branch where it first occurs). A name with no
 * participating group is bound to "" when `unset_empty` (a finished match:
 * Mathematica binds the untaken branch's name to ""), else left unbound (a
 * Condition mid-match: the name has not been reached yet).
 */
static MatchEnv* rule_bind_env(const RegexRule* r, const char* subj,
                               const size_t* ov, size_t npairs, int unset_empty) {
    MatchEnv* env = env_new();
    const char* const* gname = r->info.group_name;
    for (int g = 1; g <= r->info.ngroups; g++) {
        const char* nm = gname[g];
        if (!nm || env_get(env, nm)) continue;
        Expr* v = NULL;
        for (int h = g; h <= r->info.ngroups && !v; h++) {
            if (gname[h] == nm && (size_t)h < npairs && ov[2 * h] != REGEX_UNSET)
                v = rc_substr(subj, ov[2 * h], ov[2 * h + 1]);
        }
        if (!v && unset_empty) v = expr_new_string("");
        if (v) { env_set(env, nm, v); expr_free(v); }
    }
    return env;
}

/*
 * Callout hook for a translated pattern (see string_pattern.c): evaluate the
 * Condition / PatternTest it marks against the match in progress. Returns 0 to
 * let the match continue, 1 to make the engine backtrack.
 */
static int rule_callout(int number, const char* subj, size_t len,
                        const size_t* ov, size_t npairs, void* data) {
    (void)len;
    const RegexRule* r = (const RegexRule*)data;
    if (number < 1 || number > r->info.ncallouts) return 0;
    const SpCallout* c = &r->info.callouts[number - 1];

    Expr* test;
    if (c->kind == SP_CALLOUT_TEST) {
        Expr* text = ((size_t)c->group < npairs && ov[2 * c->group] != REGEX_UNSET)
                   ? rc_substr(subj, ov[2 * c->group], ov[2 * c->group + 1])
                   : expr_new_string("");
        Expr* targs[1] = { text };
        test = expr_new_function(expr_copy(c->expr), targs, 1);
    } else {
        MatchEnv* env = rule_bind_env(r, subj, ov, npairs, /*unset_empty=*/0);
        test = replace_bindings(c->expr, env);
        env_free(env);
    }
    Expr* v = eval_and_free(test);
    int pass = v && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True;
    if (v) expr_free(v);
    return pass ? 0 : 1;
}

/* Is e a Rule[...] or RuleDelayed[...] with two arguments? */
static int is_rule2(Expr* e, Expr** lhs, Expr** rhs) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol.name == SYM_Rule ||
         e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
        e->data.function.arg_count == 2) {
        *lhs = e->data.function.args[0];
        *rhs = e->data.function.args[1];
        return 1;
    }
    return 0;
}

/* Compile one element into `out`. Returns 0 on success, -1 on failure. */
static int build_one(Expr* e, unsigned anchored, unsigned regex_flags,
                     RegexRule* out, const char* head) {
    Expr* patt = e;
    Expr* rhs = NULL;
    Expr* l; Expr* r;
    if (is_rule2(e, &l, &r)) { patt = l; rhs = r; }

    /* Shared Wolfram-string-pattern translator (literal strings,
     * RegularExpression, character classes, StringExpression, Alternatives,
     * Repeated, Except, ...), keeping its capture-name table for the RHS. */
    char* src = wl_pattern_to_regex_wild(patt, &out->info,
                                         anchored & REGEX_RULE_WILDCARDS);
    if (!src) return -1;                    /* unsupported pattern -> unevaluated */

    if (anchored & REGEX_RULE_ANCHORED) {
        char* a = wrap_anchored(src);
        free(src);
        if (!a) { sp_info_free(&out->info); return -1; }
        src = a;
    }

    /* A translated pattern runs DOTALL (blanks match newlines); a bare
     * RegularExpression keeps PCRE's defaults. */
    unsigned flags = regex_flags;
    if (!out->info.raw_regex)  flags |= REGEX_DOTALL;

    char err[256];
    RegexProgram* prog = regex_compile_ex(src, flags,
                                          out->info.ncallouts ? rule_callout : NULL,
                                          out, err, sizeof err);
    free(src);
    if (!prog) {
        mth_message(head, "regex", "%s", err);
        sp_info_free(&out->info);
        return -1;
    }
    out->prog = prog;
    out->rhs = rhs;
    out->lhs = patt;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public: build / free a rule set                                    */
/* ------------------------------------------------------------------ */

int regex_rules_build(Expr* patt, int anchored, RegexRule** out, const char* head) {
    return regex_rules_build_ex(patt, anchored, /*caseless=*/0, out, head);
}

int regex_rules_build_ex(Expr* patt, int anchored, int caseless,
                         RegexRule** out, const char* head) {
    return regex_rules_build_opts(patt, (unsigned)anchored,
                                  caseless ? REGEX_CASELESS : 0u, out, head);
}

int regex_rules_build_opts(Expr* patt, unsigned rule_flags, unsigned regex_flags,
                           RegexRule** out, const char* head) {
    if (!regex_available()) {
        mth_message(head, "regavail",
                "regular-expression support is not available; "
                "rebuild Mathilda with PCRE2 (USE_REGEX=1).");
        return -1;
    }

    Expr** elems;
    int n;
    if (patt->type == EXPR_FUNCTION &&
        patt->data.function.head->type == EXPR_SYMBOL &&
        patt->data.function.head->data.symbol.name == SYM_List) {
        n = (int)patt->data.function.arg_count;
        elems = patt->data.function.args;
        if (n == 0) return -1;
    } else {
        n = 1;
        elems = &patt;
    }

    RegexRule* rules = calloc((size_t)n, sizeof(RegexRule));
    if (!rules) return -1;

    for (int i = 0; i < n; i++) {
        if (build_one(elems[i], rule_flags, regex_flags, &rules[i], head) != 0) {
            regex_rules_free(rules, i);     /* free the i already built */
            return -1;
        }
    }
    *out = rules;
    return n;
}

void regex_rules_free(RegexRule* rules, int n) {
    if (!rules) return;
    for (int i = 0; i < n; i++) {
        regex_free(rules[i].prog);
        sp_info_free(&rules[i].info);
    }
    free(rules);
}

/* ------------------------------------------------------------------ */
/* Shared match scanner                                               */
/* ------------------------------------------------------------------ */

/*
 * Accumulator for one scan: the growable span array plus the optional flat
 * capture pool.  Spans index the pool by offset rather than pointer so that
 * growing the pool never invalidates already-recorded spans.
 */
typedef struct {
    RegexSpan* spans;
    size_t     count, cap;
    size_t*    caps;
    size_t     caps_len, caps_cap;
    int        want_captures;
} ScanState;

/*
 * Record one match.  `ov` holds `pairs` (start, end) offset pairs as written by
 * regex_match, relative to `base` (nonzero only in the All mode, which matches
 * against subj + p); they are rebased so the pool is always subject-relative.
 * Returns 0, or -1 on allocation failure.
 */
static int scan_push(ScanState* st, size_t ms, size_t me, int rule,
                     const size_t* ov, size_t pairs, size_t base) {
    if (st->count == st->cap) {
        size_t nc = st->cap ? st->cap * 2 : 16;
        RegexSpan* ns = realloc(st->spans, nc * sizeof(RegexSpan));
        if (!ns) return -1;
        st->spans = ns;
        st->cap = nc;
    }

    size_t off = 0, np = 0;
    if (st->want_captures && ov && pairs) {
        np = pairs;
        if (st->caps_len + 2 * np > st->caps_cap) {
            size_t nc = st->caps_cap ? st->caps_cap : 32;
            while (nc < st->caps_len + 2 * np) nc *= 2;
            size_t* ncaps = realloc(st->caps, nc * sizeof(size_t));
            if (!ncaps) return -1;
            st->caps = ncaps;
            st->caps_cap = nc;
        }
        off = st->caps_len;
        for (size_t k = 0; k < 2 * np; k++)
            st->caps[off + k] =
                (ov[k] == REGEX_UNSET) ? REGEX_UNSET : ov[k] + base;
        st->caps_len += 2 * np;
    }

    st->spans[st->count].ms       = ms;
    st->spans[st->count].me       = me;
    st->spans[st->count].rule     = rule;
    st->spans[st->count].caps_off = off;
    st->spans[st->count].npairs   = np;
    st->count++;
    return 0;
}

/*
 * Offset pairs to request per rule: the whole match plus every capture group,
 * clamped to the range we expose ($0..$63).  Computed once per scan -- All-mode
 * probes O(len^2 * nr) substrings, so querying the group count inside that loop
 * would be pure overhead.  Returns a malloc'd array, or NULL on failure.
 */
static size_t* scan_pairs_table(RegexRule* rules, int nr) {
    size_t* pairs = malloc((size_t)nr * sizeof(size_t));
    if (!pairs) return NULL;
    for (int i = 0; i < nr; i++) {
        int gc = regex_group_count(rules[i].prog) + 1;
        pairs[i] = (gc > REGEX_MAX_PAIRS) ? (size_t)REGEX_MAX_PAIRS : (size_t)gc;
    }
    return pairs;
}

/*
 * Overlaps -> False.  Streaming left-to-right: at each position take the match
 * with the smallest start across all rules (ties broken by rule order), then
 * resume at its end.  Zero-width matches advance by one so the loop always
 * makes progress.
 */
static int scan_false(const char* subj, size_t len, RegexRule* rules, int nr,
                      const size_t* pairs, ScanState* st) {
    size_t pos = 0;
    while (pos <= len) {
        int found = 0, best_rule = 0;
        size_t best_ms = 0, best_me = 0, best_pairs = 0;
        size_t best_ov[REGEX_MAX_PAIRS * 2];

        for (int i = 0; i < nr; i++) {
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs[i]) != 1)
                continue;
            if (!found || ov[0] < best_ms) {
                found = 1;
                best_rule = i;
                best_ms = ov[0];
                best_me = ov[1];
                best_pairs = pairs[i];
                memcpy(best_ov, ov, sizeof(size_t) * 2 * pairs[i]);
            }
        }
        if (!found) break;

        if (scan_push(st, best_ms, best_me, best_rule, best_ov, best_pairs, 0) != 0)
            return -1;
        pos = (best_me > best_ms) ? best_me : best_ms + 1;   /* progress */
    }
    return 0;
}

/*
 * Overlaps -> True.  Per rule independently, the leftmost match at or after each
 * position, advancing by one past the match START so every distinct match start
 * is enumerated.
 */
static int scan_true(const char* subj, size_t len, RegexRule* rules, int nr,
                     const size_t* pairs, ScanState* st) {
    for (int i = 0; i < nr; i++) {
        size_t pos = 0;
        while (pos <= len) {
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs[i]) != 1) break;
            if (scan_push(st, ov[0], ov[1], i, ov, pairs[i], 0) != 0) return -1;
            pos = ov[0] + 1;   /* past the match start (overlap-friendly) */
        }
    }
    return 0;
}

/*
 * Overlaps -> All.  Every matching substring at every start: for each start p
 * ascending and end e descending, test an exact match of subj[p, e).  The rules
 * are anchored (\A(?:...)\z) by the caller, so a match means the pattern covers
 * [p, e) exactly.
 */
static int scan_all(const char* subj, size_t len, RegexRule* rules, int nr,
                    const size_t* pairs, ScanState* st) {
    for (size_t p = 0; p < len; p++) {
        for (size_t e = len; e > p; e--) {
            for (int i = 0; i < nr; i++) {
                size_t ov[REGEX_MAX_PAIRS * 2];
                if (regex_match(rules[i].prog, subj + p, e - p, 0, ov, pairs[i]) != 1)
                    continue;
                if (scan_push(st, p, e, i, ov, pairs[i], p) != 0) return -1;
            }
        }
    }
    return 0;
}

/* Stable insertion sort by start offset (ties keep discovery order: rule index,
 * then match order). qsort is not stable, and the input is nearly sorted -- it
 * is a concatenation of per-rule ascending runs -- so this is close to linear. */
static void scan_stable_sort(RegexSpan* a, size_t n) {
    for (size_t i = 1; i < n; i++) {
        RegexSpan key = a[i];
        size_t j = i;
        while (j > 0 && a[j - 1].ms > key.ms) { a[j] = a[j - 1]; j--; }
        a[j] = key;
    }
}

long regex_scan(const char* subj, size_t len, RegexRule* rules, int nr,
                RegexOverlapMode mode, int want_captures, RegexScan* out) {
    ScanState st;
    memset(&st, 0, sizeof st);
    st.want_captures = want_captures;

    out->spans = NULL;
    out->count = 0;
    out->caps  = NULL;

    size_t* pairs = scan_pairs_table(rules, nr);
    if (!pairs) return -1;

    int rc;
    if (mode == REGEX_OV_ALL)        rc = scan_all(subj, len, rules, nr, pairs, &st);
    else if (mode == REGEX_OV_FALSE) rc = scan_false(subj, len, rules, nr, pairs, &st);
    else                             rc = scan_true(subj, len, rules, nr, pairs, &st);
    free(pairs);

    if (rc != 0) {                  /* OOM mid-scan: *out stays the zeroed state */
        free(st.spans);
        free(st.caps);
        return -1;
    }

    /* The False and All scans already emit in ascending start order; only the
     * True scan concatenates one run per rule and needs merging. */
    if (mode == REGEX_OV_TRUE && nr > 1) scan_stable_sort(st.spans, st.count);

    out->spans = st.spans;
    out->count = st.count;
    out->caps  = st.caps;
    return (long)st.count;
}

void regex_scan_free(RegexScan* s) {
    if (!s) return;
    free(s->spans);
    free(s->caps);
    s->spans = NULL;
    s->caps  = NULL;
    s->count = 0;
}

int regex_match_opt(const Expr* e, const char* opt_sym, int* value, int overlaps) {
    if (e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL ||
        (e->data.function.head->data.symbol.name != SYM_Rule &&
         e->data.function.head->data.symbol.name != SYM_RuleDelayed) ||
        e->data.function.arg_count != 2 ||
        e->data.function.args[0]->type != EXPR_SYMBOL ||
        e->data.function.args[0]->data.symbol.name != opt_sym)
        return 0;

    Expr* v = e->data.function.args[1];
    if (overlaps) {
        if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_All)
            *value = REGEX_OV_ALL;
        else if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False)
            *value = REGEX_OV_FALSE;
        else
            *value = REGEX_OV_TRUE;
    } else {
        *value = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Public: $n template expansion                                      */
/* ------------------------------------------------------------------ */

char* regex_expand_template(const char* tpl, const char* subj,
                            const size_t* ov, size_t npairs) {
    RegexBuf b = {0};
    for (size_t i = 0; tpl[i];) {
        if (tpl[i] == '$') {
            char c = tpl[i + 1];
            if (c == '$') {                          /* $$ -> literal $ */
                if (regexbuf_add(&b, "$", 1)) goto oom;
                i += 2;
                continue;
            }
            if (c >= '0' && c <= '9') {              /* $n -> group n */
                size_t k = 0, j = i + 1;
                while (tpl[j] >= '0' && tpl[j] <= '9') {
                    k = k * 10 + (size_t)(tpl[j] - '0');
                    j++;
                }
                if (k < npairs && ov[2 * k] != REGEX_UNSET) {
                    size_t s = ov[2 * k], e = ov[2 * k + 1];
                    if (regexbuf_add(&b, subj + s, e - s)) goto oom;
                }
                i = j;
                continue;
            }
            if (regexbuf_add(&b, "$", 1)) goto oom;   /* lone $ */
            i += 1;
            continue;
        }
        if (regexbuf_add(&b, &tpl[i], 1)) goto oom;
        i += 1;
    }
    if (!b.p) return rc_strdup("");                   /* empty template */
    return b.p;
oom:
    free(b.p);
    return NULL;
}

Expr* regex_rule_eval_rhs(const RegexRule* r, const char* subj,
                          const size_t* ov, size_t npairs) {
    if (!r->rhs) return NULL;

    if (r->rhs->type == EXPR_STRING) {
        if (!r->info.has_regex) return expr_copy(r->rhs);   /* literal text */
        char* t = regex_expand_template(r->rhs->data.string, subj, ov, npairs);
        if (!t) return expr_copy(r->rhs);
        Expr* e = expr_new_string(t);
        free(t);
        return e;
    }

    /* No named captures: the RHS is constant per match; just evaluate it. */
    if (r->info.nnamed == 0) return evaluate(r->rhs);

    MatchEnv* env = rule_bind_env(r, subj, ov, npairs, /*unset_empty=*/1);
    Expr* inst = replace_bindings(r->rhs, env);
    env_free(env);
    return eval_and_free(inst);
}
