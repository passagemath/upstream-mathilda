/*
 * stringreplace.c - StringReplace[subject, rule | {rules...}]
 *
 * Replaces every non-overlapping match of a rule's pattern by its right-hand
 * side, scanning left to right; unmatched text is copied verbatim.  Pattern
 * names in the LHS are bound to the text they matched and substituted into the
 * RHS (regex_rule_eval_rhs); a non-string value turns the result into a
 * StringExpression.  For a RegularExpression pattern a string RHS may use
 * $0/$1... for the whole match and capture groups.
 * With several rules, at each position the leftmost match wins, ties broken by
 * rule order.  A list of subjects threads.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"

#include <string.h>

/*
 * Result assembly. Unmatched text and string replacements accumulate in `buf`;
 * a replacement that is not a string (x:DigitCharacter :> ToExpression[x])
 * cannot be spliced into text, so -- as in Mathematica -- the result becomes a
 * StringExpression of the string runs and those values:
 *   StringReplace["abc123", x:DigitCharacter :> f[x]]
 *     = "abc" ~~ f["1"] ~~ f["2"] ~~ f["3"]
 * Sequence[...] replacements splice (Sequence[] deletes).
 */
typedef struct {
    RegexBuf buf;        /* current run of text */
    Expr**   items;      /* StringExpression pieces, once any non-string appears */
    size_t   n, cap;
    int      oom;
} SrOut;

static void sr_push_item(SrOut* o, Expr* e) {
    if (o->n == o->cap) {
        size_t nc = o->cap ? o->cap * 2 : 8;
        Expr** ni = realloc(o->items, nc * sizeof(Expr*));
        if (!ni) { expr_free(e); o->oom = 1; return; }
        o->items = ni;
        o->cap = nc;
    }
    o->items[o->n++] = e;
}

/* Close the current text run as a string piece (only once pieces exist). */
static void sr_flush(SrOut* o) {
    if (o->buf.len == 0) return;
    sr_push_item(o, expr_new_string(o->buf.p));
    o->buf.len = 0;
    o->buf.p[0] = '\0';
}

static void sr_text(SrOut* o, const char* p, size_t n) {
    if (n && regexbuf_add(&o->buf, p, n)) o->oom = 1;
}

/* Append one evaluated replacement value; takes ownership of v. */
static void sr_value(SrOut* o, Expr* v) {
    if (v->type == EXPR_STRING) {
        sr_text(o, v->data.string, strlen(v->data.string));
        expr_free(v);
        return;
    }
    if (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
        v->data.function.head->data.symbol.name == SYM_Sequence) {
        for (size_t i = 0; i < v->data.function.arg_count; i++)
            sr_value(o, expr_copy(v->data.function.args[i]));
        expr_free(v);
        return;
    }
    sr_flush(o);
    sr_push_item(o, v);
}

/* Replace in one subject. Returns a new String, or StringExpression[...] when a
 * replacement value is not a string; NULL on OOM. */
static Expr* sr_scalar(const char* subj, RegexRule* rules, int nr) {
    size_t len = strlen(subj);
    SrOut o;
    memset(&o, 0, sizeof o);

    size_t pos = 0;
    while (pos <= len && !o.oom) {
        int found = 0, best_rule = 0;
        size_t best_ms = 0, best_me = 0, best_pairs = 0;
        size_t best_ov[REGEX_MAX_PAIRS * 2];

        for (int i = 0; i < nr; i++) {
            int gc = regex_group_count(rules[i].prog) + 1;
            size_t pairs = (gc > REGEX_MAX_PAIRS) ? REGEX_MAX_PAIRS : (size_t)gc;
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs) == 1) {
                if (!found || ov[0] < best_ms) {
                    found = 1;
                    best_rule = i;
                    best_ms = ov[0];
                    best_me = ov[1];
                    best_pairs = pairs;
                    memcpy(best_ov, ov, sizeof(size_t) * 2 * pairs);
                }
            }
        }
        if (!found) {
            sr_text(&o, subj + pos, len - pos);          /* tail */
            break;
        }

        sr_text(&o, subj + pos, best_ms - pos);          /* literal before match */
        sr_value(&o, regex_rule_eval_rhs(&rules[best_rule], subj,
                                         best_ov, best_pairs));

        if (best_me > best_ms) {
            pos = best_me;
        } else {                                         /* zero-width match */
            if (best_ms < len) sr_text(&o, subj + best_ms, 1);
            pos = best_ms + 1;
        }
    }

    Expr* result = NULL;
    if (!o.oom) {
        if (o.n == 0) {
            result = expr_new_string(o.buf.p ? o.buf.p : "");
        } else {
            sr_flush(&o);
            if (!o.oom) {
                result = expr_new_function(expr_new_symbol(SYM_StringExpression),
                                           o.items, o.n);
                o.n = 0;
            }
        }
    }
    for (size_t i = 0; i < o.n; i++) expr_free(o.items[i]);
    free(o.items);
    free(o.buf.p);
    return result;
}

Expr* builtin_stringreplace(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* patt = res->data.function.args[1];

    RegexRule* rules;
    int nr = regex_rules_build(patt, /*anchored=*/0, &rules, "StringReplace");
    if (nr < 0) return NULL;

    /* StringReplace requires rule(s): every element must carry an RHS. */
    for (int i = 0; i < nr; i++) {
        if (!rules[i].rhs) { regex_rules_free(rules, nr); return NULL; }
    }

    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol.name == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!out) { regex_rules_free(rules, nr); return NULL; }
        for (size_t i = 0; i < n; i++) {
            Expr* si = arg0->data.function.args[i];
            out[i] = (si->type == EXPR_STRING)
                         ? sr_scalar(si->data.string, rules, nr)
                         : expr_copy(si);
            if (!out[i]) out[i] = expr_copy(si);         /* OOM: leave as is */
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
        free(out);
        regex_rules_free(rules, nr);
        return result;
    }

    if (arg0->type != EXPR_STRING) {
        regex_rules_free(rules, nr);
        return NULL;
    }

    Expr* result = sr_scalar(arg0->data.string, rules, nr);
    regex_rules_free(rules, nr);
    return result;                           /* NULL (OOM): left unevaluated */
}
