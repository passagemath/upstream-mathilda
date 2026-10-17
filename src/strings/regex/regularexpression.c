/*
 * regularexpression.c - the RegularExpression[...] head.
 *
 * RegularExpression["re"] is an inert data head: it carries a PCRE pattern
 * string for use by StringMatchQ / StringCases / StringReplace / StringSplit
 * and does not evaluate to anything else (exactly like Wolfram Language).
 * The builtin therefore validates its argument and otherwise returns NULL so
 * the expression survives evaluation unchanged.  When the pattern does not
 * compile a RegularExpression::regex diagnostic is emitted (still inert), and
 * when Mathilda was built without PCRE2 a RegularExpression::regavail note is
 * printed the first time an invalid-but-present call is seen.
 */

#include "picostrings.h"
#include "regex_engine.h"
#include "regex_common.h"

#include <stdio.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */

Expr* builtin_regularexpression(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    /* Only the single-string form is meaningful; anything else stays inert. */
    if (res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    /* Best-effort syntax check: warn (once per evaluation) on a bad pattern,
     * but still leave the head inert so downstream code can report context.
     * The check goes through the shared rule builder, so a valid pattern is
     * compiled once and left in its cache -- under exactly the key the string
     * functions then look up, so StringCases[s, RegularExpression[re]] does
     * not compile `re` a second time. */
    if (regex_available()) {
        RegexRule* rules;
        int nr = regex_rules_build_opts(res, 0, 0, &rules, "RegularExpression");
        if (nr > 0) regex_rules_free(rules, nr);
    }

    return NULL;   /* inert: RegularExpression["re"] stays as itself */
}
