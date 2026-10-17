/*
 * stringtakedrop.c - StringTakeDrop["string", spec].
 *
 * Gives {StringTake["string", spec], StringDrop["string", spec]} for one
 * standard sequence specification (n, -n, {n}, {m, n}, {m, n, s}, UpTo[n],
 * All, None), sharing StringTake's / StringDrop's scalar kernels. A list of
 * strings gives one pair per string. Positions count characters.
 *
 * Diagnostics are StringTakeDrop's own (Mathematica): ::strse for a subject
 * that is not a string or list of strings, ::seqs for an invalid spec, ::take
 * for out-of-range positions -- all leaving the call unevaluated.
 */

#include "picostrings.h"
#include "strings_internal.h"

#include <inttypes.h>

static Expr* takedrop_pair(Expr* str, const StrSeq* q) {
    Expr* t = strings_take_seq(str, q);
    Expr* d = t ? strings_drop_seq(str, q) : NULL;
    if (!t || !d) { if (t) expr_free(t); return NULL; }
    Expr* args[2] = { t, d };
    return expr_new_function(expr_new_symbol(SYM_List), args, 2);
}

Expr* builtin_stringtakedrop(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return strings_argr_error("StringTakeDrop", argc, 2);
    Expr* subj = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    if (subj->type != EXPR_STRING &&
        !(strings_is_list(subj) && strings_all_strings(subj)))
        return strings_msg_strse("StringTakeDrop", res, 1);

    StrSeq q;
    if (!strings_seq_parse(spec, &q)) {
        char* c = expr_to_string_outputform(res);
        mth_message("StringTakeDrop", "seqs",
                    "Sequence specification (+n, -n, {+n}, {-n}, {m, n}, or {m, n, s}) "
                    "expected at position 2 in %s.", c ? c : "");
        free(c);
        return NULL;
    }

    /* Validate the positions against every string first, so an out-of-range
     * spec reports StringTakeDrop::take (not the kernels' own heads). */
    size_t n = subj->type == EXPR_STRING ? 1 : subj->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        Expr* si = subj->type == EXPR_STRING ? subj : subj->data.function.args[i];
        StrChars c;
        strings_chars_init(&c, si->data.string, strings_seq_prefix(&q));
        int64_t first, step, count;
        if (!strings_seq_resolve(&q, c.len, &first, &step, &count)) {
            char* s = expr_to_string(si);
            mth_message("StringTakeDrop", "take",
                        "Cannot take positions %" PRId64 " through %" PRId64 " in %s.",
                        q.m, q.n, s ? s : "");
            free(s);
            return NULL;
        }
    }

    if (subj->type == EXPR_STRING) return takedrop_pair(subj, &q);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
    for (size_t i = 0; out && i < n; i++) {
        out->data.function.args[i] = takedrop_pair(subj->data.function.args[i], &q);
        if (!out->data.function.args[i]) { expr_free(out); out = NULL; }
    }
    return out;
}
