/*
 * stringpartition.c - StringPartition builtin for Mathilda
 *
 * StringPartition["string", n]        - non-overlapping length-n blocks
 * StringPartition["string", n, d]     - length-n blocks starting every d chars
 * StringPartition["string", UpTo[n]]  - length-<=n blocks; final may be shorter
 * StringPartition["string", UpTo[n], d] - the same, with offset d
 *
 * Lengths and offsets count characters (code points), not UTF-8 bytes. As in
 * Mathematica, the first argument must be a single string (no threading).
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

/* StringPartition::<tag>: <what> expected at position pos in call. */
static Expr* sp_emit_pos(const char* tag, const char* what, int pos, Expr* call) {
    char* c = expr_to_string_outputform(call);
    mth_message("StringPartition", tag, "%s expected at position %d in %s.",
                what, pos, c ? c : "");
    free(c);
    return NULL;
}

/*
 * builtin_stringpartition:
 * See file header for the accepted forms. Messages (call left unevaluated):
 * argtu/argt for the wrong argument count, string for a non-string first
 * argument, intpm for a block length that is not a positive machine integer
 * (or UpTo of one) -- intnm, "non-negative", when an offset is given, since
 * n == 0 is then allowed -- intnm for UpTo[Infinity], intp for a
 * non-positive or non-integer offset. UpTo[negative] is left unevaluated silently: in
 * Mathematica the only message there is UpTo's own.
 */
Expr* builtin_stringpartition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 2 || argc > 3) return strings_argt_error("StringPartition", argc);
    if (args[0]->type != EXPR_STRING)
        return strings_msg_string("StringPartition", res, 1);

    /* Block length n, and whether a short final block is allowed. With an
     * explicit offset Mathematica accepts n == 0 (and UpTo[0]): every block
     * is then "", and the messages say "non-negative" instead of "positive". */
    bool has_d = argc == 3;
    const char* bad_tag = has_d ? "intnm" : "intpm";
    const char* bad_what = has_d ? "Non-negative machine-sized integer"
                                 : "Positive machine-sized integer";
    int64_t n;
    bool allow_short = false;
    if (args[1]->type == EXPR_INTEGER) {
        n = args[1]->data.integer;
        if (n < 0 || (n == 0 && !has_d))
            return sp_emit_pos(bad_tag, bad_what, 2, res);
    } else if (strings_has_head(args[1], SYM_UpTo) &&
               args[1]->data.function.arg_count == 1) {
        Expr* k = args[1]->data.function.args[0];
        if (k->type == EXPR_INTEGER && k->data.integer < 0) return NULL;
        if (strings_is_infinity(k))
            return sp_emit_pos("intnm", "Non-negative machine-sized integer", 2, res);
        if (k->type != EXPR_INTEGER || (k->data.integer == 0 && !has_d))
            return sp_emit_pos("intpm", "Positive machine-sized integer", 2, res);
        n = k->data.integer;
        allow_short = true;
    } else {
        return sp_emit_pos(bad_tag, bad_what, 2, res);
    }

    /* Offset d defaults to n; an explicit third argument overrides it. */
    int64_t d = n;
    if (has_d) {
        if (args[2]->type != EXPR_INTEGER || args[2]->data.integer <= 0)
            return sp_emit_pos("intp", "Positive integer", 3, res);
        d = args[2]->data.integer;
    }

    const char* s = args[0]->data.string;
    U8View v;
    if (u8_view_init(&v, s, strlen(s)) != 0) return NULL;
    int64_t len = (int64_t)v.nchars;

    /* At most one block per offset step. */
    size_t max_blocks = (size_t)(len / d) + 1;
    Expr** blocks = malloc(sizeof(Expr*) * max_blocks);
    if (!blocks) { u8_view_free(&v); return NULL; }

    /* Blocks start at 0, d, 2d, ... (characters). A fixed length n needs a
     * whole block (start + n <= len); UpTo[n] takes every start before the
     * end, each block clipped to what remains (so with d < n the tail
     * repeats shorter and shorter blocks, as in Mathematica). */
    size_t count = 0;
    for (int64_t start = 0; ; start += d) {
        int64_t avail = len - start;
        if (allow_short ? avail <= 0 : n > avail) break;   /* no start+n overflow */
        int64_t take = n <= avail ? n : avail;
        size_t b0 = STR_OFF(&v, start);
        blocks[count++] = expr_new_string_n(s + b0, STR_OFF(&v, start + take) - b0);
        if (d > avail) break;                            /* next start is past len */
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), blocks, count);
    free(blocks);
    u8_view_free(&v);
    return result;
}
