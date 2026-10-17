/*
 * stringinsert.c - StringInsert builtin for Mathilda
 *
 * StringInsert["string", "snew", spec] - inserts a string at given positions.
 * Positions count characters (code points), not UTF-8 bytes.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <inttypes.h>

/*
 * si_emit_psl:
 * StringInsert::psl. For a scalar spec Mathematica names the spec and the
 * whole call; for a list it names the 1-based index of the bad element and
 * the list.
 */
static Expr* si_emit_psl(Expr* call, Expr* spec, size_t index) {
    char* a = index ? NULL : expr_to_string_outputform(spec);
    char* b = expr_to_string_outputform(index ? spec : call);
    if (index)
        mth_message("StringInsert", "psl",
                    "Position specification %zu in %s is not a machine-sized "
                    "integer or a list of machine-sized integers.", index, b ? b : "");
    else
        mth_message("StringInsert", "psl",
                    "Position specification %s in %s is not a machine-sized "
                    "integer or a list of machine-sized integers.", a ? a : "", b ? b : "");
    free(a);
    free(b);
    return NULL;
}

/* int64 ascending, for qsort. */
static int si_cmp_i64(const void* x, const void* y) {
    int64_t a = *(const int64_t*)x, b = *(const int64_t*)y;
    return (a > b) - (a < b);
}

/*
 * insert_scalar:
 * StringInsert[str, ins, positions] where positions are already known to be
 * machine integers. StringInsert["string", "snew", n] makes the first
 * character of "snew" the nth character of the result, i.e. it inserts before
 * original character n (character offset n - 1); a negative -k makes the last
 * character of "snew" the kth from the end (offset len + 1 - k). Valid
 * positions are therefore 1..len+1 and -(len+1)..-1. All positions refer to
 * the ORIGINAL string. Returns NULL after StringInsert::ins for a position out
 * of range.
 */
static Expr* insert_scalar(Expr* str, const char* ins, Expr** pos, size_t npos) {
    const char* s = str->data.string;
    U8View v;
    if (u8_view_init(&v, s, strlen(s)) != 0) return NULL;
    int64_t len = (int64_t)v.nchars;

    /* Byte offsets of the insertion points, sorted. */
    int64_t* at = malloc(sizeof(int64_t) * (npos ? npos : 1));
    if (!at) { u8_view_free(&v); return NULL; }
    for (size_t i = 0; i < npos; i++) {
        int64_t n = pos[i]->data.integer, off;
        if (n > 0 && n <= len + 1)         off = n - 1;
        else if (n < 0 && n >= -(len + 1)) off = len + 1 + n;
        else {
            mth_message("StringInsert", "ins", "Cannot insert at position %" PRId64 " in %s.",
                        n, s);
            free(at);
            u8_view_free(&v);
            return NULL;
        }
        at[i] = (int64_t)STR_OFF(&v, off);
    }
    if (npos > 1) qsort(at, npos, sizeof(int64_t), si_cmp_i64);

    size_t inslen = strlen(ins);
    if (npos && inslen > (SIZE_MAX - v.nbytes - 1) / npos) {
        free(at); u8_view_free(&v); return NULL;
    }
    char* buf = malloc(v.nbytes + npos * inslen + 1);
    if (!buf) { free(at); u8_view_free(&v); return NULL; }
    size_t o = 0, prev = 0;
    for (size_t i = 0; i < npos; i++) {
        size_t b = (size_t)at[i];
        memcpy(buf + o, s + prev, b - prev);
        o += b - prev;
        memcpy(buf + o, ins, inslen);
        o += inslen;
        prev = b;
    }
    memcpy(buf + o, s + prev, v.nbytes - prev);
    o += v.nbytes - prev;
    buf[o] = '\0';
    free(at);
    u8_view_free(&v);
    return expr_new_string_adopt(buf);
}

/*
 * builtin_stringinsert:
 * StringInsert["string", "snew", n]           - insert "snew" before position n
 * StringInsert["string", "snew", -n]          - insert n positions from the end
 * StringInsert["string", "snew", {n1, n2, ...}] - insert a copy at each position
 * StringInsert[{s1, s2, ...}, "snew", spec]   - maps over the list of strings
 *
 * Messages (call left unevaluated): argr/argrx for the wrong argument count,
 * strse for a first argument that is not a string or list of strings, string
 * for a non-string "snew", psl for a position that is not a machine integer,
 * ins for one out of range. Over a list of strings an element with an
 * out-of-range position stays as the unevaluated StringInsert[si, snew, spec].
 */
Expr* builtin_stringinsert(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 3) return strings_argr_error("StringInsert", argc, 3);

    Expr* arg0 = res->data.function.args[0];
    Expr* snew = res->data.function.args[1];
    Expr* spec = res->data.function.args[2];

    bool threaded = strings_is_list(arg0);
    if (!(arg0->type == EXPR_STRING || (threaded && strings_all_strings(arg0))))
        return strings_msg_strse("StringInsert", res, 1);
    if (snew->type != EXPR_STRING)
        return strings_msg_string("StringInsert", res, 2);

    /* Gather the positions; each must be a machine integer. */
    Expr** pos;
    size_t npos;
    if (spec->type == EXPR_INTEGER) {
        pos = &spec;
        npos = 1;
    } else if (strings_is_list(spec)) {
        pos = spec->data.function.args;
        npos = spec->data.function.arg_count;
        for (size_t i = 0; i < npos; i++)
            if (pos[i]->type != EXPR_INTEGER) return si_emit_psl(res, spec, i + 1);
    } else {
        return si_emit_psl(res, spec, 0);
    }

    const char* ins = snew->data.string;
    if (!threaded) return insert_scalar(arg0, ins, pos, npos);

    size_t n = arg0->data.function.arg_count;
    Expr** results = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!results) return NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* si = arg0->data.function.args[i];
        results[i] = insert_scalar(si, ins, pos, npos);
        if (!results[i])
            results[i] = strings_uneval(SYM_StringInsert, si, snew, spec);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), results, n);
    free(results);
    return out;
}
