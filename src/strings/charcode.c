/*
 * charcode.c - character codes and character ranges.
 *
 *   ToCharacterCode[s]           list of the Unicode code points of s
 *   ToCharacterCode[{s1, ...}]   one list per string
 *   FromCharacterCode[n]         the one-character string for code point n
 *   FromCharacterCode[{n1, ...}] the string of those code points
 *   FromCharacterCode[{{..}, ..}] one string per sublist
 *   CharacterRange["a", "e"]     {"a", "b", "c", "d", "e"}
 *   CharacterRange[n1, n2]       the characters with codes n1 .. n2
 *   Alphabet[]                   the 26 lowercase English letters
 *
 * Array surfaces. A string's codes are a machine-integer vector, so
 * ToCharacterCode writes them straight into a packed int64 buffer (above the
 * packing threshold) instead of building one Integer node per character, and
 * FromCharacterCode -- on the packed-aware list in src/pack.c -- reads a packed
 * or NDArray int64 buffer directly instead of having the gate box every code.
 * Neither has a Compile[] lowering: one side of each is a String, which no
 * compiled register holds (recorded in tools/check_compile_coverage.py).
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"
#include "sym_names.h"
#include "message.h"
#include "print.h"
#include "pack.h"
#include "ndarray.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/* ToCharacterCode                                                    */
/* ------------------------------------------------------------------ */

/* Codes of one string: a packed int64 vector when packing applies, else a
 * List of Integers. */
static Expr* codes_of(const char* s) {
    size_t len = strlen(s);
    size_t n = u8_length(s, len);
    int64_t* buf = NULL;
    Expr* packed = ndbuild_open_i64((int64_t)n, &buf);
    if (packed) {
        size_t k = 0;
        for (size_t i = 0; i < len; ) buf[k++] = (int64_t)u8_decode(s, len, &i);
        return packed;
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
    if (!out) return NULL;
    size_t k = 0;
    for (size_t i = 0; i < len; )
        out->data.function.args[k++] = expr_new_integer((int64_t)u8_decode(s, len, &i));
    return out;
}

Expr* builtin_tocharactercode(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];
    if (a->type == EXPR_STRING) return codes_of(a->data.string);

    int ok = a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
             a->data.function.head->data.symbol.name == SYM_List;
    for (size_t i = 0; ok && i < a->data.function.arg_count; i++)
        if (a->data.function.args[i]->type != EXPR_STRING) ok = 0;
    if (!ok)
        return strings_arg_message("ToCharacterCode", "strse",
                                   "String or list of strings", 1, res);

    size_t m = a->data.function.arg_count;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, m);
    for (size_t i = 0; out && i < m; i++) {
        out->data.function.args[i] = codes_of(a->data.function.args[i]->data.string);
        if (!out->data.function.args[i]) { expr_free(out); out = NULL; }
    }
    return out;
}

/* ------------------------------------------------------------------ */
/* FromCharacterCode                                                  */
/* ------------------------------------------------------------------ */

/* A code is valid in 0 .. 0x10FFFF, surrogates excluded. */
static int code_ok(int64_t c) {
    return c >= 0 && c <= 0x10FFFF && !(c >= 0xD800 && c <= 0xDFFF);
}

/* Report an invalid code: a negative or non-integer one is ::intnm on the
 * call; an out-of-range one is ::notunicode, naming the code list. */
static Expr* fcc_bad(Expr* res, const Expr* code, size_t pos, const Expr* list) {
    if (code->type == EXPR_INTEGER && code->data.integer >= 0) {
        char* shown = expr_to_string_outputform((Expr*)(list ? list : code));
        mth_message("FromCharacterCode", "notunicode",
                    "A character code, which should be a non-negative integer "
                    "less than 1114112, is expected at position %zu in %s%s%s.",
                    pos, list ? "" : "{", shown ? shown : "", list ? "" : "}");
        free(shown);
        return NULL;
    }
    return strings_arg_message("FromCharacterCode", "intnm",
                               "Non-negative machine-sized integer", pos, res);
}

/* The string of the n codes read by get(k); NULL (message emitted) if one is
 * invalid. `res`/`list` are only used for diagnostics. */
static Expr* string_of_codes(const int64_t* codes, size_t n, Expr* res, const Expr* list) {
    char* buf = malloc(4 * n + 1);
    if (!buf) return NULL;
    size_t o = 0;
    for (size_t k = 0; k < n; k++) {
        if (!code_ok(codes[k])) {
            free(buf);
            Expr* c = expr_new_integer(codes[k]);
            fcc_bad(res, c, k + 1, list);
            expr_free(c);
            return NULL;
        }
        o += u8_encode((uint32_t)codes[k], buf + o);
    }
    buf[o] = '\0';
    return expr_new_string_adopt(buf);
}

/* A List of Integer codes -> String; NULL if not all are valid codes. */
static Expr* string_of_list(const Expr* lst, Expr* res) {
    size_t n = lst->data.function.arg_count;
    int64_t* codes = malloc((n ? n : 1) * sizeof(int64_t));
    if (!codes) return NULL;
    for (size_t k = 0; k < n; k++) {
        const Expr* c = lst->data.function.args[k];
        if (c->type != EXPR_INTEGER || !code_ok(c->data.integer)) {
            free(codes);
            return fcc_bad(res, c, k + 1, lst);
        }
        codes[k] = c->data.integer;
    }
    Expr* s = string_of_codes(codes, n, res, lst);
    free(codes);
    return s;
}

static int is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* FromCharacterCode of the value `a`; `res` is the call, for diagnostics. */
static Expr* fcc_value(Expr* a, Expr* res) {
    if (a->type == EXPR_INTEGER) {
        if (!code_ok(a->data.integer)) return fcc_bad(res, a, 1, NULL);
        char buf[5];
        size_t n = u8_encode((uint32_t)a->data.integer, buf);
        return strings_char_expr(buf, n);
    }

    /* An int64 buffer (packed List or visible NDArray): read it directly. */
    if (a->type == EXPR_NDARRAY && a->data.ndarray.dtype == NDT_INT64 &&
        (a->data.ndarray.rank == 1 || a->data.ndarray.rank == 2)) {
        const NDArrayData* nd = &a->data.ndarray;
        const int64_t* codes = (const int64_t*)nd->data;
        if (nd->rank == 1)
            return string_of_codes(codes, (size_t)nd->dims[0], res, a);
        size_t rows = (size_t)nd->dims[0], cols = (size_t)nd->dims[1];
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, rows);
        for (size_t r = 0; out && r < rows; r++) {
            out->data.function.args[r] = string_of_codes(codes + r * cols, cols, res, a);
            if (!out->data.function.args[r]) { expr_free(out); out = NULL; }
        }
        return out;
    }
    if (a->type == EXPR_NDARRAY) {                /* other dtypes: as a List */
        Expr* plain = pack_unpack(a);
        if (!plain) return NULL;
        Expr* r = fcc_value(plain, res);
        expr_free(plain);
        return r;
    }

    if (!is_list(a)) return fcc_bad(res, a, 1, NULL);

    /* {{codes}, {codes}, ...}: one string per sublist. */
    size_t m = a->data.function.arg_count;
    if (m > 0 && is_list(a->data.function.args[0])) {
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, m);
        for (size_t i = 0; out && i < m; i++) {
            const Expr* sub = a->data.function.args[i];
            out->data.function.args[i] = is_list(sub) ? string_of_list(sub, res)
                                                      : fcc_bad(res, sub, i + 1, NULL);
            if (!out->data.function.args[i]) { expr_free(out); out = NULL; }
        }
        return out;
    }
    return string_of_list(a, res);
}

Expr* builtin_fromcharactercode(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return fcc_value(res->data.function.args[0], res);
}

/* ------------------------------------------------------------------ */
/* CharacterRange / Alphabet                                          */
/* ------------------------------------------------------------------ */

/* The code point of a one-character string or a non-negative integer, or -1. */
static long range_end(const Expr* e) {
    if (e->type == EXPR_INTEGER && code_ok(e->data.integer)) return (long)e->data.integer;
    if (e->type == EXPR_STRING) {
        const char* s = e->data.string;
        size_t len = strlen(s), i = 0;
        if (len == 0) return -1;
        uint32_t cp = u8_decode(s, len, &i);
        return i == len ? (long)cp : -1;
    }
    return -1;
}

static Expr* char_list(uint32_t lo, uint32_t hi) {
    size_t n = hi >= lo ? (size_t)(hi - lo + 1) : 0;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
    for (size_t k = 0; out && k < n; k++) {
        char buf[5];
        uint32_t cp = lo + (uint32_t)k;
        size_t len = u8_encode(cp, buf);
        if (len == 0) { buf[0] = '?'; len = 1; }          /* surrogate gap */
        out->data.function.args[k] = strings_char_expr(buf, len);
    }
    return out;
}

Expr* builtin_characterrange(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    long lo = range_end(a), hi = range_end(b);
    /* Both ends must be the same kind: two characters or two codes. */
    if (lo < 0 || hi < 0 || a->type != b->type) {
        char* sa = expr_to_string_outputform(a);
        char* sb = expr_to_string_outputform(b);
        mth_message("CharacterRange", "argtype",
                    "Arguments %s and %s of CharacterRange should be either "
                    "non-negative integers or one-character strings.",
                    sa ? sa : "", sb ? sb : "");
        free(sa);
        free(sb);
        return NULL;
    }
    return char_list((uint32_t)lo, (uint32_t)hi);
}

Expr* builtin_alphabet(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) return NULL;
    return char_list('a', 'z');
}
