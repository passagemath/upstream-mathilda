/*
 * stringrepeat.c - StringRepeat builtin for Mathilda
 *
 * StringRepeat["str", n]        - "str" concatenated n times.
 * StringRepeat["str", n, max]   - up to n copies of "str", truncated so the
 *                                 total length is at most max (a partial final
 *                                 copy is allowed).
 *
 * Lengths count characters (Unicode code points): StringRepeat["a\:00e9", 2, 3]
 * is "a\:00e9a". Invalid arguments emit Mathematica's StringRepeat::string /
 * StringRepeat::intp and leave the call unevaluated.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"
#include "common.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/*
 * builtin_stringrepeat:
 * See file header for the accepted forms. Returns NULL (unevaluated) for a
 * non-string first argument, a non-integer or negative count/max, or when the
 * requested length would overflow. A wrong argument count emits
 * StringRepeat::argt (via builtin_arg_error) and leaves the call unevaluated.
 */
Expr* builtin_stringrepeat(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 2 || argc > 3)
        return builtin_arg_error("StringRepeat", argc, 2, 3);

    Expr* arg0 = args[0];
    if (arg0->type != EXPR_STRING)
        return strings_arg_message("StringRepeat", "string", "String", 1, res);

    /* Count n: non-negative integer. */
    if (args[1]->type != EXPR_INTEGER || args[1]->data.integer < 0)
        return strings_arg_message("StringRepeat", "intp", "Positive integer", 2, res);
    int64_t n = args[1]->data.integer;

    /* Optional maximum total length (characters): non-negative integer. */
    bool have_max = false;
    int64_t max = 0;
    if (argc == 3) {
        if (args[2]->type != EXPR_INTEGER || args[2]->data.integer < 0)
            return strings_arg_message("StringRepeat", "intp", "Positive integer", 3, res);
        max = args[2]->data.integer;
        have_max = true;
    }

    const char* str = arg0->data.string;
    size_t len = strlen(str);
    if (n == 0 || len == 0 || (have_max && max == 0)) return expr_new_string("");

    /* Whole copies that fit, then (when capped by max) a partial copy of the
     * first `rem` characters. Lengths count characters; copies are memcpy'd
     * whole, so only the cut needs a character offset. */
    size_t clen = u8_length(str, len);
    size_t copies = (size_t)n, rem_bytes = 0;
    if (have_max && (uint64_t)max / clen < (uint64_t)n) {
        copies = (size_t)((uint64_t)max / clen);
        rem_bytes = u8_offset(str, len, (size_t)((uint64_t)max % clen));
    }
    if (copies > (SIZE_MAX - rem_bytes - 1) / len) return NULL;   /* overflow */
    size_t out_len = copies * len + rem_bytes;

    char* buf = malloc(out_len + 1);
    if (!buf) return NULL;
    for (size_t i = 0; i < copies; i++)
        memcpy(buf + i * len, str, len);
    memcpy(buf + copies * len, str, rem_bytes);
    buf[out_len] = '\0';
    return expr_new_string_adopt(buf);
}
