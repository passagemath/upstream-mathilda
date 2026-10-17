/*
 * stringreverse.c - StringReverse builtin for Mathilda
 *
 * StringReverse["string"] - reverses the order of the characters in a string
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"
#include <string.h>
#include <stdio.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */

/*
 * sr_emit_argx:
 * Emit `StringReverse::argx: StringReverse called with N arguments; 1
 * argument is expected.` to stderr and return NULL so the call is left
 * unevaluated. Matches Mathematica's diagnostic for a wrong argument count.
 */
static Expr* sr_emit_argx(size_t argc) {
    mth_message("StringReverse", "argx", "StringReverse called with %zu argument%s; 1 argument is expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/*
 * builtin_stringreverse:
 * StringReverse["string"] reverses the order of the characters in "string",
 * returning a new string. StringReverse is Listable, so the evaluator threads
 * it element-wise over a list argument before this builtin is reached; each
 * threaded call therefore sees a single string.
 *
 * A call with a number of arguments other than one emits StringReverse::argx
 * and is left unevaluated. A single non-string argument also leaves the call
 * unevaluated (NULL) so symbolic arguments flow through unchanged.
 */
Expr* builtin_stringreverse(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 1)
        return sr_emit_argx(argc);

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return strings_arg_message("StringReverse", "string", "String", 1, res);

    const char* str = arg->data.string;
    size_t len = strlen(str);

    char* buf = malloc(len + 1);
    if (!buf)
        return NULL;

    /* Reverse the order of the CHARACTERS: each UTF-8 sequence is copied
     * whole to its mirrored position (a byte reversal would corrupt them). */
    if (u8_is_ascii(str, len)) {
        for (size_t i = 0; i < len; i++)
            buf[i] = str[len - 1 - i];
    } else {
        for (size_t i = 0; i < len; ) {
            size_t n = u8_char_len(str, len, i);
            memcpy(buf + len - i - n, str + i, n);
            i += n;
        }
    }
    buf[len] = '\0';
    return expr_new_string_adopt(buf);
}
