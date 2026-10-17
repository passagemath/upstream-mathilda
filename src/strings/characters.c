/*
 * characters.c - Characters builtin for Mathilda
 *
 * Characters["string"] - gives a list of the individual characters in a string
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"
#include "sym_names.h"
#include <string.h>

/* One shared String per ASCII character, handed out by reference (expr_copy is
 * a refcount bump): Characters / StringSplit[s, ""] / StringPart on ASCII text
 * then allocate nothing per character. Built on first use. */
static Expr* g_ascii_chars[128];

Expr* strings_char_expr(const char* s, size_t n) {
    if (n == 1 && !((unsigned char)s[0] & 0x80) && s[0] != '\0') {
        unsigned char c = (unsigned char)s[0];
        if (!g_ascii_chars[c]) {
            char buf[2] = { (char)c, '\0' };
            g_ascii_chars[c] = expr_new_string(buf);
            if (!g_ascii_chars[c]) return NULL;
        }
        return expr_copy(g_ascii_chars[c]);
    }
    return expr_new_string_n(s, n);
}

/*
 * builtin_characters:
 * Characters["string"] returns a List of one-character strings, one per
 * Unicode code point (a multi-byte UTF-8 character stays whole).
 * Listable (threads over lists). Returns NULL (unevaluated) for non-strings.
 */
Expr* builtin_characters(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    const char* str = arg->data.string;
    size_t len = strlen(str);

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), NULL,
                                     u8_length(str, len));
    if (!result) return NULL;
    Expr** out = result->data.function.args;
    size_t k = 0;
    for (size_t i = 0; i < len; ) {
        size_t n = u8_char_len(str, len, i);
        out[k] = strings_char_expr(str + i, n);
        if (!out[k]) { expr_free(result); return NULL; }   /* NULL slots are skipped */
        k++;
        i += n;
    }
    return result;
}
