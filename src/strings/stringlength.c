/*
 * stringlength.c - StringLength builtin for Mathilda
 *
 * StringLength["string"] - gives the number of characters in a string
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"
#include <string.h>

/*
 * builtin_stringlength:
 * StringLength["string"] returns the number of characters (Unicode code
 * points) in the UTF-8 string -- StringLength["h\:00e9llo"] is 5 -- skipping
 * ASCII runs eight bytes at a time. Listable (threads over lists). A
 * non-string argument emits StringLength::string and stays unevaluated.
 */
Expr* builtin_stringlength(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return strings_arg_message("StringLength", "string", "String", 1, res);

    const char* s = arg->data.string;
    return expr_new_integer((int64_t)u8_length(s, strlen(s)));
}
