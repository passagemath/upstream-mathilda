/*
 * strings_common.c - helpers shared by the String* builtins.
 */

#include "strings_internal.h"
#include "message.h"
#include "print.h"

#include <stdlib.h>

/*
 * Emit Mathematica's standard argument diagnostic
 *     head::tag: <what> expected at position <pos> in <call>.
 * with the call rendered in OutputForm (strings unquoted, as Mathematica
 * prints it), through the message funnel so Quiet / Check see it. Always
 * returns NULL, so a builtin can `return strings_arg_message(...)` to stay
 * unevaluated.
 */
Expr* strings_arg_message(const char* head, const char* tag, const char* what,
                          size_t pos, Expr* call) {
    char* shown = expr_to_string_outputform(call);
    mth_message(head, tag, "%s expected at position %zu in %s.", what, pos,
                shown ? shown : head);
    free(shown);
    return NULL;
}
