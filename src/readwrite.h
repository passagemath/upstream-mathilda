#ifndef READWRITE_H
#define READWRITE_H

#include "expr.h"

/* File I/O builtins.
 *
 *   Get["file"]                           reads and evaluates every
 *                                         expression in "file", returning
 *                                         the value of the last one.
 *   Get["ctx`"]                           the same for the file that holds
 *                                         package ctx` (found along $Path).
 *
 *   Put[expr_1, ..., expr_n, "file"]      writes the InputForm of each
 *                                         expr_i to "file", one per line,
 *                                         replacing any existing content.
 *   Put["file"]                           creates an empty "file".
 *   expr >> "file"                        parser shorthand for
 *                                         Put[expr, "file"].
 *
 *   PutAppend[expr_1, ..., expr_n, "file"]
 *                                         like Put but appends to the
 *                                         existing contents of "file".
 *   expr >>> "file"                       parser shorthand for
 *                                         PutAppend[expr, "file"].
 *
 * All three builtins return Null on success and $Failed on I/O error.
 */

void readwrite_init(void);

Expr* builtin_get(Expr* res);

/* The file Get[spec] would read, as a malloc'd path, or NULL.  `spec` is a
 * file name (tried as given, then -- when relative -- under each directory of
 * $Path) or a context name "A`B`", tried along $Path as A/B.m, A/B.wl,
 * A/B/init.m, A/B/init.wl, A/B/Kernel/init.m, A/B/Kernel/init.wl. */
char* readwrite_resolve(const char* spec);

/* Get[spec] for C callers (Needs): reads and evaluates the resolved file and
 * returns the last value (owned), or emits Get::noopen and returns $Failed. */
Expr* readwrite_get(const char* spec);
Expr* builtin_put(Expr* res);
Expr* builtin_putappend(Expr* res);

#endif /* READWRITE_H */
