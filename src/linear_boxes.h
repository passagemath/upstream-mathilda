#ifndef LINEAR_BOXES_H
#define LINEAR_BOXES_H

#include <stddef.h>
#include <stdint.h>

/* Linear box syntax: \!\( ... \) and \( ... \), with \* introducing a box
 * expression such as \*SuperscriptBox[\(mb\), \(2\)].  The control characters
 * are written either as a backslash pair in source (\!, \(, \), \*, ...) or as
 * the private-use code points a string literal decodes them to (U+F7C0 ..
 * U+F7CD, the reference's), so ToExpression of such a string parses the same
 * way. */

/* The private-use code point a string literal decodes the escape \c to, or 0
 * when \c is not a box control escape (c one of ! ( ) * @ ^ & % + / _ `). */
uint32_t lb_escape_cp(char c);

/* The escape spelling (e.g. "\\(") of a box control code point, or NULL. */
const char* lb_cp_escape(uint32_t cp);

/* Length of the \!\( or \( that opens a linear box group at p (either
 * spelling), or 0. */
size_t lb_group_open(const char* p);

/* Convert the linear box group opening at p to InputForm source text
 * (malloc'd; the caller frees), with *end set past its closing \).  The boxes
 * become the expressions they display: SuperscriptBox[a, b] is a^b,
 * FractionBox[a, b] is a/b, SqrtBox[a] is Sqrt[a], SubscriptBox[a, b] is
 * Subscript[a, b], RowBox concatenates.  NULL on malformed input. */
char* lb_group_to_input(const char* p, const char** end);

#endif
