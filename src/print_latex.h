#ifndef MATHILDA_PRINT_LATEX_H
#define MATHILDA_PRINT_LATEX_H

#include "expr.h"

/*
 * expr_to_latex(e)
 *
 * Convert an expression to a KaTeX-compatible LaTeX string in TraditionalForm
 * (textbook) style: fractions, roots, superscripts, Greek letters, trig
 * functions, etc. Equivalent to expr_to_latex_form(e, LATEX_FORM_TRADITIONAL).
 *
 * Returns a heap-allocated string — caller must free().
 * Returns NULL only on allocation failure.
 *
 * Returns the EMPTY string for an expression carrying a printer directive
 * (InputForm, FullForm, TeXForm, NumberForm anywhere in the tree): those ask
 * for a notation typesetting cannot express, and the plain-text form from
 * expr_to_string is already the answer. Every caller treats an empty result as
 * "no LaTeX" and falls back to that text, so a directive is honoured rather
 * than overridden by the typeset StandardForm the user asked not to see.
 *
 * Examples:
 *   Times[Pi, Power[E, -2]]  →  "\frac{\pi}{e^{2}}"
 *   Power[x, Rational[1,2]]  →  "\sqrt{x}"
 *   Rational[1,6]            →  "\frac{1}{6}"
 *   Sin[x]                   →  "\sin(x)"
 *   InputForm[Rational[1,6]] →  ""   (use expr_to_string's "1/6")
 */
char* expr_to_latex(const Expr* e);

/*
 * The two typeset notations (see print_latex.c):
 *   LATEX_FORM_TRADITIONAL  textbook: \sin(x), \ln(x), e^{x}, x = 1 — what
 *                           expr_to_latex produces.
 *   LATEX_FORM_STANDARD     Mathilda names, typeset: \mathrm{Sin}{\left[x\right]},
 *                           ⅇ^{x}, x == 1; fractions/powers/radicals stay 2D.
 * TraditionalForm[...] / StandardForm[...] wrappers inside the tree switch the
 * notation for their subtree in either mode.
 */
typedef enum { LATEX_FORM_TRADITIONAL = 0, LATEX_FORM_STANDARD = 1 } LatexForm;

/*
 * expr_to_latex_form(e, form)
 *
 * As expr_to_latex, in the given notation. Caller frees. Returns "" (no typeset
 * form) for a printer directive as above, and — for LATEX_FORM_STANDARD only —
 * when the OUTERMOST head is TraditionalForm: that expression asked for the
 * TraditionalForm rendering outright, so the consumer should show that one.
 *
 * Examples (StandardForm):
 *   Sin[x]                   →  "\mathrm{Sin}{\left[\mathrm{x}\right]}"
 *   Power[E, x]              →  "ⅇ^{\mathrm{x}}"
 *   TraditionalForm[Sin[x]]  →  ""   (use the TraditionalForm LaTeX)
 */
char* expr_to_latex_form(const Expr* e, LatexForm form);

#endif /* MATHILDA_PRINT_LATEX_H */
