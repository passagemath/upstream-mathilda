#ifndef MATHILDA_DIALOG_H
#define MATHILDA_DIALOG_H

/*
 * dialog.{c,h} -- the session's top-level input bookkeeping, and Dialog[].
 *
 * Two things live here because they are one thing: the $Line / In[n] / Out[n]
 * history every front end records (the interactive REPL, the NDJSON pipe, and
 * a dialog's nested input loop all go through the same three helpers, so they
 * cannot drift), and Dialog[], which is a nested input loop.
 *
 * Dialog[] reads its input lines through a reader the front end installs
 * (mth_dialog_set_io): readline in the interactive REPL, stdin under -file.
 * With no reader installed -- the NDJSON pipe, where stdin is the protocol, or
 * an embedding host -- a dialog sees end of input at once: it localises,
 * runs DialogProlog, evaluates its argument, runs its $Epilog and returns.
 */

#include "expr.h"
#include <stdbool.h>

/* ---- session history: $Line, In[n], Out[n] ------------------------------- */

/* Advance past whitespace and (* ... *) comments (nested comments allowed),
 * returning the first character that is neither -- the terminating NUL if the
 * text holds nothing else. An unterminated comment is NOT skipped: the result
 * is its opening '(', where that error starts. */
const char* mth_skip_blanks_and_comments(const char* s);

/* True if `s` holds only whitespace and (* ... *) comments: a no-op input. */
bool mth_is_blank_or_comment_only(const char* s);

/* Set $Line = line and start a fresh $MessageList. Must run BEFORE the line is
 * evaluated: `%` is Out[$Line - 1]. */
void mth_session_set_line(int line);

/* In[line] = parsed (as typed, before $Pre). Borrows `parsed`. */
void mth_session_store_in(int line, Expr* parsed);

/* Out[line] = evaluated, a top-level NumberForm stripped. Borrows `evaluated`;
 * NULL stores nothing. */
void mth_session_store_out(int line, Expr* evaluated);

/* ---- Dialog[] input and output -------------------------------------------- */

/* Return the next input line as a malloc'd string (no terminator), or NULL at
 * end of input. The reader shows `prompt` itself (readline does). */
typedef char* (*MthDialogReadFn)(const char* prompt, void* ctx);
/* Write dialog output text (Out lines, parse errors). */
typedef void (*MthDialogWriteFn)(const char* text, void* ctx);

typedef struct {
    MthDialogReadFn  read;    /* NULL: no input source (end of input at once) */
    MthDialogWriteFn write;   /* NULL: stdout */
    void*            ctx;
    /* A terminal: an end-of-file character there is a keystroke, which
     * $IgnoreEOF = True may ignore. A pipe or file at end of input always
     * ends the dialog (ignoring it would spin). */
    bool             interactive;
} MthDialogIO;

/* Install the dialog I/O (copied); NULL removes it. */
void mth_dialog_set_io(const MthDialogIO* io);

/* Number of dialogs currently open (0 at the session's own top level). */
int mth_dialog_depth(void);

/* Registers Dialog, DialogSymbols, DialogProlog, $IgnoreEOF. From core_init(). */
void dialog_init(void);

#endif /* MATHILDA_DIALOG_H */
