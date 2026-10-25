/*
 * dialog.c -- session history ($Line / In[n] / Out[n]) and Dialog[].
 *
 * Session history
 * ---------------
 * `%` is not a REPL feature, it is a KERNEL one. The parser turns it into
 * Out[-1] (parse.c) and builtin_out resolves the negative index against $Line,
 * so a front end that evaluates without recording this bookkeeping leaves `%`
 * sitting in its output unevaluated. The notebook did exactly that -- pipe mode
 * evaluated statements and recorded none of them, so `Integrate[x^5 E^x, x]`
 * followed by `% // Factor` printed a literal `Out[-1]`. Every front end, and
 * Dialog's own input loop, goes through the three mth_session_* helpers below,
 * so the history cannot drift between them again.
 *
 * Dialog
 * ------
 * Dialog[] / Dialog[expr] opens a nested session inside an evaluation: a
 * sequence of input and output lines, read through the reader the front end
 * installed (mth_dialog_set_io), until the user enters Return[] / Return[x]
 * or an end of input. In order:
 *
 *   1. localise $Line (to its current value), $MessageList, $Epilog and the
 *      DialogSymbols -- through the Block frame stack (mth_block_enter), so a
 *      TimeConstrained unwind through a dialog restores them like a Block's;
 *   2. evaluate DialogProlog;
 *   3. evaluate expr, which becomes the first dialog line's Out (the value
 *      of `%` for the first input);
 *   4. read, evaluate and print lines, numbered on from the global $Line;
 *   5. on exit, evaluate the local $Epilog (if one was set), restore, and
 *      return Return's argument -- or, for Return[] or an end of input, the
 *      value of the last line (Null if there was none).
 *
 * A dialog input line is a top-level evaluation even though the Dialog call is
 * itself inside one: an uncaught Throw is reported and made inert, Abort[]
 * gives $Aborted, Break[] is Break::nofwd -- none of them leaves the dialog
 * (eval_report_uncaught). Only an input that IS Return[] or Return[x] exits:
 * a Return produced by evaluation (`Return[2] // Identity`, `Return[4];`) is
 * an ordinary value. A trailing `;` suppresses the Out
 * line but not the value: `8;` stores Out[n] = 8, and a following Return[]
 * returns 8.
 */

#include "dialog.h"
#include "parse.h"
#include "print.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "message.h"
#include "repl_hooks.h"
#include "modular.h"
#include "core.h"
#include "quit.h"
#include "sym_names.h"
#include "show.h"
#include "render3d.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ========================================================================= */
/* Session history                                                           */
/* ========================================================================= */

const char* mth_skip_blanks_and_comments(const char* s) {
    while (*s) {
        if (isspace((unsigned char)*s)) {
            s++;
        } else if (s[0] == '(' && s[1] == '*') {
            const char* open = s;
            int depth = 1;
            s += 2;
            while (*s && depth > 0) {
                if (s[0] == '(' && s[1] == '*') { depth++; s += 2; }
                else if (s[0] == '*' && s[1] == ')') { depth--; s += 2; }
                else { s++; }
            }
            if (depth > 0) return open;  /* unterminated comment is a real error */
        } else {
            return s;
        }
    }
    return s;
}

bool mth_is_blank_or_comment_only(const char* s) {
    return *mth_skip_blanks_and_comments(s) == '\0';
}

/* A top-level NumberForm is stripped from the value stored in Out[n]/%, so
 * `%` (and arithmetic on it) sees the underlying number rather than the inert
 * format wrapper -- while the displayed Out[n]= line still uses the wrapper's
 * formatting. TOP-LEVEL only: a nested NumberForm, or one bound
 * to a variable via Set, is preserved. Returns a BORROWED pointer into `e`. */
static Expr* out_value_unwrapped(Expr* e) {
    if (e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_NumberForm
        && e->data.function.arg_count >= 1)
        return e->data.function.args[0];
    return e;
}

void mth_session_set_line(int line) {
    mth_message_list_reset();     /* $MessageList covers one input line */
    Expr* line_sym = expr_new_symbol(SYM_DollarLine);
    Expr* line_val = expr_new_integer(line);
    symtab_add_own_value("$Line", line_sym, line_val);
    expr_free(line_sym);
    expr_free(line_val);
}

/* add_down_value copies the value; expr_new_function consumes the head and arg
 * it is given, so only the pattern is ours to free. */
void mth_session_store_in(int line, Expr* parsed) {
    Expr* in_arg = expr_new_integer(line);
    Expr* in_pattern = expr_new_function(expr_new_symbol(SYM_In), &in_arg, 1);
    symtab_add_down_value("In", in_pattern, parsed);
    expr_free(in_pattern);
}

/* Stored for EVERY statement, including one ending in ';' and one that
 * evaluated to Null: a suppressed result is still the session's previous
 * result. */
void mth_session_store_out(int line, Expr* evaluated) {
    if (!evaluated) return;
    Expr* out_arg = expr_new_integer(line);
    Expr* out_pattern = expr_new_function(expr_new_symbol(SYM_Out), &out_arg, 1);
    symtab_add_down_value("Out", out_pattern, out_value_unwrapped(evaluated));
    expr_free(out_pattern);
}

/* ========================================================================= */
/* Dialog I/O                                                                */
/* ========================================================================= */

static MthDialogIO g_io;
static bool        g_io_set = false;
/* Open dialogs. Registered with the TimeConstrained unwind (dialog_init): a
 * timeout that jumps out of a dialog skips its decrement. */
static int         g_dialog_depth = 0;

void mth_dialog_set_io(const MthDialogIO* io) {
    if (io) { g_io = *io; g_io_set = true; }
    else    { memset(&g_io, 0, sizeof g_io); g_io_set = false; }
}

int mth_dialog_depth(void) { return g_dialog_depth; }

static void dlg_write(const char* s) {
    if (g_io_set && g_io.write) { g_io.write(s, g_io.ctx); return; }
    fputs(s, stdout);
    fflush(stdout);
}

static char* dlg_read(const char* prompt) {
    if (!g_io_set || !g_io.read) return NULL;
    return g_io.read(prompt, g_io.ctx);
}

/* One space per open dialog, before every prompt and Out line:
 * " In[3]:= " one level down, "  In[5]:= " two. */
static char* dlg_line_text(const char* fmt, int line, const char* body, const char* tail) {
    int depth = g_dialog_depth > 0 ? g_dialog_depth : 0;
    size_t n = (size_t)depth + strlen(fmt) + 24 + (body ? strlen(body) : 0) + strlen(tail);
    char* s = malloc(n);
    if (!s) return NULL;
    memset(s, ' ', (size_t)depth);
    int k = snprintf(s + depth, n - (size_t)depth, fmt, line);
    if (k < 0) { free(s); return NULL; }
    size_t pos = (size_t)depth + (size_t)k;
    if (body) { strcpy(s + pos, body); pos += strlen(body); }
    strcpy(s + pos, tail);
    return s;
}

/* Print a dialog line's result: nothing for Null; ?sym usage text verbatim;
 * otherwise the display form, with $PrePrint applied to the copy shown (Out[n]
 * keeps the unmodified value). A Graphics / Graphics3D result is displayed as
 * at the top level. Borrows `parsed` (may be NULL) and `value`. */
static void dlg_show(int line, const Expr* parsed, Expr* value) {
    if (value->type == EXPR_SYMBOL && value->data.symbol.name == SYM_Null) return;
    Expr* to_print = repl_apply_pre_print(expr_copy(value));
    if (!to_print) return;

    bool info = parsed && parsed->type == EXPR_FUNCTION
        && parsed->data.function.head->type == EXPR_SYMBOL
        && parsed->data.function.head->data.symbol.name == SYM_Information
        && to_print->type == EXPR_STRING;
    char* body = info ? NULL : expr_to_string_display(to_print);
    char* text = dlg_line_text("Out[%d]= ", line,
                               info ? to_print->data.string : (body ? body : ""),
                               info ? "\n" : "\n\n");
    if (text) dlg_write(text);
    free(text);
    free(body);

    if (to_print->type == EXPR_FUNCTION
        && to_print->data.function.head->type == EXPR_SYMBOL
        && to_print->data.function.arg_count >= 1
        && getenv("MATHILDA_NO_WINDOW") == NULL) {
        if (to_print->data.function.head->data.symbol.name == SYM_Graphics) graphics_show(to_print);
        else if (to_print->data.function.head->data.symbol.name == SYM_Graphics3D) graphics3d_show(to_print);
    }
    expr_free(to_print);
}

/* ========================================================================= */
/* The dialog's input loop                                                   */
/* ========================================================================= */

static bool dlg_is_head(const Expr* e, const char* head) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == head;
}

static bool dlg_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

/* The input's last significant character is ';': its Out line is suppressed.
 * (A trailing comment after the ';' hides it; the value is unaffected.) */
static bool dlg_ends_with_semicolon(const char* s) {
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) n--;
    return n > 0 && s[n - 1] == ';';
}

/* `a; b;` parses as CompoundExpression[a, b, Null]; the line's value is b,
 * not that trailing Null. Drops the Null. Takes `e`. */
static Expr* dlg_drop_trailing_null(Expr* e) {
    if (!dlg_is_head(e, SYM_CompoundExpression)) return e;
    size_t n = e->data.function.arg_count;
    if (n < 2 || !dlg_is_symbol(e->data.function.args[n - 1], SYM_Null)) return e;
    Expr* out;
    if (n == 2) {
        out = expr_copy(e->data.function.args[0]);
    } else {
        Expr** kept = malloc((n - 1) * sizeof(Expr*));
        if (!kept) return e;
        for (size_t i = 0; i + 1 < n; i++) kept[i] = expr_copy(e->data.function.args[i]);
        out = expr_new_function(expr_new_symbol(SYM_CompoundExpression), kept, n - 1);
        free(kept);
    }
    expr_free(e);
    return out;
}

static bool dlg_ignore_eof(void) {
    Rule* r = symtab_get_own_values("$IgnoreEOF");
    return r && dlg_is_symbol(r->replacement, SYM_True);
}

/* Evaluate one dialog input line as line `line`. `*last` holds the dialog's
 * current value (owned; may be NULL) and is replaced by this line's value.
 * Returns true when the dialog ends here (Return, or Quit[]). */
static bool dlg_eval_line(const char* text, int line, Expr** last) {
    mth_session_set_line(line);

    char* cooked = repl_apply_pre_read(text);
    if (!cooked) { dlg_write("Parse error\n\n"); return false; }
    bool suppressed = dlg_ends_with_semicolon(cooked);
    Expr* parsed = parse_expression(cooked);
    free(cooked);
    if (!parsed) { dlg_write("Parse error\n\n"); return false; }

    mth_session_store_in(line, parsed);

    /* Only the input itself being Return[] / Return[x] ends the dialog. */
    bool is_return = dlg_is_head(parsed, SYM_Return) && parsed->data.function.arg_count <= 1;

    Expr* input = expr_copy(parsed);
    if (suppressed) input = dlg_drop_trailing_null(input);
    input = repl_apply_pre(input);
    Expr* v = evaluate(input);
    expr_free(input);

    /* A top level of its own: finish any sentinel here, and clear what a
     * sentinel or a recursion overflow left for the enclosing frames (the
     * Dialog call is still inside the outer evaluation). */
    v = eval_report_uncaught(v);
    (void)eval_take_top_abort();
    eval_reset_recursion_depth(eval_get_recursion_depth());

    if (mth_quit_requested()) {
        expr_free(v);
        expr_free(parsed);
        return true;
    }

    if (is_return) {
        /* Return[x] at a top level stays Return[x] (no enclosing function
         * consumes it): its argument is the dialog's value. Return[] keeps
         * the value of the last line. A $Pre that rewrote the input leaves
         * some other value, which is then the result. */
        if (dlg_is_head(v, SYM_Return) && v->data.function.arg_count == 0) {
            expr_free(v);
        } else {
            Expr* val = v;
            if (dlg_is_head(v, SYM_Return) && v->data.function.arg_count == 1) {
                val = expr_copy(v->data.function.args[0]);
                expr_free(v);
            }
            expr_free(*last);
            *last = val;
        }
        expr_free(parsed);
        return true;
    }

    v = repl_apply_post(v);
    if (!v) { expr_free(parsed); return false; }

    mth_session_store_out(line, v);
    if (!suppressed) dlg_show(line, parsed, v);
    expr_free(*last);
    *last = v;
    expr_free(parsed);
    return false;
}

/* The integer value of $Line, or 0 when there is none (a -file script, an
 * embedding host): the numbering the dialog continues from. */
static int dlg_current_line(void) {
    Rule* r = symtab_get_own_values("$Line");
    if (r && r->replacement && r->replacement->type == EXPR_INTEGER
        && r->replacement->data.integer >= 0 && r->replacement->data.integer < 1000000000)
        return (int)r->replacement->data.integer;
    return 0;
}

/* Run a dialog. `arg` (may be NULL), `prolog` (may be NULL) and the entries of
 * `syms` (a List, may be NULL) are borrowed, unevaluated. Returns the dialog's
 * value, or an in-flight Throw / Abort from DialogProlog or the argument,
 * which leaves the dialog unopened and propagates to the caller. */
static Expr* dialog_run(Expr* arg, Expr* prolog, Expr* syms) {
    int line = dlg_current_line();

    size_t nsyms = syms ? syms->data.function.arg_count : 0;
    size_t n = 3 + nsyms;
    Expr** specs = malloc(n * sizeof(Expr*));
    if (!specs) return NULL;
    Expr* line_args[2] = { expr_new_symbol(SYM_DollarLine), expr_new_integer(line) };
    specs[0] = expr_new_function(expr_new_symbol(SYM_Set), line_args, 2);
    specs[1] = expr_new_symbol(SYM_DollarMessageList);
    specs[2] = expr_new_symbol(SYM_DollarEpilog);
    for (size_t i = 0; i < nsyms; i++) specs[3 + i] = syms->data.function.args[i];

    void* frame = mth_block_enter(specs, n);
    for (size_t i = 0; i < 3; i++) expr_free(specs[i]);
    free(specs);

    mth_message_list_reset();   /* the local $MessageList */
    g_dialog_depth++;

    Expr* last = NULL;
    Expr* escape = NULL;

    if (prolog && !dlg_is_symbol(prolog, SYM_Null)) {
        Expr* r = evaluate(prolog);
        if (eval_is_inflight_unwind(r)) escape = r;
        else expr_free(r);
    }

    if (!escape && arg && !mth_quit_requested()) {
        line++;
        mth_session_set_line(line);
        Expr* v = evaluate(arg);
        if (eval_is_inflight_unwind(v)) {
            escape = v;
        } else {
            mth_session_store_out(line, v);
            /* No reader means no session to show it in (and, in pipe mode,
             * stdout is the protocol). */
            if (g_io_set && g_io.read) dlg_show(line, NULL, v);
            last = v;
        }
    }

    while (!escape && !mth_quit_requested()) {
        char* prompt = dlg_line_text("In[%d]:= ", line + 1, NULL, "");
        char* text = prompt ? dlg_read(prompt) : NULL;
        free(prompt);
        if (!text) {
            /* End of input. At a terminal, $IgnoreEOF = True ignores the
             * keystroke; a pipe or file at its end always ends the dialog. */
            if (g_io_set && g_io.read) dlg_write("\n");
            if (g_io_set && g_io.read && g_io.interactive && dlg_ignore_eof()) continue;
            break;
        }
        if (mth_is_blank_or_comment_only(text)) { free(text); continue; }
        line++;
        bool done = dlg_eval_line(text, line, &last);
        free(text);
        if (done) break;
    }

    if (!escape) repl_apply_epilog();   /* the dialog's own $Epilog, if set */
    g_dialog_depth--;
    mth_block_leave(frame);

    if (escape) { expr_free(last); return escape; }
    return last ? last : expr_new_symbol(SYM_Null);
}

/* ========================================================================= */
/* Dialog[] and its options                                                  */
/* ========================================================================= */

static bool dlg_is_rule(const Expr* e) {
    return (dlg_is_head(e, SYM_Rule) || dlg_is_head(e, SYM_RuleDelayed))
        && e->data.function.arg_count == 2;
}

/* An option argument: a rule, or a non-empty list of rules. */
static bool dlg_is_option_arg(const Expr* e) {
    if (dlg_is_rule(e)) return true;
    if (!dlg_is_head(e, SYM_List) || e->data.function.arg_count == 0) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!dlg_is_rule(e->data.function.args[i])) return false;
    return true;
}

/* Record one option rule. The FIRST setting of a name wins. Returns false
 * (after Dialog::optx) on an unknown name. */
static bool dlg_take_option(Expr* rule, Expr* res, Expr** syms, Expr** prolog,
                            bool* syms_set, bool* prolog_set) {
    Expr* key = rule->data.function.args[0];
    Expr* val = rule->data.function.args[1];
    if (dlg_is_symbol(key, SYM_DialogSymbols)) {
        if (!*syms_set) { *syms = val; *syms_set = true; }
        return true;
    }
    if (dlg_is_symbol(key, SYM_DialogProlog)) {
        if (!*prolog_set) { *prolog = val; *prolog_set = true; }
        return true;
    }
    char* k = expr_to_string(key);
    char* call = expr_to_string(res);
    mth_message("Dialog", "optx", "Unknown option %s in %s.", k ? k : "?", call ? call : "?");
    free(k);
    free(call);
    return false;
}

/* DialogSymbols must be a list of symbols and sym = value assignments, each
 * symbol once. Emits Dialog::lvlist / ::lvsym / ::dup and returns false
 * otherwise. */
static bool dlg_check_symbols(Expr* syms) {
    char* spec = expr_to_string(syms);
    bool ok = true;
    if (!dlg_is_head(syms, SYM_List)) {
        mth_message("Dialog", "lvlist", "Local variable specification DialogSymbols :> %s is not a List.",
                    spec ? spec : "?");
        ok = false;
    }
    size_t n = ok ? syms->data.function.arg_count : 0;
    for (size_t i = 0; ok && i < n; i++) {
        Expr* v = syms->data.function.args[i];
        const char* name = NULL;
        if (v->type == EXPR_SYMBOL) name = v->data.symbol.name;
        else if (dlg_is_head(v, SYM_Set) && v->data.function.arg_count == 2
                 && v->data.function.args[0]->type == EXPR_SYMBOL)
            name = v->data.function.args[0]->data.symbol.name;
        if (!name) {
            char* bad = expr_to_string(v);
            mth_message("Dialog", "lvsym",
                        "Local variable specification DialogSymbols :> %s contains %s, "
                        "which is not a symbol or an assignment to a symbol.",
                        spec ? spec : "?", bad ? bad : "?");
            free(bad);
            ok = false;
            break;
        }
        for (size_t j = 0; j < i; j++) {
            Expr* w = syms->data.function.args[j];
            const char* other = w->type == EXPR_SYMBOL ? w->data.symbol.name
                                                       : w->data.function.args[0]->data.symbol.name;
            if (other == name) {   /* interned: pointer equality */
                mth_message("Dialog", "dup",
                            "Duplicate local variable %s found in local variable specification "
                            "DialogSymbols :> %s.", name, spec ? spec : "?");
                ok = false;
                break;
            }
        }
    }
    free(spec);
    return ok;
}

/* Dialog[] / Dialog[expr] / Dialog[expr, opts] -- HoldAll, Protected. */
static Expr* builtin_dialog(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    size_t npos = argc;
    while (npos > 0 && dlg_is_option_arg(args[npos - 1])) npos--;
    if (npos > 1) {
        char* bad = expr_to_string(args[npos - 1]);
        char* call = expr_to_string(res);
        mth_message("Dialog", "nonopt",
                    "Options expected (instead of %s) beyond position 1 in %s. "
                    "An option must be a rule or a list of rules.",
                    bad ? bad : "?", call ? call : "?");
        free(bad);
        free(call);
        return NULL;
    }

    /* The call's settings first (first wins), then Options[Dialog]. */
    Expr* syms = NULL;
    Expr* prolog = NULL;
    bool syms_set = false, prolog_set = false;
    for (size_t i = npos; i < argc; i++) {
        Expr* o = args[i];
        if (dlg_is_rule(o)) {
            if (!dlg_take_option(o, res, &syms, &prolog, &syms_set, &prolog_set)) return NULL;
        } else {
            for (size_t j = 0; j < o->data.function.arg_count; j++)
                if (!dlg_take_option(o->data.function.args[j], res, &syms, &prolog,
                                     &syms_set, &prolog_set)) return NULL;
        }
    }
    Expr* defaults = symtab_get_options("Dialog");
    if (dlg_is_head(defaults, SYM_List)) {
        for (size_t i = 0; i < defaults->data.function.arg_count; i++) {
            Expr* o = defaults->data.function.args[i];
            if (!dlg_is_rule(o)) continue;
            Expr* key = o->data.function.args[0];
            if (dlg_is_symbol(key, SYM_DialogSymbols) && !syms_set) { syms = o->data.function.args[1]; syms_set = true; }
            if (dlg_is_symbol(key, SYM_DialogProlog) && !prolog_set) { prolog = o->data.function.args[1]; prolog_set = true; }
        }
    }

    if (syms && !dlg_check_symbols(syms)) return NULL;
    return dialog_run(npos == 1 ? args[0] : NULL, prolog, syms);
}

/* ========================================================================= */
/* Registration                                                              */
/* ========================================================================= */

void dialog_init(void) {
    symtab_add_builtin("Dialog", builtin_dialog);
    symtab_get_def("Dialog")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("DialogSymbols")->attributes |= ATTR_PROTECTED;
    symtab_get_def("DialogProlog")->attributes |= ATTR_PROTECTED;

    Expr* opts[2];
    Expr* prolog_rule[2] = { expr_new_symbol(SYM_DialogProlog), expr_new_symbol(SYM_Null) };
    opts[0] = expr_new_function(expr_new_symbol(SYM_Rule), prolog_rule, 2);
    Expr* syms_rule[2] = { expr_new_symbol(SYM_DialogSymbols),
                           expr_new_function(expr_new_symbol(SYM_List), NULL, 0) };
    opts[1] = expr_new_function(expr_new_symbol(SYM_Rule), syms_rule, 2);
    symtab_set_options("Dialog", expr_new_function(expr_new_symbol(SYM_List), opts, 2));

    Expr* eof_sym = expr_new_symbol(SYM_DollarIgnoreEOF);
    Expr* eof_val = expr_new_symbol(SYM_False);
    symtab_add_own_value("$IgnoreEOF", eof_sym, eof_val);
    expr_free(eof_sym);
    expr_free(eof_val);

    tc_register_unwind_counter(&g_dialog_depth);

    symtab_set_docstring("Dialog",
        "Dialog[]\n"
        "\tinitiates a dialog: a nested sequence of input and output lines.\n"
        "Dialog[expr]\n"
        "\tinitiates a dialog with expr as the current value of %.\n"
        "\n"
        "Exit with Return[expr] (the value of Dialog) or Return[] or an end of\n"
        "input (the value of the last line). Line numbering continues from\n"
        "$Line and reverts on exit. $Line, $MessageList and $Epilog are\n"
        "localized; a local $Epilog is evaluated on exit. $Pre, $Post and the\n"
        "other session hooks apply to dialog lines.\n"
        "Options: DialogSymbols :> {x, y = y0, ...} localizes symbols;\n"
        "DialogProlog :> expr is evaluated after localization, before expr.\n"
        "Attributes: HoldAll, Protected.");
    symtab_set_docstring("DialogSymbols",
        "DialogSymbols\n"
        "\tis an option for Dialog: DialogSymbols :> {x, y = y0, ...} gives\n"
        "\tsymbols local values (and optional initial values) within the dialog.");
    symtab_set_docstring("DialogProlog",
        "DialogProlog\n"
        "\tis an option for Dialog: DialogProlog :> expr evaluates expr after the\n"
        "\tdialog's symbols are localized and before its argument.");
    symtab_set_docstring("$IgnoreEOF",
        "$IgnoreEOF\n"
        "\tspecifies whether an end-of-file character typed at a terminal is\n"
        "\tignored inside a Dialog (True) or ends it (False, the default).");
}
