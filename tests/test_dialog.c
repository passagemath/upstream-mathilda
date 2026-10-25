/*
 * test_dialog.c
 *
 * Dialog[]: the nested input loop, driven through a scripted reader and a
 * capturing writer (mth_dialog_set_io), so each test sees exactly the prompts
 * and Out lines a terminal would show.
 *
 * Covered: the value (Return[x], Return[], end of input, a `;` line), line
 * numbering from $Line and its restoration, `%`, which inputs exit (only one
 * that IS Return[] / Return[x]), dialog lines as top levels (Throw, Abort,
 * Break stay inside; an outer Catch does not see them), DialogSymbols /
 * DialogProlog and their order, the option-list form, $Epilog / $MessageList
 * localisation, $Pre / $Post / $PreRead / $PrePrint, nesting and indentation,
 * $IgnoreEOF, parse errors, blank lines, the argument and option messages,
 * attributes, no reader installed, and Quit[] (last: it cannot be undone).
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "parse.h"
#include "dialog.h"
#include "quit.h"
#include "message.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* scripted I/O                                                        */
/* ------------------------------------------------------------------ */

/* A NULL entry is one end of input (a Ctrl-D at a terminal); the script's
 * end is a permanent one. */
#define EOF_MARK ((const char*)1)

typedef struct {
    const char* const* lines;
    size_t n, next;
    char*  out;          /* transcript: prompts, echoed input, output */
    size_t len, cap;
} Script;

static void sc_append(Script* s, const char* t) {
    size_t k = strlen(t);
    if (s->len + k + 1 > s->cap) {
        size_t cap = (s->cap ? s->cap * 2 : 256);
        while (cap < s->len + k + 1) cap *= 2;
        s->out = realloc(s->out, cap);
        s->cap = cap;
    }
    memcpy(s->out + s->len, t, k + 1);
    s->len += k;
}

static char* sc_read(const char* prompt, void* ctx) {
    Script* s = ctx;
    sc_append(s, prompt);
    if (s->next >= s->n) return NULL;
    const char* l = s->lines[s->next++];
    if (l == EOF_MARK) return NULL;
    sc_append(s, l);
    sc_append(s, "\n");
    char* copy = malloc(strlen(l) + 1);
    strcpy(copy, l);
    return copy;
}

static void sc_write(const char* text, void* ctx) { sc_append(ctx, text); }

/* Evaluate `src` with `lines` as the dialog input. Returns the printed result
 * (caller frees); *transcript gets the transcript (caller frees). */
static char* run_dialog(const char* src, const char* const* lines, size_t n,
                        bool interactive, char** transcript) {
    Script s = { lines, n, 0, NULL, 0, 0 };
    sc_append(&s, "");
    MthDialogIO io = { sc_read, sc_write, &s, interactive };
    mth_dialog_set_io(&io);
    Expr* p = parse_expression(src);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* str = expr_to_string(r);
    expr_free(r);
    expr_free(p);
    mth_dialog_set_io(NULL);
    if (transcript) *transcript = s.out; else free(s.out);
    ASSERT(mth_dialog_depth() == 0);
    return str;
}

static void eval_is(const char* src, const char* want) {
    Expr* p = parse_expression(src);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, want) != 0) fprintf(stderr, "FAIL: %s\n  expected %s\n  got      %s\n", src, want, s);
    ASSERT(strcmp(s, want) == 0);
    free(s);
    expr_free(r);
    expr_free(p);
}

static void dialog_is(const char* src, const char* const* lines, size_t n,
                      const char* want, const char* want_transcript) {
    char* t = NULL;
    char* got = run_dialog(src, lines, n, false, &t);
    if (strcmp(got, want) != 0) fprintf(stderr, "FAIL: %s\n  expected %s\n  got      %s\n", src, want, got);
    ASSERT(strcmp(got, want) == 0);
    if (want_transcript && strcmp(t, want_transcript) != 0)
        fprintf(stderr, "FAIL transcript of %s:\n--- expected\n%s--- got\n%s---\n", src, want_transcript, t);
    if (want_transcript) ASSERT(strcmp(t, want_transcript) == 0);
    free(got);
    free(t);
}

static bool contains(const char* hay, const char* needle) { return strstr(hay, needle) != NULL; }

#define N(a) (sizeof(a) / sizeof((a)[0]))

/* Each test numbers from $Line = 4, so the expected transcripts are fixed. */
static void at_line_4(void) { mth_session_set_line(4); }

/* ------------------------------------------------------------------ */
/* registration                                                        */
/* ------------------------------------------------------------------ */

static void test_attributes_and_options(void) {
    eval_is("Attributes[Dialog]", "{HoldAll, Protected}");
    eval_is("Options[Dialog]", "{DialogProlog -> Null, DialogSymbols -> {}}");
    eval_is("Attributes[DialogSymbols]", "{Protected}");
    eval_is("Attributes[DialogProlog]", "{Protected}");
    eval_is("$IgnoreEOF", "False");
    eval_is("StringQ[Information[Dialog]]", "True");
    eval_is("Quiet[Dialog = 3]; Attributes[Dialog]", "{HoldAll, Protected}");
}

/* ------------------------------------------------------------------ */
/* the value                                                           */
/* ------------------------------------------------------------------ */

static void test_no_reader_returns_at_once(void) {
    eval_is("Dialog[]", "Null");
    eval_is("Dialog[5]", "5");
    eval_is("Dialog[t = 3; t + 1]", "4");   /* HoldAll: evaluated inside */
    eval_is("Dialog[DialogProlog :> (pr = 9)]; pr", "9");
    eval_is("Clear[t, pr]", "Null");
}

static void test_return_with_value(void) {
    at_line_4();
    const char* in[] = { "x = 7", "Return[x + 1]" };
    dialog_is("Dialog[]", in, N(in), "8",
              " In[5]:= x = 7\n Out[5]= 7\n\n In[6]:= Return[x + 1]\n");
    eval_is("Clear[x]", "Null");
}

static void test_return_without_value_gives_last_line(void) {
    at_line_4();
    const char* in[] = { "10", "20", "Return[]" };
    dialog_is("Dialog[]", in, N(in), "20", NULL);
    const char* none[] = { "Return[]" };
    dialog_is("Dialog[]", none, N(none), "Null", " In[5]:= Return[]\n");
    const char* arg[] = { "Return[]" };
    dialog_is("Dialog[7]", arg, N(arg), "7", NULL);
}

static void test_end_of_input_gives_last_line(void) {
    at_line_4();
    const char* in[] = { "10", "20" };
    dialog_is("Dialog[]", in, N(in), "20",
              " In[5]:= 10\n Out[5]= 10\n\n In[6]:= 20\n Out[6]= 20\n\n In[7]:= \n");
    dialog_is("Dialog[]", NULL, 0, "Null", " In[5]:= \n");
    dialog_is("Dialog[7]", NULL, 0, "7", " Out[5]= 7\n\n In[6]:= \n");
}

static void test_semicolon_line_is_the_value_but_not_shown(void) {
    at_line_4();
    const char* in[] = { "8;", "Return[]" };
    dialog_is("Dialog[]", in, N(in), "8", " In[5]:= 8;\n In[6]:= Return[]\n");
    const char* two[] = { "1; 2;", "%", "Return[]" };
    dialog_is("Dialog[]", two, N(two), "2", NULL);
    const char* nul[] = { "5", "Null", "Return[]" };   /* a Null line shows nothing */
    dialog_is("Dialog[]", nul, N(nul), "Null", " In[5]:= 5\n Out[5]= 5\n\n In[6]:= Null\n In[7]:= Return[]\n");
}

/* ------------------------------------------------------------------ */
/* numbering and %                                                     */
/* ------------------------------------------------------------------ */

static void test_numbering_follows_line_and_reverts(void) {
    at_line_4();
    const char* in[] = { "%", "$Line", "Return[%]" };
    dialog_is("Dialog[5]", in, N(in), "7",
              " Out[5]= 5\n\n In[6]:= %\n Out[6]= 5\n\n In[7]:= $Line\n Out[7]= 7\n\n In[8]:= Return[%]\n");
    eval_is("$Line", "4");
    eval_is("{Out[5], Out[6], Out[7]}", "{5, 5, 7}");   /* history is global */
    eval_is("In[7]", "4");                                /* In[7] is $Line, re-evaluated out here */
}

static void test_argument_line_counts_even_when_null(void) {
    at_line_4();
    const char* in[] = { "Return[$Line]" };
    dialog_is("Dialog[q = 1;]", in, N(in), "6", " In[6]:= Return[$Line]\n");
    eval_is("Clear[q]", "Null");
}

static void test_blank_and_comment_lines_take_no_number(void) {
    at_line_4();
    const char* in[] = { "", "   ", "(* note *)", "Return[$Line]" };
    dialog_is("Dialog[]", in, N(in), "5", NULL);
}

static void test_parse_error_is_reported_and_dialog_continues(void) {
    at_line_4();
    const char* in[] = { "1 +* )", "Return[$Line]" };
    char* t = NULL;
    char* got = run_dialog("Dialog[]", in, N(in), false, &t);
    ASSERT(strcmp(got, "6") == 0);
    ASSERT(contains(t, "Parse error"));
    free(got); free(t);
}

/* ------------------------------------------------------------------ */
/* what exits                                                          */
/* ------------------------------------------------------------------ */

static void test_only_an_input_return_exits(void) {
    at_line_4();
    const char* in[] = { "If[True, Return[4]]", "Return[5];", "z := Return[9]", "z",
                         "Return[2] // Identity", "(Return[3])" };
    char* t = NULL;
    char* got = run_dialog("Dialog[]", in, N(in), false, &t);
    ASSERT(strcmp(got, "3") == 0);
    ASSERT(contains(t, " Out[5]= Return[4]\n"));
    ASSERT(contains(t, " In[7]:= z := Return[9]\n In[8]:= z\n Out[8]= 9\n"));
    ASSERT(contains(t, " Out[9]= Return[2]\n"));
    free(got); free(t);
    eval_is("Clear[z]", "Null");
}

/* ------------------------------------------------------------------ */
/* dialog lines are top levels                                         */
/* ------------------------------------------------------------------ */

static void test_throw_abort_break_stay_inside(void) {
    at_line_4();
    const char* in[] = { "Quiet[Throw[3]]", "Abort[]", "Quiet[Break[]]", "Return[{1, 2}]" };
    char* t = NULL;
    char* got = run_dialog("Catch[Dialog[]; 5]", in, N(in), false, &t);
    ASSERT(strcmp(got, "5") == 0);              /* the outer Catch saw no Throw */
    ASSERT(contains(t, " Out[5]= Hold[Throw[3]]\n"));
    ASSERT(contains(t, " Out[6]= $Aborted\n"));
    ASSERT(contains(t, " Out[7]= Hold[Break[]]\n"));
    free(got); free(t);
    ASSERT(!eval_take_top_abort());             /* the dialog consumed the abort */
}

static void test_loop_around_dialog_keeps_running(void) {
    at_line_4();
    const char* in[] = { "Quiet[Break[]]", "Return[]", "Return[]" };
    char* t = NULL;
    char* got = run_dialog("n = 0; Do[n++; Dialog[], {2}]; n", in, N(in), false, &t);
    ASSERT(strcmp(got, "2") == 0);
    free(got); free(t);
    eval_is("Clear[n]", "Null");
}

static void test_throw_from_prolog_or_argument_propagates(void) {
    at_line_4();
    char* t = NULL;
    char* got = run_dialog("Catch[Dialog[DialogProlog :> Throw[7]]]", NULL, 0, false, &t);
    ASSERT(strcmp(got, "7") == 0);
    ASSERT(strcmp(t, "") == 0);                 /* never prompted */
    free(got); free(t);
    got = run_dialog("Catch[Dialog[Throw[8]]]", NULL, 0, false, &t);
    ASSERT(strcmp(got, "8") == 0);
    free(got); free(t);
    eval_is("$Line", "4");
}

/* ------------------------------------------------------------------ */
/* localisation                                                        */
/* ------------------------------------------------------------------ */

static void test_dialog_symbols_and_prolog_order(void) {
    at_line_4();
    eval_is("x = 1; y = 2;", "Null");
    const char* in[] = { "{x, y}", "x = 5", "Return[x y]" };
    char* t = NULL;
    char* got = run_dialog("Dialog[p = {x, y}; x + y, DialogSymbols :> {x, y = 10}, "
                           "DialogProlog :> (pr = {x, y})]", in, N(in), false, &t);
    ASSERT(strcmp(got, "50") == 0);
    ASSERT(contains(t, " Out[5]= 10 + x\n"));
    ASSERT(contains(t, " Out[6]= {x, 10}\n"));
    free(got); free(t);
    eval_is("{x, y, pr, p}", "{1, 2, {1, 10}, {1, 10}}");   /* localise, prolog, argument */
    eval_is("Clear[x, y, pr, p]", "Null");
}

static void test_option_list_and_rule_forms(void) {
    at_line_4();
    const char* in[] = { "Return[{w, $Line}]" };
    dialog_is("Dialog[{DialogProlog :> (pq = 1), DialogSymbols -> {w = 3}}]", in, N(in), "{3, 5}", NULL);
    eval_is("{pq, w}", "{1, w}");
    eval_is("Clear[pq]", "Null");
}

static void test_first_option_setting_wins(void) {
    at_line_4();
    const char* in[] = { "Return[{b, c}]" };
    eval_is("b = 1; c = 2;", "Null");
    dialog_is("Dialog[DialogSymbols :> {b = 10}, DialogSymbols :> {c = 20}]", in, N(in), "{10, 2}", NULL);
    eval_is("{b, c}", "{1, 2}");
    eval_is("Clear[b, c]", "Null");
}

static void test_options_follow_set_options(void) {
    at_line_4();
    eval_is("SetOptions[Dialog, DialogProlog :> (so = 1)];", "Null");
    eval_is("Dialog[]; so", "1");                /* no reader: still runs the prolog */
    eval_is("SetOptions[Dialog, DialogProlog -> Null]; Clear[so];", "Null");
    eval_is("Options[Dialog]", "{DialogProlog -> Null, DialogSymbols -> {}}");
}

static void test_epilog_is_local_and_runs_on_exit(void) {
    at_line_4();
    eval_is("g = 0; $Epilog := (g = g + 1);", "Null");
    const char* in[] = { "$Epilog", "$Epilog := (ep = \"ran\")", "Return[3]" };
    char* t = NULL;
    char* got = run_dialog("Dialog[]", in, N(in), false, &t);
    ASSERT(strcmp(got, "3") == 0);
    ASSERT(contains(t, " Out[5]= $Epilog\n"));  /* unset inside the dialog */
    free(got); free(t);
    eval_is("{ep === \"ran\", g}", "{True, 0}");             /* the local one ran, the global one did not */
    eval_is("Clear[$Epilog, g, ep]", "Null");
}

static void test_message_list_is_local(void) {
    at_line_4();
    eval_is("Quiet[$MessageList = {HoldForm[outer::tag]}]", "{outer::tag}");
    const char* in[] = { "Quiet[1/0]; $MessageList", "1/0; $MessageList", "Return[]" };
    char* t = NULL;
    char* got = run_dialog("Dialog[]", in, N(in), false, &t);
    ASSERT(strcmp(got, "{Power::infy}") == 0);
    ASSERT(contains(t, " Out[5]= {}\n"));       /* fresh per line; Quiet lists nothing */
    free(got); free(t);
    eval_is("$MessageList", "{outer::tag}");
}

/* ------------------------------------------------------------------ */
/* session hooks                                                       */
/* ------------------------------------------------------------------ */

static void test_session_hooks_apply(void) {
    at_line_4();
    eval_is("cnt = 0; $Pre = Function[e, (cnt++; e), HoldAll]; $Post = post;", "Null");
    const char* in[] = { "3", "Return[2]" };
    char* t = NULL;
    char* got = run_dialog("Dialog[]", in, N(in), false, &t);
    ASSERT(strcmp(got, "2") == 0);
    ASSERT(contains(t, " Out[5]= post[3]\n"));
    free(got); free(t);
    eval_is("Clear[$Pre, $Post]; cnt", "2");    /* $Pre saw both lines, Return's too */
    eval_is("$PreRead = StringReplace[#, \"a\" -> \"b\"] &; $PrePrint = shown;", "Null");
    const char* in2[] = { "a", "Return[%]" };
    got = run_dialog("Dialog[]", in2, N(in2), false, &t);
    ASSERT(strcmp(got, "b") == 0);              /* Out keeps the value, not shown[b] */
    ASSERT(contains(t, " Out[5]= shown[b]\n"));
    free(got); free(t);
    eval_is("Clear[$PreRead, $PrePrint, cnt]", "Null");
}

/* ------------------------------------------------------------------ */
/* nesting                                                             */
/* ------------------------------------------------------------------ */

static void test_nested_dialogs_indent_and_return(void) {
    at_line_4();
    const char* in[] = { "s = Dialog[1]", "Return[2]", "Return[{s, $Line}]" };
    dialog_is("Dialog[]", in, N(in), "{2, 6}",
              " In[5]:= s = Dialog[1]\n  Out[6]= 1\n\n  In[7]:= Return[2]\n"
              " Out[5]= 2\n\n In[6]:= Return[{s, $Line}]\n");
    eval_is("Clear[s]", "Null");
}

/* ------------------------------------------------------------------ */
/* end of input and $IgnoreEOF                                         */
/* ------------------------------------------------------------------ */

static void test_ignore_eof(void) {
    at_line_4();
    const char* in[] = { "1", EOF_MARK, "Return[6]" };
    char* t = NULL;
    char* got = run_dialog("Dialog[]", in, N(in), true, &t);    /* $IgnoreEOF = False */
    ASSERT(strcmp(got, "1") == 0);
    free(got); free(t);

    eval_is("$IgnoreEOF = True;", "Null");
    got = run_dialog("Dialog[]", in, N(in), true, &t);          /* terminal: ignored */
    ASSERT(strcmp(got, "6") == 0);
    free(got); free(t);
    got = run_dialog("Dialog[]", in, N(in), false, &t);         /* a pipe: ends anyway */
    ASSERT(strcmp(got, "1") == 0);
    free(got); free(t);
    eval_is("$IgnoreEOF = False;", "Null");
}

/* ------------------------------------------------------------------ */
/* messages                                                            */
/* ------------------------------------------------------------------ */

static void test_bad_arguments_stay_unevaluated(void) {
    eval_is("Quiet[Dialog[1, 2]]", "Dialog[1, 2]");
    eval_is("Quiet[Dialog[1, Foo -> 2]]", "Dialog[1, Foo -> 2]");
    eval_is("Quiet[Dialog[DialogSymbols -> x]]", "Dialog[DialogSymbols -> x]");
    eval_is("Quiet[Dialog[DialogSymbols :> {1}]]", "Dialog[DialogSymbols :> {1}]");
    eval_is("Quiet[Dialog[DialogSymbols :> {f[x]}]]", "Dialog[DialogSymbols :> {f[x]}]");
    eval_is("Quiet[Dialog[DialogSymbols :> {a, a = 2}]]", "Dialog[DialogSymbols :> {a, a = 2}]");
    eval_is("Check[Dialog[1, 2], bad]", "bad");
    eval_is("Check[Dialog[1, Foo -> 2], bad]", "bad");
    mth_message_list_reset();
    eval_is("1/0; Dialog[1, 2]; Dialog[1, Foo -> 2]; Dialog[DialogSymbols -> x]; "
            "Dialog[DialogSymbols :> {1}]; Dialog[DialogSymbols :> {a, a}]; $MessageList",
            "{Power::infy, Dialog::nonopt, Dialog::optx, Dialog::lvlist, Dialog::lvsym, Dialog::dup}");
    ASSERT(mth_dialog_depth() == 0);
}

/* A TimeConstrained timeout jumps out of the dialog: the Block-frame unwind
 * restores its symbols and $Line, and the registered counter its depth. */
static void test_timeout_through_a_dialog_restores_state(void) {
    at_line_4();
    eval_is("tt = 0; TimeConstrained[Dialog[DialogSymbols :> {tt = 1}, "
            "DialogProlog :> While[True]], 0.3]", "$Aborted");
    eval_is("{tt, $Line}", "{0, 4}");
    ASSERT(mth_dialog_depth() == 0);
    const char* in[] = { "Return[1]" };
    dialog_is("Dialog[]", in, N(in), "1", " In[5]:= Return[1]\n");   /* one level, not two */
    eval_is("Clear[tt]", "Null");
}

/* Last: a Quit[] request cannot be withdrawn. */
static void test_quit_ends_the_dialog(void) {
    at_line_4();
    const char* in[] = { "Quit[]", "never = 1" };
    char* t = NULL;
    char* got = run_dialog("Dialog[3]", in, N(in), false, &t);
    ASSERT(mth_quit_requested());
    ASSERT(strcmp(got, "3") == 0);
    ASSERT(!contains(t, "never"));
    free(got); free(t);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_attributes_and_options);
    TEST(test_no_reader_returns_at_once);
    TEST(test_return_with_value);
    TEST(test_return_without_value_gives_last_line);
    TEST(test_end_of_input_gives_last_line);
    TEST(test_semicolon_line_is_the_value_but_not_shown);
    TEST(test_numbering_follows_line_and_reverts);
    TEST(test_argument_line_counts_even_when_null);
    TEST(test_blank_and_comment_lines_take_no_number);
    TEST(test_parse_error_is_reported_and_dialog_continues);
    TEST(test_only_an_input_return_exits);
    TEST(test_throw_abort_break_stay_inside);
    TEST(test_loop_around_dialog_keeps_running);
    TEST(test_throw_from_prolog_or_argument_propagates);
    TEST(test_dialog_symbols_and_prolog_order);
    TEST(test_option_list_and_rule_forms);
    TEST(test_first_option_setting_wins);
    TEST(test_options_follow_set_options);
    TEST(test_epilog_is_local_and_runs_on_exit);
    TEST(test_message_list_is_local);
    TEST(test_session_hooks_apply);
    TEST(test_nested_dialogs_indent_and_return);
    TEST(test_ignore_eof);
    TEST(test_bad_arguments_stay_unevaluated);
    TEST(test_timeout_through_a_dialog_restores_state);
    TEST(test_quit_ends_the_dialog);

    printf("All Dialog tests passed.\n");
    return 0;
}
