/* test_principal_value_corpus.c
 *
 * Corpus runner for the Cauchy principal-value definite-integral stress suite
 * (tests/PrincipalValueIntegralsCorpus.m, battery
 * `extremeComplexSingularIntegralStressTest`): fractional branch cuts over real
 * poles, log branch cuts, trig poles, oscillatory-over-complex-singularity, and
 * fractional powers over higher-order poles -- all as principal values.
 *
 * Each corpus record is  {Hold[rawIntegral], Hold[expected], {param -> val, ...}}.
 * For every record the runner (in a forked child with a per-case timeout):
 *
 *   1. r   = ReleaseHold[rawIntegral]          (evaluate the definite integral)
 *   2. if r still contains an Integrate head  -> UNEVALUATED
 *   3. exp = ReleaseHold[expected]
 *   4. d   = N[(r - exp) /. instantiation]     (instantiate free params, numericize)
 *   5. |d| < 1e-6  ->  CORRECT   else  ->  WRONG
 *
 * The numeric-instantiation check is the authoritative correctness gate: it is
 * independent of whether Simplify can collapse the symbolic form (so an
 * unreduced-but-correct answer still passes).  ConditionalExpression wrappers are
 * stripped to their value before the near-zero test (as in the extreme runner).
 *
 * EXCLUDED cases (corpus_is_excluded): no finite oracle exists -- 11/14/15/22/28
 * diverge, 21/29/30 have infinitely many poles (deferred mechanism), 41/42/43 are
 * convention-dependent double-pole Hadamard finite parts.  See the corpus header
 * for the per-case justification and the full oracle-correction list.
 *
 * The hard regression bar is WRONG: a finite value that is not the (verified)
 * expected value.  The baseline ratchets DOWN as the campaign closes each case;
 * it starts at 4 (cases 5,6,49,50 currently return a wrong/non-finite value --
 * two b!=1 keyhole-PV imaginary leaks and two Log-Mellin 0/0 blow-ups).
 *
 * Modeled on test_extreme_definite_corpus.c.
 */

#include "test_utils.h"
#include "expr.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>

void symtab_init(void);
void core_init(void);

extern char* expr_to_string(struct Expr*);

/* 1-based corpus indices excluded from scoring (never CORRECT, never WRONG):
 *   11,14,15  Log / x Log over (x^n-1): DIVERGE (tail ~ Log[x]/x).
 *   22,28     Cot[ax]/x, Cos[ax]/(x(x^2-b^2)): DIVERGE at x=0.
 *   21,29,30  Tan[ax]/x, Tan/Cot[ax]/(x^2-b^2): infinitely many poles -- deferred.
 *   41,42,43  x^p/(x^n-1)^2: double axis pole -> Cauchy PV diverges (~2A/eps); only
 *             a regularization-dependent Hadamard finite part exists. */
static int corpus_is_excluded(size_t one_based_index) {
    switch (one_based_index) {
        case 11: case 14: case 15:
        case 21: case 22: case 28: case 29: case 30:
        case 41: case 42: case 43:
            return 1;
        default:
            return 0;
    }
}

static bool expr_contains_head_named(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head
            && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol.name, name) == 0) {
            return true;
        }
        if (expr_contains_head_named(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (expr_contains_head_named(e->data.function.args[i], name)) return true;
        }
    }
    return false;
}

/* True when an N[..] value is a finite numeric within tol of zero -- a real/
 * integer node near 0, or a Complex[..] with both parts near 0.  A
 * ConditionalExpression[v, cond] is unwrapped to its value v first. */
static bool n_value_near_zero(Expr* nv, double tol) {
    if (!nv) return false;
    if (nv->type == EXPR_FUNCTION
        && nv->data.function.head
        && nv->data.function.head->type == EXPR_SYMBOL
        && strcmp(nv->data.function.head->data.symbol.name, "ConditionalExpression") == 0
        && nv->data.function.arg_count >= 1) {
        return n_value_near_zero(nv->data.function.args[0], tol);
    }
    if (nv->type == EXPR_REAL) {
        double v = nv->data.real;
        return v >= -tol && v <= tol;
    }
    if (nv->type == EXPR_INTEGER) return nv->data.integer == 0;
    if (nv->type == EXPR_FUNCTION
        && nv->data.function.head
        && nv->data.function.head->type == EXPR_SYMBOL
        && strcmp(nv->data.function.head->data.symbol.name, "Complex") == 0
        && nv->data.function.arg_count == 2) {
        Expr* re = nv->data.function.args[0];
        Expr* im = nv->data.function.args[1];
        if ((re->type != EXPR_REAL && re->type != EXPR_INTEGER) ||
            (im->type != EXPR_REAL && im->type != EXPR_INTEGER)) return false;
        double rv = (re->type == EXPR_REAL) ? re->data.real : (double)re->data.integer;
        double iv = (im->type == EXPR_REAL) ? im->data.real : (double)im->data.integer;
        return rv >= -tol && rv <= tol && iv >= -tol && iv <= tol;
    }
    return false;
}

#define CORPUS_PER_CASE_TIMEOUT_SEC 20

enum {
    CORPUS_CORRECT     = 'Z',
    CORPUS_WRONG       = 'D',
    CORPUS_UNEVALUATED = 'U',
    CORPUS_EXCLUDED    = 'X',
};

/* Build ReleaseHold[copy(held)] and evaluate it. Caller owns the result. */
static Expr* release_and_eval(Expr* held) {
    Expr* rh = expr_new_function(expr_new_symbol("ReleaseHold"),
        (Expr*[]){ expr_copy(held) }, 1);
    Expr* r = evaluate(rh);
    expr_free(rh);
    return r;
}

int main(int argc, char** argv) {
    alarm(0);
    setvbuf(stderr, NULL, _IONBF, 0);

    const char* corpus_file =
        (argc > 1) ? argv[1] : "../PrincipalValueIntegralsCorpus.m";

    symtab_init();
    core_init();

    fprintf(stderr, "==> corpus file: %s\n", corpus_file);

    char get_buf[1024];
    snprintf(get_buf, sizeof(get_buf), "Get[\"%s\"]", corpus_file);
    Expr* get_call = parse_expression(get_buf);
    ASSERT(get_call != NULL);
    Expr* tests = evaluate(get_call);
    expr_free(get_call);

    if (!tests || tests->type != EXPR_FUNCTION
        || !tests->data.function.head
        || tests->data.function.head->type != EXPR_SYMBOL
        || strcmp(tests->data.function.head->data.symbol.name, "List") != 0) {
        fprintf(stderr, "FAIL: could not load %s as a List.\n", corpus_file);
        if (tests) expr_free(tests);
        ASSERT(false);
        return 1;
    }

    size_t n = tests->data.function.arg_count;
    fprintf(stderr, "==> principal-value definite-integral corpus: %zu cases\n", n);

    int correct = 0, wrong = 0, unevaluated = 0, malformed = 0,
        timed_out = 0, crashed = 0, excluded = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* rec = tests->data.function.args[i];
        if (rec->type != EXPR_FUNCTION
            || !rec->data.function.head
            || rec->data.function.head->type != EXPR_SYMBOL
            || strcmp(rec->data.function.head->data.symbol.name, "List") != 0
            || rec->data.function.arg_count < 3) {
            malformed++;
            fprintf(stderr, "  [%3zu/%zu] MALFORMED\n", i + 1, n);
            continue;
        }
        Expr* held_int = rec->data.function.args[0];
        Expr* held_exp = rec->data.function.args[1];
        Expr* inst     = rec->data.function.args[2];

        char* in_str = expr_to_string(held_int);
        fprintf(stderr, "  [%3zu/%zu] %s\n", i + 1, n, in_str);
        free(in_str);

        int pipefd[2];
        if (pipe(pipefd) != 0) { crashed++; continue; }
        pid_t pid = fork();
        if (pid < 0) { close(pipefd[0]); close(pipefd[1]); crashed++; continue; }

        if (pid == 0) {
            /* --- Child --- */
            close(pipefd[0]);
            char code = CORPUS_UNEVALUATED;
            alarm(CORPUS_PER_CASE_TIMEOUT_SEC);

            if (corpus_is_excluded(i + 1)) {
                code = CORPUS_EXCLUDED;
                ssize_t w0 = write(pipefd[1], &code, 1); (void)w0;
                close(pipefd[1]);
                _exit(0);
            }

            Expr* r = release_and_eval(held_int);

            if (expr_contains_head_named(r, "Integrate")) {
                code = CORPUS_UNEVALUATED;
            } else {
                Expr* exp = release_and_eval(held_exp);
                Expr* neg = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), exp }, 2);
                Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(r), neg }, 2);
                Expr* sub = expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ diff, expr_copy(inst) }, 2);
                Expr* strip_rule =
                    parse_expression("ConditionalExpression[cev_, _] :> cev");
                Expr* stripped = expr_new_function(expr_new_symbol("ReplaceRepeated"),
                    (Expr*[]){ sub, strip_rule }, 2);
                Expr* nv_call = expr_new_function(expr_new_symbol("N"),
                    (Expr*[]){ stripped }, 1);
                Expr* nv = evaluate(nv_call);
                expr_free(nv_call);
                code = n_value_near_zero(nv, 1e-6) ? CORPUS_CORRECT : CORPUS_WRONG;
                expr_free(nv);
            }
            alarm(0);
            if (r) expr_free(r);
            ssize_t w = write(pipefd[1], &code, 1); (void)w;
            close(pipefd[1]);
            _exit(0);
        }

        /* --- Parent --- */
        close(pipefd[1]);
        time_t deadline = time(NULL) + CORPUS_PER_CASE_TIMEOUT_SEC + 2;
        int status = 0;
        bool reaped = false;
        while (time(NULL) < deadline) {
            pid_t rr = waitpid(pid, &status, WNOHANG);
            if (rr == pid) { reaped = true; break; }
            if (rr < 0) break;
            usleep(50000);
        }
        if (!reaped) {
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            close(pipefd[0]);
            timed_out++;
            fprintf(stderr, "      -> TIMEOUT (>%ds)\n", CORPUS_PER_CASE_TIMEOUT_SEC);
            continue;
        }

        char code = 0;
        ssize_t rd = read(pipefd[0], &code, 1);
        close(pipefd[0]);

        if (rd != 1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) {
                timed_out++;
                fprintf(stderr, "      -> TIMEOUT (child alarm)\n");
                continue;
            }
            crashed++;
            fprintf(stderr, "      -> CRASH (exit %d, sig %d)\n",
                    WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                    WIFSIGNALED(status) ? WTERMSIG(status)  : 0);
            continue;
        }

        switch (code) {
            case CORPUS_CORRECT: correct++;
                fprintf(stderr, "      -> CORRECT\n");
                break;
            case CORPUS_WRONG: {
                wrong++;
                fprintf(stderr, "      -> WRONG (case %zu: finite value != expected)\n", i + 1);
                break;
            }
            case CORPUS_UNEVALUATED:
                unevaluated++;
                fprintf(stderr, "      -> UNEVALUATED\n");
                break;
            case CORPUS_EXCLUDED:
                excluded++;
                fprintf(stderr, "      -> EXCLUDED (no finite oracle)\n");
                break;
            default:
                crashed++;
                break;
        }
    }
    expr_free(tests);

    fprintf(stderr, "\n=== principal-value definite-integral corpus result ===\n");
    fprintf(stderr, "  Total cases:                  %zu\n", n);
    fprintf(stderr, "  Correct (value verified):     %d  %.1f%%\n",
            correct, n ? 100.0 * correct / n : 0.0);
    fprintf(stderr, "  WRONG (finite != expected):   %d\n", wrong);
    fprintf(stderr, "  Unevaluated:                  %d\n", unevaluated);
    fprintf(stderr, "  Timed out (>%d s/case):       %d\n",
            CORPUS_PER_CASE_TIMEOUT_SEC, timed_out);
    fprintf(stderr, "  Excluded (no finite oracle):  %d\n", excluded);
    fprintf(stderr, "  Crashed (child fault):        %d\n", crashed);
    fprintf(stderr, "  Malformed (test file issue):  %d\n", malformed);
    fprintf(stderr, "=======================================================\n");

    /* WRONG is the hard regression bar.  Ratchet DOWN as each case is fixed; never
     * raise without confirming a real new wrong answer.
     *   v0.406 start: 4 -- cases 5,6 (b!=1 keyhole-PV imaginary leak) and 49,50
     *   (Log-Mellin 0/0 -> ComplexInfinity/Indeterminate).  All four fixed: 5,6 via
     *   the keyhole-PV real projection, 49,50 via the Log-weighted Mellin handler.
     *   Now assert-empty. */
    const int CORPUS_WRONG_BASELINE = 0;
    if (wrong > CORPUS_WRONG_BASELINE) {
        fprintf(stderr,
            "FAIL: %d case(s) returned a finite WRONG value (above baseline %d).\n",
            wrong, CORPUS_WRONG_BASELINE);
        return 1;
    }

    return 0;
}
