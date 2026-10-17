/* flint_factor.c -- raw FLINT factoring driver for the FactorInteger shootout.
 *
 * Reads lines "<id> <n>" from stdin.  Mode (argv[1]):
 *   batch   time the whole input with fmpz_factor, print "#TOTAL\t<sec>" then
 *           one "<id>\t\t<p e p e ...>" line per case;
 *   single  time each case separately, "<id>\t<sec>\t<p e ...>";
 * Method (argv[2], optional): "auto" (fmpz_factor, default), "qs"
 * (qsieve_factor), "ecm:<B1>:<curves>" (fmpz_factor_ecm, one split), and for
 * word-sized n the ulong splitters "squfof" (n_factor_SQUFOF), "lehman"
 * (n_factor_lehman), "pp1" (n_factor_pp1_wrapper).
 *
 * Build: cc -O2 -o flint_factor flint_factor.c $(pkg-config --cflags --libs flint) -lgmp
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/qsieve.h>
#include <flint/ulong_extras.h>

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

typedef struct { char id[128]; fmpz_t n; fmpz_factor_t f; int ok; } Case;

static void run(Case* c, const char* method, flint_rand_t st) {
    if (strcmp(method, "qs") == 0) {
        qsieve_factor(c->f, c->n);
    } else if (strncmp(method, "ecm:", 4) == 0) {
        unsigned long B1 = 0, curves = 0;
        sscanf(method + 4, "%lu:%lu", &B1, &curves);
        fmpz_t d; fmpz_init(d);
        if (fmpz_factor_ecm(d, curves, B1, 100 * B1, st, c->n)) {
            _fmpz_factor_append(c->f, d, 1);
            fmpz_divexact(d, c->n, d);
            _fmpz_factor_append(c->f, d, 1);
        } else {
            _fmpz_factor_append(c->f, c->n, 1);
        }
        fmpz_clear(d);
    } else if (!strcmp(method, "squfof") || !strcmp(method, "lehman") || !strcmp(method, "pp1")) {
        ulong n = fmpz_get_ui(c->n), d = 0;
        if (!strcmp(method, "squfof")) d = n_factor_SQUFOF(n, 1 << 20);
        else if (!strcmp(method, "lehman")) d = n_factor_lehman(n);
        else d = n_factor_pp1_wrapper(n);
        if (d > 1 && d < n) {
            _fmpz_factor_append_ui(c->f, d, 1);
            _fmpz_factor_append_ui(c->f, n / d, 1);
        } else {
            _fmpz_factor_append(c->f, c->n, 1);
        }
    } else {
        fmpz_factor(c->f, c->n);
    }
}

static void print_factors(const Case* c) {
    for (slong i = 0; i < c->f->num; i++) {
        if (i) putchar(' ');
        fmpz_print(c->f->p + i);
        printf(" %lu", (unsigned long)c->f->exp[i]);
    }
    putchar('\n');
}

int main(int argc, char** argv) {
    const char* mode = argc > 1 ? argv[1] : "single";
    const char* method = argc > 2 ? argv[2] : "auto";
    size_t cap = 4096, nc = 0;
    Case* cs = malloc(cap * sizeof(Case));
    char id[128];
    static char num[1 << 16];
    while (scanf("%127s %65535s", id, num) == 2) {
        if (nc == cap) { cap *= 2; cs = realloc(cs, cap * sizeof(Case)); }
        strcpy(cs[nc].id, id);
        fmpz_init(cs[nc].n);
        fmpz_set_str(cs[nc].n, num, 10);
        fmpz_factor_init(cs[nc].f);
        nc++;
    }
    flint_rand_t st;
    flint_rand_init(st);
    if (strcmp(mode, "batch") == 0) {
        double t0 = now();
        for (size_t i = 0; i < nc; i++) run(&cs[i], method, st);
        printf("#TOTAL\t%.9f\n", now() - t0);
        for (size_t i = 0; i < nc; i++) { printf("%s\t\t", cs[i].id); print_factors(&cs[i]); }
    } else {
        for (size_t i = 0; i < nc; i++) {
            double t0 = now();
            run(&cs[i], method, st);
            double t = now() - t0;
            printf("%s\t%.9f\t", cs[i].id, t);
            print_factors(&cs[i]);
            fflush(stdout);
        }
    }
    flint_rand_clear(st);
    return 0;
}
