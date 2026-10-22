/*
 * simp_sf_test_util.h -- shared checks for the special-function Simplify
 * suites (tests/test_simp_sf_*.c, src/simp/simp_sf_engine.c).
 *
 * Every expectation in those suites is DERIVED from a defining identity
 * (series, integral, recurrence, ODE; see docs/design/simp_special_functions.md),
 * never transcribed from another system -- and every one is ALSO checked
 * numerically, so a wrongly derived identity fails rather than being
 * enshrined:
 *
 *   sf_zero(e)          Simplify[e] prints "0", and e vanishes numerically.
 *   sf_to(e, r)         Simplify[e] prints r, and e - r vanishes numerically.
 *   sf_sound(e)         Simplify[e] - e vanishes numerically (form free).
 *   sf_unchanged(e)     Simplify[e] === e (no identity applies).
 *   sf_identity(e)      e vanishes numerically (no Simplify involved): pins
 *                       the truth of a relation the engine encodes.
 *   *_assuming(...)     the same under Simplify[e, assum], with the numeric
 *                       check at caller-supplied points that satisfy assum.
 *   sf_fast(e, secs)    Simplify[e] finishes within secs (performance guard).
 *
 * Each of those is one CASE, counted toward the family opened by
 * sf_family(name); sf_finish() requires >= SF_MIN_CASES cases per family.
 *
 * Numeric checks run at three points: two real, one complex (SF_PTS_*),
 * to `g_sf_prec` digits, chopped at 10^-(g_sf_prec - 10).  The complex
 * point catches identities that only hold on the real line.
 *
 * Failures do NOT abort: each is reported with its input, and sf_finish()
 * returns non-zero when any check failed (immune to -DNDEBUG, unlike libc
 * assert).  sf_finish() also enforces a minimum check count per suite.
 */
#pragma once
#include "test_utils.h"
#include <time.h>

static int g_sf_checks = 0;
static int g_sf_fails  = 0;
static int g_sf_prec   = 30;

/* Per-family case counts: sf_family("Beta") opens a family; every public
 * check (one CASE: an expression and its numeric verification) counts once
 * toward the open family, and sf_finish() fails any family with fewer than
 * SF_MIN_CASES cases. */
#define SF_MIN_CASES 50
static struct { const char* name; int cases; } g_sf_fam[64];
static int g_sf_nfam = 0;
static int g_sf_cur  = -1;

/* The harness constructor arms a 120 s alarm for the whole binary (hang
 * protection).  These suites run hundreds of cases per family, so the alarm
 * is re-armed per family instead: the hang guard becomes 120 s per family.
 * stdout is line-buffered so a timeout still shows the last family run. */
static void sf_family(const char* name) {
    static int buffered = 0;
    if (!buffered) { setvbuf(stdout, NULL, _IOLBF, 0); buffered = 1; }
    alarm(120);
    printf("  [family %s]\n", name);
    for (int i = 0; i < g_sf_nfam; i++)
        if (strcmp(g_sf_fam[i].name, name) == 0) { g_sf_cur = i; return; }
    g_sf_fam[g_sf_nfam].name = name;
    g_sf_fam[g_sf_nfam].cases = 0;
    g_sf_cur = g_sf_nfam++;
}

static void sf_case(void) { if (g_sf_cur >= 0) g_sf_fam[g_sf_cur].cases++; }

/* Two real points and one complex point.  Values are chosen inside every
 * family's domain of convergence / principal branch used by the suites
 * (|z| < 1 for 2F1, 0 < m < 1 for elliptic parameters, Re > 0 for
 * Gamma-type arguments). */
#define SF_PTS_1 "{x -> 17/10, y -> 29/10, z -> 13/20, w -> 7/4, n -> 3/7, m -> 2/5, " \
                 "a -> 6/5, b -> 7/3, c -> 19/4, k -> 5/11, s -> 9/4, q -> 1/3, t -> 7/10, nu -> 4/9}"
#define SF_PTS_2 "{x -> 37/10, y -> 11/10, z -> 3/10, w -> 5/2, n -> 5/9, m -> 1/7, " \
                 "a -> 2/3, b -> 13/5, c -> 31/7, k -> 2/9, s -> 11/3, q -> 2/7, t -> 1/5, nu -> 7/8}"
#define SF_PTS_3 "{x -> 9/10 + 2/5 I, y -> 7/5 - 1/3 I, z -> 1/4 + 1/5 I, w -> 3/2 + 1/3 I, " \
                 "n -> 2/7 + 1/9 I, m -> 1/3 + 1/7 I, a -> 4/3 + 1/4 I, b -> 5/2 - 1/5 I, " \
                 "c -> 23/6 + 1/3 I, k -> 1/5 + 1/6 I, s -> 5/2 + 1/2 I, q -> 1/4 + 1/8 I, " \
                 "t -> 3/10 + 1/10 I, nu -> 2/3 + 1/4 I}"

static char* sf_eval_str(const char* input) {
    struct Expr* p = parse_expression(input);
    if (!p) {
        char* e = malloc(16);
        memcpy(e, "<parse error>", 14);
        return e;
    }
    struct Expr* r = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(r);
    expr_free(r);
    return s;
}

static void sf_fail(const char* what, const char* input, const char* expected, const char* actual) {
    g_sf_fails++;
    fprintf(stderr, "FAIL [%s]: %s\n  Expected: %s\n  Actual:   %s\n", what, input, expected, actual);
}

/* Check that the printed value of `input` is `expected`. */
static void sf_expect(const char* what, const char* input, const char* expected) {
    g_sf_checks++;
    char* got = sf_eval_str(input);
    if (strcmp(got, expected) != 0) sf_fail(what, input, expected, got);
    free(got);
}

/* `e` vanishes numerically at every point of the rule-list list `pts`
 * ("{{rules1}, {rules2}, ...}" or a single "{rules}"). */
/* The check is RELATIVE to the largest top-level term: Erfi[x^2]^2 at
 * x = 37/10 is ~1e109, and an absolute 30-digit Chop would report its
 * roundoff as a failure. */
static void sf_numeric_zero_at(const char* what, const char* e, const char** pts, size_t npts) {
    char buf[16384];
    for (size_t i = 0; i < npts; i++) {
        snprintf(buf, sizeof buf,
                 "Module[{sfE = (%s) /. %s, sfT}, sfT = If[Head[sfE] === Plus, List @@ sfE, {sfE}]; "
                 "Chop[N[sfE, %d]/Max[1, Max[Abs[N[sfT, %d]]]], 10^-%d]]",
                 e, pts[i], g_sf_prec, g_sf_prec, g_sf_prec - 10);
        sf_expect(what, buf, "0");
    }
}

static const char* SF_DEFAULT_PTS[] = { SF_PTS_1, SF_PTS_2, SF_PTS_3 };
/* real points only, for identities that hold on (part of) the real line */
static const char* SF_REAL_PTS[] = { SF_PTS_1, SF_PTS_2 };

static void sf_numeric_zero(const char* what, const char* e) {
    sf_numeric_zero_at(what, e, SF_DEFAULT_PTS, 3);
}

static void sf_identity(const char* e) { sf_case(); sf_numeric_zero("identity", e); }
static void sf_identity_real(const char* e) { sf_case(); sf_numeric_zero_at("identity-real", e, SF_REAL_PTS, 2); }

static void sf_zero(const char* e) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    sf_expect("zero", buf, "0");
    sf_numeric_zero("zero-numeric", e);
}

/* Real-line-only variant (identity valid for real arguments in the sample
 * range but not on the complex point, e.g. a branch-cut identity). */
static void sf_zero_real(const char* e) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    sf_expect("zero", buf, "0");
    sf_numeric_zero_at("zero-numeric", e, SF_REAL_PTS, 2);
}

static void sf_to(const char* e, const char* expected) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    sf_expect("to", buf, expected);
    snprintf(buf, sizeof buf, "(%s) - (%s)", e, expected);
    sf_numeric_zero("to-numeric", buf);
}

static void sf_to_real(const char* e, const char* expected) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    sf_expect("to", buf, expected);
    snprintf(buf, sizeof buf, "(%s) - (%s)", e, expected);
    sf_numeric_zero_at("to-numeric", buf, SF_REAL_PTS, 2);
}

/* Simplify[e] == 0 (no assumptions), numeric check at caller points (a
 * family whose identities were verified only on part of the plane, e.g.
 * LegendreQ's Casoratian inside (-1, 1) and off the real axis). */
static void sf_zero_at(const char* e, const char** pts, size_t npts) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    sf_expect("zero", buf, "0");
    sf_numeric_zero_at("zero-numeric", e, pts, npts);
}

static void sf_sound(const char* e) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s] - (%s)", e, e);
    sf_numeric_zero("sound", buf);
}

static void sf_unchanged(const char* e) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s] === (%s)", e, e);
    sf_expect("unchanged", buf, "True");
}

/* Under assumptions; numeric check at the given points (which must satisfy
 * the assumptions). */
static void sf_zero_assuming(const char* e, const char* assum, const char** pts, size_t npts) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s, %s]", e, assum);
    sf_expect("zero-assuming", buf, "0");
    sf_numeric_zero_at("zero-assuming-numeric", e, pts, npts);
}

static void sf_to_assuming(const char* e, const char* assum, const char* expected,
                           const char** pts, size_t npts) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s, %s]", e, assum);
    sf_expect("to-assuming", buf, expected);
    snprintf(buf, sizeof buf, "(%s) - (%s)", e, expected);
    sf_numeric_zero_at("to-assuming-numeric", buf, pts, npts);
}

static void sf_unchanged_assuming(const char* e, const char* assum) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s, %s] === (%s)", e, assum, e);
    sf_expect("unchanged-assuming", buf, "True");
}

/* Performance guard: Simplify[e] must finish within `secs` wall seconds. */
static void sf_fast(const char* e, double secs) {
    sf_case();
    char buf[16384];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    g_sf_checks++;
    struct Expr* p = parse_expression(buf);
    clock_t t0 = clock();
    struct Expr* r = evaluate(p);
    double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
    expr_free(p);
    expr_free(r);
    if (dt > secs) {
        char got[64];
        snprintf(got, sizeof got, "%.3f s", dt);
        char lim[64];
        snprintf(lim, sizeof lim, "<= %.3f s", secs);
        sf_fail("fast", buf, lim, got);
    }
}

/* Final report.  `min_checks` is the suite's floor on the number of checks
 * (so a suite that silently loses its tests also fails). */
static int sf_finish(const char* suite, int min_checks) {
    for (int i = 0; i < g_sf_nfam; i++) {
        printf("  family %-18s %4d cases\n", g_sf_fam[i].name, g_sf_fam[i].cases);
        if (g_sf_fam[i].cases < SF_MIN_CASES) {
            fprintf(stderr, "FAIL [%s]: family %s has %d cases (floor %d)\n", suite,
                    g_sf_fam[i].name, g_sf_fam[i].cases, SF_MIN_CASES);
            g_sf_fails++;
        }
    }
    if (g_sf_checks < min_checks) {
        fprintf(stderr, "FAIL [%s]: only %d checks ran (floor %d)\n", suite, g_sf_checks, min_checks);
        g_sf_fails++;
    }
    if (g_sf_fails) {
        fprintf(stderr, "%s: %d of %d checks FAILED\n", suite, g_sf_fails, g_sf_checks);
        return 1;
    }
    printf("All %s tests passed! (%d checks)\n", suite, g_sf_checks);
    return 0;
}
