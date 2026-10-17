/* Tests for MeijerG (src/special_functions/meijer_g.c), its machine Slater
 * engine (sf_machine_meijerg.c), the arb engine (flint_num_bridge.c), the
 * Compile[] / NDArray / packed surfaces, D, Series and FunctionExpand.
 *
 * Every reference value below was computed with Mathematica 14 (exact inputs
 * via SetPrecision[x, Infinity], N[..., 30+]) -- not recalled -- and the
 * machine values are compared to a relative tolerance of 2e-14.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "ndarray.h"
#include "arithmetic.h"
#include "attr.h"
#include "special_functions/sf_machine_meijerg.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- helpers --------------------------------------------------------- */

/* Machine values out of a result: Real, Complex[Real, Real], or a List /
 * NDArray of those. Returns the count found. */
static int collect(const Expr* e, double* re, double* im, int cap) {
    if (!e || cap <= 0) return 0;
    if (e->type == EXPR_REAL) { re[0] = e->data.real; im[0] = 0.0; return 1; }
    Expr *a, *b;
    if (is_complex((Expr*)e, &a, &b) && a->type == EXPR_REAL && b->type == EXPR_REAL) {
        re[0] = a->data.real; im[0] = b->data.real; return 1;
    }
    if (is_ndarray(e)) {
        Expr* l = ndarray_to_nested_list(e);
        int n = collect(l, re, im, cap);
        expr_free(l);
        return n;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "List") == 0) {
        int n = 0;
        for (size_t i = 0; i < e->data.function.arg_count && n < cap; i++)
            n += collect(e->data.function.args[i], re + n, im + n, cap - n);
        return n;
    }
    return 0;
}

static int eval_vals(const char* input, double* re, double* im, int cap) {
    Expr* p = parse_expression(input);
    ASSERT_MSG(p != NULL, "parse failed: %s", input);
    Expr* r = evaluate(p);
    expr_free(p);
    int n = collect(r, re, im, cap);
    expr_free(r);
    return n;
}

/* The printed result of evaluating `input`. Caller frees. */
static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT_MSG(p != NULL, "parse failed: %s", input);
    Expr* r = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(r);
    expr_free(r);
    return s;
}

/* `input` must evaluate to True. */
static void assert_true(const char* input) {
    char* s = eval_str(input);
    ASSERT_MSG(strcmp(s, "True") == 0, "%s -> %s", input, s);
    free(s);
}

static void mlit(char* buf, size_t sz, double x) {
    char t[64];
    snprintf(t, sizeof t, "%.17g", x);
    char* e = strchr(t, 'e');
    if (e) { *e = '\0'; snprintf(buf, sz, "%s`*^%d", t, atoi(e + 1)); }
    else   snprintf(buf, sz, "%s`", t);
}

/* ---- machine spot values vs Mathematica -------------------------------- */

typedef struct { const char* expr; double re, im; } Ref;
static const Ref REFS[] = {
    /* the documented examples */
    { "MeijerG[{{1},{2}},{{2.3,1.5},{4}},3.]", -1.8225742114495970, 0 },
    { "MeijerG[{{31,4},{2,3}},{{5,1.5},{.1}},3]", -1.2507804606047746e-28, 0 },
    { "MeijerG[{{3+I},{2.}},{{2/3,1/5},{4}},3]", -0.66457769279907138, -3.0296642161867087 },
    /* p < q: G^{1,1}_{1,2}, G^{2,1}_{2,3} small and large z */
    { "MeijerG[{{1/3},{}},{{1/2},{1}},0.7]", -0.0035085068168862805, 0 },
    { "MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},0.7]", 0.64569576912679501, 0 },
    { "MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},25.]", 0.029162489315014610, 0 },
    { "MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},0.5+1.2 I]", 0.53463815758638143, -0.53278087315953023 },
    /* p == q: both sides of the unit circle, and on it */
    { "MeijerG[{{1/3,1/2},{}},{{1/5,1/7},{}},0.4]", 3.0471803443053128, 0 },
    { "MeijerG[{{1/3,1/2},{}},{{1/5,1/7},{}},3.]", 2.0087816470898505, 0 },
    { "MeijerG[{{1/2,1/3},{}},{{1/5,1/4},{}},1.]", 2.0569165030652114, 0 },
    { "MeijerG[{{1/3,2/3},{1/2}},{{0.1},{0.2,0.3}},0.6]", 1.0805783965590123, 0 },
    { "MeijerG[{{1/3,2/3},{1/2}},{{0.1},{0.2,0.3}},1.6]", 1.1266218312991154, 0 },
    /* p > q */
    { "MeijerG[{{1/2,1},{3,4}},{{1/5},{}},2.5]", 0.61101996069060312, 0 },
    /* the log case (integer-spaced b's) */
    { "MeijerG[{{},{2}},{{1/2,3/2},{}},3.]", 0.041591159875396984, 0 },
    { "MeijerG[{{},{}},{{0,0,1},{}},1.5]", 0.13633421667032781, 0 },
    { "MeijerG[{{},{0.5}},{{0.25,0.25},{}},2.0]", 0.13209146198011134, 0 },
    { "MeijerG[{{0.5},{}},{{0.0},{0.0}},7.5]", 0.38011522630654138, 0 },
    /* negative and complex arguments, complex values */
    { "MeijerG[{{1/3},{}},{{1/2},{1}},-0.7]", 0, 1.6480274201707194 },
    { "MeijerG[{{0.3},{0.6}},{{0.1,0.2},{0.4}},-2.5]", 6.7601525852185073, -8.6356804073956984 },
    { "MeijerG[{{},{}},{{1.25},{0.5}},-3.0]", -11.976047357515326, -11.976047357515326 },
    /* cancellation between Slater terms (BesselK-like), large parameters */
    { "MeijerG[{{},{}},{{0.25,0.75},{}},10.]", 0.0056474898653363598, 0 },
    { "MeijerG[{{},{}},{{30.5,2.25},{}},4.]", 4.8795090431307002e29, 0 },
    /* the generalized form */
    { "MeijerG[{{1/3},{}},{{1/2},{1}},0.3,1/2]", 0.12611212330424155, 0 },
    { "MeijerG[{{1/2},{}},{{1/3,1/4},{1/5}},2.0,2]", 0.51716307002853507, 0 },
    { "MeijerG[{{1/2},{}},{{1/3,1/4},{1/5}},2.0,-1/3]", 1.1905566604670160, 0 },
};

void test_mg_machine_values() {
    for (size_t i = 0; i < sizeof REFS / sizeof REFS[0]; i++) {
        double re, im;
        ASSERT_MSG(eval_vals(REFS[i].expr, &re, &im, 1) == 1, "no machine value: %s", REFS[i].expr);
        double d = hypot(re - REFS[i].re, im - REFS[i].im);
        double s = hypot(REFS[i].re, REFS[i].im);
        ASSERT_MSG(d <= 2e-14 * s, "%s = %.17g%+.17gi, want %.17g%+.17gi (rel %.3g)",
                   REFS[i].expr, re, im, REFS[i].re, REFS[i].im, d / s);
    }
}

/* ---- arbitrary precision ---------------------------------------------- */

void test_mg_arbitrary_precision() {
    /* N at 50 digits; references from Mathematica at 60. */
    assert_true("Abs[N[MeijerG[{{},{2}},{{1/2,3/2},{}},3],50] - "
                "0.04159115987539698367033039074977875862652971515903948663784`60] < 10^-50");
    assert_true("Abs[N[MeijerG[{{1},{2}},{{2/3,1/5},{4}},3],50] - "
                "(-1.11627350824129439412252452587551932545307762389209666266967`60)] < 10^-48");
    /* The precision of the output tracks the input: 3`100 -> ~100 digits. */
    assert_true("Abs[MeijerG[{{1},{3}},{{2/3,1/5},{4}},3`100] - (-0.2438696021081476521086452421104677"
                "4509093659460150419408877915645122075868755865256614836627675990176208335956`110)] < 10^-97");
    assert_true("Precision[MeijerG[{{1},{3}},{{2/3,1/5},{4}},3`100]] >= 99");
    /* complex parameter, p == q outside the circle, the dual, the log case,
     * a negative z (purely imaginary value), the generalized form */
    assert_true("Abs[N[MeijerG[{{I+3},{2}},{{2/3,1/5},{4}},3],30] - "
                "(-0.66457769279907138064717710520500160030`40 - 3.02966421618670871838178752107873339197`40 I)] < 10^-28");
    assert_true("Abs[N[MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},25],30] - 0.0291624893150146101482229833989580335745`40] < 10^-31");
    assert_true("Abs[N[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},5/2],30] - 0.6110199606906031232416204130638200637143`40] < 10^-29");
    assert_true("Abs[N[MeijerG[{{},{}},{{0,0,1},{}},3/2],30] - 0.1363342166703278065085771760345031828333`40] < 10^-29");
    assert_true("Abs[N[MeijerG[{{1/3},{}},{{1/2},{1}},-7/10],30] - 1.6480274201707195339669908042365607492807`40 I] < 10^-28");
    assert_true("Abs[N[MeijerG[{{1/2},{}},{{1/3,1/4},{1/5}},2,-1/3],30] - 1.1905566604670160436786072879547898995228`40] < 10^-28");
    assert_true("Abs[N[MeijerG[{{1/3,1/2},{}},{{1/5,1/7},{}},3],30] - 2.0087816470898505035902289882105181632582`40] < 10^-28");
    /* a 1000-digit evaluation: its leading digits, and that it is fast */
    assert_true("Abs[MeijerG[{{6},{2}},{{2/3,1/5},{4}},3`1000] - "
                "0.0001083421017540477353652918178612563224185424970489029005`58] < 10^-58");
    assert_true("First[Timing[MeijerG[{{6},{2}},{{2/3,1/5},{4}},3`1000];]] < 2");
    /* a machine Real anywhere makes the answer machine precision */
    assert_true("MachineNumberQ[MeijerG[{{},{2.0`60}},{{0.5,1.5},{}},3]]");
    assert_true("Abs[MeijerG[{{},{2.0`60}},{{1/2,3/2}//N,{}},3] - 0.0415911598753969837] < 10^-15");
}

/* ---- threading and the buffer / Compile surfaces ----------------------- */

/* scalar == unpacked List == packed List == NDArray == Compile, bit for bit. */
static void check_surfaces(const char* head, const char* spec,
                           const double* xr, const double* xi, int n, int cplx) {
    char lits[16][96];
    for (int i = 0; i < n; i++) {
        if (cplx) {
            char a[48], b[48];
            mlit(a, sizeof a, xr[i]); mlit(b, sizeof b, xi[i]);
            snprintf(lits[i], sizeof lits[i], "Complex[%s, %s]", a, b);
        } else mlit(lits[i], sizeof lits[i], xr[i]);
    }
    char list[2048] = "{";
    for (int i = 0; i < n; i++) { strcat(list, lits[i]); strcat(list, i + 1 < n ? ", " : "}"); }
    double sr[16], si[16];
    for (int i = 0; i < n; i++) {
        char e[512], c[1024], body[256];
        snprintf(e, sizeof e, head, lits[i]);
        ASSERT_MSG(eval_vals(e, &sr[i], &si[i], 1) == 1, "scalar %s", e);
        snprintf(body, sizeof body, head, "x");
        snprintf(c, sizeof c, "Compile[{{x, %s}}, %s][%s]", spec, body, lits[i]);
        double cr, ci;
        ASSERT_MSG(eval_vals(c, &cr, &ci, 1) == 1, "compile %s", c);
        ASSERT_MSG(cr == sr[i] && ci == si[i], "%s: Compile %.17g%+.17gi != scalar %.17g%+.17gi",
                   c, cr, ci, sr[i], si[i]);
    }
    const char* wraps[3] = { "FromPackedArray[%s]", "ToPackedArray[%s]", "NDArray[%s]" };
    for (int w = 0; w < 3; w++) {
        char arg[2100], e[2600];
        snprintf(arg, sizeof arg, wraps[w], list);
        snprintf(e, sizeof e, head, arg);
        double vr[16], vi[16];
        int m = eval_vals(e, vr, vi, 16);
        ASSERT_MSG(m == n, "%s: %d values", e, m);
        for (int i = 0; i < n; i++)
            ASSERT_MSG(vr[i] == sr[i] && vi[i] == si[i], "%s [%d]: %.17g%+.17gi != %.17g%+.17gi",
                       e, i, vr[i], vi[i], sr[i], si[i]);
    }
    /* Compile at rank 1: the fused strip loop. */
    {
        char body[256], c[3000];
        snprintf(body, sizeof body, head, "v");
        snprintf(c, sizeof c, "Compile[{{v, %s, 1}}, %s][%s]", spec, body, list);
        double vr[16], vi[16];
        int m = eval_vals(c, vr, vi, 16);
        ASSERT_MSG(m == n, "%s: %d values", c, m);
        for (int i = 0; i < n; i++)
            ASSERT_MSG(vr[i] == sr[i] && vi[i] == si[i], "%s [%d]: %.17g%+.17gi != %.17g%+.17gi",
                       c, i, vr[i], vi[i], sr[i], si[i]);
    }
}

void test_mg_surfaces() {
    double x[] = { 0.05, 0.7, 1.9, 6.25, 31.0 };
    check_surfaces("MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},%s]", "_Real", x, NULL, 5, 0);
    check_surfaces("MeijerG[{{1/3},{}},{{1/2},{1}},%s,1/2]", "_Real", x, NULL, 5, 0);
    /* the log case: every point goes through the arb fallback */
    check_surfaces("MeijerG[{{},{2}},{{1/2,3/2},{}},%s]", "_Real", x, NULL, 5, 0);
    /* p == q, both expansions */
    double y[] = { 0.2, 0.9, 1.0, 1.4, 8.0 };
    check_surfaces("MeijerG[{{1/3,1/2},{}},{{1/5,1/7},{}},%s]", "_Real", y, NULL, 5, 0);
    /* complex z */
    double cr[] = { 0.5, -1.25, 3.0 }, ci[] = { 1.2, 0.5, -2.0 };
    check_surfaces("MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},%s]", "_Complex", cr, ci, 3, 1);
}

void test_mg_compile() {
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, MeijerG[{{1/3},{}},{{1/2},{1}},x]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, MeijerG[{{1/3},{}},{{1/2},{1}},v]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{v, _Complex, 1}}, MeijerG[{{1/3},{}},{{1/2},{1}},v]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}, {a, _Real}}, MeijerG[{{a},{}},{{1/2},{1}},x, 1/3]]) === True");
    /* A real-typed compiled function whose value is complex (z < 0) hands the
     * call to the interpreter rather than dropping the imaginary part. */
    assert_true("Abs[Compile[{{x, _Real}}, MeijerG[{{1},{2}},{{2.3,1.5},{4}},x]][-3.] - "
                "MeijerG[{{1},{2}},{{2.3,1.5},{4}},-3.]] == 0");
    /* the documented threading example */
    assert_true("Max[Abs[MeijerG[{{1/3},{}},{{1/2},{1}},{0.1, 0.3, 0.5},1/2] - "
                "{0.0511285, 0.126112, 0.133903}]] < 10^-6");
    /* the plain flat kernel declines garbage instead of crashing */
    double re[3] = { -1.0, 1.0, 0.5 }, im[3] = { 0, 0, 0 }, rr, ri;
    ASSERT(!sf_mg_eval_flat(re, im, 3, &rr, &ri));
}

/* ---- exact and symbolic evaluation -------------------------------------- */

static void assert_str(const char* input, const char* want) {
    char* s = eval_str(input);
    ASSERT_MSG(strcmp(s, want) == 0, "%s\n  got:  %s\n  want: %s", input, s, want);
    free(s);
}

void test_mg_reductions() {
    assert_str("MeijerG[{{},{}},{{v},{-v}},z]", "BesselJ[2 v, 2 Sqrt[z]]");
    assert_str("MeijerG[{{1/2,1/2,1/2},{}},{{0},{-(3/2),-(3/2)}},-x]",
               "16/9 Sqrt[Pi] HypergeometricPFQ[{1/2, 1/2, 1/2}, {5/2, 5/2}, x]");
    assert_str("MeijerG[{{},{}},{{5},{4}},5]", "625 Sqrt[5] BesselJ[1, 2 Sqrt[5]]");
    assert_str("MeijerG[{{1/2},{}},{{0},{0}},x]", "Sqrt[Pi] BesselI[0, 1/2 x] E^(-1/2 x)");
    assert_str("MeijerG[{{1/2,3/2},{}},{{0},{0}},x]", "-4 EllipticE[-x]");
    assert_str("MeijerG[{{},{}},{{b},{}},-x]", "E^x (-x)^b");
    assert_str("MeijerG[{{a},{}},{{d},{}},x]", "Gamma[1 - a + d] x^d (1 + x)^(-(1 - a + d))");
    assert_str("MeijerG[{{},{}},{{b,c},{}},z]", "2 z^(1/2 (b + c)) BesselK[b - c, 2 Sqrt[z]]");
    assert_str("MeijerG[{{1,1},{}},{{1},{0}},z]", "Log[1 + z]");
    /* stays: two 2F2 terms; p == q with m + n == p; the log case */
    assert_str("MeijerG[{{1},{2}},{{2/3,1/5},{4}},3]", "MeijerG[{{1}, {2}}, {{1/5, 2/3}, {4}}, 3]");
    assert_str("MeijerG[{{},{a}},{{b},{}},z]", "MeijerG[{{}, {a}}, {{b}, {}}, z]");
    assert_str("MeijerG[{{},{2}},{{1/2,3/2},{}},3]", "MeijerG[{{}, {2}}, {{1/2, 3/2}, {}}, 3]");
    /* the reduced forms are the same function: compare against the machine
     * numerics of an unreduced (inexact-parameter) call */
    const char* pairs[][2] = {
        { "MeijerG[{{1/2},{}},{{0},{0}},x]", "MeijerG[{{0.5},{}},{{0.},{0.}},x]" },
        { "MeijerG[{{1/2,3/2},{}},{{0},{0}},x]", "MeijerG[{{0.5,1.5},{}},{{0.},{0.}},x]" },
        { "MeijerG[{{},{}},{{1/3},{1/5}},x]", "MeijerG[{{},{}},{{1/3.},{0.2}},x]" },
        { "MeijerG[{{1/3},{}},{{1/4,1/5},{}},x]", "MeijerG[{{1/3.},{}},{{0.25,0.2},{}},x]" },
        { "MeijerG[{{1/3},{1/7}},{{1/4},{1/5}},x]", "MeijerG[{{1/3.},{1/7.}},{{0.25},{0.2}},x]" },
        { "MeijerG[{{1/3,2/3},{}},{{1/4,1/5},{}},x]", "MeijerG[{{1/3.,2/3.},{}},{{0.25,0.2},{}},x]" },
        { "MeijerG[{{},{}},{{1/3},{1/5}},x,1/3]", "MeijerG[{{},{}},{{1/3.},{0.2}},x,1/3]" },
        { "MeijerG[{{1/3,1/7},{}},{{},{1/5}},x]", "MeijerG[{{1/3.,1/7.},{}},{{},{0.2}},x]" },
    };
    double pts[] = { 0.3, 2.7, -1.6 };
    for (size_t i = 0; i < sizeof pairs / sizeof pairs[0]; i++)
        for (size_t k = 0; k < 3; k++) {
            char e[1024];
            snprintf(e, sizeof e, "Abs[N[(%s) /. x -> %g] - ((%s) /. x -> %g)] < 10^-12 Max[1, Abs[(%s) /. x -> %g]]",
                     pairs[i][0], pts[k], pairs[i][1], pts[k], pairs[i][1], pts[k]);
            assert_true(e);
        }
}

void test_mg_zero() {
    assert_str("MeijerG[{{1/2},{1/2}},{{},{}},0]", "ComplexInfinity");
    assert_str("MeijerG[{{1/2},{}},{{1},{}},0]", "0");
    assert_str("MeijerG[{{1/2},{}},{{0},{0}},0]", "Sqrt[Pi]");
    assert_str("MeijerG[{{1/2,1},{3,4}},{{1/5},{}},0]", "MeijerG[{{1/2, 1}, {3, 4}}, {{1/5}, {}}, 0]");
    assert_str("MeijerG[{{},{1}},{{0,1/3},{}},0]", "Gamma[1/3]");
    assert_str("MeijerG[{{},{1}},{{-1/4,1/3},{}},0]", "ComplexInfinity");
    assert_str("MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},0]", "0");
    assert_str("MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},0.]", "0.0");
}

void test_mg_messages() {
    assert_str("Check[MeijerG[{{},{}},{{},{0}},z], failed]", "failed");
    assert_str("Quiet[MeijerG[{{},{}},{{},{0}},z]]", "MeijerG[{{}, {}}, {{}, {0}}, z]");
    assert_str("Check[MeijerG[], failed]", "failed");
    assert_str("Quiet[MeijerG[]]", "MeijerG[]");
    assert_str("Quiet[MeijerG[1, 2, 3]]", "MeijerG[1, 2, 3]");
    /* a pole collision (a_k - b_j a positive integer) */
    assert_str("Check[MeijerG[{{3/2},{}},{{1/2},{}},z], failed]", "failed");
    assert_str("Check[MeijerG[{{1.5},{}},{{0.5,1.0},{0.25}},0.01], failed]", "failed");
    /* ... resolved by the Gamma-pair move when one exists (no message) */
    assert_str("Check[MeijerG[{{31,4},{2,3}},{{5,3/2},{1/10}},3], failed]",
               "-MeijerG[{{2, 4, 31}, {3}}, {{3/2}, {1/10, 5}}, 3]");
    /* canonical order within each group */
    assert_str("MeijerG[{{1},{2}},{{2/3,1/5},{4}},z]", "MeijerG[{{1}, {2}}, {{1/5, 2/3}, {4}}, z]");
    /* r = 1 is the ordinary function */
    assert_str("MeijerG[{{1},{2}},{{2/3,1/5},{4}},z,1]", "MeijerG[{{1}, {2}}, {{1/5, 2/3}, {4}}, z]");
}

/* ---- D, Series, FunctionExpand ----------------------------------------- */

void test_mg_derivative() {
    assert_str("D[MeijerG[{{a1,a2},{a3,a4}},{{v},{-v}},z],z]",
               "MeijerG[{{-1, -1 + a1, -1 + a2}, {-1 + a3, -1 + a4}}, {{-1 + v}, {0, -1 - v}}, z]");
    /* against a 40-digit central difference */
    const char* fs[] = {
        "MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},#]&",
        "MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},#,1/3]&",
        "MeijerG[{{1/2,1},{3,4}},{{1/5},{}},#]&",
        "MeijerG[{{},{2}},{{1/2,3/2},{}},#]&",
    };
    for (size_t i = 0; i < sizeof fs / sizeof fs[0]; i++) {
        char e[1024];
        snprintf(e, sizeof e,
                 "With[{f = %s, h = 10^-15}, Abs[N[(D[f[t], t] /. t -> 7/5) - (f[7/5 + h] - f[7/5 - h])/(2 h), 40]] < 10^-25]",
                 fs[i]);
        assert_true(e);
    }
    /* chain rule through the argument */
    assert_true("Abs[N[(D[MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},x^2],x] /. x -> 9/10) - "
                "2 (9/10) (D[MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},t],t] /. t -> 81/100), 30]] < 10^-25");
    /* a parameter depending on the variable: unevaluated derivative, no crash */
    assert_true("!FreeQ[D[MeijerG[{{x},{}},{{1/2},{1}},z],x], Derivative]");
}

void test_mg_series() {
    assert_str("Series[MeijerG[{{a},{}},{{d},{}},x],{x,0,2}]",
               "x^d (Gamma[1 - a + d] + (-1 + a - d) Gamma[1 - a + d] x + 1/2 (-2 + a - d) (-1 + a - d) Gamma[1 - a + d] x^2 + O[x]^3)");
    /* the log case: 2 z^(a+1/2) BesselK[1, 2 Sqrt[z]] */
    assert_str("Series[MeijerG[{{},{}},{{a,a+1},{}},z],{z,0,1}]",
               "(1/Sqrt[z] + (-1 + 2 EulerGamma + Log[z]) Sqrt[z] + O[z]^(3/2)) z^(1/2 (1 + 2 a))");
    /* truncated series against the value (unreduced shapes, incl. log case) */
    const char* gs[] = {
        "MeijerG[{{1/3},{1/4}},{{1/5,1/6},{1/7}},x]",
        "MeijerG[{{},{2}},{{1/2,3/2},{}},x]",
        "MeijerG[{{},{}},{{0,0,1},{}},x]",
        "MeijerG[{{1/3},{2}},{{0,1,1/2},{}},x]",
    };
    for (size_t i = 0; i < sizeof gs / sizeof gs[0]; i++) {
        char e[1024];
        snprintf(e, sizeof e, "Abs[N[Normal[Series[%s, {x, 0, 12}]] /. x -> 1/20, 30] - N[%s /. x -> 1/20, 30]] < 10^-14",
                 gs[i], gs[i]);
        assert_true(e);
    }
    /* at Infinity (p >= q): the dual expansion */
    assert_true("Abs[N[Normal[Series[MeijerG[{{1/3,1/2},{1/4}},{{1/5},{1/7}},x], {x, Infinity, 12}]] /. x -> 20, 30]"
                " - N[MeijerG[{{1/3,1/2},{1/4}},{{1/5},{1/7}},20], 30]] < 10^-14");
}

void test_mg_function_expand() {
    assert_str("FunctionExpand[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},z]]",
               "1/12 Gamma[1/5] Sqrt[Pi] HypergeometricPFQ[{1/5}, {1/2, 3, 4}, 1/z] - (256 Gamma[7/10] "
               "HypergeometricPFQ[{7/10}, {3/2, 7/2, 9/2}, 1/z])/(1575 Sqrt[Pi] Sqrt[z])");
    assert_true("Abs[N[FunctionExpand[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},z]] /. z -> 5/2, 30] - "
                "N[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},5/2], 30]] < 10^-25");
    assert_true("Abs[N[FunctionExpand[MeijerG[{{1},{2}},{{2/3,1/5},{4}},z]] /. z -> 3, 30] - "
                "N[MeijerG[{{1},{2}},{{2/3,1/5},{4}},3], 30]] < 10^-25");
    /* the log case is left alone, as is anything without a MeijerG */
    assert_str("FunctionExpand[MeijerG[{{},{2}},{{1/2,3/2},{}},z]]", "MeijerG[{{}, {2}}, {{1/2, 3/2}, {}}, z]");
    assert_str("FunctionExpand[Sin[x] + 1]", "1 + Sin[x]");
}

void test_mg_attributes() {
    SymbolDef* d = symtab_get_def("MeijerG");
    ASSERT(d && (d->attributes & ATTR_NUMERICFUNCTION) && (d->attributes & ATTR_PROTECTED));
    ASSERT(!(d->attributes & ATTR_LISTABLE));
    ASSERT(d->docstring != NULL);
    SymbolDef* f = symtab_get_def("FunctionExpand");
    ASSERT(f && (f->attributes & ATTR_PROTECTED) && f->docstring != NULL);
    /* N numericalizes an unevaluated exact call */
    assert_true("MachineNumberQ[N[MeijerG[{{1},{2}},{{2/3,1/5},{4}},3]]]");
    assert_true("NumericQ[MeijerG[{{1},{2}},{{2/3,1/5},{4}},3]]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_mg_machine_values);
    TEST(test_mg_arbitrary_precision);
    TEST(test_mg_surfaces);
    TEST(test_mg_compile);
    TEST(test_mg_reductions);
    TEST(test_mg_zero);
    TEST(test_mg_messages);
    TEST(test_mg_derivative);
    TEST(test_mg_series);
    TEST(test_mg_function_expand);
    TEST(test_mg_attributes);

    printf("All MeijerG tests passed.\n");
    return 0;
}
