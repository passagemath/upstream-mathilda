/* test_analytic_review.c -- the number-theory review, batch 7 (v0.476):
 * Zeta (IncludeSingularTerm, the symmetric convention at a <= 0 on every
 * surface, the exact shift rules, the option messages), LerchPhi (its two
 * options and their messages, the branch cut z > 1, DoublyInfinite, a <= 0),
 * DirichletL (NHoldAll, not NumericFunction), Log2 (exact powers, the
 * infinities), LogIntegral (li(0.) = 0) and AlgebraicNumberNorm (argx, nonopt,
 * optx, ToNumberField::nnfel, the first Extension setting winning).
 *
 * Every literal expectation below was verified against the reference
 * implementation, not recalled. */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("analytic review test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* |lhs - rhs| < tol, written as one input so the reference value is a literal. */
static void check_near(const char* lhs, const char* rhs, const char* tol) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Abs[(%s) - (%s)] < %s", lhs, rhs, tol);
    check(buf, "True");
}

/* ---- Zeta: exact Hurwitz values at a <= 0 and the shift rules ---- */
static void test_zeta_exact(void) {
    check("Zeta[2, -1] === 1 + Pi^2/6", "True");
    check("Zeta[2, 0] === Pi^2/6", "True");
    check("Zeta[2, -3/2] === 40/9 + Pi^2/2", "True");
    check("Zeta[2, -1/2] === 4 + Pi^2/2", "True");
    check("Zeta[2, -10] === 1968329/1270080 + Pi^2/6", "True");
    check("Zeta[2, -200] - HarmonicNumber[200, 2] === Pi^2/6", "True");
    check("Zeta[4, -1/2] === 16 + Pi^4/6", "True");
    check("Zeta[4, -3] === 1393/1296 + Pi^4/90", "True");
    check("Zeta[3, -1/2] === 8 + 7 Zeta[3]", "True");
    check("Zeta[3, -5/2] === 28216/3375 + 7 Zeta[3]", "True");
    check("Zeta[3, 5/2] === -224/27 + 7 Zeta[3]", "True");
    /* s = 0 counts the k + a = 0 term as 0^0 = 1; Re s < 0 as 0 */
    check("{Zeta[0, -1], Zeta[0, 0], Zeta[0, -2], Zeta[0, -1/2], Zeta[0, 1/3]}",
          "{3/2, 1/2, 5/2, 1, 1/6}");
    check("{Zeta[-1, -2], Zeta[-2, -1], Zeta[-1, 0], Zeta[-3, -2], Zeta[-2, -3]}",
          "{35/12, 1, -1/12, 1081/120, 14}");
    check("{Zeta[-2, -1/2], Zeta[-1, -1/2], Zeta[-1, 5/2]}", "{1/4, 13/24, -47/24}");
    check("{Zeta[1, -1], Zeta[1, 0], Zeta[1, -1/2], Zeta[1, 5/2]}",
          "{ComplexInfinity, ComplexInfinity, ComplexInfinity, ComplexInfinity}");
    /* symbolic or non-integer s: only the trivial shifts (terms 1^-s) */
    check("Zeta[x, -1] === 1 + Zeta[x]", "True");
    check("Zeta[x, 0] === Zeta[x]", "True");
    check("Zeta[x, 2] === -1 + Zeta[x]", "True");
    check("Zeta[2 + I, -1] === 1 + Zeta[2 + I]", "True");
    check("Head /@ {Zeta[x, 3], Zeta[x, -2], Zeta[x, -1/2], Zeta[2 + I, -2], "
          "Zeta[3/2, -1/2], Zeta[3/2, 3], Zeta[x, 5/2]}",
          "{Zeta, Zeta, Zeta, Zeta, Zeta, Zeta, Zeta}");
    /* odd s >= 3 has no closed form: integer a past the trivial ones stays */
    check("Head /@ {Zeta[3, 4], Zeta[3, 3], Zeta[5, -2], Zeta[3, -2]}", "{Zeta, Zeta, Zeta, Zeta}");
    check("Head /@ {Zeta[2, 1/3], Zeta[-2, 1/3], Zeta[-2, -2/3]}", "{Zeta, Zeta, Zeta}");
}

/* ---- Zeta: IncludeSingularTerm and the option messages ---- */
static void test_zeta_options(void) {
    check("Options[Zeta]", "{IncludeSingularTerm -> False}");
    check("Zeta[2, -1, IncludeSingularTerm -> True]", "ComplexInfinity");
    check("Zeta[2, 0, IncludeSingularTerm -> True]", "ComplexInfinity");
    check("Zeta[1/2, -1, IncludeSingularTerm -> True]", "ComplexInfinity");
    check("Zeta[2, -1, IncludeSingularTerm :> True]", "ComplexInfinity");
    check("Zeta[2, -1, \"IncludeSingularTerm\" -> True]", "ComplexInfinity");
    check("Zeta[2., -1., IncludeSingularTerm -> True]", "ComplexInfinity");
    check("Zeta[-1, -1, IncludeSingularTerm -> True]", "11/12");
    check("Zeta[0, 0, IncludeSingularTerm -> True]", "1/2");
    check("Zeta[0., -1., IncludeSingularTerm -> True]", "1.5");
    check("Zeta[3, -1/2, IncludeSingularTerm -> True] === 8 + 7 Zeta[3]", "True");
    check("Zeta[2, -1, IncludeSingularTerm -> False, IncludeSingularTerm -> True] === 1 + Pi^2/6",
          "True");
    check("Zeta[2, IncludeSingularTerm -> True] === Pi^2/6", "True");
    check("Zeta[2, -1, {IncludeSingularTerm -> True}]", "{ComplexInfinity}");
    check("Zeta[x, -1, IncludeSingularTerm -> True]", "Zeta[x, -1, IncludeSingularTerm -> True]");
    check("Zeta[-1/2, -1, IncludeSingularTerm -> True]",
          "Zeta[-1/2, -1, IncludeSingularTerm -> True]");
    check_near("Zeta[2., -1.5, IncludeSingularTerm -> True]", "9.3792466449891237539", "1.*^-14");
    check("Check[Zeta[2, -1, Foo -> 1], m, Zeta::optx]", "m");
    check("Check[Zeta[1, 2, 3], m, Zeta::nonopt]", "m");
    check("Check[Zeta[2, 3, IncludeSingularTerm -> x], m, Zeta::opttf]", "m");
    check("Check[Zeta[], m, Zeta::argt]", "m");
}

/* ---- Zeta: the machine kernel at a <= 0, on every surface ---- */
static void test_zeta_machine(void) {
    check_near("Zeta[2., -1.5]", "9.3792466449891237539", "1.*^-14");
    check_near("Zeta[2.5, -1.5]", "12.266852253362403251", "1.*^-14");
    check_near("Zeta[0.5, -10.25]", "5.7025731812063010903", "1.*^-14");
    check_near("Zeta[2., -100000.5]", "9.8695944011893577022", "1.*^-14");
    check_near("Zeta[2. + I, -1.5]", "6.9928629142974968305 + 4.4119344920932555624 I", "1.*^-14");
    check("{Zeta[0., -1], Zeta[0., 0.]}", "{1.5, 0.5}");
    check_near("Zeta[-1., -1.]", "11/12", "1.*^-15");
    /* HurwitzZeta keeps the principal powers at a < 0 */
    check_near("HurwitzZeta[2.5, -1.5]",
               "6.2471106345688114853 - 6.0197416187935917653 I", "1.*^-13");
    check("v = N[Range[-3, 3]/2 + 1/4]; f = Compile[{{x, _Real}}, Zeta[2.5, x]]; "
          "g = Compile[{{x, _Real, 1}}, Zeta[2.5, x]]; "
          "{Zeta[2.5, v] === (Zeta[2.5, #] & /@ v), f /@ v === Zeta[2.5, v], "
          "g[v] === Zeta[2.5, v], Normal[Zeta[2.5, NDArray[v]]] === Zeta[2.5, v]}",
          "{True, True, True, True}");
    check("Lookup[CompileDiagnostics[{{x, _Real, 1}}, Zeta[2.5, x]], \"Compiled\"]", "True");
}

/* ---- LerchPhi ---- */
static void test_lerchphi(void) {
    check("Options[LerchPhi]", "{DoublyInfinite -> False, IncludeSingularTerm -> False}");
    check("Check[LerchPhi[0.5, 2, 1, Foo -> 1], m, LerchPhi::optx]", "m");
    check("Check[LerchPhi[0.5, 2, 1/2, DoublyInfinite -> x], m, LerchPhi::opttf]", "m");
    check("Check[LerchPhi[1/2, 2, 1, 2, 3], m, LerchPhi::nonopt]", "m");
    check("Check[LerchPhi[1, 2], m, LerchPhi::argrx]", "m");
    /* the branch cut z > 1: the limit from below */
    check_near("LerchPhi[2., 2., 1.]", "Pi^2/8 - I Pi Log[2]/2", "1.*^-14");
    check_near("LerchPhi[3., 2., 0.5]",
               "4.3083059100926142947 - 1.9926622707261163742 I", "1.*^-14");
    check_near("LerchPhi[2., 2., 0.3]",
               "12.402557086272875158 - 1.7687495201211067168 I", "1.*^-13");
    check_near("LerchPhi[2., 2., -1.]",
               "5.9348022005446793094 - 4.3551721806072042610 I", "1.*^-14");
    check_near("LerchPhi[0.5, 2., 0.5, DoublyInfinite -> True]",
               "13.957728399277759068 - 3.0795717821423574190 I", "1.*^-13");
    check("LerchPhi[2., 0., 1.5]", "-1.0");
    /* z = 1 is Zeta, a <= 0 included */
    check("LerchPhi[1., 2., -1.5] == Zeta[2., -1.5]", "True");
    check_near("LerchPhi[0.5, 2., -1.5]", "3.5137308195896818897", "1.*^-14");
    check_near("LerchPhi[-0.7, 1.5, -3.2]", "-3.2183679870742323985", "1.*^-14");
}

/* ---- DirichletL, Log2, LogIntegral ---- */
static void test_dirichlet_logs(void) {
    check("Attributes[DirichletL]", "{Listable, NHoldAll, Protected}");
    check("NumericQ[DirichletL[4, 2, 3]]", "False");
    check("Head[1. + DirichletL[4, 2, 3]]", "Plus");
    check_near("N[DirichletL[4, 2, 3]]", "0.96894614625936938048", "1.*^-15");

    check("{Log2[2^1000], Log2[1/2^1000], Log2[1/8], Log2[2^63], Log10[10^40], Log10[1/1000]}",
          "{1000, -1000, -3, 63, 40, -3}");
    check("{Log2[Infinity], Log2[ComplexInfinity], Log2[Indeterminate], Log[ComplexInfinity], "
          "Log[Indeterminate]}",
          "{Infinity, Infinity, Indeterminate, Infinity, Indeterminate}");
    check("Log2[6] === Log[6]/Log[2]", "True");
    check("v = N[Range[100]]; Log2[v] === (Log2[#] & /@ v)", "True");
    check("Log2[{1., 8.}]", "{0.0, 3.0}");

    check("{LogIntegral[0.], LogIntegral[-0.], LogIntegral[0]}", "{0.0, 0.0, 0}");
    check("LogIntegral[ComplexInfinity]", "LogIntegral[ComplexInfinity]");
    check("Normal[LogIntegral[NDArray[{0., 2.}]]] === {0., LogIntegral[2.]}", "True");
    check("Compile[{{x, _Real}}, LogIntegral[x]][0.]", "0.0");
}

/* ---- AlgebraicNumberNorm: options and messages ---- */
static void test_algebraicnumbernorm(void) {
    check("Options[AlgebraicNumberNorm]", "{Extension -> None}");
    check("Check[AlgebraicNumberNorm[], m, AlgebraicNumberNorm::argx]", "m");
    check("Check[AlgebraicNumberNorm[1, 2, 3], m, AlgebraicNumberNorm::nonopt]", "m");
    check("Check[AlgebraicNumberNorm[Sqrt[2], x], m, AlgebraicNumberNorm::nonopt]", "m");
    check("Check[AlgebraicNumberNorm[Sqrt[2], Foo -> 1], m, AlgebraicNumberNorm::optx]", "m");
    check("Check[AlgebraicNumberNorm[Sqrt[2], Extension -> Sqrt[3]], m, ToNumberField::nnfel]", "m");
    check("AlgebraicNumberNorm[Sqrt[2], Extension -> Sqrt[2], Extension -> x]", "-2");
    check("AlgebraicNumberNorm[Sqrt[2], {Extension -> Sqrt[2]}]", "{-2}");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_zeta_exact);
    TEST(test_zeta_options);
    TEST(test_zeta_machine);
    TEST(test_lerchphi);
    TEST(test_dirichlet_logs);
    TEST(test_algebraicnumbernorm);

    printf("All analytic review tests passed!\n");
    return 0;
}
