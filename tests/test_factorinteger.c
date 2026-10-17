/* test_factorinteger.c — FactorInteger and every FactorInteger`Method.
 *
 * Expected outputs marked "MMA" were taken from Mathematica 14.0
 * (wolframscript) on 2026-10-17; every other expectation is known by
 * construction (products of explicit primes) or checked in-language as a
 * contract: the prime powers multiply back to n and every base is prime.
 */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "fi_internal.h"
#include <string.h>
#include <stdlib.h>

static int failures = 0;

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        failures++;
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* In-language contract: FactorInteger[n, opts] multiplies back to n and
 * every base is a prime. */
static void contract(const char* n, const char* opts) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "With[{f = Quiet[FactorInteger[%s%s%s]]}, "
             "Times @@ (#[[1]]^#[[2]] & /@ f) == %s && "
             "AllTrue[First /@ f, PrimeQ[#] || # == -1 &]]",
             n, *opts ? ", " : "", opts, n);
    check(buf, "True");
}

/* ---- Basics: signs, zero, units, rationals (MMA) ---- */
static void test_basics(void) {
    check("FactorInteger[1]", "{{1, 1}}");
    check("FactorInteger[-1]", "{{-1, 1}}");
    check("FactorInteger[0]", "{{0, 1}}");
    check("FactorInteger[2]", "{{2, 1}}");
    check("FactorInteger[-12]", "{{-1, 1}, {2, 2}, {3, 1}}");
    check("FactorInteger[2^64 + 1]", "{{274177, 1}, {67280421310721, 1}}");
    check("FactorInteger[10/21]", "{{2, 1}, {3, -1}, {5, 1}, {7, -1}}");
    check("FactorInteger[-6/35]", "{{-1, 1}, {2, 1}, {3, 1}, {5, -1}, {7, -1}}");
    check("FactorInteger[2^4096]", "{{2, 4096}}");
    check("FactorInteger[(2^61 - 1)^3 * (10^25 + 13)^2]",
          "{{2305843009213693951, 3}, {10000000000000000000000013, 2}}");
    check("FactorInteger[{12, 15}]", "{{{2, 2}, {3, 1}}, {{3, 1}, {5, 1}}}");
    check("FactorInteger[12, Automatic]", "{{2, 2}, {3, 1}}");
}

/* ---- Partial factorisation: at most k entries, the last a remainder ---- */
static void test_partial(void) {
    check("FactorInteger[100, 1]", "{{100, 1}}");                 /* MMA */
    check("FactorInteger[6/35, 1]", "{{6/35, 1}}");               /* MMA */
    check("FactorInteger[100, 2]", "{{2, 2}, {25, 1}}");
    check("FactorInteger[5040, 4]", "{{2, 4}, {3, 2}, {5, 1}, {7, 1}}");  /* MMA */
    check("Length[FactorInteger[5040, 3]] <= 3", "True");
    check("Times @@ (#[[1]]^#[[2]] & /@ FactorInteger[5040, 3]) == 5040", "True");
    check("Times @@ (#[[1]]^#[[2]] & /@ FactorInteger[-30, 2]) == -30", "True");
    check("FactorInteger[2^128 + 1, 0]", "FactorInteger[340282366920938463463374607431768211457, 0]");
}

/* ---- Gaussian integers (MMA) ---- */
static void test_gaussian(void) {
    check("FactorInteger[5, GaussianIntegers -> True]", "{{-I, 1}, {1 + 2*I, 1}, {2 + I, 1}}");
    check("FactorInteger[2, GaussianIntegers -> True]", "{{-I, 1}, {1 + I, 2}}");
    check("FactorInteger[12, GaussianIntegers -> True]", "{{-1, 1}, {1 + I, 4}, {3, 1}}");
    check("FactorInteger[-5, GaussianIntegers -> True]", "{{I, 1}, {1 + 2*I, 1}, {2 + I, 1}}");
    check("FactorInteger[-12, GaussianIntegers -> True]", "{{1 + I, 4}, {3, 1}}");
    check("FactorInteger[13*(2 + 3*I), GaussianIntegers -> True]",
          "{{-I, 1}, {2 + 3*I, 2}, {3 + 2*I, 1}}");
    check("FactorInteger[10/21, GaussianIntegers -> True]",
          "{{-1, 1}, {1 + I, 2}, {1 + 2*I, 1}, {2 + I, 1}, {3, -1}, {7, -1}}");
    check("FactorInteger[3 + 4*I]", "{{2 + I, 2}}");
    check("FactorInteger[-7 + I]", "{{1 + I, 1}, {1 + 2*I, 2}}");
    check("FactorInteger[I]", "{{I, 1}}");
    check("FactorInteger[1 + I]", "{{1 + I, 1}}");
    check("FactorInteger[2*I]", "{{1 + I, 2}}");
    check("FactorInteger[7*I]", "{{I, 1}, {7, 1}}");
    check("FactorInteger[(3 + 4*I)/(1 + I)]", "{{1 + I, -1}, {2 + I, 2}}");
    check("FactorInteger[1/2 + I/3]", "{{I, 1}, {1 + I, -2}, {3, -1}, {3 + 2*I, 1}}");
    check("FactorInteger[3*(1 + 4*I)]", "{{1 + 4*I, 1}, {3, 1}}");
}

/* ---- Every method on a showcase input ---- */
static void test_methods(void) {
    check("FactorInteger`TrialDivision[2^10*3^5*1000003*999983]",
          "{{2, 10}, {3, 5}, {999983, 1}, {1000003, 1}}");
    check("FactorInteger`PollardRho[(10^12 + 39)*(10^15 + 37)]",
          "{{1000000000039, 1}, {1000000000000037, 1}}");
    check("FactorInteger`Fermat[1000000007*1000000009]", "{{1000000007, 1}, {1000000009, 1}}");
    check("FactorInteger`Lehman[1000003*1000033]", "{{1000003, 1}, {1000033, 1}}");
    check("FactorInteger`ShanksSquareForms[1000000000000037*2000000000000021]",
          "{{1000000000000037, 1}, {2000000000000021, 1}}");
    check("FactorInteger`SQUFOF[1000003*1000033]", "{{1000003, 1}, {1000033, 1}}");
    check("FactorInteger`CFRAC[1000000000000037*2000000000000021]",
          "{{1000000000000037, 1}, {2000000000000021, 1}}");
    check("FactorInteger`Dixon[10000000019*20000000089]", "{{10000000019, 1}, {20000000089, 1}}");
    check("FactorInteger`QuadraticSieve[100000000000000000039*300000000000000000053]",
          "{{100000000000000000039, 1}, {300000000000000000053, 1}}");
    check("FactorInteger`QS[1000003*1000033]", "{{1000003, 1}, {1000033, 1}}");
    check("FactorInteger`PollardPMinus1[(2^61 - 1)*(10^40 + 121)]",
          "{{2305843009213693951, 1}, {10000000000000000000000000000000000000121, 1}}");
    check("FactorInteger`WilliamsPPlus1[91]", "{{7, 1}, {13, 1}}");
    check("FactorInteger`ECM[2^128 + 1]", "{{59649589127497217, 1}, {5704689200685129054721, 1}}");
    check("FactorInteger`ECM[2^128 + 1, \"B1\" -> 11000, \"Curves\" -> 200]",
          "{{59649589127497217, 1}, {5704689200685129054721, 1}}");
    check("FactorInteger`BlakeRationalBaseDescent[1001, \"Base\" -> 22/7]", "{{7, 1}, {11, 1}, {13, 1}}");
    /* Method -> forms and aliases agree with the FactorInteger` symbols. */
    check("FactorInteger[1001, Method -> {\"BlakeRationalBaseDescent\", \"Base\" -> 22/7}]",
          "{{7, 1}, {11, 1}, {13, 1}}");
    check("FactorInteger[1000003*1000033, Method -> \"PollardP-1\"] == "
          "FactorInteger`PollardPMinus1[1000003*1000033]", "True");
    check("FactorInteger[2^64 + 1, Method -> \"SIQS\"]", "{{274177, 1}, {67280421310721, 1}}");
    check("FactorInteger[2^64 + 1, Method -> Automatic]", "{{274177, 1}, {67280421310721, 1}}");
}

/* ---- Contract over shapes and sizes, per method ---- */
static void test_contract(void) {
    static const char* ns[] = {
        "2^31 - 1", "2^62 + 7", "(2^31 - 1)*(2^31 + 11)", "10^20 + 39",
        "(10^15 + 37)*(10^16 + 61)", "3^40 * 5^7", "2^256 + 1", "2^149 - 1",
        "Factorial[30] + 1", "10^60 + 1", "(10^12 + 39)^5", "(10^19 + 51)*(10^19 + 117)",
    };
    for (size_t i = 0; i < sizeof ns / sizeof ns[0]; i++) contract(ns[i], "");
    contract("(10^9 + 7)*(10^9 + 9)*(10^12 + 39)", "Method -> \"PollardRho\"");
    contract("(10^9 + 7)*(10^9 + 9)", "Method -> \"Fermat\"");
    contract("(10^6 + 3)*(10^6 + 33)*(10^5 + 3)", "Method -> \"Lehman\"");
    contract("(10^9 + 7)*(10^10 + 19)", "Method -> \"ShanksSquareForms\"");
    contract("(10^12 + 39)*(10^13 + 37)", "Method -> \"CFRAC\"");
    /* k = 2 with both primes == 7 (mod 8): every dependency of that expansion
     * is trivial; CFRAC must switch multiplier (it used to loop for ever). */
    contract("12791*5060983", "Method -> \"CFRAC\"");
    contract("886945583297207549892329", "Method -> \"CFRAC\"");
    contract("(10^8 + 7)*(10^9 + 7)", "Method -> \"Dixon\"");
    contract("(10^15 + 37)*(10^16 + 61)*(10^17 + 3)", "Method -> \"QuadraticSieve\"");
    contract("(10^9 + 7)*(10^12 + 39)", "Method -> \"ECM\"");
    contract("736725780226971718034605203450", "Method -> \"PollardPMinus1\"");
    contract("864041305442689345542027354876", "Method -> \"PollardPMinus1\"");
    contract("2^40 * 3^20 * 7^10", "Method -> \"ECM\"");
    /* Parallel paths: QS from 150-bit kN, ECM curves from 160-bit n. */
    contract("NextPrime[10^27]*NextPrime[3*10^27]", "Method -> \"QuadraticSieve\"");
    contract("NextPrime[10^15]*NextPrime[10^60]", "Method -> \"ECM\"");
    contract("NextPrime[10^24]*NextPrime[10^25]*NextPrime[10^26]", "");
}

/* ---- Messages route through the funnel (Quiet / Check) ---- */
static void test_messages(void) {
    check("Check[FactorInteger[10, Method -> \"NoSuchMethod\"], failed]", "failed");
    check("Quiet[FactorInteger[10, Method -> \"NoSuchMethod\"]]",
          "FactorInteger[10, Method -> \"NoSuchMethod\"]");
    check("Check[FactorInteger`ECM[10, \"Bound\" -> 3], failed]", "failed");
    check("Check[FactorInteger[2.5], failed]", "failed");
    /* A bounded method that gives up says so: nofac, never a silent composite. */
    check("Check[FactorInteger`TrialDivision[(10^20 + 39)*(10^21 + 117), \"Bound\" -> 1000], failed]",
          "failed");
    check("Quiet[FactorInteger`TrialDivision[(10^20 + 39)*(10^21 + 117), \"Bound\" -> 1000]]",
          "{{100000000000000000050700000000000000004563, 1}}");
}

/* ---- Options, attributes, docstrings ---- */
static void test_symbols(void) {
    check("Options[FactorInteger]", "{GaussianIntegers -> False, Method -> Automatic}");
    check("Attributes[FactorInteger]", "{Listable, Protected}");
    check("Attributes[FactorInteger`ECM]", "{Listable, Protected}");
    check("FactorInteger`ECM[{91, 2^64 + 1}]", "{{{7, 1}, {13, 1}}, {{274177, 1}, {67280421310721, 1}}}");
}

/* ---- C level: the word engine ---- */
static void test_u64(void) {
    /* Strong pseudoprimes to several small bases must be rejected. */
    static const uint64_t composites[] = {
        UINT64_C(3215031751), UINT64_C(2152302898747), UINT64_C(3474749660383),
        UINT64_C(341550071728321), UINT64_C(3825123056546413051),
    };
    for (size_t i = 0; i < 5; i++)
        if (fi_u64_is_prime(composites[i])) { printf("FAIL: pseudoprime %llu accepted\n", (unsigned long long)composites[i]); failures++; }
    static const uint64_t primes[] = { 2, 3, 61, UINT64_C(4294967291), UINT64_C(18446744073709551557) };
    for (size_t i = 0; i < 5; i++)
        if (!fi_u64_is_prime(primes[i])) { printf("FAIL: prime %llu rejected\n", (unsigned long long)primes[i]); failures++; }
    /* Each splitter on 1000000007 * 1000000009. */
    uint64_t n = UINT64_C(1000000007) * UINT64_C(1000000009), rng = 1, d;
    d = fi_u64_ecm(n, &rng, 512);
    if (d != UINT64_C(1000000007) && d != UINT64_C(1000000009)) { printf("FAIL: u64 ecm\n"); failures++; }
    d = fi_u64_squfof(n, 0);
    if (d != UINT64_C(1000000007) && d != UINT64_C(1000000009)) { printf("FAIL: u64 squfof\n"); failures++; }
    d = fi_u64_fermat(n, 1000);
    if (d != UINT64_C(1000000007) && d != UINT64_C(1000000009)) { printf("FAIL: u64 fermat\n"); failures++; }
    d = fi_u64_lehman(UINT64_C(1000003) * UINT64_C(1000033));
    if (d != UINT64_C(1000003) && d != UINT64_C(1000033)) { printf("FAIL: u64 lehman\n"); failures++; }
    d = 0;
    for (uint64_t c = 1; c < 8 && !d; c++) d = fi_u64_rho(n, c, UINT64_C(1) << 22);
    if (d != UINT64_C(1000000007) && d != UINT64_C(1000000009)) { printf("FAIL: u64 rho\n"); failures++; }
    /* Complete factorisation of a range of words matches trial division. */
    for (uint64_t m = 2; m < 200000; m++) {
        uint64_t p[64];
        int64_t e[64];
        int k = fi_u64_factor_array(m, p, e);
        uint64_t prod = 1;
        for (int i = 0; i < k; i++) {
            if (!fi_u64_is_prime(p[i])) { printf("FAIL: composite base %llu of %llu\n", (unsigned long long)p[i], (unsigned long long)m); failures++; break; }
            for (int64_t j = 0; j < e[i]; j++) prod *= p[i];
        }
        if (prod != m) { printf("FAIL: word factorisation of %llu\n", (unsigned long long)m); failures++; break; }
    }
}

/* ---- C level: GF(2) solvers (dense and block Lanczos) ---- */
static void test_gf2(size_t ncols) {
    /* Random sparse relations: each a random set of ~12 columns; with 80 more
     * rows than columns, dependencies must exist.  Each returned dependency
     * must XOR to zero. */
    size_t nrows = ncols + 80;
    FiSparseMat m;
    m.nrows = (uint32_t)nrows;
    m.ncols = (uint32_t)ncols;
    m.start = malloc((nrows + 1) * sizeof(uint32_t));
    m.col = malloc(nrows * 16 * sizeof(uint32_t));
    uint64_t s = 88172645463325252ULL;
    size_t q = 0;
    for (size_t r = 0; r < nrows; r++) {
        m.start[r] = (uint32_t)q;
        int w = 8 + (int)(s % 8);
        uint32_t last = UINT32_MAX;
        for (int k = 0; k < w; k++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            /* skew towards small columns, like a factor base */
            uint32_t c = (uint32_t)((s % ncols) * (s % ncols) / ncols);
            if (c != last) m.col[q++] = c;
            last = c;
        }
    }
    m.start[nrows] = (uint32_t)q;
    size_t W = fi_gf2_words(nrows);
    uint64_t* deps = calloc(64 * W, sizeof(uint64_t));
    size_t nd = fi_gf2_dependencies(&m, deps, 64);
    if (nd == 0) { printf("FAIL: gf2 (%zu cols) found no dependency\n", ncols); failures++; }
    unsigned char* acc = malloc(ncols);
    for (size_t k = 0; k < nd; k++) {
        memset(acc, 0, ncols);
        int nonzero = 0;
        for (size_t r = 0; r < nrows; r++) {
            if (!((deps[k * W + r / 64] >> (r % 64)) & 1)) continue;
            nonzero = 1;
            for (uint32_t j = m.start[r]; j < m.start[r + 1]; j++) acc[m.col[j]] ^= 1;
        }
        for (size_t c = 0; c < ncols; c++)
            if (acc[c]) { printf("FAIL: gf2 (%zu cols) dependency %zu does not vanish\n", ncols, k); failures++; break; }
        if (!nonzero) { printf("FAIL: gf2 zero dependency\n"); failures++; }
    }
    free(acc); free(deps); free(m.start); free(m.col);
}

int main(void) {
    symtab_init();
    core_init();
    test_basics();
    test_partial();
    test_gaussian();
    test_methods();
    test_contract();
    test_messages();
    test_symbols();
    test_u64();
    test_gf2(300);      /* dense */
    test_gf2(5000);     /* block Lanczos */
    if (failures) {
        printf("%d FactorInteger test(s) FAILED\n", failures);
        return 1;
    }
    printf("All FactorInteger tests passed!\n");
    return 0;
}
