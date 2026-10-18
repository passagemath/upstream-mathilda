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
#include "fi_rbd.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int failures = 0;

/* mpz == decimal string, 0 when equal. */
static int mpz_cmp_str_ok(const mpz_t z, const char* s) {
    mpz_t t;
    mpz_init_set_str(t, s, 10);
    int c = mpz_cmp(z, t);
    mpz_clear(t);
    return c;
}

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
    /* Bignum Lehman reaches multipliers k >= 2 above ~45 digits: p/q near a
     * small rational u/v is found at k = uv.  Before the (k, fraction) sweep
     * the k = 1 window alone (n^(1/6)/4 values) used the whole budget. */
    check("With[{p = NextPrime[10^25]}, FactorInteger`Lehman[p NextPrime[3 p]]]",
          "{{10000000000000000000000013, 1}, {30000000000000000000000067, 1}}");
    check("With[{p = NextPrime[10^30]}, FactorInteger`Lehman[p NextPrime[3 p]]]",
          "{{1000000000000000000000000000057, 1}, {3000000000000000000000000000181, 1}}");
    check("With[{p = NextPrime[10^30]}, FactorInteger`Lehman[p NextPrime[Floor[5 p/7]]]]",
          "{{714285714285714285714285714443, 1}, {1000000000000000000000000000057, 1}}");
    /* Near-square factors (k = 1, plain Fermat's case) are still found. */
    check("FactorInteger`Lehman[NextPrime[10^25] NextPrime[10^25 + 10^16]]",
          "{{10000000000000000000000013, 1}, {10000000010000000000000003, 1}}");
    /* "Iterations" is honoured: 100 steps cannot reach k = 3. */
    check("Quiet[Length[With[{p = NextPrime[10^25]}, "
          "FactorInteger`Lehman[p NextPrime[3 p], \"Iterations\" -> 100]]]]", "1");
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
    /* Threaded CFRAC / Dixon relation collection and p+1 seeds. */
    contract("NextPrime[3*10^14]*NextPrime[7*10^16]", "Method -> \"CFRAC\"");
    contract("NextPrime[3*10^8]*NextPrime[7*10^9]", "Method -> \"Dixon\"");
    contract("4346305828114280557845729185329137098049763216260074649658093501860550615869998510565182512056793391457142603346382853249",
             "Method -> \"WilliamsPPlus1\"");
}

/* ---- Threaded relation collection replays the serial stream ----
 * CFRAC and Dixon trial-divide on threads but feed the relation store in
 * candidate order, so the smallest "Iterations" budget that succeeds is the
 * serial method's, whatever the thread count (these are the v0.453 serial
 * thresholds).  The CFRAC case crosses multiplier restarts. */
static void test_serial_stream(void) {
    const char* ok = "Length[Quiet[FactorInteger[%s, Method -> {\"%s\", \"Iterations\" -> %d}]]]";
    static const struct { const char* n; const char* m; int thr; } T[] = {
        { "(10^12 + 39)*(10^13 + 37)", "CFRAC", 10313 },
        { "12791*5060983", "CFRAC", 423 },
        { "(10^8 + 7)*(10^9 + 7)", "Dixon", 1211465 },
    };
    char buf[512];
    for (size_t i = 0; i < sizeof T / sizeof T[0]; i++) {
        snprintf(buf, sizeof buf, ok, T[i].n, T[i].m, T[i].thr);
        check(buf, "2");
        snprintf(buf, sizeof buf, ok, T[i].n, T[i].m, T[i].thr - 1);
        check(buf, "1");
    }
    /* A timeout inside a threaded region is taken once the region ends. */
    check("TimeConstrained[FactorInteger[NextPrime[3*10^21]*NextPrime[7*10^23], Method -> \"CFRAC\"], 0.3]",
          "$Aborted");
    check("TimeConstrained[FactorInteger[NextPrime[3*10^12]*NextPrime[7*10^13], Method -> \"Dixon\"], 0.3]",
          "$Aborted");
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

/* ---- BlakeRationalBaseDescent ----
 * Instances from the two papers in factorisation_temp/ (reproduced with
 * gmpy2; the seeded lattice instance is the companion paper's Appendix A
 * driver, seed 5).  The hit records — (k, j) of every distinct gcd — were
 * checked against an independent Python implementation of the paper's
 * reference algorithm (Appendix A, recording every hit), and agree with
 * the paper's Examples 1-4 and its 533-base comb sweep. */
static const char RBD_NMIN[] = "13434917067328449643383271289062122492729438008563207396013386022436439625786814860310635092648849808707283885079425559380060603741711887741983827702536077124606526390026531236424033896417798802339681855054916007583273591721147874872971014173899334436432277372280112100845364029259802969568473368054409653616559493278947418005745364854816116209";
static const char RBD_NMAX[] = "13434917067328449643383271289062122492729438008563207396013386022436439625786814860310635946232502851058036972147222447976537871300232309110236077983704171247600333635603624662228897207255607429818163453115658688268705666669472586861511441849804920313665499850048253134034455058178126780155420375132055020929262098390932828248906338838243332477";
static const char RBD_Q255[] = "31415926535897932384626433832795028841971693993751058209749445923078164119021";
static const char RBD_NCOMB[] = "742245361114170698080188160556986382238609541182118959142360194921632978579427957105431";
static const char RBD_NLAT[] = "195539142822724777401163665279552940381006352574199957151137";   /* 197 bits */
static const char RBD_QLAT[] = "296454907933147663952202071141";
static const char RBD_NMULT[] = "155787495570451096703327230969771495548082681342268724343193395164915333146397756437241087392438717459";

/* One C-level descent; returns the hit count (0/1, or the record size). */
static int rbd_run(mpz_t g, uint64_t* k, int64_t* j, const char* N, unsigned long a,
                   unsigned long b, int exact, int batch, int all, fi_rbd_hits_t* h) {
    mpz_t n, A, B;
    mpz_init_set_str(n, N, 10);
    mpz_init_set_ui(A, a);
    mpz_init_set_ui(B, b);
    fi_rbd_opts_t o;
    memset(&o, 0, sizeof o);
    o.exact = exact;
    o.batch = batch;
    o.record_all = all;
    fi_rbd_scratch_t s;
    fi_rbd_scratch_init(&s);
    int r = fi_rbd_descend(g, k, j, n, A, B, &o, h, &s);
    fi_rbd_scratch_clear(&s);
    mpz_clears(n, A, B, NULL);
    return r;
}

#define RBD_FAIL(...) do { printf("FAIL: rbd: " __VA_ARGS__); printf("\n"); failures++; } while (0)

static void test_rbd_core(void) {
    mpz_t g;
    mpz_init(g);
    uint64_t k = 0;
    int64_t j = 0;

    /* Example 1 (paper): the 255-bit q at descent step 512, misses -2 / 0;
     * the iterated descent unbatched and batched, the exact descent. */
    for (int batch = 0; batch < 2; batch++) {
        if (!rbd_run(g, &k, &j, RBD_NMIN, 22, 7, 0, batch, 0, NULL) || k != 512 || j != -2 ||
            mpz_cmp_str_ok(g, RBD_Q255) != 0)
            RBD_FAIL("Nmin classic batch=%d: k=%llu j=%lld", batch, (unsigned long long)k, (long long)j);
        if (!rbd_run(g, &k, &j, RBD_NMAX, 22, 7, 0, batch, 0, NULL) || k != 512 || j != 0)
            RBD_FAIL("Nmax classic batch=%d: k=%llu j=%lld", batch, (unsigned long long)k, (long long)j);
    }
    if (!rbd_run(g, &k, &j, RBD_NMIN, 22, 7, 1, 1, 0, NULL) || k != 512 || mpz_cmp_str_ok(g, RBD_Q255) != 0)
        RBD_FAIL("Nmin exact: k=%llu", (unsigned long long)k);

    /* Example 2: 2^1041 + 1, base 2 — exactly four distinct gcds. */
    {
        fi_rbd_hits_t h;
        fi_rbd_hits_init(&h);
        mpz_t n;
        mpz_init(n);
        mpz_ui_pow_ui(n, 2, 1041);
        mpz_add_ui(n, n, 1);
        char* ns = mpz_get_str(NULL, 10, n);
        int r = rbd_run(g, NULL, NULL, ns, 2, 1, 0, 1, 1, &h);
        mpz_t w;
        mpz_init(w);
        mpz_ui_pow_ui(w, 2, 347);
        mpz_add_ui(w, w, 1);
        static const struct { unsigned long g; uint64_t k; int64_t j; } want[] = {
            { 3, 1, 1 }, { 9, 1, 4 }, { 2083, 35, -3 }, { 0, 345, 4 } };
        if (r != 4 || h.n != 4) RBD_FAIL("2^1041+1: %zu hits", h.n);
        else
            for (int i = 0; i < 4; i++) {
                int geq = want[i].g ? mpz_cmp_ui(h.v[i].g, want[i].g) == 0 : mpz_cmp(h.v[i].g, w) == 0;
                if (!geq || h.v[i].k != want[i].k || h.v[i].j != want[i].j)
                    RBD_FAIL("2^1041+1 hit %d at (%llu, %lld)", i, (unsigned long long)h.v[i].k, (long long)h.v[i].j);
            }
        /* The coprime refinement of {n, hits}: pairwise coprime divisors. */
        size_t np;
        mpz_t* pc = fi_rbd_coprime_base(&h, n, &np);
        for (size_t a = 0; a < np; a++) {
            if (!mpz_divisible_p(n, pc[a])) RBD_FAIL("coprime piece does not divide n");
            for (size_t b = a + 1; b < np; b++) {
                mpz_gcd(g, pc[a], pc[b]);
                if (mpz_cmp_ui(g, 1) != 0) RBD_FAIL("coprime pieces share a factor");
            }
        }
        fi_rbd_free_mpz_array(pc, np);
        free(ns);
        mpz_clears(n, w, NULL);
        fi_rbd_hits_clear(&h);
    }

    /* Example 3: R_2224, base 10 — 14 distinct gcds, R_1112 at (1111, 1). */
    {
        fi_rbd_hits_t h;
        fi_rbd_hits_init(&h);
        mpz_t n, r1112;
        mpz_inits(n, r1112, NULL);
        mpz_ui_pow_ui(n, 10, 2224); mpz_sub_ui(n, n, 1); mpz_divexact_ui(n, n, 9);
        mpz_ui_pow_ui(r1112, 10, 1112); mpz_sub_ui(r1112, r1112, 1); mpz_divexact_ui(r1112, r1112, 9);
        char* ns = mpz_get_str(NULL, 10, n);
        rbd_run(g, NULL, NULL, ns, 10, 1, 0, 1, 1, &h);
        int at = 0;
        for (size_t i = 0; i < h.n; i++)
            if (mpz_cmp(h.v[i].g, r1112) == 0 && h.v[i].k == 1111 && h.v[i].j == 1) at = 1;
        if (h.n != 14 || !at) RBD_FAIL("R_2224: %zu hits, R_1112 at (1111, 1): %d", h.n, at);
        free(ns);
        mpz_clears(n, r1112, NULL);
        fi_rbd_hits_clear(&h);
    }

    /* Comb sweep (lattice paper, Section 4): the 533 primitive bases with
     * a + b <= 60 against the planted c = 3, n = 150, 22/7 instance; only
     * 22/7 and its harmonic 11/7 fire, both at (150, -1). */
    {
        fi_rbd_sum_t it;
        fi_rbd_sum_init(&it, 0);
        fi_rbd_base_t bs[64];
        size_t nb = 0, got;
        int fired = 0, ok = 1;
        while ((got = fi_rbd_sum_next(&it, bs, 64)) > 0) {
            for (size_t i = 0; i < got; i++) {
                if (bs[i].a + bs[i].b > 60) { got = 0; break; }
                nb++;
                if (rbd_run(g, &k, &j, RBD_NCOMB, (unsigned long)bs[i].a, (unsigned long)bs[i].b, 0, 1, 0, NULL)) {
                    fired++;
                    if (!((bs[i].a == 22 || bs[i].a == 11) && bs[i].b == 7 && k == 150 && j == -1)) ok = 0;
                }
            }
            if (!got) break;
        }
        if (nb != 533 || fired != 2 || !ok) RBD_FAIL("comb sweep: %zu bases, %d fired", nb, fired);
    }

    /* Stern-Brocot order and the perfect-power filter. */
    {
        fi_rbd_sb_t sb;
        fi_rbd_sb_init(&sb, 5);
        fi_rbd_base_t bs[16];
        size_t got = fi_rbd_sb_next(&sb, bs, 16);
        static const uint64_t want[][2] = { {2,1}, {3,2}, {3,1}, {4,3}, {5,3}, {5,2}, {5,4}, {5,1} };
        int ok = got == 8;
        for (size_t i = 0; ok && i < 8; i++) ok = bs[i].a == want[i][0] && bs[i].b == want[i][1];
        if (!ok) RBD_FAIL("Stern-Brocot order (%zu bases)", got);
        fi_rbd_sb_clear(&sb);
        if (fi_rbd_is_primitive(4, 1) || fi_rbd_is_primitive(9, 4) || fi_rbd_is_primitive(27, 8) ||
            !fi_rbd_is_primitive(12, 1) || !fi_rbd_is_primitive(22, 7) || !fi_rbd_is_primitive(8, 3))
            RBD_FAIL("perfect-power filter");
    }

    /* Simplest rational in an interval (brute-force oracle). */
    {
        static const char* T[][3] = {
            { "157/50", "3927/1250", "157/50" }, { "1/3", "1/2", "1/2" }, { "2/5", "3/5", "1/2" },
            { "22/7", "22/7", "22/7" }, { "7/3", "5/2", "5/2" }, { "10001/10000", "5001/5000", "5001/5000" },
            { "3", "7/2", "3" } };
        mpq_t lo, hi, r;
        mpz_t p, q;
        mpq_inits(lo, hi, r, NULL);
        mpz_inits(p, q, NULL);
        for (size_t i = 0; i < sizeof T / sizeof T[0]; i++) {
            mpq_set_str(lo, T[i][0], 10); mpq_canonicalize(lo);
            mpq_set_str(hi, T[i][1], 10); mpq_canonicalize(hi);
            mpq_set_str(r, T[i][2], 10); mpq_canonicalize(r);
            fi_rbd_simplest_between(p, q, lo, hi);
            if (mpz_cmp(p, mpq_numref(r)) || mpz_cmp(q, mpq_denref(r))) RBD_FAIL("simplest in [%s, %s]", T[i][0], T[i][1]);
        }
        mpq_clears(lo, hi, r, NULL);
        mpz_clears(p, q, NULL);
    }

    /* Howgrave-Graham detector at the aligned depth of the planted-miss
     * instance (miss -268435481, 29 bits), on both LLLs. */
    {
        mpz_t N, Q, X, x0, t;
        mpz_inits(N, Q, X, x0, t, NULL);
        mpz_set_str(N, RBD_NLAT, 10);
        mpz_ui_pow_ui(Q, 7, 59); mpz_mul(Q, Q, N); mpz_ui_pow_ui(t, 22, 59); mpz_fdiv_q(Q, Q, t);
        mpz_ui_pow_ui(X, 2, 30);
        if (!fi_rbd_hg_detect(g, x0, N, Q, X, 6, 12) || mpz_cmp_si(x0, -268435481) || mpz_cmp_str_ok(g, RBD_QLAT))
            RBD_FAIL("HG detector (6, 12)");
        if (!fi_rbd_hg_detect_native(g, x0, N, Q, X, 6, 12) || mpz_cmp_si(x0, -268435481))
            RBD_FAIL("HG detector, native LLL");
        mpz_ui_pow_ui(X, 2, 41);
        if (!fi_rbd_hg_detect(g, x0, N, Q, X, 5, 12)) RBD_FAIL("HG detector at 2^41");
        int tt = 0, dd = 0;
        double lx = 0;
        if (!fi_rbd_hg_params(&tt, &dd, &lx, 197, 0.497, -1.0) || dd != 12 || lx < 43 || lx > 44)
            RBD_FAIL("hg_params: t=%d d=%d lx=%g", tt, dd, lx);
        mpz_clears(N, Q, X, x0, t, NULL);
    }
    mpz_clear(g);
}

static void test_rbd_language(void) {
    char buf[8192];
    /* Example 1 through every entry point. */
    snprintf(buf, sizeof buf, "MemberQ[FactorInteger`BlakeRationalBaseDescent[%s, \"Base\" -> 22/7][[All, 1]], %s]", RBD_NMIN, RBD_Q255);
    check(buf, "True");
    snprintf(buf, sizeof buf, "MemberQ[FactorInteger`BlakeRationalBaseDescentClassic[%s, \"Base\" -> 22/7][[All, 1]], %s]", RBD_NMAX, RBD_Q255);
    check(buf, "True");
    snprintf(buf, sizeof buf, "MemberQ[FactorInteger[%s, Method -> {\"BlakeRationalBaseDescent\", \"Variant\" -> \"Classic\", \"Base\" -> 7/22}][[All, 1]], %s]", RBD_NMAX, RBD_Q255);
    check(buf, "True");
    /* The unknown-base search finds the planted 22/7 instance; both orders. */
    snprintf(buf, sizeof buf, "FactorInteger`BlakeRationalBaseDescent[%s]", RBD_NCOMB);
    check(buf, "{{623347347967, 1}, {1190741187132579502456731851112421085210615890827643793135496647841845583593, 1}}");
    snprintf(buf, sizeof buf, "Length[FactorInteger`BlakeRationalBaseDescent[%s, \"Order\" -> \"Sum\"]]", RBD_NCOMB);
    check(buf, "2");
    /* Open Question 2's multiplier; the composite cofactor is split by the
     * cofactor fallback. */
    snprintf(buf, sizeof buf, "FactorInteger`BlakeRationalBaseDescent[%s, \"Base\" -> 3, \"Multiplier\" -> 1443663791]", RBD_NMULT);
    check(buf, "{{24473, 1}, {195311, 1}, {32592575621351777380295131014550050576823494298654980010178247189670100796213387298934358053, 1}}");
    /* Truncation divisibility: the algebraic factors of 2^1041 + 1. */
    check("Select[Quiet[FactorInteger`BlakeRationalBaseDescentAllHits[2^1041 + 1, \"Base\" -> 2]][[All, 1]], # < 10^120 &] == {3, 2083, (2^347 + 1)/3}",
          "True");
    /* Lattice-augmented descent: the planted 2^28 miss defeats the window
     * (nofac) and not the lattice; auto parameters, explicit ones, the
     * fixed-base Coppersmith case, and TwoPass. */
    snprintf(buf, sizeof buf, "Length[Quiet[FactorInteger`BlakeRationalBaseDescent[%s, \"Base\" -> 22/7]]]", RBD_NLAT);
    check(buf, "1");
    const char* lat[] = {
        "FactorInteger`BlakeRationalBaseDescentLattice[%s, \"Base\" -> 22/7]",
        "FactorInteger`BlakeRationalBaseDescentLattice[%s, \"Base\" -> 22/7, \"Tolerance\" -> 2^30, \"Multiplicity\" -> 6, \"Dimension\" -> 12]",
        "FactorInteger`BlakeRationalBaseDescentLattice[%s, \"Base\" -> 22/7, \"Tolerance\" -> 2^36]",
        "FactorInteger`BlakeRationalBaseDescentCoppersmith[%s, \"Base\" -> 22/7, \"Coefficient\" -> 3, \"Exponent\" -> 59]",
        "FactorInteger`BlakeRationalBaseDescentTwoPass[%s, \"Base\" -> {3/2, 22/7}]",
    };
    for (size_t i = 0; i < sizeof lat / sizeof lat[0]; i++) {
        char call[4096];
        snprintf(call, sizeof call, lat[i], RBD_NLAT);
        snprintf(buf, sizeof buf, "First /@ %s", call);
        check(buf, "{296454907933147663952202071141, 659591518271709676788151702157}");
    }
    /* Small inputs, primes, powers. */
    check("FactorInteger`BlakeRationalBaseDescent[15]", "{{3, 1}, {5, 1}}");
    check("FactorInteger`BlakeRationalBaseDescent[2^89 - 1]", "{{618970019642690137449562111, 1}}");
    check("FactorInteger`BlakeRationalBaseDescent[3^10 * 7^4]", "{{3, 10}, {7, 4}}");
    check("FactorInteger`BlakeRationalBaseDescent[1001, \"Base\" -> 44/14]", "{{7, 1}, {11, 1}, {13, 1}}");
    /* Malformed suboptions and the Coppersmith requirements: messages. */
    check("Check[FactorInteger`BlakeRationalBaseDescent[1001, \"Base\" -> 1], failed]", "failed");
    check("Check[FactorInteger`BlakeRationalBaseDescent[1001, \"Variant\" -> \"Nope\"], failed]", "failed");
    check("Check[FactorInteger`BlakeRationalBaseDescent[1001, \"Beta\" -> 2], failed]", "failed");
    check("Check[FactorInteger`BlakeRationalBaseDescentCoppersmith[1001, \"Base\" -> 2], failed]", "failed");
    check("Quiet[FactorInteger`BlakeRationalBaseDescentCoppersmith[1001, \"Base\" -> 2]]",
          "FactorInteger`BlakeRationalBaseDescentCoppersmith[1001, \"Base\" -> 2]");
    check("Attributes[FactorInteger`BlakeRationalBaseDescentLattice]", "{Listable, Protected}");
    /* The parent docstring explains the method (it is new), and every
     * variant head has its own. */
    const char* doc = symtab_get_docstring("FactorInteger`BlakeRationalBaseDescent");
    if (!doc || !strstr(doc, "truncation divisibility") || !strstr(doc, "Howgrave-Graham") ||
        !strstr(doc, "\"MaxCoefficient\"") || strlen(doc) < 3000)
        RBD_FAIL("parent docstring");
    static const char* vh[] = { "Exact", "Classic", "AllHits", "Lattice", "TwoPass", "Coppersmith" };
    for (size_t i = 0; i < 6; i++) {
        char nm[96];
        snprintf(nm, sizeof nm, "FactorInteger`BlakeRationalBaseDescent%s", vh[i]);
        if (!symtab_get_docstring(nm)) RBD_FAIL("no docstring for %s", nm);
    }
}

int main(void) {
    symtab_init();
    core_init();
    test_basics();
    test_partial();
    test_gaussian();
    test_methods();
    test_contract();
    test_serial_stream();
    test_messages();
    test_symbols();
    test_u64();
    test_gf2(300);      /* dense */
    test_gf2(5000);     /* block Lanczos */
    test_rbd_core();
    test_rbd_language();
    if (failures) {
        printf("%d FactorInteger test(s) FAILED\n", failures);
        return 1;
    }
    printf("All FactorInteger tests passed!\n");
    return 0;
}
