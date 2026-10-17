/* Machine-precision (raw double re/im) Meijer G-function.
 *
 *   MeijerG[{{a1..an},{an+1..ap}}, {{b1..bm},{bm+1..bq}}, z, r]
 *
 * This is the ONE implementation behind every machine surface of MeijerG: the
 * scalar builtin (src/special_functions/meijer_g.c), its NDArray / packed
 * buffer path, and the n-ary kernel registered in ndkernels.c that Compile[]
 * and auto-compilation reach. Sharing it is what makes scalar == List ==
 * packed == NDArray == Compile bit-for-bit (CLAUDE.md).
 *
 * Method (Slater's theorem). With w = z^(1/r) = exp(log z / r):
 *
 *   p < q, or p == q and |w| < 1:  the sum over h = 1..m of
 *
 *     C_h  w^{b_h}  pFq~(1+b_h-a_*; 1+b_h-b_{*!=h}; (-1)^{p-m-n} w)
 *
 *     C_h = prod_{j<=m, j!=h} pi / sin(pi (b_j - b_h))
 *           prod_{j<=n} Gamma(1+b_h-a_j) / prod_{j>n} Gamma(a_j-b_h)
 *
 *   where pFq~ is the REGULARIZED hypergeometric function (so lower parameters
 *   that are non-positive integers are finite, not poles). The pi/sin factor is
 *   Gamma(b_j-b_h) Gamma(1+b_h-b_j) -- the regularization of the m-1 lower
 *   parameters that come from the b_j, j <= m -- folded into one closed form.
 *
 *   p > q, or p == q and |w| > 1: the same sum applied to the dual function,
 *     G^{m,n}_{p,q}(a; b; w) = G^{n,m}_{q,p}(1-b; 1-a; 1/w).
 *
 * Everything that does not depend on z (the C_h, the shifted parameter lists,
 * the regularization shifts) is computed ONCE into an SfMgPlan, so a buffer or
 * a compiled loop pays only the per-point pFq evaluations.
 *
 * The log case (two of b_1..b_m an integer apart, so the pi/sin factors of C_h
 * have poles that cancel between terms) is not a Slater sum of pFq's at all;
 * it is summed directly as the residue series of the Mellin-Barnes integrand
 *
 *   F(s) = prod_{j<=m} Gamma(b_j - s) prod_{j<=n} Gamma(1 - a_j + s)
 *          / ( prod_{j>m} Gamma(1 - b_j + s) prod_{j>n} Gamma(a_j - s) ),
 *
 *   G = - sum over the poles s* of Res F(s) w^s,
 *
 * grouping b_1..b_m into classes of integer-spaced members with base beta (the
 * member of least real part). The poles of a class sit at s* = beta + k, of
 * order up to the class size T, so the residue there is
 *
 *   w^{beta+k} sum_{j<T} c[k][j] (log w)^j,   c[k][j] = -f^{(k)}_{-1-j} / j!
 *
 * with f^{(k)} the Laurent coefficients of F about beta + k. Those are
 * z-independent, so the plan tabulates them once (SfMgLogTab): the expansion
 * about beta comes from arb (flint_mg_laurent: Gamma / 1/Gamma series), and
 * each next one from the exact step F(s+1) = F(s) R(s),
 *
 *   R(s) = prod_{j<=n} (1-a_j+s) prod_{j>n} (a_j-s-1)
 *          / ( prod_{j<=m} (b_j-s-1) prod_{j>m} (1-b_j+s) ),
 *
 * as truncated power series in eps carried in double-double (a linear factor
 * whose constant term is an exact zero only moves the pole order). Per point
 * the series is summed in double while measuring its own cancellation (the
 * largest-term mass against the sum); past SF_MG_LOG_D_LOST_BITS it is re-summed
 * with log w and w in double-double (sf_dd.h), which absorbs up to
 * SF_MG_DD_MAX_LOST_BITS -- the e^{+|w|} terms of a decaying e^{-|w|} G at
 * large |w|. Integer-valued parameter differences are snapped to the exact
 * integer by the same is_int(double difference) rule that detected the case.
 *
 * Decline contract: sf_mg_eval_machine returns false when it cannot produce a
 * machine answer it trusts -- more than SF_MG_MAX_LOST_BITS of cancellation
 * between the Slater terms (or between the log-case classes), a pFq kernel
 * declining, a log-case series that loses more than SF_MG_DD_MAX_LOST_BITS or
 * runs past its SF_MG_LOG_KMAX tabulated terms (|w| near 1 when p == q), a
 * build without FLINT (no log table), or overflow. sf_mg_eval then retries the
 * same point through the arbitrary-precision arb engine (flint_mg_eval_d) at 53
 * good bits, so the surfaces never disagree and a declined point never becomes
 * a NaN.
 *
 * Complex arithmetic is spelled out on (re, im) pairs: no <complex.h>.
 */
#ifndef MATHILDA_SF_MACHINE_MEIJERG_H
#define MATHILDA_SF_MACHINE_MEIJERG_H

#include <stdbool.h>
#include <stddef.h>

/* Largest p and q handled on the machine path (and by the n-ary kernel). */
#define SF_MG_MAXP 12

/* Cancellation (bits) between the Slater terms beyond which the machine sum is
 * not trusted and the point is re-evaluated in arb. */
#define SF_MG_MAX_LOST_BITS 6.0

/* Log case: the cancellation a double-double re-sum still certifies to double
 * precision (106 - 53 - 8 guard bits), and the number of residue terms the
 * plan tabulates per class. */
#define SF_MG_DD_MAX_LOST_BITS 45.0
#define SF_MG_LOG_D_LOST_BITS 3.0
#define SF_MG_LOG_KMAX 256

/* One integer-spaced class of b_1..b_m: residues at w^{beta + k}, k < K, each
 * a polynomial of degree < T in log w. c holds K * T complex double-doubles,
 * c[4 (k T + j) + 0..3] = (re.hi, re.lo, im.hi, im.lo) of c[k][j]. */
typedef struct {
    double  br, bi;             /* beta */
    size_t  T, K;
    double* c;
} SfMgLogClass;

typedef struct {
    size_t       nc;
    SfMgLogClass cl[SF_MG_MAXP];
} SfMgLogTab;

/* One Slater expansion (in w, or in 1/w for the dual): nt terms of the shape
 * coef * w^{e} * x^{k} * pFq(u; l; x), x = sgn * w. */
typedef struct {
    bool   valid;               /* parameters admit this expansion */
    bool   log_case;            /* integer-differing b's: summed from lt instead */
    SfMgLogTab* lt;             /* the log-case residue table (owned; NULL: decline) */
    size_t nt;                  /* number of terms (m, or n for the dual) */
    size_t pu, ql;              /* pFq shape of every term */
    double sgn;                 /* (-1)^{p-m-n} */
    double cr[SF_MG_MAXP], ci[SF_MG_MAXP];   /* C_h (regularization folded in) */
    double er[SF_MG_MAXP], ei[SF_MG_MAXP];   /* exponent of w */
    int    kx[SF_MG_MAXP];                   /* extra integer power of x (shift) */
    double ur[SF_MG_MAXP][SF_MG_MAXP], ui[SF_MG_MAXP][SF_MG_MAXP];  /* upper */
    double lr[SF_MG_MAXP][SF_MG_MAXP], li[SF_MG_MAXP][SF_MG_MAXP];  /* lower */
    bool   real;                /* every coefficient and parameter is real */
} SfMgSeries;

typedef struct {
    size_t n, p, m, q;          /* MeijerG shape */
    double r;                   /* the generalized form's r (1 by default) */
    bool   ok;                  /* the plan was built */
    bool   real_params;         /* every a_j, b_j is real */
    /* The raw parameters, kept for the arb fallback. a[0..p), b[0..q). */
    double ar[SF_MG_MAXP], ai[SF_MG_MAXP], br[SF_MG_MAXP], bi[SF_MG_MAXP];
    SfMgSeries zs;              /* expansion in w   */
    SfMgSeries ds;              /* expansion in 1/w */
} SfMgPlan;

/* Build the plan for MeijerG[{{a[0..n)},{a[n..p)}}, {{b[0..m)},{b[m..q)}}, ., r].
 * ai / bi may be NULL for real parameters. Returns false (plan unusable) when
 * p or q exceeds SF_MG_MAXP, m > q, n > p, a parameter is not finite, or r is
 * zero / not finite. A plan that builds may still decline individual points.
 * A log-case plan owns heap tables: release every plan with sf_mg_plan_free
 * (also after a failed build, which is harmless). */
bool sf_mg_plan(const double* ar, const double* ai, size_t n, size_t p,
                const double* br, const double* bi, size_t m, size_t q,
                double r, SfMgPlan* P);
void sf_mg_plan_free(SfMgPlan* P);

/* The plan for these parameters from the process-wide memo, building it on a
 * miss. A built plan is published to the memo when `memo_all` is set or when
 * it carries a log table (the expensive build); otherwise, or once the memo is
 * full, it is built into *scratch. Returns the plan to use, or NULL when it
 * does not build; when the result is `scratch` the caller frees it with
 * sf_mg_plan_free. Memoised plans are immutable and live for the process, so
 * the result may be shared across threads. */
const SfMgPlan* sf_mg_plan_get(const double* ar, const double* ai, size_t n, size_t p,
                               const double* br, const double* bi, size_t m, size_t q,
                               double r, bool memo_all, SfMgPlan* scratch);

/* The pure machine evaluation at z = zr + i zi. *is_c reports whether the
 * value is complex (false only for the all-real lane: real parameters and
 * z > 0). Returns false to decline (see the contract above). */
bool sf_mg_eval_machine(const SfMgPlan* P, double zr, double zi,
                        double* rr, double* ri, bool* is_c);

/* The full machine-precision answer: the machine evaluation, and on decline
 * the arb engine rounded to double. Returns false only when the value does not
 * exist (a genuine pole / ComplexInfinity / non-finite input) or, in a build
 * without FLINT, when the machine path declined. */
bool sf_mg_eval(const SfMgPlan* P, double zr, double zi,
                double* rr, double* ri, bool* is_c);

/* Flattened n-ary layout shared by the ndkernels.c kernel and Compile[]:
 *   v[0] = n + 64 p + 4096 m   (q is implied by the operand count)
 *   v[1] = r
 *   v[2 .. 2+p)       a_1 .. a_p
 *   v[2+p .. 2+p+q)   b_1 .. b_q
 *   v[2+p+q]          z
 * so count = p + q + 3. */
double sf_mg_shape_code(size_t n, size_t p, size_t m);
bool   sf_mg_eval_flat(const double* re, const double* im, size_t count,
                       double* rr, double* ri);

#endif /* MATHILDA_SF_MACHINE_MEIJERG_H */
