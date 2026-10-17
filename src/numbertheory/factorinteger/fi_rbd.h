/*
 * fi_rbd.h — Blake rational base descent (RBD): the reentrant core.
 *
 * References (factorisation_temp/ in the source tree):
 *   [RBD]  S. Blake, "Rational Base Descent: A Deterministic Algorithm for
 *          Factoring Structured Semiprimes".
 *   [LRBD] S. Blake, "Lattice-Augmented Rational Base Descent".
 *
 * The method factors N = P M when the larger factor is P ~ c (a/b)^n for a
 * rational base a/b > 1.  Descending Q <- floor(b Q / a) from N, after n
 * steps Q lies within a constant of c M, and gcd(Q - j, N) for a small j
 * exposes M ([RBD] Theorem 1).  Everything in this header is plain GMP (the
 * lattice reduction uses FLINT's fmpz_lll when built with USE_FLINT, else a
 * native integral LLL): no global state, no I/O, no allocation shared
 * between calls, so every function is reentrant and safe to run on worker
 * threads with per-thread scratch.  The Mathilda glue (option parsing,
 * threads, the FactorInteger driver) lives in fi_rbd.c.
 *
 * Conventions.  A base is a coprime pair a > b >= 1.  A hit is reported as
 * (k, j, g): at descent depth k the candidate Q_k - j had g = gcd(Q_k - j, N)
 * with 1 < g < N ([RBD]'s (k, j) notation; for the exact descent Q_k is the
 * scaled value Qhat_k of [RBD] Section 5.1).
 */
#ifndef MATHILDA_FI_RBD_H
#define MATHILDA_FI_RBD_H

#include <gmp.h>
#include <stddef.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/* Hit records                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t    g;     /* 1 < g < N */
    uint64_t k;     /* descent depth */
    int64_t  j;     /* window offset: the candidate was Q_k - j */
} fi_rbd_hit_t;

typedef struct {
    fi_rbd_hit_t* v;
    size_t        n, cap;
} fi_rbd_hits_t;

void fi_rbd_hits_init(fi_rbd_hits_t* h);
void fi_rbd_hits_clear(fi_rbd_hits_t* h);
/* Append (g, k, j) unless an equal g is already recorded.  Returns 1 when
 * the hit was new. */
int  fi_rbd_hits_add(fi_rbd_hits_t* h, const mpz_t g, uint64_t k, int64_t j);

/* Refine the recorded hits into pairwise-coprime divisors of N whose
 * product is N (factor refinement: every hit, and N divided by it, is a
 * product of the returned pieces).  Pieces equal to 1 are dropped; the
 * array is malloc'd, *count set.  Returns NULL when there are no hits. */
mpz_t* fi_rbd_coprime_base(const fi_rbd_hits_t* h, const mpz_t N, size_t* count);
void   fi_rbd_free_mpz_array(mpz_t* v, size_t n);

/* ------------------------------------------------------------------ */
/* Descent                                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    /* 0: the iterated descent of [RBD] R1-R4 with the symmetric window
     *    |j| <= J = ceil(2 + a/(a-b)) (the paper's reference algorithm);
     * 1: the scaled exact descent of [RBD] Section 5.1: the descent runs on
     *    S N with S = 2^s >= 2a/(a-b), and the one-sided base-independent
     *    window {Qhat, ..., Qhat + 3} is tested (j = 0, -1, -2, -3). */
    int        exact;
    /* Accumulate window products modulo N and take one gcd per block of
     * steps ([RBD] "batched gcds"); a block whose gcd is non-trivial is
     * replayed candidate by candidate, so the hits found are exactly those
     * of the unbatched search. */
    int        batch;
    /* Run to the end of the descent and record every distinct hit in the
     * caller's hit list (the truncation-divisibility mode of [RBD] Section
     * 6); otherwise stop at the first hit. */
    int        record_all;
    /* Descend from mult * N instead of N ([RBD] Open Question 2); gcds are
     * still taken with N.  NULL means 1. */
    mpz_srcptr mult;
    /* Deferred testing ([RBD] Section 5.1): no candidate is tested while
     * Q > ceiling.  The exact descent also jumps straight to the first
     * depth below the ceiling with one division.  NULL tests every depth. */
    mpz_srcptr ceiling;
} fi_rbd_opts_t;

/* Per-thread scratch, reusable across bases. */
typedef struct {
    mpz_t  Q, Qh, t, u, v, acc, g, Neff, Nr;
    mpz_t* blk;              /* replay buffer: the tested Q of each step */
    uint64_t* blk_k;
    size_t blk_cap;
} fi_rbd_scratch_t;

void fi_rbd_scratch_init(fi_rbd_scratch_t* s);
void fi_rbd_scratch_clear(fi_rbd_scratch_t* s);

/* The search radius J = ceil(2 + a/(a-b)) of [RBD] R3, or 0 when it does
 * not fit in 2^24 (the iterated descent then declines). */
uint64_t fi_rbd_window(const mpz_t a, const mpz_t b);

/* The deferred-testing ceiling floor(sqrt(cmax * mult * N)) + slack: no
 * target multiple c M with c <= cmax can lie above it, since P > c M
 * forces (c M)^2 < c N ([RBD] Section 5.1). */
void fi_rbd_ceiling(mpz_t T, const mpz_t N, mpz_srcptr mult, const mpz_t cmax, uint64_t slack);

/* Run one descent of N (odd, composite, > 3) with base a/b (coprime,
 * a > b >= 1).  On a hit sets g (and *k, *j when non-NULL) and returns 1;
 * with o->record_all it returns the number of new hits added to `hits` (g
 * is then the first of them).  Returns 0 when the descent finds nothing. */
int fi_rbd_descend(mpz_t g, uint64_t* k, int64_t* j, const mpz_t N, const mpz_t a,
                   const mpz_t b, const fi_rbd_opts_t* o, fi_rbd_hits_t* hits,
                   fi_rbd_scratch_t* s);

/* The exact descent values Q'_k = floor(mult N b^k / a^k) (to within one:
 * the scaled value floor(Qtilde_k / S)) for every depth k with
 * lo <= Q'_k <= hi, in descending order of Q (ascending k).  The arrays are
 * malloc'd; returns the count.  Used by the lattice descent, which tests
 * only the band X <~ Q <~ sqrt(cmax N) + X of [LRBD] Section 6. */
size_t fi_rbd_band(mpz_t** Qs, uint64_t** ks, const mpz_t N, mpz_srcptr mult,
                   const mpz_t a, const mpz_t b, const mpz_t lo, const mpz_t hi);

/* ------------------------------------------------------------------ */
/* Base enumeration (fi_rbd_bases.c)                                   */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t a, b; } fi_rbd_base_t;

/* 0 when a = x^k and b = y^k for some k >= 2 (the base is then redundant
 * with x/y: (a/b)^n = (x/y)^(kn), [RBD] Section 4.1). */
int fi_rbd_is_primitive(uint64_t a, uint64_t b);

/* Stern-Brocot enumeration of the primitive bases a/b > 1 with a <= H, in
 * breadth-first (coarse-to-fine) order: every rational appears once, in
 * lowest terms ([RBD] Section 4.1, [LRBD] Section 6). */
typedef struct {
    uint64_t* q;            /* ring of intervals (pl, ql, pr, qr) */
    size_t    head, len, cap;
    uint64_t  H;
} fi_rbd_sb_t;
void   fi_rbd_sb_init(fi_rbd_sb_t* it, uint64_t H);
void   fi_rbd_sb_clear(fi_rbd_sb_t* it);
/* Up to max bases into out; returns the count (0 once exhausted). */
size_t fi_rbd_sb_next(fi_rbd_sb_t* it, fi_rbd_base_t* out, size_t max);

/* The reference order of [RBD] Appendix A: a + b = 3, 4, ..., b ascending,
 * coprime and primitive pairs only, a <= H (H = 0: unbounded). */
typedef struct { uint64_t s, b, H; } fi_rbd_sum_t;
void   fi_rbd_sum_init(fi_rbd_sum_t* it, uint64_t H);
size_t fi_rbd_sum_next(fi_rbd_sum_t* it, fi_rbd_base_t* out, size_t max);

/* The simplest rational p/q (least q, then least p) in the closed interval
 * [lo, hi], 0 < lo <= hi: the Stern-Brocot convergent that realises a grid
 * point of [LRBD] pass P2 with the lowest possible height. */
void fi_rbd_simplest_between(mpz_t p, mpz_t q, mpq_srcptr lo, mpq_srcptr hi);

/* ------------------------------------------------------------------ */
/* Lattice detector (fi_rbd_lattice.c, fi_rbd_lll.c)                   */
/* ------------------------------------------------------------------ */

/* LLL-reduce the rows of the rows x cols integer matrix B (row-major,
 * initialised mpz_t) in place, delta = 0.99.  FLINT's fmpz_lll when built
 * with USE_FLINT, else a native integral (fraction-free) LLL.  Returns 0 on
 * success, 1 when the rows are linearly dependent (B is then unspecified). */
int fi_rbd_lll(mpz_t* B, int rows, int cols);
/* The native reduction, always available (the FLINT build uses it only as
 * a cross-check in the tests). */
int fi_rbd_lll_native(mpz_t* B, int rows, int cols);

/* Choose Howgrave-Graham parameters.  With log2X >= 0, the least dimension d
 * (and its multiplicity t) whose LLL guarantee covers tolerance 2^log2X for
 * a divisor M >= N^beta; with log2X < 0, d = *d (0: 12) and the t that
 * maximises the tolerance, which is returned in *log2X_out.  Returns 0 when
 * no d <= 48 suffices. */
int fi_rbd_hg_params(int* t, int* d, double* log2X_out, double log2N, double beta,
                     double log2X);

/* Howgrave-Graham's partially approximate common divisor detector
 * ([LRBD] Theorem 2): find x0 with |x0| <= X such that Q + x0 shares a
 * divisor M >= N^beta with N, using the lattice of N^(t-i) (Q + x)^i
 * (0 <= i <= t) and x^j (Q + x)^t (1 <= j < d - t).  On success sets
 * g = gcd(Q + x0, N) (1 < g < N) and x0, and returns 1. */
int fi_rbd_hg_detect(mpz_t g, mpz_t x0, const mpz_t N, const mpz_t Q, const mpz_t X,
                     int t, int d);
/* The same detector, always on the native LLL (a cross-check for builds
 * where fi_rbd_lll dispatches to FLINT). */
int fi_rbd_hg_detect_native(mpz_t g, mpz_t x0, const mpz_t N, const mpz_t Q, const mpz_t X,
                            int t, int d);

/* Integer roots x of the polynomial h[0] + h[1] x + ... + h[deg] x^deg
 * (h[deg] != 0) with |x| <= bound, found p-adically (roots modulo a small
 * prime, Hensel-lifted past 2 bound).  Writes at most deg roots to roots[]
 * (initialised mpz_t) and returns their count. */
int fi_rbd_int_roots(mpz_t* roots, const mpz_t* h, int deg, const mpz_t bound);

#endif /* MATHILDA_FI_RBD_H */
