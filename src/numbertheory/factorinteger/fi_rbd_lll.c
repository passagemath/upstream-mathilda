/*
 * fi_rbd_lll.c — lattice reduction for the Howgrave-Graham detector of
 * lattice-augmented rational base descent.
 *
 * With USE_FLINT the reduction is FLINT's fmpz_lll (floating-point L^2 with
 * exact fallbacks): the fast path.  The native fallback, always compiled, is
 * the integral (fraction-free) LLL of de Weger / Cohen, "A Course in
 * Computational Algebraic Number Theory", Algorithm 2.6.7, generalised to the
 * Lovasz constant delta = p/q:
 *
 *   d_0 = 1, d_i = det(Gram(b_1..b_i)), lambda_{k,j} = d_j mu_{k,j}
 *
 * are all integers, every update divides exactly, and no rational or
 * floating-point quantity appears, so the result is exact for any entry
 * size.  It is slower than L^2 (its arithmetic is on numbers of about twice
 * the entry size) but the detector's lattices are small (d <= 48, typically
 * 12), so it is a practical fallback rather than a token one.
 */
#include "fi_rbd.h"

#include <stdlib.h>

#ifdef USE_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#endif

/* Lovasz constant delta = LLL_P / LLL_Q = 0.99. */
#define LLL_P 99
#define LLL_Q 100

/* Row-major access (1-indexed rows, 0-indexed columns, as in Cohen). */
#define ROWV(i) (B + (size_t)((i) - 1) * (size_t)cols)

static void dot(mpz_t r, mpz_t* x, mpz_t* y, int cols) {
    mpz_set_ui(r, 0);
    for (int c = 0; c < cols; c++) mpz_addmul(r, x[c], y[c]);
}

/* lambda is an (n+1) x (n+1) lower-triangular table, lam(k, j) for j < k. */
#define LAM(k, j) lam[(size_t)(k) * (size_t)(n + 1) + (size_t)(j)]

/* REDI(k, l): size-reduce b_k against b_l. */
static void redi(mpz_t* B, int cols, mpz_t* lam, mpz_t* d, int n, int k, int l, mpz_t q, mpz_t t) {
    /* |2 lambda_{k,l}| > d_l ? */
    mpz_mul_2exp(t, LAM(k, l), 1);
    mpz_abs(t, t);
    if (mpz_cmp(t, d[l]) <= 0) return;
    /* q = round(lambda_{k,l} / d_l) = floor((2 lambda + d) / (2 d)). */
    mpz_mul_2exp(t, LAM(k, l), 1);
    mpz_add(t, t, d[l]);
    mpz_mul_2exp(q, d[l], 1);
    mpz_fdiv_q(q, t, q);
    mpz_t* bk = ROWV(k);
    mpz_t* bl = ROWV(l);
    for (int c = 0; c < cols; c++) mpz_submul(bk[c], q, bl[c]);
    mpz_submul(LAM(k, l), q, d[l]);
    for (int i = 1; i < l; i++) mpz_submul(LAM(k, i), q, LAM(l, i));
}

/* SWAPI(k): exchange b_k and b_{k-1}, updating lambda and d exactly. */
static void swapi(mpz_t* B, int cols, mpz_t* lam, mpz_t* d, int n, int k, int kmax,
                  mpz_t lm, mpz_t Bv, mpz_t t, mpz_t u) {
    mpz_t* bk = ROWV(k);
    mpz_t* bk1 = ROWV(k - 1);
    for (int c = 0; c < cols; c++) mpz_swap(bk[c], bk1[c]);
    for (int j = 1; j <= k - 2; j++) mpz_swap(LAM(k, j), LAM(k - 1, j));
    mpz_set(lm, LAM(k, k - 1));
    /* B = (d_{k-2} d_k + lambda^2) / d_{k-1} */
    mpz_mul(Bv, d[k - 2], d[k]);
    mpz_addmul(Bv, lm, lm);
    mpz_divexact(Bv, Bv, d[k - 1]);
    for (int i = k + 1; i <= kmax; i++) {
        mpz_set(t, LAM(i, k));
        /* lambda_{i,k} = (d_k lambda_{i,k-1} - lambda t) / d_{k-1} */
        mpz_mul(u, d[k], LAM(i, k - 1));
        mpz_submul(u, lm, t);
        mpz_divexact(LAM(i, k), u, d[k - 1]);
        /* lambda_{i,k-1} = (B t + lambda lambda_{i,k}) / d_k */
        mpz_mul(u, Bv, t);
        mpz_addmul(u, lm, LAM(i, k));
        mpz_divexact(LAM(i, k - 1), u, d[k]);
    }
    mpz_set(d[k - 1], Bv);
}

int fi_rbd_lll_native(mpz_t* B, int rows, int cols) {
    int n = rows;
    if (n <= 1) return 0;
    mpz_t* lam = malloc((size_t)(n + 1) * (size_t)(n + 1) * sizeof(mpz_t));
    mpz_t* d = malloc((size_t)(n + 1) * sizeof(mpz_t));
    if (!lam || !d) { free(lam); free(d); return 1; }
    for (size_t i = 0; i < (size_t)(n + 1) * (size_t)(n + 1); i++) mpz_init(lam[i]);
    for (int i = 0; i <= n; i++) mpz_init(d[i]);
    mpz_t u, q, t, lm, Bv, w;
    mpz_inits(u, q, t, lm, Bv, w, NULL);

    int dependent = 0;
    mpz_set_ui(d[0], 1);
    dot(d[1], ROWV(1), ROWV(1), cols);
    if (mpz_sgn(d[1]) == 0) dependent = 1;
    int k = 2, kmax = 1;
    while (!dependent && k <= n) {
        if (k > kmax) {
            /* Incremental Gram-Schmidt for the new row k. */
            kmax = k;
            for (int j = 1; j <= k; j++) {
                dot(u, ROWV(k), ROWV(j), cols);
                for (int i = 1; i <= j - 1; i++) {
                    mpz_mul(u, u, d[i]);
                    mpz_submul(u, LAM(k, i), LAM(j, i));
                    mpz_divexact(u, u, d[i - 1]);
                }
                if (j < k) mpz_set(LAM(k, j), u);
                else       mpz_set(d[k], u);
            }
            if (mpz_sgn(d[k]) == 0) { dependent = 1; break; }
        }
        for (;;) {
            redi(B, cols, lam, d, n, k, k - 1, q, t);
            /* Lovasz: swap when q d_k d_{k-2} < p d_{k-1}^2 - q lambda^2. */
            mpz_mul(t, d[k], d[k - 2]);
            mpz_mul_ui(t, t, LLL_Q);
            mpz_mul(w, d[k - 1], d[k - 1]);
            mpz_mul_ui(w, w, LLL_P);
            mpz_mul(u, LAM(k, k - 1), LAM(k, k - 1));
            mpz_submul_ui(w, u, LLL_Q);
            if (mpz_cmp(t, w) < 0) {
                swapi(B, cols, lam, d, n, k, kmax, lm, Bv, t, u);
                if (k > 2) k--;
                continue;
            }
            for (int l = k - 2; l >= 1; l--) redi(B, cols, lam, d, n, k, l, q, t);
            k++;
            break;
        }
    }

    mpz_clears(u, q, t, lm, Bv, w, NULL);
    for (size_t i = 0; i < (size_t)(n + 1) * (size_t)(n + 1); i++) mpz_clear(lam[i]);
    for (int i = 0; i <= n; i++) mpz_clear(d[i]);
    free(lam);
    free(d);
    return dependent;
}

int fi_rbd_lll(mpz_t* B, int rows, int cols) {
#ifdef USE_FLINT
    fmpz_mat_t M;
    fmpz_mat_init(M, rows, cols);
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
            fmpz_set_mpz(fmpz_mat_entry(M, r, c), B[(size_t)r * (size_t)cols + (size_t)c]);
    fmpz_lll_t fl;
    fmpz_lll_context_init(fl, 0.99, 0.51, Z_BASIS, APPROX);
    fmpz_lll(M, NULL, fl);
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
            fmpz_get_mpz(B[(size_t)r * (size_t)cols + (size_t)c], fmpz_mat_entry(M, r, c));
    fmpz_mat_clear(M);
    return 0;
#else
    return fi_rbd_lll_native(B, rows, cols);
#endif
}
