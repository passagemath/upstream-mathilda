/*
 * fi_lanczos.c — block Lanczos over GF(2) (Montgomery, 1995).
 *
 * Finds vectors x over the relations with B x = 0, where B is the sparse
 * (primes x relations) exponent-parity matrix, by running Lanczos on the
 * symmetric A = B^T B with 64-bit blocks: each iteration costs two sparse
 * passes over B (one for B v, one for B^T (B v)) plus a few n x 64 by
 * 64 x 64 products, and about n / 63 iterations are needed — O(w n / 64)
 * word operations for w nonzeros, against O(n^3 / 64) for elimination.
 *
 * The recurrence follows Montgomery's paper as implemented in msieve:
 *   v_{i+1} = A v_i S_i S_i^T + v_i D_{i+1} + v_{i-1} E_{i+1} + v_{i-2} F_{i+1}
 * with S_i the columns of v_i^T A v_i chosen so that its restriction is
 * invertible (find_nonsingular_sub), and the solution accumulated as
 * x += v_i W_i^inv v_i^T v_0 from a random x with v_0 = A x.
 *
 * Post-processing: the iteration leaves x and the final v whose 128 columns
 * span (mostly) nullspace vectors of A, not necessarily of B.  Gaussian
 * elimination on the 128 vectors B x_j, B v_j finds the combinations that B
 * maps to zero; each nonzero one is a dependency among the relations.
 */
#include "fi_internal.h"

#include <stdlib.h>
#include <string.h>

typedef uint64_t M64[64];          /* 64 x 64 GF(2) matrix: row i, bit j */

static void mul_64x64(const uint64_t* a, const uint64_t* b, uint64_t* c) {
    uint64_t t[64];
    for (int i = 0; i < 64; i++) {
        uint64_t r = 0, ai = a[i];
        for (int j = 0; ai; j++, ai >>= 1) if (ai & 1) r ^= b[j];
        t[i] = r;
    }
    memcpy(c, t, sizeof t);
}

/* c = x^T y for n x 64 blocks x, y. */
static void tmul_64xN_Nx64(const uint64_t* x, const uint64_t* y, uint64_t* c, size_t n) {
    /* Byte-sliced: for each of the 8 bytes of x[k], accumulate y[k] into a
     * 256-entry table, then fold the tables into rows. */
    static uint64_t tab[8][256];
    memset(tab, 0, sizeof tab);
    for (size_t k = 0; k < n; k++) {
        uint64_t xk = x[k], yk = y[k];
        for (int b = 0; b < 8; b++) tab[b][(xk >> (8 * b)) & 0xff] ^= yk;
    }
    for (int i = 0; i < 64; i++) {
        int b = i / 8, bit = i % 8;
        uint64_t r = 0;
        for (int v = 0; v < 256; v++) if ((v >> bit) & 1) r ^= tab[b][v];
        c[i] = r;
    }
}

/* y ^= v m for an n x 64 block v and 64 x 64 m (8-bit table lookups). */
static void mul_Nx64_64x64_acc(const uint64_t* v, const uint64_t* m, uint64_t* y, size_t n) {
    uint64_t tab[8][256];
    for (int b = 0; b < 8; b++) {
        tab[b][0] = 0;
        for (int vv = 1; vv < 256; vv++) {
            int low = vv & -vv, bit = 0;
            while (!((low >> bit) & 1)) bit++;
            tab[b][vv] = tab[b][vv ^ low] ^ m[8 * b + bit];
        }
    }
    for (size_t k = 0; k < n; k++) {
        uint64_t vk = v[k], r = 0;
        for (int b = 0; b < 8; b++) r ^= tab[b][(vk >> (8 * b)) & 0xff];
        y[k] ^= r;
    }
}

/* y = B x (length ncols), z = B^T y (length nrows): z = A x. */
static void mul_A(const FiSparseMat* m, const uint64_t* x, uint64_t* z, uint64_t* y) {
    memset(y, 0, (size_t)m->ncols * sizeof(uint64_t));
    for (uint32_t r = 0; r < m->nrows; r++) {
        uint64_t xr = x[r];
        for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) y[m->col[k]] ^= xr;
    }
    for (uint32_t r = 0; r < m->nrows; r++) {
        uint64_t acc = 0;
        for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) acc ^= y[m->col[k]];
        z[r] = acc;
    }
}

static void mul_B(const FiSparseMat* m, const uint64_t* x, uint64_t* y) {
    memset(y, 0, (size_t)m->ncols * sizeof(uint64_t));
    for (uint32_t r = 0; r < m->nrows; r++) {
        uint64_t xr = x[r];
        for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) y[m->col[k]] ^= xr;
    }
}

static uint32_t find_nonsingular_sub(const uint64_t* t, uint32_t* s, const uint32_t* last_s,
                                     uint32_t last_dim, uint64_t* w) {
    uint64_t M[64][2];
    for (int i = 0; i < 64; i++) { M[i][0] = t[i]; M[i][1] = UINT64_C(1) << i; }
    uint64_t mask = 0;
    for (uint32_t i = 0; i < last_dim; i++) mask |= UINT64_C(1) << last_s[i];
    uint32_t j = 0;
    for (uint32_t i = 0; i < 64; i++) if (!(mask & (UINT64_C(1) << i))) s[j++] = i;
    for (uint32_t i = 0; i < last_dim; i++) s[j++] = last_s[i];
    uint32_t dim = 0;
    for (uint32_t i = 0; i < 64; i++) {
        uint64_t bm = UINT64_C(1) << s[i];
        uint64_t* ri = M[s[i]];
        for (j = i; j < 64; j++) {
            uint64_t* rj = M[s[j]];
            if (rj[0] & bm) {
                uint64_t m0 = rj[0], m1 = rj[1];
                rj[0] = ri[0]; rj[1] = ri[1]; ri[0] = m0; ri[1] = m1;
                break;
            }
        }
        if (j < 64) {
            for (j = 0; j < 64; j++) {
                uint64_t* rj = M[s[j]];
                if (rj != ri && (rj[0] & bm)) { rj[0] ^= ri[0]; rj[1] ^= ri[1]; }
            }
            s[dim++] = s[i];
            continue;
        }
        for (j = i; j < 64; j++) {
            uint64_t* rj = M[s[j]];
            if (rj[1] & bm) {
                uint64_t m0 = rj[0], m1 = rj[1];
                rj[0] = ri[0]; rj[1] = ri[1]; ri[0] = m0; ri[1] = m1;
                break;
            }
        }
        if (j == 64) return 0;                  /* submatrix not invertible */
        for (j = 0; j < 64; j++) {
            uint64_t* rj = M[s[j]];
            if (rj != ri && (rj[1] & bm)) { rj[0] ^= ri[0]; rj[1] ^= ri[1]; }
        }
        ri[0] = ri[1] = 0;
    }
    for (int i = 0; i < 64; i++) w[i] = M[i][1];
    mask = 0;
    for (uint32_t i = 0; i < dim; i++) mask |= UINT64_C(1) << s[i];
    for (uint32_t i = 0; i < last_dim; i++) mask |= UINT64_C(1) << last_s[i];
    if (mask != ~UINT64_C(0)) return 0;         /* recurrence needs every column */
    return dim;
}

static uint64_t lrand(uint64_t* s) {
    uint64_t x = *s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    *s = x;
    return x * UINT64_C(0x2545F4914F6CDD1D);
}

/* Combine the 128 columns of x and v into vectors with B y = 0. */
static size_t combine(const FiSparseMat* m, const uint64_t* x, const uint64_t* v,
                      uint64_t* deps, size_t max_deps) {
    size_t n = m->nrows, nc = m->ncols, W = (nc + 63) / 64;
    uint64_t* bx = malloc(nc * sizeof(uint64_t));
    uint64_t* bv = malloc(nc * sizeof(uint64_t));
    mul_B(m, x, bx);
    mul_B(m, v, bv);
    /* 128 column vectors over the nc primes, as bitsets; plus 128-bit
     * combination trackers. */
    uint64_t* col = calloc(128 * W, sizeof(uint64_t));
    uint64_t comb[128][2];
    for (size_t k = 0; k < nc; k++)
        for (int j = 0; j < 64; j++) {
            if ((bx[k] >> j) & 1) col[(size_t)j * W + k / 64] |= UINT64_C(1) << (k % 64);
            if ((bv[k] >> j) & 1) col[(size_t)(64 + j) * W + k / 64] |= UINT64_C(1) << (k % 64);
        }
    for (int j = 0; j < 128; j++) { comb[j][0] = j < 64 ? UINT64_C(1) << j : 0; comb[j][1] = j >= 64 ? UINT64_C(1) << (j - 64) : 0; }
    size_t piv_row[128];
    int npiv = 0;
    int piv_idx[128];
    size_t nd = 0;
    size_t Wn = (n + 63) / 64;
    for (int j = 0; j < 128 && nd < max_deps; j++) {
        uint64_t* cj = col + (size_t)j * W;
        for (int p = 0; p < npiv; p++) {
            size_t r = piv_row[p];
            if ((cj[r / 64] >> (r % 64)) & 1) {
                const uint64_t* cp = col + (size_t)piv_idx[p] * W;
                for (size_t q = 0; q < W; q++) cj[q] ^= cp[q];
                comb[j][0] ^= comb[piv_idx[p]][0];
                comb[j][1] ^= comb[piv_idx[p]][1];
            }
        }
        size_t first = 0;
        int zero = 1;
        for (size_t q = 0; q < W; q++) if (cj[q]) { zero = 0; first = q * 64 + (size_t)fi_ctz64(cj[q]); break; }
        if (!zero) { piv_row[npiv] = first; piv_idx[npiv] = j; npiv++; continue; }
        /* B (combination) = 0: materialise it over the relations. */
        uint64_t* out = deps + nd * Wn;
        memset(out, 0, Wn * sizeof(uint64_t));
        int any = 0;
        for (size_t k = 0; k < n; k++) {
            uint64_t a = (x[k] & comb[j][0]) ^ (v[k] & comb[j][1]);
            if (fi_parity64(a)) { out[k / 64] |= UINT64_C(1) << (k % 64); any = 1; }
        }
        if (any) nd++;
    }
    free(bx); free(bv); free(col);
    return nd;
}

size_t fi_gf2_lanczos(const FiSparseMat* m, uint64_t* deps, size_t max_deps);

size_t fi_gf2_lanczos(const FiSparseMat* m, uint64_t* deps, size_t max_deps) {
    size_t n = m->nrows;
    uint64_t* buf = calloc(9 * n + m->ncols + 8, sizeof(uint64_t));
    uint64_t* x = buf;
    uint64_t* v[3] = { buf + n, buf + 2 * n, buf + 3 * n };
    uint64_t* vnext = buf + 4 * n;
    uint64_t* v0 = buf + 5 * n;
    uint64_t* tmp = buf + 6 * n;
    uint64_t* y = buf + 9 * n;
    uint64_t seed = UINT64_C(0x1234567887654321) ^ n;
    size_t nd = 0;
    for (int attempt = 0; attempt < 4 && nd == 0; attempt++) {
        M64 winv[3], vt_a_v[2], vt_a2_v[2];
        uint64_t d[64], e[64], f[64], f2[64];
        uint32_t s[2][64], dim0 = 0, dim1 = 64;
        uint64_t mask0 = 0, mask1 = ~UINT64_C(0);
        for (size_t k = 0; k < n; k++) x[k] = lrand(&seed);
        mul_A(m, x, v[0], y);
        memcpy(v0, v[0], n * sizeof(uint64_t));
        memset(v[1], 0, n * sizeof(uint64_t));
        memset(v[2], 0, n * sizeof(uint64_t));
        memset(winv, 0, sizeof winv);
        memset(vt_a_v, 0, sizeof vt_a_v);
        memset(vt_a2_v, 0, sizeof vt_a2_v);
        for (int i = 0; i < 64; i++) s[1][i] = (uint32_t)i;
        int ok = 1;
        size_t max_iter = n / 60 + 100;
        for (size_t iter = 0;; iter++) {
            if (iter > max_iter) { ok = 0; break; }
            mul_A(m, v[0], vnext, y);
            tmul_64xN_Nx64(v[0], vnext, vt_a_v[0], n);
            tmul_64xN_Nx64(vnext, vnext, vt_a2_v[0], n);
            int allzero = 1;
            for (int i = 0; i < 64; i++) if (vt_a_v[0][i]) { allzero = 0; break; }
            if (allzero) break;
            dim0 = find_nonsingular_sub(vt_a_v[0], s[0], s[1], dim1, winv[0]);
            if (dim0 == 0) { ok = 0; break; }
            mask0 = 0;
            for (uint32_t i = 0; i < dim0; i++) mask0 |= UINT64_C(1) << s[0][i];
            if (mask0 != ~UINT64_C(0)) for (size_t k = 0; k < n; k++) vnext[k] &= mask0;
            /* d = I - Winv0 (vt_a2_v0 mask0 + vt_a_v0) */
            for (int i = 0; i < 64; i++) d[i] = (vt_a2_v[0][i] & mask0) ^ vt_a_v[0][i];
            mul_64x64(winv[0], d, d);
            for (int i = 0; i < 64; i++) d[i] ^= UINT64_C(1) << i;
            /* e = - Winv1 vt_a_v0 mask0 */
            mul_64x64(winv[1], vt_a_v[0], e);
            for (int i = 0; i < 64; i++) e[i] &= mask0;
            /* f = - Winv2 (I - vt_a_v1 Winv1) (vt_a2_v1 mask1 + vt_a_v1) mask0 */
            mul_64x64(vt_a_v[1], winv[1], f);
            for (int i = 0; i < 64; i++) f[i] ^= UINT64_C(1) << i;
            mul_64x64(winv[2], f, f);
            for (int i = 0; i < 64; i++) f2[i] = ((vt_a2_v[1][i] & mask1) ^ vt_a_v[1][i]) & mask0;
            mul_64x64(f, f2, f);
            mul_Nx64_64x64_acc(v[0], d, vnext, n);
            mul_Nx64_64x64_acc(v[1], e, vnext, n);
            mul_Nx64_64x64_acc(v[2], f, vnext, n);
            /* x += v0 Winv0 (v0^T v_initial) */
            tmul_64xN_Nx64(v[0], v0, d, n);
            mul_64x64(winv[0], d, d);
            mul_Nx64_64x64_acc(v[0], d, x, n);
            /* rotate */
            uint64_t* t = v[2]; v[2] = v[1]; v[1] = v[0]; v[0] = vnext; vnext = t;
            memcpy(winv[2], winv[1], sizeof(M64));
            memcpy(winv[1], winv[0], sizeof(M64));
            memcpy(vt_a_v[1], vt_a_v[0], sizeof(M64));
            memcpy(vt_a2_v[1], vt_a2_v[0], sizeof(M64));
            memcpy(s[1], s[0], sizeof s[0]);
            mask1 = mask0;
            dim1 = dim0;
        }
        if (ok) nd = combine(m, x, v[0], deps, max_deps);
        (void)tmp;
    }
    free(buf);
    return nd;
}
