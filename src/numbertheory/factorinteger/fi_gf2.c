/*
 * fi_gf2.c — relations, GF(2) linear algebra and the square-root step shared
 * by the congruence-of-squares methods (CFRAC, Dixon, the quadratic sieve).
 *
 * A relation is X^2 == sign * prod p_i^e_i * E^2 (mod n) over a factor base,
 * where E collects squared large primes from combined partial relations.
 * Partial relations (one large prime L outside the factor base) are kept in
 * a hash table keyed by L; a second partial with the same L combines with
 * the first into a full relation with E = L.
 *
 * Linear algebra finds subsets of relations whose exponent vectors sum to
 * zero mod 2.  The matrix is first pruned (singleton columns removed
 * repeatedly, since a relation holding the only odd occurrence of a prime
 * can never be in a dependency), then the survivors are solved:
 *   - dense Gaussian elimination on bit-packed columns below DENSE_MAX
 *     columns: simple, cache friendly, cost ~ c^3/64 word ops;
 *   - block Lanczos (Montgomery 1995, 64-bit blocks) above, O(w c^2 / 64)
 *     for w nonzeros, the only practical solver for the 10^4..10^5-column
 *     matrices of a 70..100 digit quadratic sieve.
 * Each dependency gives X = prod X_r, Y = sign-free sqrt of the product of
 * the right-hand sides (exponents halved) times prod E_r, and gcd(X - Y, n)
 * splits n with probability 1/2 per prime pair.
 */
#include "fi_internal.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Relation store                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t     X, E;
    uint32_t  off, cnt;      /* into the shared (col, exp) pool */
} Rel;

typedef struct {
    uint64_t  L;             /* large prime, 0 = empty slot */
    mpz_t     X;
    uint32_t  off, cnt;
} Partial;

struct FiRels {
    size_t    ncols;         /* factor-base columns, column 0 = sign */
    const uint32_t* fb;      /* fb[c] = prime for column c (fb[0] unused) */
    Rel*      rel;
    size_t    nrel, caprel;
    uint32_t* pcol;          /* pool of column indices */
    uint32_t* pexp;          /* matching exponents */
    size_t    npool, cappool;
    Partial*  part;          /* open-addressed hash of partials */
    size_t    partcap, npart;
    uint64_t  ncombined;
    unsigned char* colseen;  /* column ever odd in a full relation */
    size_t    nactive;       /* number of such columns */
};

static uint32_t pool_push(FiRels* R, const uint32_t* col, const uint32_t* exps, size_t cnt) {
    if (R->npool + cnt > R->cappool) {
        while (R->npool + cnt > R->cappool) R->cappool = R->cappool ? 2 * R->cappool : 4096;
        R->pcol = realloc(R->pcol, R->cappool * sizeof(uint32_t));
        R->pexp = realloc(R->pexp, R->cappool * sizeof(uint32_t));
    }
    uint32_t off = (uint32_t)R->npool;
    memcpy(R->pcol + off, col, cnt * sizeof(uint32_t));
    memcpy(R->pexp + off, exps, cnt * sizeof(uint32_t));
    R->npool += cnt;
    return off;
}

FiRels* fi_rels_new(const uint32_t* fb, size_t ncols) {
    FiRels* R = calloc(1, sizeof *R);
    R->fb = fb;
    R->ncols = ncols;
    R->partcap = 1 << 12;
    R->part = calloc(R->partcap, sizeof(Partial));
    R->colseen = calloc(ncols + 1, 1);
    return R;
}

void fi_rels_free(FiRels* R) {
    if (!R) return;
    for (size_t i = 0; i < R->nrel; i++) mpz_clears(R->rel[i].X, R->rel[i].E, NULL);
    for (size_t i = 0; i < R->partcap; i++) if (R->part[i].L) mpz_clear(R->part[i].X);
    free(R->rel); free(R->pcol); free(R->pexp); free(R->part); free(R->colseen);
    free(R);
}

size_t fi_rels_count(const FiRels* R) { return R->nrel; }
size_t fi_rels_active(const FiRels* R) { return R->nactive; }
size_t fi_rels_partials(const FiRels* R) { return R->npart; }

static void add_full(FiRels* R, const mpz_t X, const mpz_t E, const uint32_t* col,
                     const uint32_t* exps, size_t cnt) {
    if (R->nrel == R->caprel) {
        R->caprel = R->caprel ? 2 * R->caprel : 1024;
        R->rel = realloc(R->rel, R->caprel * sizeof(Rel));
    }
    Rel* r = &R->rel[R->nrel++];
    mpz_init_set(r->X, X);
    mpz_init_set(r->E, E);
    r->cnt = (uint32_t)cnt;
    r->off = pool_push(R, col, exps, cnt);
    for (size_t k = 0; k < cnt; k++)
        if ((exps[k] & 1) && !R->colseen[col[k]]) { R->colseen[col[k]] = 1; R->nactive++; }
}

static void part_grow(FiRels* R);

/* Merge two sparse (col, exp) lists into tmp buffers; returns count. */
static size_t merge_lists(const uint32_t* c1, const uint32_t* e1, size_t n1,
                          const uint32_t* c2, const uint32_t* e2, size_t n2,
                          uint32_t* oc, uint32_t* oe) {
    size_t i = 0, j = 0, k = 0;
    while (i < n1 || j < n2) {
        if (j >= n2 || (i < n1 && c1[i] < c2[j])) { oc[k] = c1[i]; oe[k++] = e1[i++]; }
        else if (i >= n1 || c2[j] < c1[i])        { oc[k] = c2[j]; oe[k++] = e2[j++]; }
        else { oc[k] = c1[i]; oe[k++] = e1[i++] + e2[j++]; }
    }
    return k;
}

/* Add a relation X^2 == prod fb[col]^exp * L (mod n).  col must be sorted
 * ascending (column 0 is the sign, -1).  L = 1 for a full relation, else a
 * prime outside the factor base.  Returns 1 if a full relation was added. */
int fi_rels_add(FiRels* R, const mpz_t X, const uint32_t* col, const uint32_t* exps,
                size_t cnt, uint64_t L, const mpz_t n) {
    mpz_t one;
    if (L <= 1) {
        mpz_init_set_ui(one, 1);
        add_full(R, X, one, col, exps, cnt);
        mpz_clear(one);
        return 1;
    }
    size_t mask = R->partcap - 1;
    size_t h = (size_t)((L * UINT64_C(0x9E3779B97F4A7C15)) >> 20) & mask;
    while (R->part[h].L && R->part[h].L != L) h = (h + 1) & mask;
    if (!R->part[h].L) {
        Partial* p = &R->part[h];
        p->L = L;
        mpz_init_set(p->X, X);
        p->cnt = (uint32_t)cnt;
        p->off = pool_push(R, col, exps, cnt);
        if (++R->npart * 2 > R->partcap) part_grow(R);
        return 0;
    }
    /* Combine with the stored partial: X = X1 X2, E = L. */
    Partial* p = &R->part[h];
    uint32_t* oc = malloc((cnt + p->cnt) * sizeof(uint32_t));
    uint32_t* oe = malloc((cnt + p->cnt) * sizeof(uint32_t));
    size_t k = merge_lists(R->pcol + p->off, R->pexp + p->off, p->cnt, col, exps, cnt, oc, oe);
    mpz_t XX, E;
    mpz_inits(XX, E, NULL);
    mpz_mul(XX, X, p->X);
    mpz_mod(XX, XX, n);
    mpz_set_ui(E, (unsigned long)L);
    add_full(R, XX, E, oc, oe, k);
    mpz_clears(XX, E, NULL);
    free(oc); free(oe);
    R->ncombined++;
    return 1;
}

static void part_grow(FiRels* R) {
    size_t oldcap = R->partcap;
    Partial* old = R->part;
    R->partcap *= 2;
    R->part = calloc(R->partcap, sizeof(Partial));
    size_t mask = R->partcap - 1;
    for (size_t i = 0; i < oldcap; i++) {
        if (!old[i].L) continue;
        size_t h = (size_t)((old[i].L * UINT64_C(0x9E3779B97F4A7C15)) >> 20) & mask;
        while (R->part[h].L) h = (h + 1) & mask;
        R->part[h] = old[i];            /* moves the mpz_t struct by value */
    }
    free(old);
}

/* ------------------------------------------------------------------ */
/* Dense Gaussian elimination                                          */
/* ------------------------------------------------------------------ */

/* Columns of the GF(2) system are the relations (nr of them); rows are the
 * factor-base primes (nc).  Returns dependencies as row-bit vectors over
 * relations. */
static size_t dense_deps(uint64_t* rows, size_t nc, size_t nr, uint64_t* deps, size_t max_deps) {
    size_t W = fi_gf2_words(nr);
    uint32_t* pivcol = malloc(nc * sizeof(uint32_t));
    unsigned char* is_piv = calloc(nr, 1);
    size_t rank = 0;
    for (size_t j = 0; j < nr && rank < nc; j++) {
        size_t w = j / 64;
        uint64_t bit = UINT64_C(1) << (j % 64);
        size_t r = rank;
        while (r < nc && !(rows[r * W + w] & bit)) r++;
        if (r == nc) continue;
        if (r != rank)
            for (size_t k = 0; k < W; k++) {
                uint64_t t = rows[r * W + k]; rows[r * W + k] = rows[rank * W + k]; rows[rank * W + k] = t;
            }
        uint64_t* pr = rows + rank * W;
        for (size_t i = 0; i < nc; i++) {
            if (i == rank) continue;
            uint64_t* ri = rows + i * W;
            if (ri[w] & bit) for (size_t k = w; k < W; k++) ri[k] ^= pr[k];
        }
        pivcol[rank++] = (uint32_t)j;
        is_piv[j] = 1;
    }
    size_t nd = 0;
    for (size_t f = 0; f < nr && nd < max_deps; f++) {
        if (is_piv[f]) continue;
        uint64_t* x = deps + nd * W;
        memset(x, 0, W * sizeof(uint64_t));
        x[f / 64] |= UINT64_C(1) << (f % 64);
        for (size_t i = 0; i < rank; i++)
            if ((rows[i * W + f / 64] >> (f % 64)) & 1)
                x[pivcol[i] / 64] |= UINT64_C(1) << (pivcol[i] % 64);
        nd++;
    }
    free(pivcol);
    free(is_piv);
    return nd;
}

/* ------------------------------------------------------------------ */
/* Sparse front end: prune singletons, dispatch to a solver            */
/* ------------------------------------------------------------------ */

size_t fi_gf2_lanczos(const FiSparseMat* m, uint64_t* deps, size_t max_deps);

/* Lanczos needs room for 64-wide blocks; tiny systems go dense. */
#define DENSE_MAX 1000

size_t fi_gf2_dependencies(const FiSparseMat* m, uint64_t* deps, size_t max_deps) {
    size_t nr = m->nrows, nc = m->ncols;
    /* Prune: drop rows containing a column of weight 1, repeatedly. */
    unsigned char* alive = malloc(nr);
    memset(alive, 1, nr);
    uint32_t* w = calloc(nc, sizeof(uint32_t));
    for (size_t r = 0; r < nr; r++)
        for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) w[m->col[k]]++;
    for (int changed = 1; changed; ) {
        changed = 0;
        for (size_t r = 0; r < nr; r++) {
            if (!alive[r]) continue;
            for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) {
                if (w[m->col[k]] == 1) {
                    alive[r] = 0;
                    for (uint32_t q = m->start[r]; q < m->start[r + 1]; q++) w[m->col[q]]--;
                    changed = 1;
                    break;
                }
            }
        }
    }
    /* Compact columns and rows. */
    uint32_t* cmap = malloc(nc * sizeof(uint32_t));
    size_t nc2 = 0;
    for (size_t c = 0; c < nc; c++) cmap[c] = w[c] ? (uint32_t)nc2++ : UINT32_MAX;
    uint32_t* rmap = malloc(nr * sizeof(uint32_t));
    size_t nr2 = 0;
    for (size_t r = 0; r < nr; r++) if (alive[r]) rmap[nr2++] = (uint32_t)r;
    size_t nd = 0;
    if (nr2 > nc2 || (nr2 > 0 && nc2 == 0)) {
        size_t W2 = fi_gf2_words(nr2);
        uint64_t* d2 = calloc(max_deps * W2 + 1, sizeof(uint64_t));
        size_t dense_max = DENSE_MAX;
        FI_TUNE_U("FI_DENSE_MAX", dense_max);
        if (nc2 <= dense_max) {
            uint64_t* rows = calloc(nc2 * W2 + 1, sizeof(uint64_t));
            for (size_t i = 0; i < nr2; i++) {
                size_t r = rmap[i];
                for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) {
                    uint32_t c = cmap[m->col[k]];
                    rows[c * W2 + i / 64] ^= UINT64_C(1) << (i % 64);
                }
            }
            nd = dense_deps(rows, nc2, nr2, d2, max_deps);
            free(rows);
        } else {
            FiSparseMat sm;
            sm.nrows = (uint32_t)nr2;
            sm.ncols = (uint32_t)nc2;
            sm.start = malloc((nr2 + 1) * sizeof(uint32_t));
            size_t tot = 0;
            for (size_t i = 0; i < nr2; i++) tot += m->start[rmap[i] + 1] - m->start[rmap[i]];
            sm.col = malloc((tot + 1) * sizeof(uint32_t));
            size_t q = 0;
            for (size_t i = 0; i < nr2; i++) {
                sm.start[i] = (uint32_t)q;
                size_t r = rmap[i];
                for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) sm.col[q++] = cmap[m->col[k]];
            }
            sm.start[nr2] = (uint32_t)q;
            nd = fi_gf2_lanczos(&sm, d2, max_deps);
            free(sm.start);
            free(sm.col);
            if (nd == 0 && nc2 <= 4 * DENSE_MAX) {
                /* Lanczos can fail on small or degenerate systems: fall back. */
                uint64_t* rows = calloc(nc2 * W2 + 1, sizeof(uint64_t));
                for (size_t i = 0; i < nr2; i++) {
                    size_t r = rmap[i];
                    for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) {
                        uint32_t c = cmap[m->col[k]];
                        rows[c * W2 + i / 64] ^= UINT64_C(1) << (i % 64);
                    }
                }
                nd = dense_deps(rows, nc2, nr2, d2, max_deps);
                free(rows);
            }
        }
        /* Map dependencies back to original row numbering. */
        size_t W = fi_gf2_words(nr);
        for (size_t k = 0; k < nd; k++) {
            uint64_t* out = deps + k * W;
            memset(out, 0, W * sizeof(uint64_t));
            for (size_t i = 0; i < nr2; i++)
                if ((d2[k * W2 + i / 64] >> (i % 64)) & 1)
                    out[rmap[i] / 64] |= UINT64_C(1) << (rmap[i] % 64);
        }
        free(d2);
    }
    free(alive); free(w); free(cmap); free(rmap);
    return nd;
}

/* ------------------------------------------------------------------ */
/* Square root step                                                    */
/* ------------------------------------------------------------------ */

/* Try every dependency until one splits n.  Returns 1 with d set. */
int fi_rels_solve(FiRels* R, const mpz_t n, mpz_t d) {
    size_t nr = R->nrel;
    if (nr == 0) return 0;
    FiSparseMat m;
    m.nrows = (uint32_t)nr;
    m.ncols = (uint32_t)R->ncols;
    m.start = malloc((nr + 1) * sizeof(uint32_t));
    m.col = malloc((R->npool + 1) * sizeof(uint32_t));
    size_t q = 0;
    for (size_t r = 0; r < nr; r++) {
        m.start[r] = (uint32_t)q;
        const Rel* x = &R->rel[r];
        for (uint32_t k = 0; k < x->cnt; k++)
            if (R->pexp[x->off + k] & 1) m.col[q++] = R->pcol[x->off + k];
    }
    m.start[nr] = (uint32_t)q;
    const size_t MAXD = 64;
    size_t W = fi_gf2_words(nr);
    uint64_t* deps = calloc(MAXD * W, sizeof(uint64_t));
    size_t nd = fi_gf2_dependencies(&m, deps, MAXD);
    free(m.start);
    free(m.col);

    int found = 0;
    uint64_t* ex = calloc(R->ncols, sizeof(uint64_t));
    mpz_t X, Y, t;
    mpz_inits(X, Y, t, NULL);
    for (size_t k = 0; k < nd && !found; k++) {
        memset(ex, 0, R->ncols * sizeof(uint64_t));
        mpz_set_ui(X, 1);
        mpz_set_ui(Y, 1);
        for (size_t r = 0; r < nr; r++) {
            if (!((deps[k * W + r / 64] >> (r % 64)) & 1)) continue;
            const Rel* x = &R->rel[r];
            mpz_mul(X, X, x->X); mpz_mod(X, X, n);
            mpz_mul(Y, Y, x->E); mpz_mod(Y, Y, n);
            for (uint32_t j = 0; j < x->cnt; j++) ex[R->pcol[x->off + j]] += R->pexp[x->off + j];
        }
        int bad = 0;
        for (size_t c = 1; c < R->ncols; c++) {
            if (ex[c] & 1) { bad = 1; break; }
            if (ex[c]) {
                mpz_set_ui(t, R->fb[c]);
                mpz_powm_ui(t, t, (unsigned long)(ex[c] / 2), n);
                mpz_mul(Y, Y, t); mpz_mod(Y, Y, n);
            }
        }
        if (bad || (ex[0] & 1)) continue;
        mpz_sub(t, X, Y);
        mpz_gcd(t, t, n);
        if (mpz_cmp_ui(t, 1) > 0 && mpz_cmp(t, n) < 0) { mpz_set(d, t); found = 1; }
    }
    mpz_clears(X, Y, t, NULL);
    free(ex);
    free(deps);
    return found;
}
