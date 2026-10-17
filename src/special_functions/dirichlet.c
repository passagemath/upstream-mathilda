/* dirichlet.c -- the Dirichlet character engine. See dirichlet.h for the
 * labelling convention and the table layout.
 *
 * Discrete logarithms: dense tables are built by walking the powers of each
 * component generator (O(f), no logs at all); a single value for a large
 * modulus uses Pohlig-Hellman over the component order, with baby-step
 * giant-step in each prime-order subgroup.
 */
#include "dirichlet.h"
#include "facint.h"
#include "eval.h"
#include "sym_names.h"
#include "arithmetic.h"   /* make_rational, make_complex, is_rational_like */
#include "message.h"
#include "print.h"
#include <gmp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef MATHILDA_THREADS
#include <pthread.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/*  Modular arithmetic on uint64 (moduli < 2^63)                       */
/* ------------------------------------------------------------------ */

uint64_t dch_gcd(uint64_t a, uint64_t b) {
    while (b) { uint64_t r = a % b; a = b; b = r; }
    return a;
}

uint64_t dch_mulmod(uint64_t a, uint64_t b, uint64_t m) {
    a %= m; b %= m;
    if (a < UINT64_C(4294967296) && b < UINT64_C(4294967296)) return (a * b) % m;
    /* Russian-peasant multiplication: every intermediate stays below 2m < 2^64. */
    uint64_t r = 0;
    while (b) {
        if (b & 1) r = (r >= m - a) ? r - (m - a) : r + a;
        b >>= 1;
        a = (a >= m - a) ? a - (m - a) : a + a;
    }
    return r;
}

uint64_t dch_powmod(uint64_t b, uint64_t e, uint64_t m) {
    if (m == 1) return 0;
    uint64_t r = 1; b %= m;
    while (e) {
        if (e & 1) r = dch_mulmod(r, b, m);
        e >>= 1;
        if (e) b = dch_mulmod(b, b, m);
    }
    return r;
}

/* Inverse of a mod m (gcd(a, m) = 1, m < 2^63), by the extended Euclid
 * recurrence on signed coefficients bounded by m. */
static uint64_t inv_mod(uint64_t a, uint64_t m) {
    int64_t t = 0, nt = 1;
    uint64_t r = m, nr = a % m;
    while (nr) {
        uint64_t q = r / nr;
        int64_t tmp = t - (int64_t)q * nt; t = nt; nt = tmp;
        uint64_t rr = r - q * nr; r = nr; nr = rr;
    }
    if (t < 0) t += (int64_t)m;
    return (uint64_t)t;
}

/* ------------------------------------------------------------------ */
/*  Factorisation of a uint64                                          */
/* ------------------------------------------------------------------ */

/* Ascending prime factorisation of n >= 1 (n < 2^63). Trial division removes
 * every prime below 2^16; a cofactor below 2^32 is then prime, and a larger one
 * goes to the proven factoriser. Returns false when that cannot complete. */
static bool factor_u64(uint64_t n, uint64_t* p, int* e, int cap, int* count) {
    int c = 0;
    if (n % 2 == 0) {
        int k = 0; while (n % 2 == 0) { n /= 2; k++; }
        if (c >= cap) return false;
        p[c] = 2; e[c] = k; c++;
    }
    uint64_t d = 3;
    for (; d < 65536 && d * d <= n; d += 2) {
        if (n % d) continue;
        int k = 0; while (n % d == 0) { n /= d; k++; }
        if (c >= cap) return false;
        p[c] = d; e[c] = k; c++;
    }
    if (n > 1) {
        /* Prime when trial division ran past its square root (every factor
         * below 2^16 is gone, so a cofactor below 2^32 is prime too). */
        if (d * d > n || n < UINT64_C(4294967296)) {
            if (c >= cap) return false;
            p[c] = n; e[c] = 1; c++;
        } else {
            mpz_t N; mpz_init(N);
            mpz_import(N, 1, 1, sizeof(uint64_t), 0, 0, &n);
            mpz_t fp[DCH_MAX_FAC]; int64_t fe[DCH_MAX_FAC]; int fc = 0;
            bool ok = facint_factor_complete(N, fp, fe, DCH_MAX_FAC, &fc);
            mpz_clear(N);
            if (!ok) return false;
            for (int i = 0; i < fc; i++) {
                if (ok && c < cap) {
                    uint64_t v = 0;
                    mpz_export(&v, NULL, 1, sizeof(uint64_t), 0, 0, fp[i]);
                    p[c] = v; e[c] = (int)fe[i]; c++;
                } else ok = false;
                mpz_clear(fp[i]);
            }
            if (!ok) return false;
            /* trial division left only primes > 2^16, all ascending already from
             * facint, and all larger than the trial primes: order is preserved. */
        }
    }
    *count = c;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Discrete logarithm (Pohlig-Hellman + BSGS)                         */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t v, i; } BsgsEnt;

static int bsgs_cmp(const void* a, const void* b) {
    uint64_t x = ((const BsgsEnt*)a)->v, y = ((const BsgsEnt*)b)->v;
    return (x > y) - (x < y);
}

/* log base gam of h in the subgroup of prime order l, modulo q. */
static bool dlog_prime(uint64_t h, uint64_t gam, uint64_t l, uint64_t q, uint64_t* out) {
    if (h == 1) { *out = 0; return true; }
    if (l <= 4096) {
        uint64_t cur = 1;
        for (uint64_t i = 0; i < l; i++) {
            if (cur == h) { *out = i; return true; }
            cur = dch_mulmod(cur, gam, q);
        }
        return false;
    }
    if (l > (UINT64_C(1) << 44)) return false;           /* > 2^22 baby steps */
    uint64_t m = (uint64_t)ceil(sqrt((double)l));
    BsgsEnt* tab = (BsgsEnt*)malloc(m * sizeof(BsgsEnt));
    if (!tab) return false;
    uint64_t cur = 1;
    for (uint64_t i = 0; i < m; i++) { tab[i].v = cur; tab[i].i = i; cur = dch_mulmod(cur, gam, q); }
    qsort(tab, m, sizeof(BsgsEnt), bsgs_cmp);
    uint64_t step = dch_powmod(gam, l - (m % l), q);     /* gam^-m */
    uint64_t y = h;
    bool found = false;
    for (uint64_t g = 0; g <= m && !found; g++) {
        BsgsEnt key; key.v = y; key.i = 0;
        BsgsEnt* hit = (BsgsEnt*)bsearch(&key, tab, m, sizeof(BsgsEnt), bsgs_cmp);
        if (hit) { *out = (g * m + hit->i) % l; found = true; }
        y = dch_mulmod(y, step, q);
    }
    free(tab);
    return found;
}

/* log base g of x modulo q, where g generates a cyclic group of order N with
 * the factorisation in c. */
static bool dlog_comp(uint64_t x, uint64_t g, uint64_t q, const DchComp* c, uint64_t* out) {
    uint64_t N = c->ord;
    if (N == 1) { *out = 0; return true; }
    uint64_t t = 0;
    for (int i = 0; i < c->nfac; i++) {
        uint64_t l = c->fac_p[i];
        int e = c->fac_e[i];
        uint64_t pe = 1; for (int k = 0; k < e; k++) pe *= l;
        uint64_t gg = dch_powmod(g, N / pe, q);            /* order pe */
        uint64_t hh = dch_powmod(x, N / pe, q);
        uint64_t gam = dch_powmod(gg, pe / l, q);          /* order l  */
        uint64_t y = 0, lk = 1;
        for (int k = 0; k < e; k++) {
            /* (hh * gg^-y)^(pe / l^(k+1)) */
            uint64_t ginv = dch_powmod(gg, (pe - (y % pe)) % pe, q);
            uint64_t z = dch_mulmod(hh, ginv, q);
            z = dch_powmod(z, pe / (lk * l), q);
            uint64_t dig;
            if (!dlog_prime(z, gam, l, q, &dig)) return false;
            y += dig * lk;
            lk *= l;
        }
        /* CRT: t = t (mod N/..) combined with y (mod pe). */
        uint64_t M = N / pe;
        uint64_t coef = dch_mulmod(M, inv_mod(M % pe, pe), N);
        t = (t + dch_mulmod(y % pe, coef, N)) % N;
    }
    *out = t;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Table construction                                                 */
/* ------------------------------------------------------------------ */

/* Reduced digit of a component: d/ord = dr/o in lowest terms. */
static void comp_reduced(const DchComp* c, uint64_t* dr, uint64_t* o) {
    if (c->d == 0) { *dr = 0; *o = 1; return; }
    uint64_t g = dch_gcd(c->d, c->ord);
    *dr = c->d / g; *o = c->ord / g;
}

/* Merge the factorisation of p-1 with p^(a-1) into c->fac. */
static bool comp_set_fac_cyc(DchComp* c, uint64_t p, int a) {
    uint64_t fp[DCH_MAX_FAC]; int fe[DCH_MAX_FAC]; int n = 0;
    if (p > 2 && !factor_u64(p - 1, fp, fe, DCH_MAX_FAC, &n)) return false;
    c->nfac = 0;
    bool placed = (a <= 1);
    for (int i = 0; i < n; i++) {
        if (!placed && fp[i] > p) {
            c->fac_p[c->nfac] = p; c->fac_e[c->nfac] = a - 1; c->nfac++; placed = true;
        }
        if (c->nfac >= DCH_MAX_FAC) return false;
        c->fac_p[c->nfac] = fp[i]; c->fac_e[c->nfac] = fe[i]; c->nfac++;
    }
    if (!placed) {
        if (c->nfac >= DCH_MAX_FAC) return false;
        c->fac_p[c->nfac] = p; c->fac_e[c->nfac] = a - 1; c->nfac++;
    }
    return true;
}

/* Smallest primitive root mod q = p^a (p odd). */
static uint64_t smallest_primroot(uint64_t p, uint64_t q, const DchComp* c) {
    for (uint64_t g = 2; g < q; g++) {
        if (g % p == 0) continue;
        bool ok = true;
        for (int i = 0; i < c->nfac && ok; i++)
            if (dch_powmod(g, c->ord / c->fac_p[i], q) == 1) ok = false;
        if (ok) return g;
    }
    return 0;
}

/* Log of x (a unit mod part->q) in component c of that part. */
static bool comp_log(const DchPart* P, const DchComp* c, uint64_t x, uint64_t* t) {
    switch (c->kind) {
        case DCH_SIGN: *t = (x % 4 == 3) ? 1 : 0; return true;
        case DCH_FIVE: {
            uint64_t m = x % P->q;
            if (m % 4 == 3) m = P->q - m;
            return dlog_comp(m, 5, P->q, c, t);
        }
        default:
            return dlog_comp(x % P->q, c->g, P->q, c, t);
    }
}

/* Contribution of component c at log t, as a numerator over `order`. */
static uint64_t comp_contrib(const DchComp* c, uint64_t t, uint64_t order) {
    uint64_t dr, o;
    comp_reduced(c, &dr, &o);
    if (o == 1) return 0;
    uint64_t v = dch_mulmod(dr, t % o, o);
    return dch_mulmod(v, order / o, order);
}

bool dch_star_exp(const DchTable* T, uint64_t n, uint64_t* e) {
    uint64_t acc = 0;
    for (int i = 0; i < T->nparts; i++) {
        const DchPart* P = &T->part[i];
        if (P->b == 0) continue;
        if (n % P->p == 0) return false;
        for (int k = 0; k < P->ncomp; k++) {
            const DchComp* c = &P->comp[k];
            if (c->d == 0) continue;
            uint64_t t;
            if (!comp_log(P, c, n, &t)) return false;
            acc = (acc + comp_contrib(c, t, T->order)) % T->order;
        }
    }
    *e = acc;
    return true;
}

bool dch_exp(const DchTable* T, uint64_t n, uint64_t* e) {
    if (T->k > 1 && dch_gcd(n % (uint64_t)T->k, (uint64_t)T->k) != 1) return false;
    if (T->principal) { *e = 0; return true; }
    return dch_star_exp(T, n, e);
}

/* Dense chi* tables over 0..f-1 by walking the component generators. */
static bool build_tables(DchTable* T) {
    uint64_t f = T->f;
    T->estar = (int32_t*)malloc(f * sizeof(int32_t));
    T->cre = (double*)malloc(f * sizeof(double));
    T->cim = (double*)malloc(f * sizeof(double));
    if (!T->estar || !T->cre || !T->cim) return false;
    for (uint64_t a = 0; a < f; a++) T->estar[a] = 0;
    if (f == 1) T->estar[0] = 0;
    uint64_t O = T->order;
    for (int i = 0; i < T->nparts; i++) {
        const DchPart* P = &T->part[i];
        if (P->b == 0) continue;
        uint64_t qb = P->qb;
        int32_t* E = (int32_t*)malloc(qb * sizeof(int32_t));
        if (!E) return false;
        for (uint64_t x = 0; x < qb; x++) E[x] = -1;
        if (P->p != 2) {
            const DchComp* c = &P->comp[0];
            uint64_t phib = qb / P->p * (P->p - 1);
            uint64_t g = c->g % qb, x = 1;
            for (uint64_t t = 0; t < phib; t++) {
                E[x] = (int32_t)comp_contrib(c, t, O);
                x = dch_mulmod(x, g, qb);
            }
        } else {
            const DchComp* cs = &P->comp[0];
            const DchComp* c5 = (P->ncomp > 1) ? &P->comp[1] : NULL;
            uint64_t n5 = (P->b >= 3) ? (qb >> 2) : 1;      /* |<5> mod 2^b| */
            for (int sg = 0; sg < 2; sg++) {
                uint64_t x = 1;
                for (uint64_t t = 0; t < n5; t++) {
                    uint64_t y = sg ? (qb - x) % qb : x;
                    uint64_t v = comp_contrib(cs, (uint64_t)sg, O);
                    if (c5) v = (v + comp_contrib(c5, t, O)) % O;
                    E[y] = (int32_t)v;
                    x = dch_mulmod(x, 5, qb);
                }
            }
        }
        for (uint64_t a = 0; a < f; a++) {
            if (T->estar[a] < 0) continue;
            int32_t v = E[a % qb];
            T->estar[a] = (v < 0) ? -1 : (int32_t)(((uint64_t)T->estar[a] + (uint64_t)v) % O);
        }
        free(E);
    }
    if (f > 1) T->estar[0] = -1;
    /* Machine values: exact on the axes, cos/sin of the reduced angle elsewhere. */
    for (uint64_t a = 0; a < f; a++) {
        int32_t e = T->estar[a];
        if (e < 0) { T->cre[a] = 0.0; T->cim[a] = 0.0; continue; }
        uint64_t e4 = (uint64_t)e * 4;
        if (e4 % O == 0) {
            static const double ur[4] = { 1.0, 0.0, -1.0, 0.0 };
            static const double ui[4] = { 0.0, 1.0, 0.0, -1.0 };
            int qd = (int)(e4 / O);
            T->cre[a] = ur[qd]; T->cim[a] = ui[qd];
        } else {
            double num = (double)e;
            if (2 * (uint64_t)e > O) num -= (double)O;     /* angle in (-pi, pi] */
            double th = 2.0 * M_PI * num / (double)O;
            T->cre[a] = cos(th); T->cim[a] = sin(th);
        }
    }
    /* Root number eps = tau(chi*)/(I^odd Sqrt[f]). Exactly 1 for a real
     * primitive character (Gauss); otherwise summed with exact phases. */
    if (T->real || f == 1) {
        T->eps_re = 1.0; T->eps_im = 0.0;
    } else {
        uint64_t L = O / dch_gcd(O, f) * f;              /* lcm(O, f) <= 2^44 */
        double sr = 0.0, si = 0.0;
        for (uint64_t a = 1; a < f; a++) {
            int32_t e = T->estar[a];
            if (e < 0) continue;
            uint64_t num = ((uint64_t)e * (L / O) + a * (L / f)) % L;
            double fr = (double)num;
            if (2 * num > L) fr -= (double)L;
            double th = 2.0 * M_PI * fr / (double)L;
            sr += cos(th); si += sin(th);
        }
        /* divide by I^odd */
        if (T->odd) { double t = sr; sr = si; si = -t; }
        double m = hypot(sr, si);
        if (!(m > 0.0)) return false;
        T->eps_re = sr / m; T->eps_im = si / m;
    }
    return true;
}

void dch_free(DchTable* T) {
    if (!T) return;
    free(T->estar); free(T->cre); free(T->cim);
    free(T);
}

/* Build the table for (k, j); *phi_out = EulerPhi[k] when k is handled. */
static DchTable* dch_build(int64_t k, int64_t j, bool tables, uint64_t* phi_out, bool* range_err) {
    *range_err = false;
    if (k < 1 || j < 1) return NULL;
    DchTable* T = (DchTable*)calloc(1, sizeof(DchTable));
    if (!T) return NULL;
    T->k = k; T->j = j;
    uint64_t pr[DCH_MAX_PRIMES]; int ex[DCH_MAX_PRIMES]; int np = 0;
    if (!factor_u64((uint64_t)k, pr, ex, DCH_MAX_PRIMES, &np)) { free(T); return NULL; }
    T->nparts = np;
    uint64_t phi = 1;
    for (int i = 0; i < np; i++) {
        DchPart* P = &T->part[i];
        P->p = pr[i]; P->a = ex[i];
        P->q = 1; for (int k2 = 0; k2 < P->a; k2++) P->q *= P->p;
        if (P->p == 2) {
            if (P->a >= 2) {
                DchComp* s = &P->comp[P->ncomp++];
                s->kind = DCH_SIGN; s->g = P->q - 1; s->ord = 2;
                s->nfac = 1; s->fac_p[0] = 2; s->fac_e[0] = 1;
            }
            if (P->a >= 3) {
                DchComp* f5 = &P->comp[P->ncomp++];
                f5->kind = DCH_FIVE; f5->g = 5; f5->ord = P->q >> 2;
                f5->nfac = 1; f5->fac_p[0] = 2; f5->fac_e[0] = P->a - 2;
            }
        } else {
            DchComp* c = &P->comp[P->ncomp++];
            c->kind = DCH_CYC;
            c->ord = P->q / P->p * (P->p - 1);
            if (!comp_set_fac_cyc(c, P->p, P->a)) { free(T); return NULL; }
            c->g = smallest_primroot(P->p, P->q, c);
            if (c->g == 0) { free(T); return NULL; }
        }
        for (int c = 0; c < P->ncomp; c++) phi *= P->comp[c].ord;
    }
    if (phi_out) *phi_out = phi;
    if ((uint64_t)j > phi) { *range_err = true; free(T); return NULL; }

    /* Mixed-radix digits: the last component is least significant. */
    uint64_t r = (uint64_t)(j - 1);
    for (int i = np - 1; i >= 0; i--)
        for (int c = T->part[i].ncomp - 1; c >= 0; c--) {
            DchComp* C = &T->part[i].comp[c];
            C->d = r % C->ord; r /= C->ord;
        }

    /* Order, parity, conductor. */
    uint64_t O = 1;
    int odd = 0;
    T->f = 1;
    for (int i = 0; i < np; i++) {
        DchPart* P = &T->part[i];
        for (int c = 0; c < P->ncomp; c++) {
            uint64_t dr, o;
            comp_reduced(&P->comp[c], &dr, &o);
            O = O / dch_gcd(O, o) * o;
            if (P->comp[c].kind != DCH_FIVE) odd ^= (int)(P->comp[c].d & 1);
        }
        P->b = 0;
        if (P->p != 2) {
            uint64_t dr, o;
            comp_reduced(&P->comp[0], &dr, &o);
            if (o > 1) { int v = 0; while (o % P->p == 0) { o /= P->p; v++; } P->b = v + 1; }
        } else if (P->ncomp > 0) {
            uint64_t ds = P->comp[0].d;
            uint64_t o5 = 1, dr;
            if (P->ncomp > 1) comp_reduced(&P->comp[1], &dr, &o5);
            if (o5 > 1) { int c2 = 0; while (o5 > 1) { o5 >>= 1; c2++; } P->b = c2 + 2; }
            else if (ds) P->b = 2;
        }
        P->qb = 1; for (int k2 = 0; k2 < P->b; k2++) P->qb *= P->p;
        T->f *= P->qb;
    }
    T->order = O;
    T->principal = (O == 1);
    T->real = (O <= 2);
    T->odd = odd;

    /* Euler factors. */
    for (int i = 0; i < np; i++) {
        const DchPart* P = &T->part[i];
        if (P->b != 0) continue;
        uint64_t e = 0;
        if (!dch_star_exp(T, P->p, &e)) { dch_free(T); return NULL; }
        T->eul_p[T->neul] = P->p; T->eul_e[T->neul] = e; T->neul++;
    }

    if (tables && T->f <= DCH_TABLE_MAX && !build_tables(T)) { dch_free(T); return NULL; }
    return T;
}

/* ------------------------------------------------------------------ */
/*  Cache                                                              */
/* ------------------------------------------------------------------ */

/* Append-only, immutable-once-published: a reader never sees a table being
 * freed, so the threaded fused-map loops may hold the pointer without the lock.
 * Only tables with small conductor are cached (bounded memory). */
#define DCH_CACHE_SLOTS 64
#define DCH_CACHE_FMAX  ((uint64_t)1 << 16)
static DchTable* g_cache[DCH_CACHE_SLOTS];
static size_t g_cache_n = 0;
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_cache_lock = PTHREAD_MUTEX_INITIALIZER;
#define DCH_LOCK()   pthread_mutex_lock(&g_cache_lock)
#define DCH_UNLOCK() pthread_mutex_unlock(&g_cache_lock)
#else
#define DCH_LOCK()   ((void)0)
#define DCH_UNLOCK() ((void)0)
#endif

/* A cached entry serves a request when it carries the dense tables, or the
 * request does not need them, or they can never exist (f too large). */
static bool dch_serves(const DchTable* T, bool tables) {
    return !tables || T->estar || T->f > DCH_TABLE_MAX;
}

const DchTable* dch_get(int64_t k, int64_t j, bool tables, DchTable** owned, uint64_t* phi_out) {
    *owned = NULL;
    if (phi_out) *phi_out = 0;
    if (k < 1 || j < 1) return NULL;
    DchTable* hit = NULL;
    DCH_LOCK();
    for (size_t i = 0; i < g_cache_n; i++)
        if (g_cache[i]->k == k && g_cache[i]->j == j && dch_serves(g_cache[i], tables)) {
            hit = g_cache[i]; break;
        }
    DCH_UNLOCK();
    if (hit) {
        if (phi_out) {
            uint64_t phi = 1;
            for (int i = 0; i < hit->nparts; i++)
                for (int c = 0; c < hit->part[i].ncomp; c++) phi *= hit->part[i].comp[c].ord;
            *phi_out = phi;
        }
        return hit;
    }
    bool range_err;
    DchTable* T = dch_build(k, j, tables, phi_out, &range_err);
    if (!T) return NULL;
    if (T->f <= DCH_CACHE_FMAX) {
        DchTable* pub = NULL;
        DCH_LOCK();
        for (size_t i = 0; i < g_cache_n; i++)          /* lost a race? */
            if (g_cache[i]->k == k && g_cache[i]->j == j && dch_serves(g_cache[i], tables)) {
                pub = g_cache[i]; break;
            }
        if (!pub && g_cache_n < DCH_CACHE_SLOTS) { g_cache[g_cache_n++] = T; pub = T; T = NULL; }
        DCH_UNLOCK();
        if (pub) { dch_free(T); return pub; }
    }
    *owned = T;
    return T;
}

/* ------------------------------------------------------------------ */
/*  Exact values                                                       */
/* ------------------------------------------------------------------ */

Expr* dch_root_of_unity_expr(int quarter, uint64_t num, uint64_t den) {
    if (den == 0) den = 1;
    num %= den;
    uint64_t g = dch_gcd(num, den);
    if (g == 0) g = den;
    if (num) { num /= g; den /= g; }
    if (num && (4 * num) % den == 0) { quarter += (int)(4 * num / den); num = 0; }
    quarter = ((quarter % 4) + 4) % 4;
    Expr* unit;
    switch (quarter) {
        case 0:  unit = expr_new_integer(1); break;
        case 1:  unit = make_complex(expr_new_integer(0), expr_new_integer(1)); break;
        case 2:  unit = expr_new_integer(-1); break;
        default: unit = make_complex(expr_new_integer(0), expr_new_integer(-1)); break;
    }
    if (num == 0) return unit;
    /* exp(2 Pi I num/den) = E^(I Pi r), r = 2 num/den folded into (-1, 1). */
    int64_t rn = (int64_t)(2 * num), rd = (int64_t)den;
    if (2 * num > den) rn -= 2 * rd;
    Expr* ex = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ make_complex(expr_new_integer(0), make_rational(rn, rd)),
                              expr_new_symbol(SYM_Pi) }, 2);
    Expr* pw = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_new_symbol(SYM_E), ex }, 2);
    if (quarter == 0) { expr_free(unit); return eval_and_free(pw); }
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ unit, pw }, 2));
}

Expr* dch_value_expr_mpz_mod(const DchTable* T, uint64_t n) {
    uint64_t k = (uint64_t)T->k;
    if (k > 1) n %= k;
    if (k > 1 && dch_gcd(n, k) != 1) return expr_new_integer(0);
    int quarter = 0;
    uint64_t R = 0, O = T->order;
    for (int i = 0; i < T->nparts; i++) {
        const DchPart* P = &T->part[i];
        for (int c = 0; c < P->ncomp; c++) {
            const DchComp* C = &P->comp[c];
            if (C->d == 0) continue;
            uint64_t t, dr, o;
            if (!comp_log(P, C, n, &t)) return NULL;
            comp_reduced(C, &dr, &o);
            uint64_t v = dch_mulmod(dr, t % o, o);
            if ((4 * v) % o == 0) quarter += (int)(4 * v / o);
            else R = (R + dch_mulmod(v, O / o, O)) % O;
        }
    }
    return dch_root_of_unity_expr(quarter, R, O);
}

/* ------------------------------------------------------------------ */
/*  Argument validation shared by DirichletCharacter and DirichletL    */
/* ------------------------------------------------------------------ */

/* A number that is not a positive machine-sized integer earns `intp`; a
 * symbolic argument earns nothing (the call just stays unevaluated). */
static bool dch_is_number(Expr* e) {
    Expr *re, *im;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL
        || dch_is_mpfr(e)) return true;
    if (is_rational_like(e)) return true;
    if (is_complex(e, &re, &im)) return true;
    return false;
}

static void dch_msg_pos(const char* head, const char* tag, const char* what,
                        int pos, Expr* res) {
    char* s = expr_to_string(res);
    mth_message(head, tag, "%s expected at position %d in %s.", what, pos, s ? s : "?");
    free(s);
}

const DchTable* dch_validate_kj(const char* head, Expr* res, bool tables, DchTable** owned) {
    *owned = NULL;
    Expr* K = res->data.function.args[0];
    Expr* J = res->data.function.args[1];
    if (!(K->type == EXPR_INTEGER && K->data.integer >= 1)) {
        if (dch_is_number(K) && !(K->type == EXPR_BIGINT && mpz_sgn(K->data.bigint) > 0))
            dch_msg_pos(head, "intp", "Positive integer", 1, res);
        return NULL;   /* symbolic, or a modulus beyond 2^63: unevaluated */
    }
    bool j_big = (J->type == EXPR_BIGINT && mpz_sgn(J->data.bigint) > 0);
    if (!(J->type == EXPR_INTEGER && J->data.integer >= 1) && !j_big) {
        if (dch_is_number(J)) dch_msg_pos(head, "intp", "Positive integer", 2, res);
        return NULL;
    }
    uint64_t phi = 0;
    const DchTable* T = j_big ? NULL : dch_get(K->data.integer, J->data.integer, tables, owned, &phi);
    if (T) return T;
    if (j_big) {                              /* j >= 2^63 > EulerPhi[k] */
        DchTable* tmp = NULL;
        if (dch_get(K->data.integer, 1, false, &tmp, &phi)) dch_free(tmp);
    }
    if (phi > 0) {                            /* k fine, j out of range */
        char* js = expr_to_string(J);
        char* s = expr_to_string(res);
        mth_message(head, "invl",
                    "Argument %s at position 2 in %s should be a positive integer less than or equal to %llu.",
                    js ? js : "?", s ? s : "?", (unsigned long long)phi);
        free(js); free(s);
    }
    return NULL;
}

bool dch_labels_real(const Expr* k, const Expr* j) {
    if (k->type != EXPR_INTEGER || j->type != EXPR_INTEGER) return false;
    DchTable* owned = NULL;
    const DchTable* T = dch_get(k->data.integer, j->data.integer, false, &owned, NULL);
    bool r = T && T->real;
    dch_free(owned);
    return r;
}
