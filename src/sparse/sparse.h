#ifndef MATHILDA_SPARSE_H
#define MATHILDA_SPARSE_H

/*
 * SparseArray -- native compressed sparse storage (EXPR_SPARSE).
 *
 * A SparseArray represents a full rank-r array of which only the entries that
 * differ from one implicit BACKGROUND value are stored. It is a raw object: an
 * atom to AtomQ and to the pattern matcher (only `_SparseArray` and a literal
 * SameQ array match it), self-evaluating, and never re-evaluated -- after
 * `s = SparseArray[{1 -> x}]; x = 7`, `s` still holds `x`.
 *
 * STORAGE (compressed sparse row, generalised to any rank):
 *
 *   rank >= 2   one "row" per index of the first dimension; each stored entry
 *               carries the remaining rank-1 coordinates as a tuple.
 *   rank == 1   a single pseudo-row; each entry carries its one coordinate.
 *
 *   rowptr[r] .. rowptr[r+1]-1   the entries of row r (0-based, nrows+1 slots)
 *   colidx[k*width + j]          0-based coordinate j of entry k's tuple
 *   vals[k]                      entry k's value
 *
 * Within a row the tuples are strictly increasing in row-major order, so the
 * representation of a given array is UNIQUE. expr_eq / expr_hash rely on that,
 * and so does the value-kind rule below.
 *
 * VALUE KINDS. Values are kept in one typed buffer whenever every value has the
 * same machine type, chosen deterministically from the values themselves:
 *
 *   SPV_INT64  every value an EXPR_INTEGER                int64_t[nnz]
 *   SPV_F64    every value an EXPR_REAL                   double[nnz]
 *   SPV_C128   every value Complex[Real, Real]            double[2*nnz] (re,im)
 *   SPV_EXPR   anything else (exact, symbolic, mixed)     Expr*[nnz], owned
 *
 * sp_new() normalises a boxed buffer to a typed one when it can, so the kind is
 * a function of the value list and two SameQ arrays always share a kind. The
 * background is an ordinary owned Expr, independent of the kind: an F64 array
 * may have the exact background 0 (SparseArray[{1 -> 1.}, 3]) or 0.
 *
 * SHARING. The index structure lives in a separate, immutable, reference-counted
 * SparsePattern. An operation that only transforms values -- 2 s, s + 1, Sin[s],
 * N[s] -- reuses its operand's pattern and allocates only the value buffer, and
 * a binary operation whose operands share one pattern is a pure dense-vector
 * operation. The SparseArray payload itself is also immutable once built and
 * reference-counted, so expr_copy and expr_unshare share it.
 *
 * INDEX WIDTH. colidx is int32 when every coordinate it stores fits (every
 * inner dimension < 2^31), which halves index bandwidth; otherwise int64. Use
 * sp_col() to read it, or the SP_COLIDX_DISPATCH macro in a hot loop.
 *
 * SEMANTICS (one invariant): for every element-wise operation f,
 *     Normal[f[A, B, ...]] === f[Normal[A], Normal[B], ...]
 * The background of the result is f of the operands' backgrounds; entries are
 * the union of the operands' explicit positions. Times with two zero
 * backgrounds is the one exception: its explicit set is the intersection,
 * since the product is the background there anyway. Arithmetic does not
 * re-minimise: an explicit entry that cancels to the background stays explicit
 * until SparseArray[s] minimises it.
 */

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

/* Largest rank a SparseArray may have. */
#define SP_MAX_RANK 32

typedef enum {
    SPV_INT64 = 0,
    SPV_F64,
    SPV_C128,
    SPV_EXPR
} SpValKind;

typedef struct SparsePattern {
    unsigned refcount;
    int64_t  nrows;      /* dims[0] for rank >= 2, 1 for rank 1 */
    int64_t  nnz;
    int      width;      /* coordinates per tuple: rank - 1, or 1 for rank 1 */
    bool     idx32;      /* colidx is int32_t (else int64_t) */
    int64_t* rowptr;     /* nrows + 1 */
    void*    colidx;     /* nnz * width, 0-based */
} SparsePattern;

typedef struct SparseArray {
    unsigned       refcount;
    int            rank;
    int64_t*       dims;       /* rank entries, owned */
    Expr*          bg;         /* background value, owned */
    SparsePattern* pat;        /* one reference held */
    SpValKind      kind;
    void*          vals;       /* see SpValKind; owned */
} SparseArray;

/* ---- lifecycle ----------------------------------------------------------- */

/* A pattern with uninitialised rowptr/colidx storage for `nnz` tuples of
 * `width` coordinates over `nrows` rows (refcount 1). idx32 chooses the colidx
 * element width. NULL on allocation failure. */
SparsePattern* sp_pattern_new(int64_t nrows, int64_t nnz, int width, bool idx32);
SparsePattern* sp_pattern_retain(SparsePattern* p);
void           sp_pattern_release(SparsePattern* p);
/* True when `a` and `b` describe the same index structure (pointer-equal, or
 * equal contents). */
bool           sp_pattern_same(const SparsePattern* a, const SparsePattern* b);

/* Should the colidx of an array with these dims use 32-bit indices? */
bool sp_dims_idx32(int rank, const int64_t* dims);

/* Assemble an array (refcount 1). Takes ownership of `bg`, of one reference to
 * `pat` and of `vals` (whose layout `kind` describes; for SPV_EXPR the Expr*
 * elements are owned too). `dims` is copied. A boxed buffer whose values all
 * share one machine type is converted to that typed buffer here, so the kind
 * is canonical. NULL on allocation failure, in which case everything passed
 * in has been released. */
SparseArray* sp_new(int rank, const int64_t* dims, Expr* bg, SparsePattern* pat,
                    SpValKind kind, void* vals);
SparseArray* sp_retain(SparseArray* sa);
void         sp_release(SparseArray* sa);

/* Wrap a payload in an EXPR_SPARSE node, adopting the one reference. */
Expr* sp_to_expr(SparseArray* sa);

/* The payload of an EXPR_SPARSE node (NULL for anything else). Borrowed. */
static inline const SparseArray* sp_of(const Expr* e) {
    return (e && e->type == EXPR_SPARSE) ? e->data.sparse : NULL;
}
static inline bool is_sparse(const Expr* e) {
    return e && e->type == EXPR_SPARSE;
}

/* ---- element access ------------------------------------------------------ */

static inline int64_t sp_col(const SparsePattern* p, int64_t k, int j) {
    return p->idx32 ? (int64_t)((const int32_t*)p->colidx)[k * p->width + j]
                    : ((const int64_t*)p->colidx)[k * p->width + j];
}

static inline int64_t sp_nnz(const SparseArray* sa) { return sa->pat->nnz; }

/* Bytes per value for a typed kind (0 for SPV_EXPR). */
size_t sp_kind_size(SpValKind kind);

/* The value of explicit entry k as a fresh Expr (caller frees). */
Expr* sp_value_expr(const SparseArray* sa, int64_t k);

/* The element at 0-based coordinates `c` (rank entries), as a fresh Expr --
 * the stored value, or a copy of the background. Binary search in the row. */
Expr* sp_get(const SparseArray* sa, const int64_t* c);

/* Total element count prod(dims), or -1 when it overflows int64. */
int64_t sp_total(int rank, const int64_t* dims);

/* ---- identity ------------------------------------------------------------ */

bool     sp_eq(const SparseArray* a, const SparseArray* b);
uint64_t sp_hash(const SparseArray* sa);
int      sp_compare(const SparseArray* a, const SparseArray* b);

/* ---- conversion ---------------------------------------------------------- */

/* SparseArray[Automatic, dims, bg, {1, {rowptr, colidx}, vals}] as an ordinary
 * (unevaluated) EXPR_FUNCTION, with 1-based column indices. This is the form
 * InputForm and FullForm print, and the constructor reads it back. */
Expr* sp_canonical_expr(const SparseArray* sa);

/* The dense array Normal[s] denotes. A packed List when the values and the
 * background share one machine type, nested boxed Lists otherwise. NULL when
 * the dense array would be too large to allocate. */
Expr* sp_normal(const SparseArray* sa);

/* Does `v` equal the background `bg` for the purpose of leaving it implicit?
 * Numerically equal machine/exact numbers count (0. against 0), otherwise
 * SameQ. */
bool sp_value_is_background(const Expr* v, const Expr* bg);

/* ---- construction from COO ----------------------------------------------- */

/*
 * A coordinate-list builder. Entries are appended in specification order;
 * sp_coo_finish sorts them into canonical order (stably, so when a position
 * repeats the FIRST entry wins), drops entries equal to the background, and
 * assembles the array.
 *
 * Values are either boxed (SPV_EXPR: owned Expr*) or one typed machine buffer,
 * fixed at sp_coo_init.
 */
typedef struct {
    int       rank;
    int64_t   dims[SP_MAX_RANK];
    SpValKind kind;
    int64_t   n, cap;
    int64_t*  coords;     /* n * rank, 0-based */
    void*     vals;       /* n values of `kind` */
    bool      ordered;    /* caller guarantees canonical order, no repeats */
    /* Alternatively (a producer that reads a position buffer directly), the
     * entries' composite keys (row << kbits) | inner, as sp_coo_composite
     * defines them, with coords left NULL. Owned. */
    uint64_t* ck;
} SpCoo;

void sp_coo_init(SpCoo* c, int rank, const int64_t* dims, SpValKind kind);
/* Reserve room for `extra` more entries. False on allocation failure. */
bool sp_coo_reserve(SpCoo* c, int64_t extra);
/* Append one entry (0-based coordinates, already range-checked). The boxed
 * form takes ownership of `v`. */
bool sp_coo_push_expr(SpCoo* c, const int64_t* coords, Expr* v);
bool sp_coo_push_f64(SpCoo* c, const int64_t* coords, double v);
bool sp_coo_push_i64(SpCoo* c, const int64_t* coords, int64_t v);
bool sp_coo_push_c128(SpCoo* c, const int64_t* coords, double re, double im);
void sp_coo_free(SpCoo* c);
/* The composite key layout for an array of these dims: entry key
 * (row << *kbits) | inner, with inner the row-major index of the inner
 * coordinates. False when it does not fit 63 bits (then coords must be used). */
bool sp_coo_composite(int rank, const int64_t* dims, int* kbits);

/* Build the array and release the builder. Takes ownership of `bg`. When
 * `drop_bg` is set, entries equal to the background (sp_value_is_background)
 * are not stored. NULL on allocation failure. */
SparseArray* sp_coo_finish(SpCoo* c, Expr* bg, bool drop_bg);

/* ---- operations (sparse_arith.c / sparse_map.c / sparse_part.c) ---------- */

/* Thread a Listable call with at least one EXPR_SPARSE argument and no List
 * argument. Returns the result (an EXPR_SPARSE, or NULL to leave the call
 * alone, e.g. on a dimension mismatch, after the message). Borrows `call`. */
Expr* sp_thread_listable(const Expr* call);

/* Part[s, i1, ..., ik]. NULL with *err set when a part does not exist (the
 * caller reports Part::partw); NULL with *err clear when the spec is not one
 * this path handles. Borrows everything. */
Expr* sp_part(const Expr* s, Expr* const* idx, size_t n, bool* err);

/* s with Part[s, idx...] = rhs applied, or NULL when the assignment is not one
 * this path handles (*err set when a part does not exist). Borrows
 * everything. */
Expr* sp_part_assign(const Expr* s, Expr* const* idx, size_t n, const Expr* rhs, bool* err);

/* A copy of `s` with every explicit value and the background mapped through
 * `fn` (which borrows its argument and returns a fresh Expr, or NULL to give
 * up). The pattern is shared. When `int_to_real` is set an SPV_INT64 value
 * buffer is converted to SPV_F64 directly and SPV_F64 / SPV_C128 buffers are
 * reused untouched -- the machine-precision N[] fast path; only the
 * background then goes through `fn`. NULL when `fn` gave up. Borrows `s`. */
typedef Expr* (*SpValueFn)(const Expr* v, void* ctx);
Expr* sp_map_values(const Expr* s, SpValueFn fn, void* ctx, bool int_to_real);

/* Does `test` give True on every element of the array (the background once,
 * standing for every implicit element)? Borrows both. */
bool sp_elements_pass(const Expr* s, const Expr* test);

/* s["prop"]: a property of the array (sparse_builtins.c). NULL, after
 * SparseArray::nomthd for an unknown property name, leaves s[...] alone.
 * Borrows everything. */
Expr* sp_apply(const Expr* s, Expr* const* args, size_t n);

/* ---- transparency gate (sparse_gate.c) ----------------------------------- */

/* SymbolDef.sparse_mode: how a head treats a SparseArray argument. 0 (the
 * default) leaves it an atom. See sparse_gate.c for the modes. */
enum {
    SPG_NONE     = 0,
    SPG_DENSE    = 1,     /* answer as on Normal[s] */
    SPG_RESPARSE = 2,     /* ... and an array answer comes back sparse */
    SPG_ROWS     = 3,     /* replace the array by the List of its level-1 parts */
    SPG_MAP      = 4,     /* SPG_ROWS for rank >= 2; SPG_RESPARSE for a vector */
    SPG_MODE_MASK = 0x0f,
    SPG_NATIVE   = 0x10,  /* the builtin handles sparse itself; gate on decline */
    SPG_DEEP     = 0x20,  /* also convert arrays inside List arguments (2 levels) */
    SPG_ANYPOS   = 0x40,  /* SPG_RESPARSE: any array argument makes the answer sparse */
    SPG_ALLARR   = 0x80   /* SPG_RESPARSE: ... only if EVERY array argument is sparse */
};

/* Apply the gate to `call` (borrowed), whose head has sparse_mode `mode_bits`.
 * Returns the answer, or NULL: with *rest set the call must be left as it is
 * (the head declined, or densifying was refused), otherwise the gate had
 * nothing to do. */
Expr* sp_gate(const Expr* call, int mode_bits, bool* rest);

/* A native SparseArray kernel (SymbolDef.sparse_native): borrows the whole
 * call, which has at least one SparseArray argument, and returns the answer,
 * or NULL to decline (the head's sparse_mode gate, if any, then answers).
 * Register with symtab_set_sparse_native. */
typedef Expr* (*SpNativeFn)(const Expr* call);

/* Register the gate modes of the heads that opt in (called by sparse_init). */
void sparse_gate_init(void);

/* Has any SparseArray been built in this session? The evaluator's argument
 * scan for the gate is skipped until one has. */
extern bool sp_g_any_created;
static inline bool sp_any_created(void) { return sp_g_any_created; }

/* ---- module -------------------------------------------------------------- */

void sparse_init(void);

#endif /* MATHILDA_SPARSE_H */
