#ifndef MATHILDA_SPARSE_INTERNAL_H
#define MATHILDA_SPARSE_INTERNAL_H

/* Shared internals of the SparseArray module (not part of its public API). */

#include "sparse.h"
#include "ndarray.h"

/* A rank-1 NDArray node that BORROWS `n` values of `kind` from `vals` -- the
 * value buffer of a SparseArray seen as a dense vector, so the element-wise ND
 * kernels run over it unchanged. Release with sp_view_free (never expr_free:
 * the buffer is not the view's). NULL for SPV_EXPR. */
Expr* sp_view(SpValKind kind, void* vals, int64_t n);
void  sp_view_free(Expr* view);

/* Take the buffer out of an NDArray result (rank 1) as a value buffer of the
 * matching kind, freeing the node. False (node freed) for a dtype that has no
 * kind (float32, bool). */
bool sp_adopt_ndarray(Expr* nd, SpValKind* kind, void** vals);

/* The union of several arrays' explicit positions (equal dims), in canonical
 * order. src[k * m + j] is the entry index of output position k in operand j,
 * or -1 where operand j holds its background there. When `need_all` is
 * non-NULL, a position is kept only if every operand j with need_all[j] set has
 * an entry there (Times over zero backgrounds: the intersection). */
typedef struct {
    SparsePattern* pat;      /* one reference */
    int64_t*       src;      /* nnz * m */
    int            m;
} SpUnion;

bool sp_union(const SparseArray* const* ops, int m, const bool* need_all, SpUnion* u);
void sp_union_free(SpUnion* u);

/* A machine complex with a zero imaginary part is a Real (the scalar
 * arithmetic collapses Complex[x, 0.] to x), so a C128 buffer holding one is
 * not canonical: convert it to a boxed buffer of Reals and Complexes (which
 * sp_new then re-types). No-op for any other kind or a buffer with no such
 * element. */
void sp_c128_canonicalise(SpValKind* kind, void** vals, int64_t n);

/* Is `e` an exact or machine zero (0, 0., 0. + 0. I)? */
bool sp_is_zero(const Expr* e);

/* Typed Plus / Times / Power fast paths (sparse_arith.c). Each returns the
 * result, or NULL to let the generic boxed path answer. Borrow everything. */
Expr* sp_arith_scalar(bool is_plus, const SparseArray* sa, const Expr* c);
Expr* sp_arith_binary(bool is_plus, const SparseArray* a, const SparseArray* b);
Expr* sp_power_scalar(const SparseArray* sa, const Expr* c, bool sparse_is_base);

#endif /* MATHILDA_SPARSE_INTERNAL_H */
