/* End-to-end tests for the NDArray[...] dense N-dimensional array type.
 *
 * Each assert_eval_eq drives the full pipeline: parse -> evaluate -> print,
 * so these exercise the parser, the evaluator/builtins, and the printer
 * together. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ---------- Construction, round-trip, printing ---------- */

void test_ndarray_construct_rank2() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}]",
                   "NDArray[{{1.0, 2.0}, {3.0, 4.0}}]", 0);
}

void test_ndarray_construct_rank1() {
    assert_eval_eq("NDArray[{1, 2, 3}]", "NDArray[{1.0, 2.0, 3.0}]", 0);
}

void test_ndarray_construct_rank3() {
    assert_eval_eq("NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]",
                   "NDArray[{{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}}]", 0);
}

void test_ndarray_construct_rank4() {
    assert_eval_eq("Dimensions[NDArray[{{{{1,2,3}}}, {{{4,5,6}}}}]]", "{2, 1, 1, 3}", 0);
}

void test_ndarray_construct_reals() {
    assert_eval_eq("NDArray[{{1.5, 2.5}, {3.5, 4.5}}]",
                   "NDArray[{{1.5, 2.5}, {3.5, 4.5}}]", 0);
}

void test_ndarray_reject_jagged() {
    /* Ragged nested lists can't form a rectangular array: a NDArray::ragged
     * warning is printed (to stderr, not checked here) and the expression is
     * left unevaluated. Covers unequal sublist lengths and mixed scalar/list
     * siblings. */
    assert_eval_eq("NDArray[{{1, 2}, {3}}]", "NDArray[{{1, 2}, {3}}]", 0);
    assert_eval_eq("NDArray[{1, {2}}]", "NDArray[{1, {2}}]", 0);
    assert_eval_eq("NDArray[{1, {2, 3}}]", "NDArray[{1, {2, 3}}]", 0);
}

void test_ndarray_reject_symbolic() {
    /* A symbolic entry can't be packed to machine precision -- unevaluated. */
    assert_eval_eq("NDArray[{{1, x}, {3, 4}}]", "NDArray[{{1, x}, {3, 4}}]", 0);
}

void test_ndarray_reject_non_list() {
    assert_eval_eq("NDArray[5]", "NDArray[5]", 0);
}

void test_ndarray_normal_roundtrip() {
    /* Normal[NDArray[list]] gives back the original list (as reals). */
    assert_eval_eq("Normal[NDArray[{{1, 2}, {3, 4}}]]",
                   "{{1.0, 2.0}, {3.0, 4.0}}", 0);
    assert_eval_eq("Normal[NDArray[{1, 2, 3}]]", "{1.0, 2.0, 3.0}", 0);
    /* Rank-3 round-trip. */
    assert_eval_eq("Normal[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]",
                   "{{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}}", 0);
}

/* ---------- Head / Dimensions / Depth / Length / predicates ---------- */

void test_ndarray_head() {
    assert_eval_eq("Head[NDArray[{{1, 2}, {3, 4}}]]", "NDArray", 0);
}

void test_ndarray_dimensions() {
    assert_eval_eq("Dimensions[NDArray[{{1, 2}, {3, 4}}]]", "{2, 2}", 0);
    assert_eval_eq("Dimensions[NDArray[{1, 2, 3}]]", "{3}", 0);
    assert_eval_eq("Dimensions[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]",
                   "{2, 2, 2}", 0);
}

void test_ndarray_depth() {
    /* Depth of a rank-r NDArray matches the equivalent nested list: r + 1. */
    assert_eval_eq("Depth[NDArray[{1, 2, 3}]]", "2", 0);
    assert_eval_eq("Depth[NDArray[{{1, 2}, {3, 4}}]]", "3", 0);
    assert_eval_eq("Depth[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]", "4", 0);
}

void test_ndarray_length() {
    /* Length is the leading-axis length (numpy len). */
    assert_eval_eq("Length[NDArray[{1, 2, 3}]]", "3", 0);
    assert_eval_eq("Length[NDArray[{{1, 2}, {3, 4}}]]", "2", 0);
    assert_eval_eq("Length[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]", "2", 0);
}

void test_ndarray_ndarrayq() {
    assert_eval_eq("NDArrayQ[NDArray[{{1, 2}, {3, 4}}]]", "True", 0);
    assert_eval_eq("NDArrayQ[{{1, 2}, {3, 4}}]", "False", 0);
    assert_eval_eq("NDArrayQ[5]", "False", 0);
}

void test_ndarray_matrixq() {
    assert_eval_eq("MatrixQ[NDArray[{{1, 2}, {3, 4}}]]", "True", 0);
    assert_eval_eq("MatrixQ[NDArray[{1, 2, 3}]]", "False", 0);
}

void test_ndarray_vectorq() {
    assert_eval_eq("VectorQ[NDArray[{1, 2, 3}]]", "True", 0);
    assert_eval_eq("VectorQ[NDArray[{{1, 2}, {3, 4}}]]", "False", 0);
}

void test_ndarray_listq_always_false() {
    /* Visibly distinct from List: ListQ must never call an NDArray a List,
     * even though a plain nested-list equivalent answers True. */
    assert_eval_eq("ListQ[NDArray[{{1, 2}, {3, 4}}]]", "False", 0);
    assert_eval_eq("ListQ[{{1, 2}, {3, 4}}]", "True", 0);
}

/* ---------- Equality / structural comparison ---------- */

void test_ndarray_equality() {
    /* NDArray values compare as their Lists do: equal ones are True and
     * ones with a differing element are False, as {1, 2} == {1, 3} is. */
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] == NDArray[{{1, 2}, {3, 4}}]", "True", 0);
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] == NDArray[{{1, 2}, {3, 5}}]", "False", 0);
    assert_eval_eq("{{1, 2}, {3, 4}} == {{1, 2}, {3, 5}}", "False", 0);
}

void test_ndarray_sameq() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] === NDArray[{{1, 2}, {3, 4}}]", "True", 0);
    /* Not structurally identical to the equivalent nested List -- distinct types. */
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] === {{1.0, 2.0}, {3.0, 4.0}}", "False", 0);
}

/* ---------- Dot fast path ---------- */

void test_ndarray_dot_matrix_matrix() {
    assert_eval_eq("Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{5, 6}, {7, 8}}]]",
                   "NDArray[{{19.0, 22.0}, {43.0, 50.0}}]", 0);
}

void test_ndarray_dot_vector_vector() {
    assert_eval_eq("Dot[NDArray[{1, 2, 3}], NDArray[{4, 5, 6}]]", "32.0", 0);
}

void test_ndarray_dot_matrix_vector() {
    assert_eval_eq("Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{5, 6}]]",
                   "NDArray[{17.0, 39.0}]", 0);
}

void test_ndarray_dot_matches_list_path() {
    /* The NDArray fast path and the generic nested-list path must agree
     * (modulo Integer vs. machine-Real representation). */
    assert_eval_eq("Normal[Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{5, 6}, {7, 8}}]]]",
                   "{{19.0, 22.0}, {43.0, 50.0}}", 0);
    assert_eval_eq("Dot[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]", "{{19, 22}, {43, 50}}", 0);
}

void test_ndarray_dot_shape_mismatch() {
    /* Incompatible inner dimensions: Dot stays unevaluated (error printed
     * to stderr, not checked here). */
    assert_eval_eq("Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]]",
                   "NDArray[{{1.0, 2.0}, {3.0, 4.0}}] . "
                   "NDArray[{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}]", 0);
}

/* ---------- Elementwise Plus / Times fast paths ---------- */

void test_ndarray_plus_elementwise() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] + NDArray[{{5, 6}, {7, 8}}]",
                   "NDArray[{{6.0, 8.0}, {10.0, 12.0}}]", 0);
}

void test_ndarray_times_elementwise() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] * NDArray[{{5, 6}, {7, 8}}]",
                   "NDArray[{{5.0, 12.0}, {21.0, 32.0}}]", 0);
}

void test_ndarray_plus_matches_list_path() {
    assert_eval_eq("{{1, 2}, {3, 4}} + {{5, 6}, {7, 8}}", "{{6, 8}, {10, 12}}", 0);
    assert_eval_eq("Normal[NDArray[{{1, 2}, {3, 4}}] + NDArray[{{5, 6}, {7, 8}}]]",
                   "{{6.0, 8.0}, {10.0, 12.0}}", 0);
}

void test_ndarray_times_matches_list_path() {
    assert_eval_eq("{{1, 2}, {3, 4}} * {{5, 6}, {7, 8}}", "{{5, 12}, {21, 32}}", 0);
    assert_eval_eq("Normal[NDArray[{{1, 2}, {3, 4}}] * NDArray[{{5, 6}, {7, 8}}]]",
                   "{{5.0, 12.0}, {21.0, 32.0}}", 0);
}

void test_ndarray_plus_shape_mismatch_degrades() {
    /* Mismatched shapes: falls through and stays symbolic/unevaluated. */
    assert_eval_eq("NDArray[{1, 2}] + NDArray[{1, 2, 3}]",
                   "NDArray[{1.0, 2.0}] + NDArray[{1.0, 2.0, 3.0}]", 0);
}

void test_ndarray_symbolic_combine_degrades() {
    /* NDArray is purely numeric: combining with a symbolic term prints an
     * NDArray::sym warning (to stderr, not checked here) and stays unevaluated.
     * A numeric scalar still broadcasts (no warning), which is verified by the
     * scalar-broadcast tests. */
    assert_eval_eq("NDArray[{1., 3.}] + a",
                   "a + NDArray[{1.0, 3.0}]", 0);
    assert_eval_eq("c * NDArray[{1., 3.}]",
                   "c NDArray[{1.0, 3.0}]", 0);
    assert_eval_eq("NDArray[{1., 3.}] ^ n",
                   "NDArray[{1.0, 3.0}]^n", 0);
    assert_eval_eq("b ^ NDArray[{1., 3.}]",
                   "b^NDArray[{1.0, 3.0}]", 0);
}

/* ---------- DataType (dtype) ---------- */

void test_ndarray_datatype_default() {
    /* Default dtype is float64; DataType[] reports it. */
    assert_eval_eq("DataType[NDArray[{1, 2, 3}]] === \"float64\"", "True", 0);
    /* DataType stays symbolic on a non-array. */
    assert_eval_eq("DataType[5]", "DataType[5]", 0);
}

void test_ndarray_datatype_infer_complex() {
    /* With no explicit DataType, complex leaves auto-infer complex64 rather than
     * failing to pack and leaving NDArray[...] unevaluated. */
    assert_eval_eq("NDArrayQ[NDArray[{Complex[1, 2], Complex[-1, 3]}]]", "True", 0);
    assert_eval_eq("DataType[NDArray[{1.0 + 2.0*I, -0.5 + 0.3*I}]] === \"complex64\"",
                   "True", 0);
    /* An element-wise map over the inferred-complex array matches the List path. */
    assert_eval_eq("Normal[Sin[NDArray[{1.0 + 2.0*I, -0.5 + 0.3*I}]]] == "
                   "Sin[{1.0 + 2.0*I, -0.5 + 0.3*I}]", "True", 0);
    /* Inference also fires for nested (matrix) complex data. */
    assert_eval_eq("DataType[NDArray[{{1.0 + I, 2.0}, {3.0, 4.0*I}}]] === \"complex64\"",
                   "True", 0);
    /* Purely real data still defaults to float64. */
    assert_eval_eq("DataType[NDArray[{1.0, 2.0, 3.0}]] === \"float64\"", "True", 0);
}

void test_ndarray_datatype_options() {
    /* The option is surfaced by Options[NDArray]. */
    assert_eval_eq("Options[NDArray] === {DataType -> \"float64\"}", "True", 0);
}

void test_ndarray_construct_dtypes() {
    assert_eval_eq("NDArrayQ[NDArray[{1, 2, 3}, DataType -> \"float32\"]]", "True", 0);
    assert_eval_eq("DataType[NDArray[{1, 2, 3}, DataType -> \"float32\"]] === \"float32\"", "True", 0);
    assert_eval_eq("NDArrayQ[NDArray[{Complex[1, 2]}, DataType -> \"complex64\"]]", "True", 0);
    assert_eval_eq("DataType[NDArray[{Complex[1, 2]}, DataType -> \"complex32\"]] === \"complex32\"", "True", 0);
    /* A complex dtype accepts bare real leaves too (im = 0). */
    assert_eval_eq("NDArrayQ[NDArray[{1, 2}, DataType -> \"complex64\"]]", "True", 0);
}

void test_ndarray_dtype_roundtrip() {
    /* float32 values that are exactly representable round-trip cleanly. */
    assert_eval_eq("Normal[NDArray[{1.5, 2.5}, DataType -> \"float32\"]]", "{1.5, 2.5}", 0);
    /* Complex dtypes rebuild Complex[] leaves on Normal. */
    assert_eval_eq("Normal[NDArray[{Complex[1, 2], Complex[3, 4]}, DataType -> \"complex32\"]]",
                   "{1.0 + 2.0*I, 3.0 + 4.0*I}", 0);
}

void test_ndarray_dtype_narrowing() {
    /* 0.1 isn't representable in float32, so the stored value differs from the
     * float64 literal. */
    assert_eval_eq("First[Normal[NDArray[{0.1}, DataType -> \"float32\"]]] == 0.1", "False", 0);
}

void test_ndarray_dtype_equality_distinguishes() {
    /* Same values, different dtype -> not SameQ (dtype is part of identity). */
    assert_eval_eq("NDArray[{1, 2}, DataType -> \"float32\"] === NDArray[{1, 2}, DataType -> \"float64\"]",
                   "False", 0);
    assert_eval_eq("NDArray[{1, 2}, DataType -> \"float32\"] === NDArray[{1, 2}, DataType -> \"float32\"]",
                   "True", 0);
}

void test_ndarray_complex_plus() {
    assert_eval_eq("NDArray[{Complex[1, 1]}, DataType -> \"complex64\"] + "
                   "NDArray[{Complex[2, 3]}, DataType -> \"complex64\"]",
                   "NDArray[{3.0 + 4.0*I}]", 0);
}

void test_ndarray_mixed_dtype_promotion() {
    assert_eval_eq("DataType[NDArray[{1, 2}, DataType -> \"float32\"] + NDArray[{3, 4}]] === \"float64\"",
                   "True", 0);
    assert_eval_eq("DataType[NDArray[{1, 2}, DataType -> \"float32\"] + "
                   "NDArray[{3, 4}, DataType -> \"float32\"]] === \"float32\"", "True", 0);
}

/* ---------- Boolean dtype ---------- */

void test_ndarray_bool_construct() {
    /* Constructs from True/False leaves and round-trips through DataType. */
    assert_eval_eq("NDArrayQ[NDArray[{True, False, True}, DataType -> \"bool\"]]", "True", 0);
    assert_eval_eq("DataType[NDArray[{True, False}, DataType -> \"bool\"]] === \"bool\"", "True", 0);
    /* The "Boolean" alias is accepted on input and canonicalises to "bool". */
    assert_eval_eq("DataType[NDArray[{True, False}, DataType -> \"Boolean\"]] === \"bool\"", "True", 0);
    /* A non-True/False leaf declines the bool dtype (stays unevaluated). */
    assert_eval_eq("NDArrayQ[NDArray[{1, 0}, DataType -> \"bool\"]]", "False", 0);
}

void test_ndarray_bool_roundtrip() {
    assert_eval_eq("Normal[NDArray[{True, False, True}, DataType -> \"bool\"]]",
                   "{True, False, True}", 0);
    /* Rank-2 bool arrays round-trip too. */
    assert_eval_eq("Normal[NDArray[{{True, False}, {False, True}}, DataType -> \"bool\"]]",
                   "{{True, False}, {False, True}}", 0);
    /* A packed bool list is indistinguishable from the plain list. */
    assert_eval_eq("ToNDArray[{True, False, True}] === {True, False, True}", "True", 0);
    assert_eval_eq("DataType[ToNDArray[{True, False, True}]] === \"bool\"", "True", 0);
    /* Mixed True/number does not pack to bool. */
    assert_eval_eq("NDArrayQ[ToNDArray[{True, 1, False}]]", "False", 0);
}

void test_ndarray_bool_not_numeric() {
    /* Arithmetic and transcendentals over a bool array delist to the symbolic
     * List, exactly as Mathematica leaves Sin[True] / True + True unevaluated. */
    assert_eval_eq("Sin[ToNDArray[{True, False}]]", "{Sin[True], Sin[False]}", 0);
    assert_eval_eq("Sin[ToNDArray[{True, False}]] === Sin[{True, False}]", "True", 0);
    assert_eval_eq("ToNDArray[{True, False}] + ToNDArray[{True, True}] === "
                   "{True, False} + {True, True}", "True", 0);
}

void test_ndarray_bool_producers_consumers() {
    /* Sign predicates now answer with a packed bool array (gap C.1 closed). */
    assert_eval_eq("DataType[Positive[ToNDArray[{-1, 0, 2}]]] === \"bool\"", "True", 0);
    assert_eval_eq("Positive[ToNDArray[{-1, 0, 2}]] === Positive[{-1, 0, 2}]", "True", 0);
    assert_eval_eq("Negative[ToNDArray[{-1., 0., 2.}]]", "{True, False, False}", 0);
    /* Boole bridges a bool array back to an int64 one. */
    assert_eval_eq("DataType[Boole[NDArray[{True, False, True}, DataType -> \"bool\"]]] === \"int64\"",
                   "True", 0);
    assert_eval_eq("Normal[Boole[NDArray[{True, False, True}, DataType -> \"bool\"]]]",
                   "{1, 0, 1}", 0);
    /* AllTrue/AnyTrue/NoneTrue scan the byte buffer; the whole pipeline stays packed. */
    assert_eval_eq("AllTrue[Positive[ToNDArray[{1, 2, 3}]], TrueQ]", "True", 0);
    assert_eval_eq("AllTrue[ToNDArray[{True, False}], TrueQ]", "False", 0);
    assert_eval_eq("AnyTrue[ToNDArray[{False, False, True}], TrueQ]", "True", 0);
    assert_eval_eq("NoneTrue[ToNDArray[{False, False}], TrueQ]", "True", 0);
}

/* ---------- Scalar broadcasting (numpy-style) ---------- */

void test_ndarray_scalar_broadcast_plus() {
    assert_eval_eq("1 + NDArray[{1, 2, 3}]", "NDArray[{2.0, 3.0, 4.0}]", 0);
    /* Weak scalar: an integer/real scalar keeps the array's float width. */
    assert_eval_eq("DataType[1 + NDArray[{1, 2}, DataType -> \"float32\"]] === \"float32\"", "True", 0);
    /* A complex scalar moves the array onto the complex axis (same width). */
    assert_eval_eq("DataType[I + NDArray[{1, 2}, DataType -> \"float32\"]] === \"complex32\"", "True", 0);
}

void test_ndarray_scalar_broadcast_times() {
    assert_eval_eq("3 * NDArray[{1., 2., 3.}]", "NDArray[{3.0, 6.0, 9.0}]", 0);
    /* Subtraction of arrays reduces via Plus + broadcast Times[-1, .]. */
    assert_eval_eq("NDArray[{1., 2., 3.}] - NDArray[{1., 1., 1.}]", "NDArray[{0.0, 1.0, 2.0}]", 0);
}

/* ---------- Power fast path ---------- */

void test_ndarray_power_scalar() {
    assert_eval_eq("NDArray[{2., 3.}]^2", "NDArray[{4.0, 9.0}]", 0);
    assert_eval_eq("NDArray[{4., 9.}]^0.5", "NDArray[{2.0, 3.0}]", 0);
    /* A negative real base with a fractional exponent promotes to complex. */
    assert_eval_eq("DataType[NDArray[{-1.0}]^0.5] === \"complex64\"", "True", 0);
}

void test_ndarray_power_complex_exact() {
    /* I^2 == -1 exactly (integer exponent uses exact complex multiplication). */
    assert_eval_eq("NDArray[{Complex[0, 1]}, DataType -> \"complex64\"]^2",
                   "NDArray[{-1.0 + 0.0*I}]", 0);
}

void test_ndarray_power_elementwise_and_scalar_base() {
    assert_eval_eq("NDArray[{2., 3.}]^NDArray[{3., 2.}]", "NDArray[{8.0, 9.0}]", 0);
    assert_eval_eq("2^NDArray[{1., 2., 3.}]", "NDArray[{2.0, 4.0, 8.0}]", 0);
}

/* ---------- Complex Dot / BLAS closed system ---------- */

void test_ndarray_complex_dot() {
    /* (1+i, 2) . (1, i) = (1+i) + 2 i = 1 + 3 i. */
    assert_eval_eq("NDArray[{Complex[1, 1], 2}, DataType -> \"complex64\"] . "
                   "NDArray[{1, Complex[0, 1]}, DataType -> \"complex64\"]",
                   "1.0 + 3.0*I", 0);
}

void test_ndarray_float32_dot() {
    assert_eval_eq("DataType[NDArray[{{1., 2.}, {3., 4.}}, DataType -> \"float32\"] . "
                   "NDArray[{1., 1.}, DataType -> \"float32\"]] === \"float32\"", "True", 0);
    assert_eval_eq("Normal[NDArray[{{1., 2.}, {3., 4.}}, DataType -> \"float32\"] . "
                   "NDArray[{1., 1.}, DataType -> \"float32\"]]", "{3.0, 7.0}", 0);
}

/* ---------- Part (indexing) ---------- */

void test_ndarray_part_scalar() {
    /* Full integer indexing yields a scalar leaf. */
    assert_eval_eq("NDArray[{10., 20., 30., 40.}][[2]]", "20.0", 0);
    assert_eval_eq("NDArray[{10., 20., 30., 40.}][[-1]]", "40.0", 0);
    assert_eval_eq("NDArray[{{1., 2., 3.}, {4., 5., 6.}}][[2, 3]]", "6.0", 0);
    /* Complex leaf rebuilds Complex[]. */
    assert_eval_eq("NDArray[{Complex[1, 2], Complex[3, -4]}][[2]]", "3.0 - 4.0*I", 0);
    /* [[0]] gives the head, and it is NDArray (not List). */
    assert_eval_eq("NDArray[{1., 2., 3.}][[0]]", "NDArray", 0);
}

void test_ndarray_part_subarray() {
    /* A partial integer index yields a sub-NDArray of the trailing shape. */
    assert_eval_eq("NDArrayQ[NDArray[{{1., 2., 3.}, {4., 5., 6.}}][[1]]]", "True", 0);
    assert_eval_eq("Normal[NDArray[{{1., 2., 3.}, {4., 5., 6.}}][[1]]]", "{1.0, 2.0, 3.0}", 0);
    assert_eval_eq("Dimensions[NDArray[Table[1.0*i + j, {i, 2}, {j, 3}]][[2]]]", "{3}", 0);
}

void test_ndarray_part_matches_list_path() {
    /* Native integer indexing agrees with the delisted reference, including a
     * rank-3 array and negative subscripts. */
    assert_eval_eq("Module[{t = NDArray[Table[100.0 i + 10.0 j + k, {i, 2}, {j, 3}, {k, 4}]], "
                   "l = Table[100.0 i + 10.0 j + k, {i, 2}, {j, 3}, {k, 4}]}, "
                   "{Normal[t[[2]]] == l[[2]], Normal[t[[1, 3]]] == l[[1, 3]], "
                   "t[[2, 3, 4]] == l[[2, 3, 4]], t[[-1, -1, -1]] == l[[-1, -1, -1]]}]",
                   "{True, True, True, True}", 0);
}

void test_ndarray_part_slices_stay_packed() {
    /* Span / All / a List of positions produce a packed sub-NDArray natively. */
    assert_eval_eq("NDArrayQ[NDArray[{10., 20., 30., 40.}][[2 ;; 3]]]", "True", 0);
    assert_eval_eq("Normal[NDArray[{10., 20., 30., 40.}][[2 ;; 3]]]", "{20.0, 30.0}", 0);
    /* Stepped and reversed spans. */
    assert_eval_eq("Normal[NDArray[{10., 20., 30., 40.}][[1 ;; 4 ;; 2]]]", "{10.0, 30.0}", 0);
    assert_eval_eq("Normal[NDArray[{10., 20., 30., 40.}][[-1 ;; 1 ;; -1]]]",
                   "{40.0, 30.0, 20.0, 10.0}", 0);
    /* Fancy (list of positions) gather stays packed. */
    assert_eval_eq("NDArrayQ[NDArray[{10., 20., 30., 40.}][[{1, 3}]]]", "True", 0);
    assert_eval_eq("Normal[NDArray[{10., 20., 30., 40.}][[{1, 3}]]]", "{10.0, 30.0}", 0);
    /* A column of a matrix (All on axis 1, integer on axis 2) is a packed vector. */
    assert_eval_eq("NDArrayQ[NDArray[{{1., 2.}, {3., 4.}}][[All, 2]]]", "True", 0);
    assert_eval_eq("Normal[NDArray[{{1., 2.}, {3., 4.}}][[All, 2]]]", "{2.0, 4.0}", 0);
    /* A sub-matrix (Span on both axes) keeps rank 2. */
    assert_eval_eq("Dimensions[NDArray[{{1., 2., 3.}, {4., 5., 6.}}][[All, 2 ;; 3]]]",
                   "{2, 2}", 0);
    assert_eval_eq("Normal[NDArray[{{1., 2., 3.}, {4., 5., 6.}}][[All, 2 ;; 3]]]",
                   "{{2.0, 3.0}, {5.0, 6.0}}", 0);
    /* dtype is preserved through a slice. */
    assert_eval_eq("DataType[NDArray[{1., 2., 3., 4.}, DataType -> \"float32\"][[2 ;; 3]]] "
                   "=== \"float32\"", "True", 0);
    assert_eval_eq("DataType[NDArray[{Complex[1, 1], 2, 3}][[2 ;; 3]]] === \"complex64\"",
                   "True", 0);
}

void test_ndarray_part_slices_match_list_path() {
    /* Mixed integer / All / Span / fancy subscripts over a rank-3 array agree
     * with the delisted reference (Normal[t][[...]]). */
    assert_eval_eq("Module[{t = NDArray[Table[100.0 i + 10.0 j + k, {i, 3}, {j, 4}, {k, 5}]], "
                   "l = Table[100.0 i + 10.0 j + k, {i, 3}, {j, 4}, {k, 5}]}, "
                   "{Normal[t[[2 ;; 3]]] == l[[2 ;; 3]], "
                   "Normal[t[[All, 2 ;; 4 ;; 2]]] == l[[All, 2 ;; 4 ;; 2]], "
                   "Normal[t[[2, All, 3]]] == l[[2, All, 3]], "
                   "Normal[t[[{1, 3}, 2, {1, 4, 5}]]] == l[[{1, 3}, 2, {1, 4, 5}]], "
                   "Normal[t[[-1 ;; 1 ;; -1, 2, All]]] == l[[-1 ;; 1 ;; -1, 2, All]]}]",
                   "{True, True, True, True, True}", 0);
}

/* ---- Empty arrays ------------------------------------------------------
 * NDArray[{}] used to stay the symbolic NDArray[List[]] (Length 1). A
 * rectangular list with a zero-length axis now builds an empty array of that
 * shape, and the heads that reduce a buffer answer as on the empty List
 * instead of reading past it (Median crashed) or returning their seed
 * (Max gave Overflow[]). */
static void test_ndarray_empty_construct(void) {
    assert_eval_eq("NDArray[{}]", "NDArray[{}]", 0);
    assert_eval_eq("NDArrayQ[NDArray[{}]]", "True", 0);
    assert_eval_eq("Dimensions[NDArray[{}]]", "{0}", 0);
    assert_eval_eq("Length[NDArray[{}]]", "0", 0);
    assert_eval_eq("Normal[NDArray[{}]]", "{}", 0);
    assert_eval_eq("DataType[NDArray[{}]]", "\"float64\"", 0);
    assert_eval_eq("DataType[NDArray[{}, DataType -> \"int64\"]]", "\"int64\"", 0);
    assert_eval_eq("Dimensions[NDArray[{{}, {}}]]", "{2, 0}", 0);
    assert_eval_eq("Normal[NDArray[{{}, {}}]]", "{{}, {}}", 0);
    assert_eval_eq("Dimensions[NDArray[{{{}}}]]", "{1, 1, 0}", 0);
    assert_eval_eq("NDArray[NDArray[{}]]", "NDArray[{}]", 0);
    /* Ragged or symbolic shapes still do not form an array. */
    assert_eval_eq("NDArrayQ[NDArray[{{}, {1}}]]", "False", 0);
    assert_eval_eq("NDArrayQ[NDArray[{{}, x}]]", "False", 0);
    /* Automatic packing is unchanged: an empty result stays an ordinary List. */
    assert_eval_eq("PackedArrayQ[Range[0]]", "False", 0);
    assert_eval_eq("NDArrayQ[Select[NDArray[{1., 2.}], # > 5 &]]", "False", 0);
}

static void test_ndarray_empty_operations(void) {
    assert_eval_eq("e = NDArray[{}]; {e + 1, 2 e, Sin[e], Sort[e], Reverse[e]}",
                   "{NDArray[{}], NDArray[{}], NDArray[{}], NDArray[{}], NDArray[{}]}", 0);
    assert_eval_eq("e = NDArray[{}]; {Total[e], Map[f, e], Apply[f, e], Count[e, _]}",
                   "{0.0, {}, f[], 0}", 0);
    assert_eval_eq("e = NDArray[{}]; Join[e, NDArray[{1., 2.}]]", "NDArray[{1.0, 2.0}]", 0);
    assert_eval_eq("Dimensions[Transpose[NDArray[{{}, {}}]]]", "{0, 2}", 0);
    assert_eval_eq("e = NDArray[{}]; {Max[e], Min[e], MinMax[e]}",
                   "{-Infinity, Infinity, {Infinity, -Infinity}}", 0);
    assert_eval_eq("e = NDArray[{}, DataType -> \"int64\"]; {Max[e], Min[e]}",
                   "{-Infinity, Infinity}", 0);
    assert_eval_eq("{Max[NDArray[{{}, {}}]], MinMax[NDArray[{{}, {}}]]}",
                   "{-Infinity, {Infinity, -Infinity}}", 0);
    assert_eval_eq("e = NDArray[{}]; {Max[e, 3], Min[e, 3], Max[e, e]}", "{3, 3, -Infinity}", 0);
    assert_eval_eq("Median[NDArray[{}]]", "Median[{}]", 0);
    assert_eval_eq("Quiet[Quartiles[NDArray[{}]]]", "Quartiles[{}]", 0);
    assert_eval_eq("Rescale[NDArray[{}]]", "{}", 0);
    assert_eval_eq("Quiet[ListInterpolation[NDArray[{}]]]", "ListInterpolation[{}]", 0);
}

/* Max / Min over several arguments, one of them an array: each array is
 * reduced through its buffer path and joins the other arguments. A visible
 * NDArray used to leave the call unevaluated, and a packed one was
 * materialised element by element. */
static void test_ndarray_minmax_multi_arg(void) {
    assert_eval_eq("Max[NDArray[{1., 5.}], 2]", "5.0", 0);
    assert_eval_eq("Min[NDArray[{1., 5.}], 2]", "1.0", 0);
    assert_eval_eq("Max[NDArray[{1., 5.}], NDArray[{7., 0.}]]", "7.0", 0);
    assert_eval_eq("Max[Range[10], {20}]", "20", 0);
    assert_eval_eq("Max[Range[10], 3, {4, 50}]", "50", 0);
    assert_eval_eq("Min[Range[10], -3]", "-3", 0);
    assert_eval_eq("Max[Range[1., 1000.], x]", "Max[1000.0, x]", 0);
    assert_eval_eq("Max[NDArray[{3, 9}, DataType -> \"int64\"], 4]", "9", 0);
}

/* A Listable head with no array kernel for its arity left a visible array as
 * an atom while the List threaded. Boole of a numeric array and the
 * one-argument IntegerExponent now answer as on the List. */
static void test_ndarray_boole_integerexponent(void) {
    assert_eval_eq("Boole[NDArray[{1., 0.}]]", "{Boole[1.0], Boole[0.0]}", 0);
    assert_eval_eq("Boole[NDArray[{1, 0}, DataType -> \"int64\"]]", "{Boole[1], Boole[0]}", 0);
    assert_eval_eq("Boole[NDArray[{True, False}, DataType -> \"bool\"]]", "{1, 0}", 0);
    assert_eval_eq("Boole[NDArray[{}]]", "{}", 0);
    assert_eval_eq("IntegerExponent[NDArray[{10, 20, 300}, DataType -> \"int64\"]]",
                   "NDArray[{1, 1, 2}]", 0);
    assert_eval_eq("IntegerExponent[NDArray[{0, 10}, DataType -> \"int64\"]]", "{Infinity, 1}", 0);
    assert_eval_eq("IntegerExponent[Range[10, 100, 10]]", "{1, 1, 1, 1, 1, 1, 1, 1, 1, 2}", 0);
    assert_eval_eq("PackedArrayQ[IntegerExponent[Range[10, 100, 10]]]", "True", 0);
    assert_eval_eq("Quiet[IntegerExponent[NDArray[{10., 20.}]]]",
                   "{IntegerExponent[10.0], IntegerExponent[20.0]}", 0);
}

/* A visible array beside a plain List the lift could not absorb (unequal
 * length, symbolic entries) was threaded as a scalar. It now threads as its
 * List: a length mismatch is Thread::tdlen, a symbolic List pairs elementwise.
 * A same-shape numeric List still takes the buffer path. */
static void test_ndarray_plus_plain_list(void) {
    assert_eval_eq("Quiet[NDArray[{1., 2., 3.}] + {1, 2}]", "{1, 2} + {1.0, 2.0, 3.0}", 0);
    assert_eval_eq("Check[NDArray[{1., 2., 3.}] + {1, 2}, failed]", "failed", 0);
    assert_eval_eq("NDArray[{1., 2.}] + {x, y}", "{1.0 + x, 2.0 + y}", 0);
    assert_eval_eq("NDArray[{2., 3.}] {x, 2}", "{2.0 x, 6.0}", 0);
    assert_eval_eq("Quiet[BitShiftLeft[NDArray[{1, 2, 3}, DataType -> \"int64\"], {1, 2}]]",
                   "BitShiftLeft[{1, 2, 3}, {1, 2}]", 0);
    assert_eval_eq("NDArray[{1., 2.}] + {10, 20}", "NDArray[{11.0, 22.0}]", 0);
    assert_eval_eq("NDArray[{{1., 2.}, {3., 4.}}] + {10, 20}",
                   "NDArray[{{11.0, 12.0}, {23.0, 24.0}}]", 0);
    assert_eval_eq("NDArray[{1., 2.}] + x", "x + NDArray[{1.0, 2.0}]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_ndarray_construct_rank2);
    TEST(test_ndarray_construct_rank1);
    TEST(test_ndarray_construct_rank3);
    TEST(test_ndarray_construct_rank4);
    TEST(test_ndarray_construct_reals);
    TEST(test_ndarray_reject_jagged);
    TEST(test_ndarray_reject_symbolic);
    TEST(test_ndarray_reject_non_list);
    TEST(test_ndarray_normal_roundtrip);

    TEST(test_ndarray_head);
    TEST(test_ndarray_dimensions);
    TEST(test_ndarray_depth);
    TEST(test_ndarray_length);
    TEST(test_ndarray_ndarrayq);
    TEST(test_ndarray_matrixq);
    TEST(test_ndarray_vectorq);
    TEST(test_ndarray_listq_always_false);

    TEST(test_ndarray_equality);
    TEST(test_ndarray_sameq);

    TEST(test_ndarray_dot_matrix_matrix);
    TEST(test_ndarray_dot_vector_vector);
    TEST(test_ndarray_dot_matrix_vector);
    TEST(test_ndarray_dot_matches_list_path);
    TEST(test_ndarray_dot_shape_mismatch);

    TEST(test_ndarray_plus_elementwise);
    TEST(test_ndarray_times_elementwise);
    TEST(test_ndarray_plus_matches_list_path);
    TEST(test_ndarray_times_matches_list_path);
    TEST(test_ndarray_plus_shape_mismatch_degrades);
    TEST(test_ndarray_symbolic_combine_degrades);

    TEST(test_ndarray_datatype_default);
    TEST(test_ndarray_datatype_infer_complex);
    TEST(test_ndarray_datatype_options);
    TEST(test_ndarray_construct_dtypes);
    TEST(test_ndarray_dtype_roundtrip);
    TEST(test_ndarray_dtype_narrowing);
    TEST(test_ndarray_dtype_equality_distinguishes);
    TEST(test_ndarray_bool_construct);
    TEST(test_ndarray_bool_roundtrip);
    TEST(test_ndarray_bool_not_numeric);
    TEST(test_ndarray_bool_producers_consumers);
    TEST(test_ndarray_complex_plus);
    TEST(test_ndarray_mixed_dtype_promotion);

    TEST(test_ndarray_scalar_broadcast_plus);
    TEST(test_ndarray_scalar_broadcast_times);

    TEST(test_ndarray_power_scalar);
    TEST(test_ndarray_power_complex_exact);
    TEST(test_ndarray_power_elementwise_and_scalar_base);

    TEST(test_ndarray_complex_dot);
    TEST(test_ndarray_float32_dot);

    TEST(test_ndarray_part_scalar);
    TEST(test_ndarray_part_subarray);
    TEST(test_ndarray_part_matches_list_path);
    TEST(test_ndarray_part_slices_stay_packed);
    TEST(test_ndarray_part_slices_match_list_path);
    TEST(test_ndarray_empty_construct);
    TEST(test_ndarray_empty_operations);
    TEST(test_ndarray_minmax_multi_arg);
    TEST(test_ndarray_boole_integerexponent);
    TEST(test_ndarray_plus_plain_list);

    printf("All NDArray tests passed.\n");
    return 0;
}
