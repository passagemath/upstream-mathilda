/*
 * bitshiftleft.c
 *
 * BitShiftLeft[n, k] -- the bits of the integer n shifted k places to the
 * left, padding with zeros on the right: n * 2^k. BitShiftLeft[n] shifts by
 * one place. A negative k shifts RIGHT, and a right shift truncates toward
 * zero: BitShiftLeft[n, -k] is IntegerPart[n / 2^k], so BitShiftLeft[-5, -1]
 * is -2 (not the arithmetic shift's -3). Works on integers of any size.
 *
 * This file also holds bitshift_eval, the core BitShiftRight
 * (bitshiftright.c) shares: the two heads differ only in the sign of k.
 *
 * Attributes Listable, Protected. The evaluator threads plain Lists before the
 * builtin runs, so it only ever sees one (n, k) pair.
 *
 * Domain:
 *   - n an Integer: the exact shift (int64 or GMP).
 *   - n another exact or inexact number: n * 2^k. Machine reals scale exactly
 *     with ldexp and MPFR with mpfr_mul_2si, keeping the precision; Rational
 *     goes through Times. A Complex n shifts each part, so integer parts
 *     still truncate (BitShiftRight[3 + 5 I, 1] is 1 + 2 I) while a rational
 *     part halves exactly.
 *   - Indeterminate and every infinity are unchanged, as n * 2^k leaves them.
 *   - k not a machine-sized Integer but a number (1.5, 2., 2^70) gives
 *     ::intm and stays unevaluated; 0 or 3+ arguments give ::argt.
 *   - A symbolic n or k (x, Pi) stays unevaluated with no message.
 *   - A left shift too large for GMP to represent gives General::ovfl and
 *     Overflow[], instead of aborting the process; a non-integer shifted
 *     right that far gives General::unfl and Underflow[] (a machine real,
 *     a signed 0.). Zero is 0 for every count.
 *
 * Surfaces:
 *   - packed / NDArray: the int64 kernels ndk_BitShiftLeft_* and
 *     ndk_BitShiftRight_* in src/ndinteger.c answer before this builtin runs.
 *   - Compile[]: OP_BSHL_I / OP_BSHR_I in src/compile/compile_emit_int.c.
 * All of them, and the machine path below, call ndint_bitshiftleft_i64, so
 * every surface agrees bit for bit; where it declines (a bignum answer) each
 * surface falls back to this builtin's GMP path.
 */

#include "bitwise.h"
#include "ndinteger.h"   /* ndint_bitshiftleft_i64: the shared machine helper */
#include "common.h"      /* head_is, builtin_arg_error */
#include "arithmetic.h"  /* is_infinity_sym and friends */
#include "eval.h"
#include "message.h"
#include "print.h"       /* expr_to_string */
#include "sym_names.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>

/* The largest shift GMP can honour: an mpz holds at most INT_MAX limbs, and
 * mp_bitcnt_t is an unsigned long, only 32 bits on LLP64 Windows. A larger
 * left shift has no representation at all (the BitSet bound). */
static bool bitshift_count_representable(uint64_t s) {
    if (s > (uint64_t)ULONG_MAX) return false;
    return s < (uint64_t)INT_MAX * GMP_NUMB_BITS;
}

/* General::ovfl and Overflow[], for a left shift past bitshift_count_representable. */
static Expr* bitshift_overflow(void) {
    mth_message("General", "ovfl", "Overflow occurred in computation.");
    return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
}

/* General::unfl and Underflow[]: a non-integer shifted right past every
 * representable exponent. (An integer there is simply 0.) */
static Expr* bitshift_underflow(void) {
    mth_message("General", "unfl", "Underflow occurred in computation.");
    return expr_new_function(expr_new_symbol(SYM_Underflow), NULL, 0);
}

/* `Head::intm: Machine-sized integer expected at position 2 in <call>.` */
static Expr* bitshift_emit_intm(const char* head, Expr* res) {
    char* call_str = expr_to_string(res);
    mth_message(head, "intm", "Machine-sized integer expected at position 2 in %s.",
                call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* An exact integer (int64 or bignum) shifted LEFT by s, where a negative s
 * shifts right and truncates toward zero. The operand may be shared and is
 * never written. A left shift GMP cannot represent gives Overflow[]. */
static Expr* bitshift_integer(const Expr* n, int64_t s) {
    if (n->type == EXPR_INTEGER) {
        int64_t r;
        if (ndint_bitshiftleft_i64(n->data.integer, s, &r)) return expr_new_integer(r);
        /* Declined: a left shift whose answer is a bignum (so s > 0). */
    }

    if (s >= 0) {
        if (!bitshift_count_representable((uint64_t)s)) return bitshift_overflow();
        /* Shift into a fresh node; never mpz_set_si, whose long is 32 bits
         * on LLP64. */
        Expr* r;
        if (n->type == EXPR_INTEGER) {
            r = expr_new_bigint_from_int64(n->data.integer);
            mpz_mul_2exp(r->data.bigint, r->data.bigint, (mp_bitcnt_t)s);
        } else {
            r = expr_new_bigint_from_int64(0);
            mpz_mul_2exp(r->data.bigint, n->data.bigint, (mp_bitcnt_t)s);
        }
        return expr_bigint_normalize(r);
    }

    /* A bignum shifted right, truncating toward zero as IntegerPart does. A
     * shift past every bit (or past mp_bitcnt_t) leaves 0. The result often
     * fits in int64 again, and expr_bigint_normalize demotes it. */
    uint64_t m = (uint64_t)0 - (uint64_t)s;
    if (m > (uint64_t)ULONG_MAX || m >= mpz_sizeinbase(n->data.bigint, 2))
        return expr_new_integer(0);
    Expr* r = expr_new_bigint_from_int64(0);
    mpz_tdiv_q_2exp(r->data.bigint, n->data.bigint, (mp_bitcnt_t)m);
    return expr_bigint_normalize(r);
}

/* n * 2^s for a number n (Integer, Real, Rational, MPFR, Complex), or NULL to
 * leave the call unevaluated. Indeterminate and the infinities are handled by
 * the caller. */
static Expr* bitshift_number(const Expr* n, int64_t s) {
    if (n->type == EXPR_INTEGER || n->type == EXPR_BIGINT)
        return bitshift_integer(n, s);

    if (n->type == EXPR_REAL) {
        /* ldexp is exact while the result stays a normal double. Clamp the
         * count to a range that already overflows or underflows every double,
         * so it fits an int. On overflow or underflow, fall through to Times,
         * whose machine arithmetic decides what an out-of-range real becomes. */
        double x = n->data.real;
        int e = s > 4096 ? 4096 : s < -4096 ? -4096 : (int)s;
        double r = ldexp(x, e);
        if (x == 0.0 || (isfinite(r) && fabs(r) >= DBL_MIN)) return expr_new_real(r);
        /* Past every exponent Times could represent, a machine real shifted
         * right is a signed zero. */
        if (!bitshift_count_representable(s >= 0 ? (uint64_t)s : (uint64_t)0 - (uint64_t)s)
            && s < 0)
            return expr_new_real(copysign(0.0, x));
    }

#ifdef USE_MPFR
    if (n->type == EXPR_MPFR) {
        /* Scaling by 2^s is exact in binary floating point, so the result
         * keeps n's precision bit for bit (Times would treat 2^s as another
         * operand and re-derive it). The count bound keeps the cast to long
         * safe where long is 32 bits; past it, or if the exponent leaves
         * MPFR's range, Times decides. */
        if (s > -(INT64_C(1) << 30) && s < (INT64_C(1) << 30)) {
            Expr* r = expr_new_mpfr_copy(n->data.mpfr);
            mpfr_mul_2si(r->data.mpfr, r->data.mpfr, (long)s, MPFR_RNDN);
            if (mpfr_regular_p(r->data.mpfr) || mpfr_zero_p(n->data.mpfr)) return r;
            expr_free(r);
        }
    }
#endif

    if (head_is(n, SYM_Complex) && n->data.function.arg_count == 2) {
        Expr* re = bitshift_number(n->data.function.args[0], s);
        if (!re) return NULL;
        Expr* im = bitshift_number(n->data.function.args[1], s);
        if (!im) { expr_free(re); return NULL; }
        Expr* parts[2] = { re, im };
        Expr* call = expr_new_function(expr_new_symbol(SYM_Complex), parts, 2);
        Expr* r = evaluate(call);
        expr_free(call);
        return r;
    }

    if (!expr_is_numeric_like(n)) return NULL;

    /* Rational, MPFR, or a machine real that left the double range: n * 2^s.
     * Power[2, s] is an exact integer or rational, so precision is that of n. */
    if (!bitshift_count_representable(s >= 0 ? (uint64_t)s : (uint64_t)0 - (uint64_t)s))
        return s >= 0 ? bitshift_overflow() : bitshift_underflow();
    Expr* pw_args[2] = { expr_new_integer(2), expr_new_integer(s) };
    Expr* t_args[2] = { expr_copy((Expr*)n),
                        expr_new_function(expr_new_symbol(SYM_Power), pw_args, 2) };
    Expr* call = expr_new_function(expr_new_symbol(SYM_Times), t_args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

Expr* bitshift_eval(Expr* res, const char* head, int dir) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error(head, argc, 1, 2);

    int64_t k = 1;
    if (argc == 2) {
        Expr* kx = res->data.function.args[1];
        if (kx->type != EXPR_INTEGER) {
            /* A number that is not a machine integer (1.5, 2., 2^70, 1/2) is
             * an error; anything symbolic simply waits. */
            if (expr_is_numeric_like(kx)) return bitshift_emit_intm(head, res);
            return NULL;
        }
        k = kx->data.integer;
    }

    /* The left-shift count. -INT64_MIN saturates to INT64_MAX: both are far
     * past any representable shift, so the answer (0 stays 0, anything else
     * overflows) is the same. */
    int64_t s = dir > 0 ? k : (k == INT64_MIN ? INT64_MAX : -k);

    Expr* n = res->data.function.args[0];
    /* n * 2^k leaves Indeterminate and every infinity as it is. */
    if ((n->type == EXPR_SYMBOL && n->data.symbol.name == SYM_Indeterminate)
        || is_infinity_sym(n) || is_complex_infinity_sym(n)
        || is_neg_infinity_form(n) || head_is(n, SYM_DirectedInfinity))
        return expr_copy(n);
    return bitshift_number(n, s);
}

Expr* builtin_bitshiftleft(Expr* res) {
    return bitshift_eval(res, "BitShiftLeft", +1);
}
