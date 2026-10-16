#!/usr/bin/env python3
"""Generate the coefficient tables of src/special_functions/sf_machine_expint_coef.h.

Every table is a near-minimax relative-error fit computed here from mpmath at 50
digits (fitlib.py: Sanathanan-Koerner + Lawson), so nothing is copied from an
external library. The STRUCTURE follows the classical double-precision designs:

  * Si/Ci  -- Cephes `sici` (S. L. Moshier, Methods and Programs for Mathematical
    Functions, 1989): a power series in x^2 for small x, and for large x the
    auxiliary functions f, g of DLMF 6.2.17-6.2.18 as rationals in 1/x^2.
  * Shi/Chi -- Cephes `shichi`: a power series in x^2 for small x.
  * Ei / E1 -- Cody & Thacher, "Rational Chebyshev approximations for the
    exponential integral E1(x)", Math. Comp. 22 (1968) 641-649, and "Chebyshev
    approximations for the exponential integral Ei(x)", Math. Comp. 23 (1969)
    289-303: a power series near the origin and e^{+-x}/x times a rational in 1/x
    beyond, on a few sub-intervals.

The large-argument factors are written as 1 -/+ u R(u) with R = P/Q, where u is
1/x or 1/x^2. Each underlying function (DLMF 6.7.13-6.7.14 for f, g; the
Stieltjes integral x e^x E1(x) = Int_0^oo e^{-s}/(1 + s/x) ds) is a Stieltjes
transform in u, so P and Q come out with coefficients of ONE sign and the Horner
evaluation has no cancellation, and the correction u R is small, which scales
its rounding error down.

Run:  benchmarks/.venv/bin/python -I tools/sf_accuracy/fit_expint.py > \
          src/special_functions/sf_machine_expint_coef.h
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fitlib import cheb_nodes, fit_rational, max_rel_error, mp  # noqa: E402

TARGET = mp.mpf("2e-18")   # approximation error, << 1 ulp (1.1e-16)


def emit(name, coeffs):
    vals = ", ".join(repr(float(c)) for c in coeffs)
    return f"static const double {name}[{len(coeffs)}] = {{ {vals} }};"


def poly_in_t(f, a, b, name, degrees):
    """Polynomial in t = (2v - a - b)/(b - a) on v in [a, b]."""
    g = lambda t: f(mp.mpf(a) + (mp.mpf(b) - a) * (t + 1) / 2)
    for d in degrees:
        p, q = fit_rational(g, d, 0, iters=1, lawson=12)
        e, _ = max_rel_error(g, p, q, M=400)
        if e < TARGET or d == degrees[-1]:
            ed, _ = max_rel_error(g, p, q, M=400, as_double=True)
            out = [f"/* {name}: degree {d} in t on [{a}, {b}]; approx err {float(e):.2e}, "
                   f"double eval {float(ed):.2e} (relative) */", emit(name, p)]
            return "\n".join(out)


def poly_in_v(f, a, b, name, degrees):
    """Polynomial in v = (x - a)/(b - a) in [0, 1].  For an ABSOLUTELY MONOTONIC
    function (every derivative positive: the Ei, Shi, Chi series) the fitted
    coefficients are all positive in this variable, so Horner never cancels --
    the centred variable would alternate them at the left end and amplify the
    rounding by f(b)/f(a)."""
    g = lambda v: f(mp.mpf(a) + (mp.mpf(b) - a) * v)
    for d in degrees:
        p, q = fit_rational(g, d, 0, iters=1, lawson=12, a=0, b=1)
        e, _ = max_rel_error(g, p, q, M=400, a=0, b=1)
        if e < TARGET or d == degrees[-1]:
            ed, _ = max_rel_error(g, p, q, M=400, as_double=True, a=0, b=1)
            out = [f"/* {name}: degree {d} in v on [{a}, {b}]; approx err {float(e):.2e}, "
                   f"double eval {float(ed):.2e} (relative) */", emit(name, p)]
            return "\n".join(out)


def rational_in_u(f, ua, ub, name, degrees):
    """R(u) = P(u)/Q(u) on u in [ua, ub], fitted in v = u/ub and rescaled to u."""
    g = lambda v: f(v * ub)
    for (m, n) in degrees:
        p, q = fit_rational(g, m, n, iters=10, lawson=12, a=ua / ub, b=1)
        e, _ = max_rel_error(g, p, q, M=400, a=ua / ub, b=1)
        if e < TARGET or (m, n) == degrees[-1]:
            ed, _ = max_rel_error(g, p, q, M=400, as_double=True, a=ua / ub, b=1)
            pu = [c / ub ** k for k, c in enumerate(p)]
            qu = [c / ub ** k for k, c in enumerate(q)]
            out = [f"/* {name}: [{m}/{n}] in u on [{float(ua):.6g}, {float(ub):.6g}]; "
                   f"approx err {float(e):.2e}, double eval {float(ed):.2e} (relative) */",
                   emit(name + "_P", pu), emit(name + "_Q", qu)]
            return "\n".join(out)


def rational_in_t(f, ua, ub, name, degrees):
    """R = P(t)/Q(t), t = (2u - ua - ub)/(ub - ua) in [-1, 1].  For a function
    whose coefficients in u alternate (x e^-x Ei(x) is NOT a Stieltjes function:
    its asymptotic series sum k!/x^k does not alternate), the centred variable is
    the well-conditioned basis."""
    g = lambda t: f(ua + (ub - ua) * (t + 1) / 2)
    for (m, n) in degrees:
        p, q = fit_rational(g, m, n, iters=10, lawson=12)
        e, _ = max_rel_error(g, p, q, M=400)
        if e < TARGET or (m, n) == degrees[-1]:
            ed, _ = max_rel_error(g, p, q, M=400, as_double=True)
            out = [f"/* {name}: [{m}/{n}] in t on u in [{float(ua):.6g}, {float(ub):.6g}]; "
                   f"approx err {float(e):.2e}, double eval {float(ed):.2e} (relative) */",
                   emit(name + "_P", p), emit(name + "_Q", q)]
            return "\n".join(out)


def series_tables():
    """Power-series coefficients for the complex kernels' Horner sums, each the
    correctly rounded double of an exact rational.  Stops at the first entry that
    underflows to zero."""
    out = []
    def table(name, gen, comment):
        vals = []
        k = 1
        while True:
            v = gen(k)
            if v is None or float(v) == 0.0:
                break
            vals.append(float(v))
            k += 1
        out.append(f"/* {comment} */")
        out.append(f"#define {name}_N {len(vals)}")
        body = ", ".join(repr(x) for x in vals)
        out.append(f"static const double {name}[{len(vals)}] = {{ {body} }};")
    table("EIN_C", lambda k: 1 / (mp.mpf(k) * mp.factorial(k)) if k <= 170 else None,
          "EIN_C[k-1] = 1/(k k!), k = 1..: Ei(z) = gamma + log z + sum EIN_C z^k")
    table("SIN_C", lambda k: 1 / (mp.mpf(2 * k - 1) * mp.factorial(2 * k - 1)) if 2 * k - 1 <= 170 else None,
          "SIN_C[k] = 1/((2k+1)(2k+1)!): Si(z) = z sum_k SIN_C[k] (-z^2)^k, Shi with +z^2")
    table("COS_C", lambda k: 1 / (mp.mpf(2 * k) * mp.factorial(2 * k)) if 2 * k <= 170 else None,
          "COS_C[k-1] = 1/(2k (2k)!), k = 1..: Ci(z) = gamma + log z + sum COS_C (-z^2)^k")
    return "\n".join(out)


def fg(x):
    s, c = mp.sin(x), mp.cos(x)
    si, ci = mp.si(x), mp.ci(x)
    return ci * s + (mp.pi / 2 - si) * c, -ci * c + (mp.pi / 2 - si) * s


def main():
    print("/* GENERATED by tools/sf_accuracy/fit_expint.py -- do not edit by hand.")
    print(" * Near-minimax relative-error fits computed with mpmath at 50 digits. */")
    print("#ifndef MATHILDA_SF_MACHINE_EXPINT_COEF_H")
    print("#define MATHILDA_SF_MACHINE_EXPINT_COEF_H\n")

    # ---- small argument: Si/Ci, Shi/Chi as series in y = x^2 on [0, 16] ---------
    def S(y):  return mp.mpf(1) if y == 0 else mp.si(mp.sqrt(y)) / mp.sqrt(y)
    def C(y):  return mp.mpf(-1) / 4 if y == 0 else (mp.ci(mp.sqrt(y)) - mp.euler - mp.log(mp.sqrt(y))) / y
    def Sh(y): return mp.mpf(1) if y == 0 else mp.shi(mp.sqrt(y)) / mp.sqrt(y)
    def Ch(y): return mp.mpf(1) / 4 if y == 0 else (mp.chi(mp.sqrt(y)) - mp.euler - mp.log(mp.sqrt(y))) / y
    print("/* Si(x) = x SI_S(t), Ci(x) = gamma + log x + x^2 SI_C(t), t = x^2/8 - 1, 0 <= x <= 4 */")
    print(poly_in_t(S, 0, 16, "SI_S", list(range(9, 16))))
    print(poly_in_t(C, 0, 16, "SI_C", list(range(9, 16))))
    print("/* Shi(x) = x SH_S(v), Chi(x) = gamma + log x + x^2 SH_C(v), v = x^2/16, 0 <= x <= 4 */")
    print(poly_in_v(Sh, 0, 16, "SH_S", list(range(9, 18))))
    print(poly_in_v(Ch, 0, 16, "SH_C", list(range(9, 18))))

    # ---- Ei near the origin: Ei(x) = gamma + log x + x EI_R(t), x in [0, 2] ---------
    def R(x): return mp.mpf(1) if x == 0 else (mp.ei(x) - mp.euler - mp.log(abs(x))) / x
    print("/* Ei(x) = gamma + log x + x EI_Rk(v), v = (x - a)/(b - a), on [a, b] = "
          "[0,2] [2,4] [4,6] [6,8] [8,12] [12,16] */")
    for k, (a, b) in enumerate([(0, 2), (2, 4), (4, 6), (6, 8), (8, 12), (12, 16)]):
        print(poly_in_v(R, a, b, f"EI_R{k}", list(range(10, 26))))
    print("/* E1(x) = x E1_R(t) - gamma - log x, t = (2x - 0.6)/0.6, 0 < x <= 0.6;"
          "  E1_R(t) = (E1(x) + gamma + log x)/x, the Ei series at -x */")
    print(poly_in_t(lambda x: R(-x), 0, mp.mpf("0.6"), "E1_R", list(range(7, 16))))

    # ---- Si/Ci auxiliary f, g:  x f(x) = 1 - w RF(w),  x^2 g(x) = 1 - w RG(w), w = 1/x^2
    degs = [(m, m) for m in range(4, 11)]
    for (xa, xb, tag) in [(4, 8, "A"), (8, mp.inf, "B")]:
        wb = 1 / mp.mpf(xa) ** 2
        wa = mp.mpf(0) if xb == mp.inf else 1 / mp.mpf(xb) ** 2

        def RF(w):
            if w == 0:
                return mp.mpf(2)
            x = 1 / mp.sqrt(w)
            return (1 - fg(x)[0] * x) / w

        def RG(w):
            if w == 0:
                return mp.mpf(6)
            x = 1 / mp.sqrt(w)
            return (1 - fg(x)[1] * x * x) / w
        print(f"/* x f(x) = 1 - w SI_RF{tag}(w), x^2 g(x) = 1 - w SI_RG{tag}(w), w = 1/x^2, x in [{xa}, {xb}] */")
        print(rational_in_u(RF, wa, wb, "SI_RF" + tag, degs))
        print(rational_in_u(RG, wa, wb, "SI_RG" + tag, degs))

    # ---- E1: x e^x E1(x) = 1 - u RJ(u), u = 1/x ------------------------------------
    for (xa, xb, tag) in [(mp.mpf("0.6"), 2, "A"), (2, 8, "B"), (8, mp.inf, "C")]:
        ub = 1 / mp.mpf(xa)
        ua = mp.mpf(0) if xb == mp.inf else 1 / mp.mpf(xb)

        def RJ(u):
            if u == 0:
                return mp.mpf(1)
            x = 1 / u
            return (1 - x * mp.exp(x) * mp.e1(x)) / u
        print(f"/* x e^x E1(x) = 1 - u E1_J{tag}(u), u = 1/x, x in [{xa}, {xb}] */")
        print(rational_in_u(RJ, ua, ub, "E1_J" + tag, degs))

    # ---- Ei: x e^-x Ei(x) = 1 + u RH(u), u = 1/x -------------------------------------
    for (xa, xb, tag) in [(16, 32, "D"), (32, 64, "E")]:
        ub = 1 / mp.mpf(xa)
        ua = 1 / mp.mpf(xb)

        def RH(u):
            x = 1 / u
            return (x * mp.exp(-x) * mp.ei(x) - 1) / u
        print(f"/* x e^-x Ei(x) = 1 + u EI_H{tag}(t), u = 1/x, x in [{xa}, {xb}], t centred */")
        print(rational_in_t(RH, ua, ub, "EI_H" + tag, degs))

    print(series_tables())

    print("\n#endif /* MATHILDA_SF_MACHINE_EXPINT_COEF_H */")


if __name__ == "__main__":
    main()
