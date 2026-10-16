"""Small rational/polynomial fitting library for the machine special-function kernels.

The coefficients in src/special_functions/sf_machine_expint.c are produced by
tools/sf_accuracy/fit_expint.py with these routines and mpmath at 50 digits, so
they are reproducible from first principles (no vendored tables).

Method: near-minimax RELATIVE-error rational approximation P(s)/Q(s) on a set of
Chebyshev-distributed nodes, by the linearised Sanathanan-Koerner iteration
(Sanathanan & Koerner, IEEE Trans. Autom. Control 8 (1963) 56-58) followed by
Lawson reweighting (Lawson 1961; see Trefethen, Approximation Theory and
Approximation Practice, ch. 24/26 for the AAA/Lawson family):

    minimise  sum_i  w_i | (f_i Q(s_i) - P(s_i)) / (f_i Q_prev(s_i)) |^2 ,  q_0 = 1

A polynomial (n = 0) fit is the same problem with Q == 1.
"""
from __future__ import annotations

import math

import mpmath as mp

mp.mp.dps = 50


def cheb_nodes(N):
    return [mp.cos(mp.pi * (2 * i + 1) / (2 * N)) for i in range(N)]


def horner(c, s):
    r = mp.mpf(0)
    for a in reversed(c):
        r = r * s + a
    return r


def fit_rational(f, m, n, N=None, iters=12, lawson=30, nodes=None, a=-1, b=1):
    """f: callable on s in [a, b] (mp numbers) -> mp value; the fit variable IS s
    (choose [a, b] so the monomial basis is well conditioned there, e.g. [0, 1]
    for a Stieltjes-type function, whose P and Q then have one sign).
    Returns (p, q) coefficient lists (low order first), q[0] = 1."""
    if N is None:
        N = max(4 * (m + n + 2), 60)
    S = nodes if nodes is not None else [mp.mpf(a) + (mp.mpf(b) - a) * (c + 1) / 2
                                         for c in cheb_nodes(N)]
    F = [f(s) for s in S]
    w = [mp.mpf(1)] * len(S)
    qprev = [mp.mpf(1)] * len(S)
    best = None
    for it in range(iters + lawson):
        rows, rhs = [], []
        for s, fv, wt, qp in zip(S, F, w, qprev):
            sc = mp.sqrt(wt) / (fv * qp)
            row = [-(s ** k) * sc for k in range(m + 1)] + [fv * s ** k * sc for k in range(1, n + 1)]
            rows.append(row)
            rhs.append(-fv * sc)
        A = mp.matrix(rows)
        b = mp.matrix(rhs)
        x, _ = mp.qr_solve(A, b)
        p = [x[k] for k in range(m + 1)]
        q = [mp.mpf(1)] + [x[m + 1 + k] for k in range(n)]
        errs = []
        for s, fv in zip(S, F):
            qs = horner(q, s)
            errs.append(abs((horner(p, s) / qs - fv) / fv))
        emax = max(errs)
        if best is None or emax < best[0]:
            best = (emax, p, q)
        if it < iters:
            qprev = [horner(q, s) for s in S]
        else:
            # Lawson: emphasise the nodes where the error is largest
            tot = mp.mpf(0)
            for i, e in enumerate(errs):
                w[i] = w[i] * (e + mp.mpf(10) ** -40)
                tot += w[i]
            w = [wi / tot * len(w) for wi in w]
            qprev = [horner(q, s) for s in S]
    return best[1], best[2]


def max_rel_error(f, p, q, M=2000, a=-1, b=1, as_double=False):
    """Max relative error of P/Q against f on M points in [a, b].
    as_double: round the coefficients to double and evaluate in double Horner."""
    worst = 0.0
    wx = None
    pd = [float(c) for c in p]
    qd = [float(c) for c in q]
    for i in range(M + 1):
        s = mp.mpf(a) + (mp.mpf(b) - a) * i / M
        if as_double:
            s = mp.mpf(float(s))     # judge the evaluation, not the input rounding
        fv = f(s)
        if as_double:
            sd = float(s)
            num = 0.0
            for c in reversed(pd):
                num = num * sd + c
            den = 0.0
            for c in reversed(qd):
                den = den * sd + c
            v = mp.mpf(num / den)
        else:
            v = horner(p, s) / horner(q, s)
        e = float(abs((v - fv) / fv))
        if e > worst:
            worst, wx = e, s
    return worst, wx


def c_array(name, coeffs, comment=""):
    body = ",\n".join("    " + mp.nstr(c, 20, min_fixed=-1, max_fixed=-1, strip_zeros=False)
                      .replace("e", "e") for c in coeffs)
    hdr = f"/* {comment} */\n" if comment else ""
    return f"{hdr}static const double {name}[{len(coeffs)}] = {{\n{body}\n}};\n"


def dbl(c):
    """Repr of the nearest double, round-trip exact."""
    return repr(float(c))
