#!/usr/bin/env python3
"""Fit the machine approximations for the complete elliptic integrals.

The form is Cephes' (S. L. Moshier, ellpk.c / ellpe.c), which in turn is the
structure of Abramowitz & Stegun 17.3.34 / 17.3.36 (Hastings):

    K(m) = P(t) - log(t) Q(t)          t = 1 - m in (0, 1]
    E(m) = P(t) - t log(t) Q(t)

P and Q are polynomials.  The pair is linear in its coefficients, so the
relative-error minimax fit is a LINEAR minimax problem over the basis
{t^k} u {-log(t) t^k}; it is solved here by Lawson's iteratively-reweighted
least squares (Lawson 1961), in mpmath at 50 digits, on a point set that is
dense near both ends: log-spaced down to t = 1e-300 so that Q(0) is pinned (the
log term reaches |log t| ~ 690 there, and any error in Q(0) is multiplied by
it), and Chebyshev-clustered near t = 1 (m = 0, where the log singularity of
the individual pieces sits).

The coefficients are NOT recalled from Cephes: they are re-derived by this
script, and printed as C arrays (highest degree first, for Horner).

Usage:
  benchmarks/.venv/bin/python -I tools/sf_accuracy/fit_ellk.py [--deg 10] [--iters 60]
"""
from __future__ import annotations

import argparse
import math

import mpmath as mp

mp.mp.dps = 50


def points(n):
    pts = set()
    for i in range(n):                       # Chebyshev-clustered on (0, 1]
        x = (1 - mp.cos(mp.pi * (i + 0.5) / n)) / 2
        pts.add(x)
    for i in range(n // 2):                  # log-spaced toward 0
        e = -300 * (i / (n // 2 - 1)) ** 2
        pts.add(mp.mpf(10) ** e)
    pts.add(mp.mpf(1))
    return sorted(pts)


def target(kind, t):
    # m = 1 - t must be formed exactly: t reaches 1e-300.
    with mp.workdps(400):
        m = 1 - mp.mpf(t)
        v = mp.ellipk(m) if kind == "K" else mp.ellipe(m)
    return +v


def basis(kind, t, deg):
    lt = mp.log(t)
    row = [t ** k for k in range(deg + 1)]
    if kind == "K":
        row += [-lt * t ** k for k in range(deg + 1)]
    else:
        row += [-t * lt * t ** k for k in range(deg + 1)]
    return row


def fit(kind, deg, iters, npts):
    ts = points(npts)
    F = [target(kind, t) for t in ts]
    A = [basis(kind, t, deg) for t in ts]
    nc = len(A[0])
    w = [mp.mpf(1) / len(ts)] * len(ts)
    coef = None
    for it in range(iters):
        # weighted LS on relative error: minimise sum w_i ((A c - F)/F)^2
        M = mp.zeros(nc, nc)
        b = mp.zeros(nc, 1)
        for Ai, Fi, wi in zip(A, F, w):
            s = wi / (Fi * Fi)
            for j in range(nc):
                aj = Ai[j] * s
                b[j] += aj * Fi
                for k in range(j, nc):
                    M[j, k] += aj * Ai[k]
        for j in range(nc):
            for k in range(j):
                M[j, k] = M[k, j]
        coef = mp.lu_solve(M, b)
        errs = [abs(sum(Ai[j] * coef[j] for j in range(nc)) - Fi) / Fi for Ai, Fi in zip(A, F)]
        emax = max(errs)
        # Lawson update
        w = [wi * e for wi, e in zip(w, errs)]
        sw = sum(w)
        w = [wi / sw for wi in w]
        if it % 10 == 0 or it == iters - 1:
            print(f"  {kind} iter {it}: max rel err {mp.nstr(emax, 5)}")
    return [coef[j] for j in range(nc)]


def check_double(kind, coef, deg):
    P = [float(c) for c in coef[:deg + 1]]
    Q = [float(c) for c in coef[deg + 1:]]

    def horner(c, t):
        r = 0.0
        for a in reversed(c):
            r = r * t + a
        return r

    worst = 0.0
    import random
    rng = random.Random(1)
    for i in range(4000):
        if i % 2:
            t = rng.uniform(0, 1) or 0.5
        else:
            t = 10.0 ** rng.uniform(-300, 0)
        lt = math.log(t)
        if kind == "K":
            v = horner(P, t) - lt * horner(Q, t)
        else:
            v = horner(P, t) - t * lt * horner(Q, t)
        ref = target(kind, mp.mpf(t))
        e = abs(mp.mpf(v) - ref) / ref / mp.mpf(2) ** -53
        worst = max(worst, float(e))
    print(f"  {kind}: max error in double evaluation {worst:.2f} ulp(ish, units of 2^-53 rel)")


def emit(kind, coef, deg):
    for name, part in (("P", coef[:deg + 1]), ("Q", coef[deg + 1:])):
        print(f"static const double ELLK_{kind}{name}[{deg + 1}] = {{   /* highest degree first */")
        for c in reversed(part):
            print(f"    {float(c)!r},")
        print("};")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg", type=int, default=10)
    ap.add_argument("--iters", type=int, default=60)
    ap.add_argument("--npts", type=int, default=400)
    ap.add_argument("--kinds", default="K,E")
    a = ap.parse_args()
    for kind in a.kinds.split(","):
        c = fit(kind, a.deg, a.iters, a.npts)
        check_double(kind, c, a.deg)
        emit(kind, c, a.deg)


if __name__ == "__main__":
    main()
