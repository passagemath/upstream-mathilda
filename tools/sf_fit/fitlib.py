"""Piecewise polynomial fitting for Mathilda's machine special-function kernels.

Each piece approximates g on [lo, hi] by a polynomial in the CENTRED variable
u = (t - mid)/half in [-1, 1], from interpolation at Chebyshev nodes
(mpmath.chebyfit), which is within a small factor of the minimax polynomial.
The coefficients of u^1..u^d are rounded to double and the constant term is
kept as a double-double (c0 = hi + lo), so the piece's own rounding floor is
far below an ulp; the C evaluation is  p = c0h + (c0l + u*(c1 + u*(...))).
The approximation error of the ROUNDED table is re-measured in high precision,
so the figure reported is what the C table actually carries (Horner rounding
is measured separately, end to end, by tools/sf_accuracy/error.py).

Run with the benchmark venv:  benchmarks/.venv/bin/python -I <script>.py
"""
import mpmath as mp

mp.mp.dps = 40
EPS = mp.mpf(2) ** -53   # unit roundoff; errors below are in units of it


def round_coeffs(poly):
    """poly: mp coefficients, highest degree first. Returns (hi-first list of
    doubles for degrees d..1, c0 hi, c0 lo)."""
    cs = [float(c) for c in poly[:-1]]
    c0 = poly[-1]
    c0h = float(c0)
    c0l = float(c0 - mp.mpf(c0h))
    return cs, c0h, c0l


def eval_rounded(cs, c0h, c0l, u):
    p = mp.mpf(0)
    for c in cs:
        p = p * u + mp.mpf(c)
    return mp.mpf(c0h) + (mp.mpf(c0l) + p * u)


def fit_piece(g, lo, hi, deg, rel=True, npts=60):
    lo = mp.mpf(lo); hi = mp.mpf(hi)
    mid = (lo + hi) / 2; half = (hi - lo) / 2
    gu = lambda u: g(mid + half * u)
    poly, _ = mp.chebyfit(gu, [-1, 1], deg + 1, error=True)
    cs, c0h, c0l = round_coeffs(poly)
    worst = mp.mpf(0)
    for i in range(npts + 1):
        u = mp.cos(mp.pi * (i + mp.mpf(0.5)) / (npts + 1)) if i < npts else mp.mpf(1)
        for uu in (u, -u):
            ref = gu(uu)
            e = abs(eval_rounded(cs, c0h, c0l, uu) - ref)
            if rel and ref != 0:
                e /= abs(ref)
            worst = max(worst, e)
    return (cs, c0h, c0l), float(mid), float(half), float(worst / EPS)


def fit_auto(g, lo, hi, target=0.15, rel=True, dmin=3, dmax=30, npts=60):
    """Smallest degree whose rounded-table error is <= target (units of 2^-53);
    if the target is out of reach, the degree after which the error stops
    improving."""
    best = None
    for d in range(dmin, dmax + 1):
        r = fit_piece(g, lo, hi, d, rel, npts)
        if best is None or r[3] < best[3] * 0.9:
            best = r + (d,)
        if r[3] <= target:
            return r + (d,)
        if d - best[4] >= 3:          # three degrees without a 10% gain
            return best
    return best


def c_table(name, pieces):
    """pieces: list of ((cs, c0h, c0l), ...) all of the same degree. Emits a
    [n][deg+2] array: c_d .. c_1, c0h, c0l."""
    deg = len(pieces[0][0])
    rows = []
    for (cs, c0h, c0l) in pieces:
        assert len(cs) == deg
        rows.append(list(cs) + [c0h, c0l])
    out = ["static const double %s[%d][%d] = {" % (name, len(rows), deg + 2)]
    for r in rows:
        out.append("    { " + ", ".join("%.17g" % c for c in r) + " },")
    out.append("};")
    return "\n".join(out)


def erfcinv_mp(q):
    """erfc^{-1}(q) in working precision: double seed from SciPy, then Newton
    on erfc(z) - q (mpmath's own erfinv runs a generic root finder and is far
    too slow to sample thousands of nodes)."""
    import scipy.special as sc
    q = mp.mpf(q)
    qf = float(q)
    if qf == 0.0 or qf < 1e-300:
        z = mp.sqrt(-mp.log(q))
    else:
        z = mp.mpf(float(sc.erfcinv(qf)))
    for _ in range(60):
        f = mp.erfc(z) - q
        dz = f / (-2 / mp.sqrt(mp.pi) * mp.exp(-z * z))
        z -= dz
        if abs(dz) < abs(z) * mp.mpf(10) ** (-mp.mp.dps + 3):
            break
    return z


def erfinv_mp(p):
    p = mp.mpf(p)
    if p == 0:
        return mp.mpf(0)
    if p < 0:
        return -erfinv_mp(-p)
    if p > 0.5:
        return erfcinv_mp(1 - p)
    import scipy.special as sc
    z = mp.mpf(float(sc.erfinv(float(p))))
    for _ in range(60):
        dz = (mp.erf(z) - p) / (2 / mp.sqrt(mp.pi) * mp.exp(-z * z))
        z -= dz
        if abs(dz) < abs(z) * mp.mpf(10) ** (-mp.mp.dps + 3):
            break
    return z
