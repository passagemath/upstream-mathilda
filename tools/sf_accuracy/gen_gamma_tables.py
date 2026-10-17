#!/usr/bin/env python3
"""Generate the coefficient tables of src/special_functions/sf_machine_gamma.c.

    benchmarks/.venv/bin/python -I tools/sf_accuracy/gen_gamma_tables.py > tables.inc

Prints the block between the "Coefficient tables" banner and the complex helpers
of sf_machine_gamma.c. Everything is computed with mpmath at 60-120 digits:

  SFG_DIG12        degree-22 Chebyshev interpolant (nodes on [1, 2]) of
                   psi(x) / (x - x0), as monomial coefficients in s = x - 3/2
  SFG_B2K_FACT     B_2k / (2k)!                     (Hurwitz-zeta / polygamma tail)
  SFG_DIGAMMA_ASY  B_2k / (2k)                      (digamma asymptotic series)
  SFG_STIRLING     B_2k / (2k (2k-1))               (Stirling's series)
  SFG_FACT         k!, k <= 22
  SFG_TILE_*       Taylor series of lnGamma and psi about 30 centres
                   {1.25, 1.75, 2.25} + {0, 1/2, ..., 9/2} i, with per-tile term
                   counts for 2^-57 over the tile |u| <= sqrt(2)/4
  SFG_PSI_ROOT     Taylor series of psi about its positive zero x0 = 1.4616...
"""
import mpmath as mp

def arr(name, vals, comment):
    s = "/* %s */\nstatic const double %s[%d] = {\n" % (comment, name, len(vals))
    line = "   "
    for v in vals:
        t = " %.17e," % float(v)
        if len(line) + len(t) > 80:
            s += line + "\n"; line = "   "
        line += t
    return s + line + "\n};\n"

mp.mp.dps = 60
x0 = mp.findroot(mp.digamma, 1.46)
deg = 22; N = deg + 1
nodes = [mp.cos(mp.pi*(k + mp.mpf(1)/2)/N) for k in range(N)]
xs = [mp.mpf(1.5) + n/2 for n in nodes]
ys = [mp.digamma(x)/(x - x0) for x in xs]
A = mp.matrix([[(x - mp.mpf(1.5))**j for j in range(N)] for x in xs])
dig = mp.lu_solve(A, mp.matrix(ys))
a = x0 - 1
out = "/* x0 - 1, x0 = 1.46163... the positive zero of psi, as a double-double */\n"
out += "#define SFG_X0M1_HI %.17e\n#define SFG_X0M1_LO %.17e\n" % (float(a), float(a - mp.mpf(float(a))))
out += "#define SFG_LN_SQRT_2PI %.17e\n" % float(mp.log(2*mp.pi)/2)
out += arr("SFG_DIG12", [dig[j] for j in range(N)] + [0],
           "psi(x) / (x - x0) on [1,2] as a polynomial in s = x - 3/2 (degree-22 Chebyshev interpolant; zero-padded to 24)")
out += arr("SFG_B2K_FACT", [mp.bernoulli(2*k)/mp.factorial(2*k) for k in range(1, 25)], "B_2k / (2k)!, k = 1..24")
out += arr("SFG_DIGAMMA_ASY", [mp.bernoulli(2*k)/(2*k) for k in range(1, 13)], "B_2k / (2k), k = 1..12")
out += arr("SFG_STIRLING", [mp.bernoulli(2*k)/(2*k*(2*k-1)) for k in range(1, 14)], "B_2k / (2k (2k-1)), k = 1..13")
out += arr("SFG_FACT", [mp.factorial(k) for k in range(23)], "k!, k = 0..22 (exact in double)")
out += "/* x0 = 1.4616..., the positive zero of psi, as a double-double */\n"
out += "#define SFG_X0_HI %.17e\n#define SFG_X0_LO %.17e\n" % (float(x0), float(x0 - mp.mpf(float(x0))))

mp.mp.dps = 120
KS = [8, 12, 16, 20, 24, 28, 32, 36]
IMS = [mp.mpf(b)/2 for b in range(10)]
RE = [mp.mpf('1.25'), mp.mpf('1.75'), mp.mpf('2.25')]   # the strip 1 <= Re < 5/2
RE_LG = RE; RE_PS = RE
NC = 40   # coefficients stored (K <= 36, plus the a_0-split offset)
EXTRA = 16
def lg_coefs(c, n):
    return [mp.loggamma(c), mp.digamma(c)] + [mp.polygamma(k-1, c)/mp.factorial(k) for k in range(2, n)]
def psi_coefs(c, n):
    return [mp.polygamma(k, c)/mp.factorial(k) for k in range(0, n)]
def kfix(co, scale):
    for K in KS:
        r = mp.sqrt(mp.mpf(2))/4
        tail = sum(abs(co[k])*r**k for k in range(K, len(co)))
        if tail <= mp.mpf(2)**-57 * scale: return K
    raise SystemExit("tile needs more than %d terms" % KS[-1])
lgt=[]; pst=[]; lgk=[]; psk=[]; a0lo=[]
for b in IMS:
    for a in RE_LG:
        c = mp.mpc(a, b)
        A = lg_coefs(c, NC + EXTRA)
        if b == 0:
            # lnGamma(c+u) = (c+u-1)(c+u-2) h(u): both real zeros carried exactly
            q0 = (c.real - 1)*(c.real - 2); q1 = 2*c.real - 3
            H = []
            for k in range(len(A)):
                v = A[k] - (q1*H[k-1] if k >= 1 else 0) - (H[k-2] if k >= 2 else 0)
                H.append(v/q0)
            A = H
            scale = abs(A[0]) - mp.sqrt(2)/4*abs(A[1])
            lo = 0
        else:
            scale = max(abs(A[0]) - mp.sqrt(2)/4*abs(A[1]), mp.mpf('0.05'))
            lo = A[0] - mp.mpc(float(A[0].real), float(A[0].imag))
        lgk.append(kfix(A, scale)); lgt.append(A[:NC]); a0lo.append(lo)
    for a in RE_PS:
        c = mp.mpc(a, b)
        D = psi_coefs(c, NC + EXTRA)
        # psi's zero x0 = 1.4616 lies in the strip, but |z - x0| <= 0.35 is taken by
        # the root series, and outside that disk |psi| >= 0.3 on the tiles
        scale = max(abs(D[0]) - mp.sqrt(2)/4*abs(D[1]), mp.mpf('0.25'))
        psk.append(kfix(D, scale)); pst.append(D[:NC])
def emit3(name, tabs, res):
    s = "static const double %s[30][%d][2] = {\n" % (name, NC)
    for i, t in enumerate(tabs):
        s += "  { /* centre %s + %si */\n" % (mp.nstr(res[i % 3], 3), mp.nstr(IMS[i // 3], 2))
        for v in t:
            v = mp.mpc(v); s += "    { %.17e, %.17e },\n" % (float(v.real), float(v.imag))
        s += "  },\n"
    return s + "};\n"
tout  = "/* Taylor tiles (tools/sf_accuracy/gen_gamma_tables.py): centres c = {1.25, 1.75, 2.25}\n"
tout += " * + {0, 1/2, ..., 9/2} i covering the strip 1 <= Re < 5/2, tile index 3*(2 Im c) + j.\n"
tout += " *   lnGamma: lnGamma(c+u) = sum_k a_k u^k, a_0 = lnGamma(c) as a double-double (low\n"
tout += " *            part in SFG_TILE_LG_A0LO), except on the real axis (Im c = 0), where\n"
tout += " *            the table holds h_k with lnGamma(c+u) = (c+u-1)(c+u-2) sum_k h_k u^k,\n"
tout += " *            so the zeros at z = 1 and 2 are carried exactly.\n"
tout += " *   psi:     psi(c+u) = sum_k d_k u^k. */\n"
tout += emit3("SFG_TILE_LG", lgt, RE_LG) + emit3("SFG_TILE_PS", pst, RE_PS)
tout += "static const double SFG_TILE_LG_A0LO[30][2] = {\n" + "".join("    { %.17e, %.17e },\n" % (float(mp.mpc(v).real), float(mp.mpc(v).imag)) for v in a0lo) + "};\n"
tout += "/* Terms needed over the whole tile (|u| <= sqrt(2)/4) for 2^-57, per centre. */\n"
tout += "static const unsigned char SFG_TILE_LG_K[30] = { %s };\n" % ", ".join(map(str, lgk))
tout += "static const unsigned char SFG_TILE_PS_K[30] = { %s };\n" % ", ".join(map(str, psk))
TILES = tout

mp.mp.dps = 60
x0 = mp.findroot(mp.digamma, 1.46)
e = [mp.polygamma(k, x0)/mp.factorial(k) for k in range(1, 31)]
out += TILES
out += arr("SFG_PSI_ROOT", e, "psi(x0 + u) = sum_{k>=1} e_k u^k, x0 = 1.4616...")
print(out)
