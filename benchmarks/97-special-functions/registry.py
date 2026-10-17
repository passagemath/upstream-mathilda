"""registry.py -- the function<->engine spine of the special-functions shootout.

One record per special-function case.  A record carries, for a given Mathilda head:

  * the WL/Mathilda call template (shared verbatim by Mathilda and wolframscript,
    which speak the same surface syntax -- this is exactly why run_all.py feeds one
    .m to both);
  * the Maxima call template (or None where Maxima has no numeric counterpart);
  * a Python callable (scipy.special) for the machine-real / complex regimes, and an
    mpmath callable for arbitrary precision (and as the fallback when scipy declines
    or lacks the function);
  * the regimes the case is expected to support, and the sampling domain per slot.

Templates use ``{0}``, ``{1}``, ... placeholders for the *varying* argument slots;
fixed parameters (an integer order, the list indices of a pFq) are baked into the
template string.  The orchestrator (run_shootout.py) renders the placeholders into
each engine's/regime's numeric literal.

Conventions baked in here, with the inventory's gotchas:
  * All Elliptic heads take the PARAMETER m = k^2 (not the modulus) -- scipy ellipk(m),
    mpmath ellipk(m), Maxima elliptic_kc(m) all agree on m.
  * Mathematica Sinc[z] = sin z / z (NOT the normalised sinc); rendered explicitly so
    no engine's normalisation convention leaks in.
  * scipy.special.fresnel returns (S, C); sici returns (Si, Ci); shichi returns (Shi, Chi).
  * Pochhammer -> rf (rising factorial); InverseErf -> erfinv; ProductLog -> lambertw.

Domains are chosen to sit inside each head's real principal branch so timing measures
the kernel, not a decline-to-fallback path (same care as benchmarks/11-special-functions).
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass, field
from typing import Callable, Optional

import mpmath as mp
import scipy.special as sc

# Regime tags
MR, C, AP = "MR", "C", "AP"
ALL = (MR, C, AP)


@dataclass
class SF:
    head: str                      # Mathilda / WL head (also the display name)
    family: str
    nslots: int                    # number of varying argument slots
    wl: str                        # WL/Mathilda call template (shared with wolframscript)
    mx: Optional[str]              # Maxima call template, or None if Maxima lacks it
    mpf: Callable                  # mpmath callable(*vals) -> mpf/mpc (AP + fallback)
    scf: Optional[Callable] = None # scipy callable(*vals) for MR/C; None => use mpf
    sc_complex: bool = True        # does scf accept complex args?
    regimes: tuple = ALL           # which of (MR, C, AP) this case supports
    # per-regime slot domains.  MR/AP slot = (lo, hi) real interval.
    # C slot = (re_lo, re_hi, im_lo, im_hi).
    dom_r: tuple = ((0.3, 6.0),)   # real interval(s), one per slot (MR and AP)
    dom_c: tuple = ((0.5, 4.0, -3.0, 3.0),)  # complex box(es), one per slot
    note: str = ""
    pylib: str = "scipy"           # which library backs the Python MR/C column ("scipy" or "mpmath")

    @property
    def label(self) -> str:
        return self.head


def _airy(idx):
    return lambda z: sc.airy(z)[idx]


def _mp_airyprime(which):
    # which: 'ai' or 'bi'
    return (lambda z: mp.airyai(z, 1)) if which == "ai" else (lambda z: mp.airybi(z, 1))


# mpmath lacks erfcinv: erfc^{-1}(x) = erf^{-1}(1 - x)
def _mp_erfcinv(x):
    return mp.erfinv(1 - x)


def _sinc_mp(z):
    z = mp.mpf(z) if not isinstance(z, (mp.mpf, mp.mpc)) else z
    return mp.sin(z) / z


def _sinc_sc(z):
    return cmath.sin(z) / z if isinstance(z, complex) else math.sin(z) / z


# =====================================================================================
# NUMERIC heads.  Grouped by family for the write-up.
# =====================================================================================
NUMERIC: list[SF] = [
    # ---- Gamma family ---------------------------------------------------------------
    SF("Gamma", "Gamma", 1, "Gamma[{0}]", "gamma({0})",
       mp.gamma, sc.gamma, True, ALL, ((0.3, 6.0),), ((0.5, 4.0, -3.0, 3.0),)),
    SF("LogGamma", "Gamma", 1, "LogGamma[{0}]", "log_gamma({0})",
       mp.loggamma, sc.loggamma, True, ALL, ((0.3, 8.0),), ((0.5, 4.0, -3.0, 3.0),)),
    SF("PolyGamma0", "Gamma", 1, "PolyGamma[0, {0}]", "psi[0]({0})",
       mp.digamma, sc.digamma, True, ALL, ((0.3, 8.0),), ((0.5, 4.0, -3.0, 3.0),),
       note="digamma (PolyGamma[0, z])"),
    SF("PolyGamma2", "Gamma", 1, "PolyGamma[2, {0}]", "psi[2]({0})",
       (lambda z: mp.polygamma(2, z)), (lambda x: sc.polygamma(2, x)), False, ALL,
       ((0.3, 8.0),), ((0.5, 4.0, -3.0, 3.0),), note="scipy polygamma is real-only"),
    SF("Beta", "Gamma", 2, "Beta[{0}, {1}]", "beta({0}, {1})",
       mp.beta, sc.beta, False, ALL, ((0.3, 5.0), (0.3, 5.0)),
       ((0.5, 3.0, -2.0, 2.0), (0.5, 3.0, -2.0, 2.0)), note="scipy beta real-only"),
    SF("Pochhammer", "Gamma", 1, "Pochhammer[{0}, 3]", "pochhammer({0}, 3)",
       (lambda a: mp.rf(a, 3)), (lambda a: sc.poch(a, 3)), False, ALL,
       ((0.3, 5.0),), ((0.5, 3.0, -2.0, 2.0),), note="rising factorial, n = 3; scipy poch real-only"),

    # ---- Zeta family ----------------------------------------------------------------
    SF("Zeta", "Zeta", 1, "Zeta[{0}]", "zeta({0})",
       mp.zeta, (lambda s: sc.zeta(s, 1)), False, ALL, ((1.5, 9.0),),
       ((1.5, 5.0, -3.0, 3.0),), note="s > 1; scipy zeta real-only"),
    SF("HurwitzZeta", "Zeta", 2, "HurwitzZeta[{0}, {1}]", None,
       mp.zeta, (lambda s, a: sc.zeta(s, a)), False, ALL,
       ((1.5, 6.0), (0.5, 4.0)), ((1.5, 4.0, -2.0, 2.0), (0.5, 3.0, -1.0, 1.0)),
       note="Maxima leaves hurwitz_zeta symbolic"),
    SF("LerchPhi", "Zeta", 3, "LerchPhi[{0}, {1}, {2}]", None,
       mp.lerchphi, None, False, ALL,
       ((0.0, 0.9), (2.0, 4.0), (0.5, 2.0)),
       ((-0.6, 0.6, -0.6, 0.6), (2.0, 4.0, 0.0, 0.0), (0.5, 2.0, 0.0, 0.0)),
       note="Python column is mpmath (scipy has no LerchPhi)", pylib="mpmath"),
    # DirichletL: mpmath.dirichlet takes the character as a periodic value list,
    # built here in Mathematica's DirichletCharacter labelling (the values were
    # checked against Mathematica 14) at the working precision.
    SF("DirichletL", "Zeta", 1, "DirichletL[5, 3, {0}]", None,
       (lambda s: mp.dirichlet(s, [0, 1, -1, -1, 1])), None, False, ALL,
       ((-4.0, 6.0),), ((0.2, 4.0, -10.0, 10.0),),
       note="real character mod 5 (k = 5, j = 3); Python is mpmath; scipy/Maxima lack it",
       pylib="mpmath"),
    SF("DirichletLc", "Zeta", 1, "DirichletL[7, 2, {0}]", None,
       (lambda s: mp.dirichlet(s, [0, 1, mp.expjpi(mp.mpf(2) / 3), mp.expjpi(mp.mpf(1) / 3),
                                   mp.expjpi(mp.mpf(-2) / 3), mp.expjpi(mp.mpf(-1) / 3), -1])),
       None, False, ALL,
       ((-4.0, 6.0),), ((0.2, 4.0, -10.0, 10.0),),
       note="complex character mod 7 (k = 7, j = 2, order 6); Python is mpmath",
       pylib="mpmath"),
    SF("PolyLog", "Zeta", 1, "PolyLog[2, {0}]", None,
       (lambda z: mp.polylog(2, z)), None, False, ALL,
       ((-0.9, 0.9),), ((-0.8, 0.8, -0.8, 0.8),),
       note="dilog (s = 2); Python is mpmath; Maxima symbolic", pylib="mpmath"),

    # ---- Bessel / Airy --------------------------------------------------------------
    SF("BesselJ", "Bessel/Airy", 1, "BesselJ[2, {0}]", "bessel_j(2, {0})",
       (lambda z: mp.besselj(2, z)), (lambda z: sc.jv(2, z)), True, ALL,
       ((0.5, 15.0),), ((0.5, 8.0, -4.0, 4.0),), note="order n = 2"),
    SF("BesselY", "Bessel/Airy", 1, "BesselY[1, {0}]", "bessel_y(1, {0})",
       (lambda z: mp.bessely(1, z)), (lambda z: sc.yv(1, z)), True, ALL,
       ((0.5, 15.0),), ((0.5, 8.0, 0.1, 4.0),), note="order n = 1"),
    SF("BesselI", "Bessel/Airy", 1, "BesselI[0, {0}]", "bessel_i(0, {0})",
       (lambda z: mp.besseli(0, z)), (lambda z: sc.iv(0, z)), True, ALL,
       ((0.1, 8.0),), ((0.5, 4.0, -3.0, 3.0),), note="order n = 0"),
    SF("BesselK", "Bessel/Airy", 1, "BesselK[0, {0}]", "bessel_k(0, {0})",
       (lambda z: mp.besselk(0, z)), (lambda z: sc.kv(0, z)), True, ALL,
       ((0.1, 8.0),), ((0.5, 4.0, 0.1, 3.0),), note="order n = 0"),
    SF("AiryAi", "Bessel/Airy", 1, "AiryAi[{0}]", "airy_ai({0})",
       mp.airyai, _airy(0), True, ALL, ((-5.0, 5.0),), ((-3.0, 3.0, -3.0, 3.0),)),
    SF("AiryBi", "Bessel/Airy", 1, "AiryBi[{0}]", "airy_bi({0})",
       mp.airybi, _airy(2), True, ALL, ((-5.0, 4.0),), ((-3.0, 2.0, -3.0, 3.0),)),
    SF("AiryAiPrime", "Bessel/Airy", 1, "AiryAiPrime[{0}]", "airy_dai({0})",
       _mp_airyprime("ai"), _airy(1), True, ALL, ((-5.0, 5.0),), ((-3.0, 3.0, -3.0, 3.0),)),
    SF("AiryBiPrime", "Bessel/Airy", 1, "AiryBiPrime[{0}]", "airy_dbi({0})",
       _mp_airyprime("bi"), _airy(3), True, ALL, ((-5.0, 4.0),), ((-3.0, 2.0, -3.0, 3.0),)),
    SF("StruveH", "Bessel/Airy", 1, "StruveH[0, {0}]", "struve_h(0, {0})",
       (lambda z: mp.struveh(0, z)), (lambda z: sc.struve(0, z)), False, ALL,
       ((0.5, 10.0),), ((0.5, 4.0, -2.0, 2.0),),
       note="order 0; Mathilda forces the Arb path for complex (no machine-complex kernel)"),

    # ---- Error functions ------------------------------------------------------------
    SF("Erf", "Error", 1, "Erf[{0}]", "erf({0})",
       mp.erf, sc.erf, True, ALL, ((-3.0, 3.0),), ((-2.0, 2.0, -2.0, 2.0),)),
    SF("Erfc", "Error", 1, "Erfc[{0}]", "erfc({0})",
       mp.erfc, sc.erfc, True, ALL, ((-3.0, 3.0),), ((-2.0, 2.0, -2.0, 2.0),)),
    SF("Erfi", "Error", 1, "Erfi[{0}]", "erfi({0})",
       mp.erfi, sc.erfi, True, ALL, ((-3.0, 3.0),), ((-2.0, 2.0, -2.0, 2.0),)),
    SF("InverseErf", "Error", 1, "InverseErf[{0}]", "inverse_erf({0})",
       mp.erfinv, sc.erfinv, False, (MR, AP), ((-0.95, 0.95),), (),
       note="real-only in Mathilda; domain (-1, 1)"),
    SF("InverseErfc", "Error", 1, "InverseErfc[{0}]", "inverse_erfc({0})",
       _mp_erfcinv, sc.erfcinv, False, (MR, AP), ((0.05, 1.95),), (),
       note="real-only; domain (0, 2); mpmath via erfinv(1-x)"),
    SF("FresnelC", "Error", 1, "FresnelC[{0}]", "fresnel_c({0})",
       mp.fresnelc, (lambda z: sc.fresnel(z)[1]), True, ALL,
       ((-3.0, 3.0),), ((-2.0, 2.0, -2.0, 2.0),), note="scipy fresnel returns (S, C)"),
    SF("FresnelS", "Error", 1, "FresnelS[{0}]", "fresnel_s({0})",
       mp.fresnels, (lambda z: sc.fresnel(z)[0]), True, ALL,
       ((-3.0, 3.0),), ((-2.0, 2.0, -2.0, 2.0),)),

    # ---- Elliptic integrals (parameter m = k^2) -------------------------------------
    SF("EllipticK", "Elliptic", 1, "EllipticK[{0}]", "elliptic_kc({0})",
       mp.ellipk, sc.ellipk, False, ALL, ((0.0, 0.95),), ((-0.9, 0.9, -0.9, 0.9),),
       note="parameter m; scipy ellipk real-only"),
    SF("EllipticE", "Elliptic", 1, "EllipticE[{0}]", "elliptic_ec({0})",
       mp.ellipe, sc.ellipe, False, ALL, ((0.0, 0.95),), ((-0.9, 0.9, -0.9, 0.9),),
       note="complete, parameter m"),
    SF("EllipticF", "Elliptic", 2, "EllipticF[{0}, {1}]", "elliptic_f({0}, {1})",
       mp.ellipf, (lambda phi, m: sc.ellipkinc(phi, m)), False, (MR, AP),
       ((0.0, 1.5), (0.0, 0.95)), (), note="incomplete F(phi | m)"),
    SF("EllipticEinc", "Elliptic", 2, "EllipticE[{0}, {1}]", "elliptic_e({0}, {1})",
       (lambda phi, m: mp.ellipe(phi, m)), (lambda phi, m: sc.ellipeinc(phi, m)), False, (MR, AP),
       ((0.0, 1.5), (0.0, 0.95)), (), note="incomplete E(phi | m)"),
    SF("EllipticPi", "Elliptic", 2, "EllipticPi[{0}, {1}]", "elliptic_pi({0}, %pi/2, {1})",
       mp.ellippi, None, False, (MR, AP), ((0.0, 0.8), (0.0, 0.9)), (),
       note="complete Pi(n | m); Python is mpmath (scipy has none)", pylib="mpmath"),

    # ---- Orthogonal polynomials (integer degree) ------------------------------------
    SF("LegendreP", "Orthogonal", 1, "LegendreP[3, {0}]", "legendre_p(3, {0})",
       (lambda x: mp.legendre(3, x)), (lambda x: sc.eval_legendre(3, x)), True, ALL,
       ((-0.95, 0.95),), ((-0.9, 0.9, -0.9, 0.9),), note="degree 3; Maxima float wraps in interval()"),
    SF("LegendreQ", "Orthogonal", 1, "LegendreQ[2, {0}]", "legendre_q(2, {0})",
       (lambda x: mp.legenq(2, 0, x)), None, False, ALL,
       ((-0.9, 0.9),), ((-0.8, 0.8, -0.8, 0.8),),
       note="degree 2; Python is mpmath (scipy has no Legendre Q)", pylib="mpmath"),
    SF("HermiteH", "Orthogonal", 1, "HermiteH[4, {0}]", "hermite(4, {0})",
       (lambda x: mp.hermite(4, x)), (lambda x: sc.eval_hermite(4, x)), True, ALL,
       ((-3.0, 3.0),), ((-2.0, 2.0, -2.0, 2.0),), note="physicists' H, degree 4"),
    SF("LaguerreL", "Orthogonal", 1, "LaguerreL[3, {0}]", "laguerre(3, {0})",
       (lambda x: mp.laguerre(3, 0, x)), (lambda x: sc.eval_laguerre(3, x)), True, ALL,
       ((-2.0, 6.0),), ((-2.0, 4.0, -2.0, 2.0),), note="degree 3"),
    SF("ChebyshevT", "Orthogonal", 1, "ChebyshevT[5, {0}]", "chebyshev_t(5, {0})",
       (lambda x: mp.chebyt(5, x)), (lambda x: sc.eval_chebyt(5, x)), True, ALL,
       ((-1.0, 1.0),), ((-1.0, 1.0, -1.0, 1.0),), note="degree 5"),
    SF("ChebyshevU", "Orthogonal", 1, "ChebyshevU[5, {0}]", "chebyshev_u(5, {0})",
       (lambda x: mp.chebyu(5, x)), (lambda x: sc.eval_chebyu(5, x)), True, ALL,
       ((-0.99, 0.99),), ((-0.99, 0.99, -0.99, 0.99),), note="degree 5"),

    # ---- Integral functions ---------------------------------------------------------
    SF("SinIntegral", "Integral", 1, "SinIntegral[{0}]", "expintegral_si({0})",
       mp.si, (lambda z: sc.sici(z)[0]), False, ALL,
       ((0.1, 15.0),), ((0.5, 6.0, -3.0, 3.0),), note="scipy sici returns (Si, Ci)"),
    SF("CosIntegral", "Integral", 1, "CosIntegral[{0}]", "expintegral_ci({0})",
       mp.ci, (lambda z: sc.sici(z)[1]), False, ALL,
       ((0.1, 15.0),), ((0.5, 6.0, 0.1, 3.0),), note="x > 0 on the real axis"),
    SF("SinhIntegral", "Integral", 1, "SinhIntegral[{0}]", "expintegral_shi({0})",
       mp.shi, (lambda z: sc.shichi(z)[0]), False, ALL,
       ((0.1, 8.0),), ((0.5, 4.0, -2.0, 2.0),), note="scipy shichi returns (Shi, Chi)"),
    SF("CoshIntegral", "Integral", 1, "CoshIntegral[{0}]", "expintegral_chi({0})",
       mp.chi, (lambda z: sc.shichi(z)[1]), False, ALL,
       ((0.1, 8.0),), ((0.5, 4.0, 0.1, 2.0),), note="x > 0; Maxima adds a spurious i*pi branch"),
    SF("ExpIntegralEi", "Integral", 1, "ExpIntegralEi[{0}]", "expintegral_ei({0})",
       mp.ei, sc.expi, False, ALL, ((0.1, 8.0),), ((0.5, 4.0, -3.0, 3.0),)),
    SF("ExpIntegralE", "Integral", 1, "ExpIntegralE[2, {0}]", "expintegral_e(2, {0})",
       (lambda z: mp.expint(2, z)), (lambda z: sc.expn(2, z)), False, (MR,),
       ((0.2, 8.0),), (), note="order n = 2; Mathilda machine-real only (no complex/AP path)"),
    SF("LogIntegral", "Integral", 1, "LogIntegral[{0}]", "expintegral_li({0})",
       mp.li, (lambda z: sc.expi(math.log(z))), False, ALL,
       ((2.0, 20.0),), ((2.0, 10.0, -3.0, 3.0),), note="z > 1; li(z) = Ei(ln z) for scipy"),

    # ---- Hypergeometric / confluent -------------------------------------------------
    SF("Hypergeometric0F1", "Hypergeometric", 1, "Hypergeometric0F1[2, {0}]",
       "hypergeometric([], [2], {0})",
       (lambda z: mp.hyp0f1(2, z)), (lambda z: sc.hyp0f1(2, z)), True, ALL,
       ((-5.0, 5.0),), ((-3.0, 3.0, -3.0, 3.0),), note="0F1(; 2; z)"),
    SF("Hypergeometric1F1", "Hypergeometric", 1, "Hypergeometric1F1[1, 2, {0}]",
       "hypergeometric([1], [2], {0})",
       (lambda z: mp.hyp1f1(1, 2, z)), (lambda z: sc.hyp1f1(1, 2, z)), True, ALL,
       ((-5.0, 5.0),), ((-3.0, 3.0, -3.0, 3.0),), note="1F1(1; 2; z)"),
    SF("Hypergeometric2F1", "Hypergeometric", 1, "Hypergeometric2F1[1, 2, 3, {0}]",
       "hypergeometric([1, 2], [3], {0})",
       (lambda z: mp.hyp2f1(1, 2, 3, z)), (lambda z: sc.hyp2f1(1, 2, 3, z)), True, ALL,
       ((-0.9, 0.9),), ((-0.8, 0.8, -0.8, 0.8),), note="2F1(1,2;3;z), |z| < 1"),
    SF("HypergeometricPFQ", "Hypergeometric", 1, "HypergeometricPFQ[{{1, 2}}, {{3, 4}}, {0}]",
       "hypergeometric([1, 2], [3, 4], {0})",
       (lambda z: mp.hyper([1, 2], [3, 4], z)), None, False, ALL,
       ((-3.0, 3.0),), ((-2.0, 2.0, -2.0, 2.0),),
       note="2F2; Python is mpmath (scipy has no general pFq)", pylib="mpmath"),
    SF("HypergeometricU", "Hypergeometric", 1, "HypergeometricU[1, 2, {0}]", None,
       (lambda z: mp.hyperu(1, 2, z)), (lambda z: sc.hyperu(1, 2, z)), False, ALL,
       ((0.2, 8.0),), ((0.3, 4.0, -2.0, 2.0),), note="Tricomi U(1, 2, z), z > 0; Maxima symbolic"),
    SF("WhittakerM", "Hypergeometric", 1, "WhittakerM[1, 1/2, {0}]", None,
       (lambda z: mp.whitm(1, mp.mpf(1) / 2, z)), None, False, ALL,
       ((0.5, 8.0),), ((0.5, 4.0, -2.0, 2.0),),
       note="k=1, m=1/2; Python is mpmath", pylib="mpmath"),
    SF("WhittakerW", "Hypergeometric", 1, "WhittakerW[1, 1/2, {0}]", None,
       (lambda z: mp.whitw(1, mp.mpf(1) / 2, z)), None, False, ALL,
       ((0.5, 8.0),), ((0.5, 4.0, -2.0, 2.0),),
       note="k=1, m=1/2; Python is mpmath", pylib="mpmath"),

    # MeijerG: no scipy and no numeric Maxima G-function, so the Python column is
    # mpmath (meijerg) and Maxima is absent. Three shapes: a two-term p < q
    # Slater sum, a p == q function (both sides of the unit circle, where the
    # expansion switches from z to 1/z), and the documented log case of
    # integer-spaced b's (mpmath and Mathilda both perturb; Mathematica too).
    SF("MeijerG", "Hypergeometric", 1, "MeijerG[{{{{1/3}}, {{1/4}}}}, {{{{1/5, 1/6}}, {{1/7}}}}, {0}]", None,
       (lambda z: mp.meijerg([[mp.mpf(1) / 3], [mp.mpf(1) / 4]],
                             [[mp.mpf(1) / 5, mp.mpf(1) / 6], [mp.mpf(1) / 7]], z)), None, False, ALL,
       ((0.3, 6.0),), ((0.5, 4.0, -3.0, 3.0),),
       note="G^{2,1}_{2,3}; Python is mpmath (scipy has none); Maxima has no numeric G",
       pylib="mpmath"),
    SF("MeijerG22", "Hypergeometric", 1, "MeijerG[{{{{1/3, 1/2}}, {{}}}}, {{{{1/5, 1/7}}, {{}}}}, {0}]", None,
       (lambda z: mp.meijerg([[mp.mpf(1) / 3, mp.mpf(1) / 2], []],
                             [[mp.mpf(1) / 5, mp.mpf(1) / 7], []], z)), None, False, ALL,
       ((0.2, 4.0),), ((0.3, 3.0, -2.0, 2.0),),
       note="G^{2,2}_{2,2}: expansion in z inside |z| = 1, in 1/z outside", pylib="mpmath"),
    SF("MeijerGlog", "Hypergeometric", 1, "MeijerG[{{{{}}, {{2}}}}, {{{{1/2, 3/2}}, {{}}}}, {0}]", None,
       (lambda z: mp.meijerg([[], [2]], [[mp.mpf(1) / 2, mp.mpf(3) / 2], []], z)), None, False, ALL,
       ((0.3, 6.0),), ((0.5, 4.0, -3.0, 3.0),),
       note="log case: b's an integer apart (the documented example)", pylib="mpmath"),

    # ---- Misc -----------------------------------------------------------------------
    SF("Sinc", "Misc", 1, "Sinc[{0}]", "sin({0})/({0})",
       _sinc_mp, _sinc_sc, True, ALL, ((0.2, 10.0),), ((0.5, 5.0, -3.0, 3.0),),
       note="Mathematica Sinc = sin z / z"),
    SF("ProductLog", "Misc", 1, "ProductLog[{0}]", "lambert_w({0})",
       mp.lambertw, sc.lambertw, True, ALL, ((-0.3, 10.0),), ((-0.3, 5.0, -2.0, 2.0),),
       note="principal branch W0; z > -1/e on the real axis"),
]


# =====================================================================================
# EXACT heads -- integer/rational results, compared by exact string equality.
# scf/mpf unused; carries a size list and per-engine exact-call templates.
# =====================================================================================
@dataclass
class ExactSF:
    head: str
    wl: str                 # WL template, {0} = size n
    mx: Optional[str]       # Maxima template
    py: Optional[Callable]  # Python exact callable(n) -> int/str, or None
    sizes: tuple
    note: str = ""


def _harmonic_py(n):
    from fractions import Fraction
    s = Fraction(0)
    for k in range(1, n + 1):
        s += Fraction(1, k)
    return f"{s.numerator}/{s.denominator}"


EXACT: list[ExactSF] = [
    ExactSF("Fibonacci", "Fibonacci[{0}]", "fib({0})",
            (lambda n: str(__import__("sympy").fibonacci(n))), (500, 5000, 50000),
            note="exact bigint, fast doubling"),
    ExactSF("LucasL", "LucasL[{0}]", "lucas({0})",
            (lambda n: str(__import__("sympy").lucas(n))), (500, 5000, 50000),
            note="exact bigint"),
    ExactSF("BernoulliB", "BernoulliB[{0}]", "bern({0})",
            (lambda n: (lambda f: f"{f[0]}/{f[1]}")(mp.bernfrac(n))), (50, 200, 500),
            note="exact rational"),
    ExactSF("EulerE", "EulerE[{0}]", "euler({0})",
            (lambda n: str(__import__("sympy").euler(n))), (20, 50, 100),
            note="exact integer"),
    ExactSF("HarmonicNumber", "HarmonicNumber[{0}]", None, _harmonic_py, (100, 1000),
            note="exact rational sum; Maxima leaves it as a nounform"),
    ExactSF("BarnesG", "BarnesG[{0}]", None, None, (20, 60),
            note="exact superfactorial (GMP); only Mathilda + Wolfram"),
    ExactSF("Hyperfactorial", "Hyperfactorial[{0}]", None, None, (20, 60),
            note="exact (GMP); only Mathilda + Wolfram"),
    ExactSF("DirichletL", "DirichletL[12, 3, -{0}]", None, None, (10, 50, 150),
            note="exact generalized Bernoulli value (real character mod 12); only Mathilda + Wolfram"),
]


# Symbolic-only heads, excluded from the numeric sweep (documented in the write-up).
EXCLUDED = {
    "BesselJZero": "symbolic placeholder (Hadamard-product recogniser), no numeric path",
    "StieltjesGamma": "inert except StieltjesGamma[0] -> EulerGamma",
}
