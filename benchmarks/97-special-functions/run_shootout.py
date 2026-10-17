#!/usr/bin/env python3
"""run_shootout.py -- the special-functions shootout orchestrator.

Runs every registry case through four engines -- Mathilda, wolframscript (Mathematica),
Maxima, and Python (scipy for machine real/complex, mpmath for arbitrary precision) --
over three regimes (MR = machine real, C = complex, AP = arbitrary precision), and
records, per (head, regime, engine):

  * coverage     -- did the engine return a number at all?
  * scalar time  -- per-call latency (ms), timed in-engine over DISTINCT points so no
                    engine can serve a memoised result; min over R blocks; the whole
                    sweep is one process per engine so startup is never charged per call.
  * value        -- parsed to mpmath; used only for cross-engine AGREEMENT (no oracle).

Accuracy is agreement-only: for each point we report the number of leading significant
digits the present engines mutually share, and flag any engine that falls out of the
consensus (component-wise median of the others).  This doubles as the mapping-bug
detector -- a mis-mapped convention shows up as a disagreement.

Also runs (optional) a machine-real vectorized-throughput comparison (Mathilda packed
buffer vs scipy array) and an exact-arithmetic comparison (Fibonacci, BernoulliB, ...).

Usage:
    benchmarks/.venv/bin/python benchmarks/97-special-functions/run_shootout.py \
        [--only Gamma,BesselJ] [--regimes MR,C,AP] [--engines mathilda,wolfram,maxima,python] \
        [--ap-dps 50] [--no-throughput] [--no-exact] [--out results/<date>.json]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import random
import re
import subprocess
import sys
import tempfile
import time
import zlib
from fractions import Fraction

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

import registry as R  # noqa: E402

# ----------------------------------------------------------------------------- knobs
KVAL = 8           # value points per (case, regime), for agreement
KTIME = 512        # timing points per block (MR/C); large enough that the null-loop
                   # subtraction sits well above AbsoluteTiming resolution
KTIME_AP = 16      # timing points per block (AP; each call is costly)
RBLK = 5           # timing blocks (take the min per-call across blocks)
RBLK_AP = 2
AP_DPS = 50        # requested digits for the AP regime
DENOM = 100000     # denominator for AP rational sample points
SEED0 = 20261016
THRU_N = 10 ** 6   # vectorized-throughput buffer size

MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(REPO, "Mathilda"))
WOLFRAM = os.environ.get(
    "WOLFRAMSCRIPT",
    "/Applications/Mathematica.app/Contents/MacOS/wolframscript",
)
MAXIMA = os.environ.get("MAXIMA_BIN", "maxima")

ENGINES = ("mathilda", "wolfram", "maxima", "python")
MR, C, AP = R.MR, R.C, R.AP


# ============================================================ numeric-literal rendering
def wl_real(x):
    """A WL literal for the double x that parses as a MACHINE number.

    ``repr`` gives the shortest round-trip string, which is up to 17 significant
    digits; both Mathilda and Mathematica read a literal with more digits than
    $MachinePrecision as ARBITRARY precision (Precision[1.2345678901234567] ~ 16.26),
    which silently routed ~40% of the "machine" sample through the MPFR paths.  The
    trailing backtick (``1.2345678901234567```) pins it to machine precision while
    keeping the value bit-identical to the double the Python engine sees.
    """
    s = repr(float(x))
    if "e" in s or "E" in s:
        s = s.lower()
        m, _, e = s.partition("e")
        return f"{m}`*^{int(e)}"
    return s + "`"


def wl_complex(z):
    return f"Complex[{wl_real(z.real)}, {wl_real(z.imag)}]"


def wl_rat(fr: Fraction):
    return f"{fr.numerator}/{fr.denominator}" if fr.denominator != 1 else f"{fr.numerator}"


def mx_real(x):
    x = float(x)
    s = repr(x)
    if "e" in s:
        m, _, e = s.partition("e")
        if "." not in m:
            m += ".0"
        return f"{m}e{int(e)}"
    if "." not in s:
        s += ".0"
    return s


def mx_complex(z):
    return f"{mx_real(z.real)} + ({mx_real(z.imag)})*%i"


def mx_rat(fr: Fraction):
    return f"{fr.numerator}/{fr.denominator}" if fr.denominator != 1 else f"{fr.numerator}"


def render_slot(regime, val):
    """Render one argument value as the engine literal (WL, MX) for a regime."""
    if regime == MR:
        return wl_real(val), mx_real(val)
    if regime == C:
        return wl_complex(val), mx_complex(val)
    # AP: val is a Fraction
    return wl_rat(val), mx_rat(val)


# ============================================================ value-string parsing
def parse_wl(s):
    """Parse a Mathilda/wolframscript InputForm value string to mpmath (mpf/mpc).

    wolframscript tags every inexact number with a backtick precision suffix
    (``1.329...`50.``) and uses ``*^`` for scientific notation; Mathilda prints neither.
    We strip the backtick suffix, normalise ``*^`` -> ``e`` and ``*I`` -> ``j``.
    Parsing happens with mp.mp.dps already raised (see main), so 50-digit AP strings
    are not silently truncated to machine precision.
    """
    s = s.strip()
    if not s or any(b in s for b in ("NONNUMERIC", "Indeterminate", "Infinity", "Failed", "[")):
        return None
    s = s.replace(" ", "")
    s = re.sub(r"`[0-9.]+", "", s)   # drop Mathematica's `precision suffix
    s = s.replace("*^", "e")
    try:
        if "I" in s:
            # AP complex is not sampled, so machine-precision float() round-trip is fine here
            return mp.mpc(complex(s.replace("*I", "j").replace("I", "j")))
        return mp.mpf(s)
    except (ValueError, TypeError):
        return None


def _mx_scalar(s):
    s = s.replace(" ", "").replace("b", "e")
    try:
        return mp.mpf(s)
    except (ValueError, TypeError):
        return None


def parse_mx(s):
    """Parse a Maxima value string: real, interval(c,r), complex 'a %i + b', bigfloat '..bN'.

    Maxima sometimes returns a nounform (``'realpart(psi[0](..))+%i*'imagpart(..)``) or a
    partially-evaluated expression still carrying ``%pi`` (e.g. legendre_q at a complex
    argument): those are NOT numbers and must read as absent, not as a bogus value that
    would show up as a spurious disagreement.
    """
    s = s.strip()
    if not s or "'" in s:            # Maxima nounform quote -> symbolic
        return None
    if "interval(" in s:
        inner = s[s.index("interval(") + 9:]
        c = inner.split(",")[0]
        return _mx_scalar(c)
    if "%i" in s:
        probe = s.replace("%i", "").replace("j", "").replace("e", "").replace("E", "")
        if re.search(r"[A-Za-z%]", probe):   # leftover symbol (%pi, function name) -> not numeric
            return None
        # forms seen: "C*%i+R", "R-C*%i", "-(C*%i)-R", "C*%i-R", bare "%i".
        t = (s.replace(" ", "").replace("(", "").replace(")", "")
             .replace("*%i", "j").replace("%i", "j"))
        re_sum, im_sum = mp.mpf(0), mp.mpf(0)
        for p in re.findall(r"[+-]?[0-9.]+j?|[+-]?j", t):
            if p.endswith("j"):
                coeff = p[:-1]
                coeff = "1" if coeff in ("", "+") else ("-1" if coeff == "-" else coeff)
                im_sum += mp.mpf(coeff)
            else:
                re_sum += mp.mpf(p)
        return mp.mpc(re_sum, im_sum)
    # symbolic result -> reject (allow only digits, sign, dot, and b/e exponents)
    if re.search(r"[a-zA-Z]", s.replace("b", "").replace("e", "").replace("E", "")):
        return None
    return _mx_scalar(s)


# ============================================================ point generation
def gen_points(sf: R.SF, regime, seed):
    rng = random.Random(seed)

    def one_real(lo, hi):
        return lo + (hi - lo) * rng.random()

    def one_rat(lo, hi):
        x = lo + (hi - lo) * rng.random()
        return Fraction(round(x * DENOM), DENOM)

    def one_cplx(box):
        rl, rh, il, ih = box
        return complex(rl + (rh - rl) * rng.random(), il + (ih - il) * rng.random())

    n = sf.nslots

    def one():
        if regime == MR:
            return tuple(one_real(*sf.dom_r[j]) for j in range(n))
        if regime == AP:
            return tuple(one_rat(*sf.dom_r[j]) for j in range(n))
        return tuple(one_cplx(sf.dom_c[j]) for j in range(n))

    kt = KTIME_AP if regime == AP else KTIME
    rb = RBLK_AP if regime == AP else RBLK
    vpts = [one() for _ in range(KVAL)]
    tblocks = [[one() for _ in range(kt)] for _ in range(rb)]
    return vpts, tblocks


def case_regimes(sf: R.SF, want_regimes):
    return [r for r in (MR, C, AP) if r in sf.regimes and r in want_regimes]


def _seed(head, regime):
    """Stable per-(head,regime) seed so every engine samples the SAME points."""
    return SEED0 ^ (zlib.crc32(f"{head}|{regime}".encode()) & 0x7FFFFFFF)


# ============================================================ WL / Mathilda script gen
def _wl_timed(K, body):
    """One timed WL loop: ``ClearSystemCache`` first (Mathematica memoises some special
    functions; on Mathilda the symbol is undefined and the call is an inert no-op), then
    a bare ``Do`` -- no ``Quiet``: timing points are in-domain, and messages are the
    accuracy pass's concern, not the timer's."""
    return f"(ClearSystemCache[]; First[AbsoluteTiming[Do[{body}, {{i, 1, {K}}}]]])"


def gen_wl_script(cases, regimes, ap_dps, mode="values"):
    """Emit one .m consumed verbatim by both Mathilda and wolframscript.

    ``mode="values"`` evaluates the KVAL agreement points; ``mode="timing"`` times the
    blocks.  The two run as SEPARATE engine processes so no accuracy-pass evaluation can
    warm a cache the timer then reads.

    Timing protocol, per block b of K distinct points:
      * ``null``   -- the same ``Do`` + ``Part`` loop with no call (the loop floor);
      * ``raw``    -- the loop calling the head.
    The driver reports ``net = min_b raw - min_b null``: the head's own cost, with the
    ``Do``/``Part`` overhead removed.  After the blocks, block 0 is re-timed verbatim on
    its already-seen points (the CACHE CANARY): an engine that memoised argument-level
    results answers that repeat far faster than the first pass, and the cell is flagged.
    """
    L = []
    L.append(r'Print["VERSION\t", ToString[InputForm[$VersionNumber]], "\t", $Version];')
    L.append("Off[General::munfl]; Off[General::ovfl]; Off[General::infy];")
    for sf in cases:
        for regime in case_regimes(sf, regimes):
            seed = _seed(sf.head, regime)
            vpts, tblocks = gen_points(sf, regime, seed)
            lab = f"{sf.head}"
            if mode == "values":
                for idx, tup in enumerate(vpts):
                    lits = [render_slot(regime, v)[0] for v in tup]
                    body = sf.wl.format(*lits)
                    if regime == AP:
                        expr = f"N[{body}, {ap_dps}]"
                        show = "ToString[InputForm[r]]"
                    else:
                        expr = body
                        show = "ToString[InputForm[SetPrecision[r, 17]]]"
                    L.append(
                        f"With[{{r = Quiet[{expr}]}}, "
                        f'Print["RESULT\\t{regime}\\t{lab}\\t{idx}\\t", '
                        f'If[NumericQ[r], {show}, "NONNUMERIC"]]];'
                    )
                continue
            # ---- timing
            tag = abs(seed) % 100000
            blk_vars = []
            body0 = None
            for b, block in enumerate(tblocks):
                slot_lists = []
                for j in range(sf.nslots):
                    lits = ", ".join(render_slot(regime, tup[j])[0] for tup in block)
                    nm = f"w{tag}b{b}s{j}"
                    L.append(f"{nm} = {{{lits}}};")
                    slot_lists.append(nm)
                accessor = [f"{nm}[[i]]" for nm in slot_lists]
                body = sf.wl.format(*accessor)
                inner = f"N[{body}, {ap_dps}]" if regime == AP else body
                null = "; ".join(accessor) if len(accessor) == 1 else "(" + "; ".join(accessor) + ")"
                K = len(block)
                tn, tr = f"tn{tag}b{b}", f"tr{tag}b{b}"
                L.append(f"{tn} = {_wl_timed(K, null)};")
                L.append(f"{tr} = {_wl_timed(K, inner)};")
                blk_vars.append((tn, tr, K))
                if b == 0:
                    body0 = (inner, K)
            if blk_vars:
                tc = f"tc{tag}"
                L.append(f"{tc} = {_wl_timed(body0[1], body0[0])};")
                fmt = lambda v, K: f"ToString[InputForm[1000.0 {v}/{K}]]"
                parts = ' ,"\\t", '.join(
                    f'{fmt(tr, K)}, ":", {fmt(tn, K)}' for tn, tr, K in blk_vars
                )
                L.append(f'Print["TIMEBLK\\t{regime}\\t{lab}\\t", {parts}];')
                L.append(f'Print["CANARY\\t{regime}\\t{lab}\\t", {fmt(tc, body0[1])}];')
    return "\n".join(L) + "\n"


# ============================================================ Maxima script gen
def gen_mx_script(cases, regimes, ap_dps, mode="values"):
    """Maxima twin of gen_wl_script (same values/timing split, null loop and canary)."""
    L = ["tb: ascii(9)$", "load(orthopoly)$", "load(bffac)$", "ratprint: false$"]
    for sf in cases:
        if sf.mx is None:
            continue
        for regime in case_regimes(sf, regimes):
            seed = _seed(sf.head, regime)
            vpts, tblocks = gen_points(sf, regime, seed)
            lab = f"{sf.head}"
            if regime == AP:
                L.append(f"fpprec: {ap_dps + 12}$")
            # ---- values
            for idx, tup in (enumerate(vpts) if mode == "values" else ()):
                lits = [render_slot(regime, v)[1] for v in tup]
                body = sf.mx.format(*lits)
                if regime == AP:
                    conv = f"bfloat({body})"
                elif regime == C:
                    conv = f"rectform(float({body}))"
                else:
                    conv = f"float({body})"
                L.append(
                    f"block([r], r: errcatch({conv}), "
                    f'if r = [] then printf(true, "RESULT~a{regime}~a{lab}~a{idx}~aABSENT~%", tb, tb, tb, tb) '
                    f'else printf(true, "RESULT~a{regime}~a{lab}~a{idx}~a~a~%", tb, tb, tb, tb, string(r[1])))$'
                )
            if mode != "timing":
                continue
            # ---- timing (per block: raw loop and null loop; then the cache canary)
            parts = []
            for b, block in enumerate(tblocks):
                for j in range(sf.nslots):
                    lits = ", ".join(render_slot(regime, tup[j])[1] for tup in block)
                    L.append(f"mb{b}s{j}: [{lits}]$")
                accessor = [f"mb{b}s{j}[i]" for j in range(sf.nslots)]
                body = sf.mx.format(*accessor)
                conv = f"bfloat({body})" if regime == AP else f"float({body})"
                null = f"({', '.join(accessor)})"
                K = len(block)
                L.append(
                    f"mn{b}: elapsed_real_time()$ "
                    f"for i thru {K} do {null}$ "
                    f"mn{b}: 1000.0 * (elapsed_real_time() - mn{b}) / {K}$"
                )
                L.append(
                    f"mt{b}: elapsed_real_time()$ "
                    f"for i thru {K} do errcatch({conv})$ "
                    f"mt{b}: 1000.0 * (elapsed_real_time() - mt{b}) / {K}$"
                )
                parts.append((f"mt{b}", f"mn{b}"))
                if b == 0:
                    canary = (conv, K)
            if parts:
                L.append(
                    f"mc: elapsed_real_time()$ "
                    f"for i thru {canary[1]} do errcatch({canary[0]})$ "
                    f"mc: 1000.0 * (elapsed_real_time() - mc) / {canary[1]}$"
                )
                pieces = ['"TIMEBLK"', "tb", f'"{regime}"', "tb", f'"{lab}"']
                for (pt, pn) in parts:
                    pieces += ["tb", f"string({pt})", '":"', f"string({pn})"]
                L.append(f'printf(true, "~a~%", sconcat({", ".join(pieces)}))$')
                L.append(f'printf(true, "~a~%", sconcat("CANARY", tb, "{regime}", tb, "{lab}", tb, string(mc)))$')
    return "\n".join(L) + "\n"


# ============================================================ subprocess runners
def run_proc(cmd, script_text, suffix, timeout):
    with tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False, dir=os.environ.get("CLAUDE_JOB_DIR", None) and os.path.join(os.environ["CLAUDE_JOB_DIR"], "tmp") or None) as f:
        f.write(script_text)
        path = f.name
    full = cmd + [path]
    try:
        p = subprocess.run(full, capture_output=True, text=True, timeout=timeout)
        return p.stdout, p.stderr, p.returncode
    except subprocess.TimeoutExpired as e:
        return (e.stdout or ""), (e.stderr or "") + "\n[TIMEOUT]", -9
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def parse_engine_stdout(stdout, parser):
    """Return (values, times, version).

    values[(regime,head,idx)] -> mpval|None;
    times[(regime,head)] -> {"raw": [ms/call per block], "null": [...], "canary": ms/call},
    reduced to a single net figure by summarize_timing().
    """
    values, times, version = {}, {}, None
    for line in stdout.splitlines():
        if "\t" not in line:
            continue
        parts = line.split("\t")
        tag = parts[0]
        if tag == "VERSION":
            version = parts[1] if len(parts) > 1 else None
            if len(parts) > 2:
                version = f"{parts[1]} ({parts[2].strip()})"
        elif tag == "RESULT" and len(parts) >= 5:
            regime, head, idx, val = parts[1], parts[2], parts[3], parts[4]
            values[(regime, head, int(idx))] = None if val.strip() in ("ABSENT", "NONNUMERIC", "") else parser(val)
        elif tag == "TIMEBLK" and len(parts) >= 4:
            regime, head = parts[1], parts[2]
            raw, null = [], []
            for x in parts[3:]:
                r, _, n = x.partition(":")
                try:
                    raw.append(float(r))
                    null.append(float(n) if n else 0.0)
                except ValueError:
                    pass
            if raw:
                times.setdefault((regime, head), {}).update(raw=raw, null=null)
        elif tag == "CANARY" and len(parts) >= 4:
            try:
                times.setdefault((parts[1], parts[2]), {})["canary"] = float(parts[3])
            except ValueError:
                pass
    return values, times, version


# ============================================================ Python engine (in-process)
def run_python(cases, regimes, ap_dps):
    values, times, libs = {}, {}, {}
    for sf in cases:
        for regime in case_regimes(sf, regimes):
            seed = _seed(sf.head, regime)
            vpts, tblocks = gen_points(sf, regime, seed)
            use_mp = (sf.scf is None) or (regime == AP) or (regime == C and not sf.sc_complex)

            if regime == AP:
                mp.mp.dps = ap_dps + 12
                fn = lambda *vs, _f=sf.mpf: _f(*[mp.mpf(v.numerator) / v.denominator for v in vs])
            elif use_mp:
                mp.mp.dps = 20
                fn = sf.mpf
            else:
                fn = sf.scf
                # probe: scipy declines some complex args (e.g. eval_hermite) -> fall back to mpmath
                try:
                    fn(*vpts[0])
                except Exception:
                    use_mp, fn = True, sf.mpf
                    mp.mp.dps = 20
            libs[(regime, sf.head)] = "mpmath" if use_mp else "scipy"

            # values (coerce scipy's numpy scalar/0-d array via complex(); keep mpmath exact for AP)
            for idx, tup in enumerate(vpts):
                try:
                    r = fn(*tup)
                    values[(regime, sf.head, idx)] = mp.mpc(r) if regime == AP else mp.mpc(complex(r))
                except Exception:
                    values[(regime, sf.head, idx)] = None
            # timing
            wrng = random.Random(seed ^ 0x5A5A5A)
            warm = gen_points(sf, regime, wrng.randrange(1 << 30))[0][:3]
            t = time_py(fn, tblocks, regime, ap_dps, warm)
            if t is not None:
                times[(regime, sf.head)] = t
    return values, times, libs


def time_py(fn, tblocks, regime, ap_dps, warm):
    """Python twin of the WL timing protocol: per block a raw loop and a null loop
    (same iteration, no call), then the block-0 cache canary.  ``warm`` points come from a
    separate RNG stream, never from a timed block -- mpmath caches, so warming on a timed
    point would hand the timer a hot value."""
    mp.mp.dps = (ap_dps + 12) if regime == AP else 20
    try:
        for tup in warm:
            fn(*tup)
    except Exception:
        return None
    raw, null = [], []
    pc = time.perf_counter
    for blk in tblocks:
        t0 = pc()
        for tup in blk:
            pass
        null.append((pc() - t0) / len(blk) * 1000.0)
        t0 = pc()
        try:
            for tup in blk:
                fn(*tup)
        except Exception:
            return None
        raw.append((pc() - t0) / len(blk) * 1000.0)
    t0 = pc()
    for tup in tblocks[0]:
        fn(*tup)
    canary = (pc() - t0) / len(tblocks[0]) * 1000.0
    return {"raw": raw, "null": null, "canary": canary}


CANARY_RATIO = 0.5   # repeat-pass faster than this x first pass => argument-level memo


def summarize_timing(t):
    """Reduce {"raw","null","canary"} to the figures the report uses.

    net = min(raw) - min(null): each min is the best estimate of its own true cost, so
    the difference is the head's cost with the loop/Part overhead removed (floored at a
    nanosecond so a sub-resolution head still ranks as "fastest" rather than vanishing).
    cache_suspect: the canary re-ran block 0 on already-seen points; if that repeat is
    much faster than block 0's first pass, the engine is serving memoised results and
    its number measures the cache, not the kernel.
    """
    if not t or not t.get("raw"):
        return None
    raw_min = min(t["raw"])
    null_min = min(t["null"]) if t.get("null") else 0.0
    net = max(raw_min - null_min, 1e-6)
    canary = t.get("canary")
    first = t["raw"][0]
    suspect = bool(canary is not None and first > 0 and canary < CANARY_RATIO * first)
    return {"net": net, "raw": raw_min, "null": null_min, "canary": canary,
            "first": first, "cache_suspect": suspect}


# ============================================================ agreement / accuracy
def consensus(vals):
    """Component-wise median of a list of mpc values."""
    re_s = sorted(v.real for v in vals)
    im_s = sorted(v.imag for v in vals)
    n = len(vals)
    med = lambda s: s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    return mp.mpc(med(re_s), med(im_s))


def rel_err(v, ref):
    d = abs(v - ref)
    scale = abs(ref)
    if scale < mp.mpf(10) ** -300:
        scale = mp.mpf(1)
    return d / scale


def digits(relerr):
    if relerr <= 0:
        return 99.0
    return float(-mp.log10(relerr))


def analyze_agreement(eng_values, cases, regimes, ap_dps):
    """Return per (head,regime): {agree_digits, disagree:[{engine,digits}], coverage:{...}}."""
    out = {}
    for sf in cases:
        for regime in case_regimes(sf, regimes):
            mp.mp.dps = (ap_dps + 15) if regime == AP else 30
            present = {e: eng_values[e] for e in eng_values}
            # coverage: engine produced >=1 numeric value for this case/regime
            cov = {}
            for e in ENGINES:
                vs = present.get(e, {})
                cov[e] = any(vs.get((regime, sf.head, i)) is not None for i in range(KVAL))
            # per-point agreement
            per_eng_worst = {e: 99.0 for e in ENGINES if cov[e]}
            pt_digits = []
            for i in range(KVAL):
                vals = {e: present[e].get((regime, sf.head, i)) for e in ENGINES if cov[e]}
                vals = {e: v for e, v in vals.items() if v is not None}
                if len(vals) < 2:
                    continue
                cons = consensus(list(vals.values()))
                errs = {e: digits(rel_err(v, cons)) for e, v in vals.items()}
                pt_digits.append(min(errs.values()))
                for e, d in errs.items():
                    per_eng_worst[e] = min(per_eng_worst[e], d)
            agree = min(pt_digits) if pt_digits else None
            thr = (ap_dps - 6) if regime == AP else 9.0
            disagree = []
            if len([e for e in per_eng_worst if cov[e]]) >= 2:
                best = max(per_eng_worst.values())
                for e, d in per_eng_worst.items():
                    if d < thr and best - d > 3:
                        disagree.append({"engine": e, "digits": round(d, 2)})
            out[(sf.head, regime)] = {
                "agree_digits": None if agree is None else round(agree, 2),
                "disagree": disagree,
                "coverage": cov,
            }
    return out


# ============================================================ throughput (machine real)
def gen_wl_throughput(cases):
    L = [r'Print["THRU_VERSION\t", $Version];', "SeedRandom[42];"]
    for sf in cases:
        if sf.nslots != 1 or MR not in sf.regimes:
            continue
        lo, hi = sf.dom_r[0]
        L.append(f"vv = RandomReal[{{{lo}, {hi}}}, {THRU_N}];")
        body = sf.wl.format("vv")
        L.append(f"{sf.head};")  # touch symbol
        L.append(f"tt = First[AbsoluteTiming[{body};]];")
        L.append(f'Print["THRU\\t{sf.head}\\t", ToString[InputForm[1000.0 tt]]];')
    return "\n".join(L) + "\n"


def run_python_throughput(cases):
    import numpy as np
    rng = np.random.default_rng(42)
    out = {}
    for sf in cases:
        if sf.nslots != 1 or MR not in sf.regimes or sf.scf is None:
            continue
        lo, hi = sf.dom_r[0]
        v = rng.uniform(lo, hi, THRU_N)
        try:
            sf.scf(v[:16])  # warm + check it vectorizes
            best = None
            for _ in range(3):
                t0 = time.perf_counter()
                sf.scf(v)
                dt = (time.perf_counter() - t0) * 1000.0
                best = dt if best is None else min(best, dt)
            out[sf.head] = best
        except Exception:
            pass
    return out


# ============================================================ exact arithmetic
def gen_wl_exact(exacts):
    L = [r'Print["EXACT_VERSION\t", $Version];']
    for ex in exacts:
        for n in ex.sizes:
            body = ex.wl.format(n)
            L.append(f"te = First[AbsoluteTiming[r = {body};]];")
            if n == ex.sizes[0]:
                L.append(f'Print["EXACTVAL\\t{ex.head}\\t{n}\\t", ToString[InputForm[r]]];')
            L.append(f'Print["EXACTT\\t{ex.head}\\t{n}\\t", ToString[InputForm[1000.0 te]]];')
    return "\n".join(L) + "\n"


def gen_mx_exact(exacts):
    L = ["tb: ascii(9)$"]
    for ex in exacts:
        if ex.mx is None:
            continue
        for n in ex.sizes:
            body = ex.mx.format(n)
            L.append(f"et: elapsed_real_time()$ r: {body}$ et: 1000.0*(elapsed_real_time()-et)$")
            if n == ex.sizes[0]:
                L.append(f'printf(true, "EXACTVAL~a{ex.head}~a{n}~a~a~%", tb, tb, tb, string(r))$')
            L.append(f'printf(true, "EXACTT~a{ex.head}~a{n}~a~a~%", tb, tb, tb, string(et))$')
    return "\n".join(L) + "\n"


def run_python_exact(exacts):
    import sympy  # noqa: F401  -- pay the import cost once, up front, not inside a timed call
    vals, times = {}, {}
    for ex in exacts:
        if ex.py is None:
            continue
        try:
            ex.py(ex.sizes[0])  # warm-up: absorb lazy init / per-head caching
        except Exception:
            continue
        for n in ex.sizes:
            try:
                t0 = time.perf_counter()
                r = ex.py(n)
                dt = (time.perf_counter() - t0) * 1000.0
            except Exception:
                continue
            times[(ex.head, n)] = dt
            if n == ex.sizes[0]:
                vals[(ex.head, n)] = str(r)
    return vals, times


# ============================================================ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--regimes", default="MR,C,AP")
    ap.add_argument("--engines", default=",".join(ENGINES))
    ap.add_argument("--ap-dps", type=int, default=AP_DPS)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--no-throughput", action="store_true")
    ap.add_argument("--no-exact", action="store_true")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    want_regimes = set(args.regimes.split(","))
    want_engines = set(args.engines.split(","))
    only = set(h.strip() for h in args.only.split(",") if h.strip())
    cases = [sf for sf in R.NUMERIC if not only or sf.head in only]
    exacts = [ex for ex in R.EXACT if not only or ex.head in only]
    ap_dps = args.ap_dps
    # Raise working precision BEFORE any engine stdout is parsed: parse_wl/parse_mx call
    # mp.mpf() on the CAS value strings, and at the default dps=15 a 50-digit AP string
    # would be silently truncated to machine precision (making the full-precision Python
    # column look like the outlier).  run_python manages its own per-regime dps.
    mp.mp.dps = ap_dps + 25

    print(f"[shootout] {len(cases)} numeric cases, regimes={sorted(want_regimes)}, "
          f"engines={sorted(want_engines)}, ap_dps={ap_dps}", flush=True)

    eng_values, eng_times, versions = {}, {}, {}

    # ---- Mathilda + wolframscript share the generated .m.  Values and timing run as two
    # separate processes per engine, so the accuracy pass can never warm the timer.
    def run_two_pass(cmd, gen, suffix, parser):
        vals, tms, ver = {}, {}, None
        for mode in ("values", "timing"):
            so, se, rc = run_proc(cmd, gen(cases, want_regimes, ap_dps, mode=mode), suffix, args.timeout)
            v, t, vr = parse_engine_stdout(so, parser)
            vals.update(v)
            tms.update(t)
            ver = ver or vr
            if se.strip() and cmd[0] == MATHILDA:
                print(f"  [{mode} stderr]", se.strip()[:400], flush=True)
        return vals, tms, ver

    if "mathilda" in want_engines:
        print("[shootout] running Mathilda ...", flush=True)
        v, t, ver = run_two_pass([MATHILDA, "-file"], gen_wl_script, ".m", parse_wl)
        eng_values["mathilda"], eng_times["mathilda"], versions["mathilda"] = v, t, ver
    if "wolfram" in want_engines and os.path.exists(WOLFRAM):
        print("[shootout] running wolframscript ...", flush=True)
        v, t, ver = run_two_pass([WOLFRAM, "-file"], gen_wl_script, ".m", parse_wl)
        eng_values["wolfram"], eng_times["wolfram"], versions["wolfram"] = v, t, ver

    # ---- Maxima
    if "maxima" in want_engines:
        print("[shootout] running Maxima ...", flush=True)
        v, t, _ = run_two_pass([MAXIMA, "--very-quiet", "-b"], gen_mx_script, ".mac", parse_mx)
        eng_values["maxima"], eng_times["maxima"] = v, t
        try:
            versions["maxima"] = subprocess.run([MAXIMA, "--version"], capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:
            versions["maxima"] = "Maxima"

    # ---- Python
    if "python" in want_engines:
        print("[shootout] running Python (scipy / mpmath) ...", flush=True)
        import scipy, numpy
        v, t, libs = run_python(cases, want_regimes, ap_dps)
        eng_values["python"], eng_times["python"] = v, t
        versions["python"] = f"scipy {scipy.__version__}, mpmath {mp.__version__}, numpy {numpy.__version__}"
    else:
        libs = {}

    # ---- agreement
    print("[shootout] analyzing agreement ...", flush=True)
    agree = analyze_agreement(eng_values, cases, want_regimes, ap_dps)

    # ---- assemble per-case JSON
    cases_json = {}
    for sf in cases:
        rj = {}
        for regime in case_regimes(sf, want_regimes):
            a = agree.get((sf.head, regime), {})
            tj, tdet, suspect = {}, {}, []
            for e in ENGINES:
                st = summarize_timing(eng_times.get(e, {}).get((regime, sf.head)))
                tj[e] = round(st["net"], 7) if st else None
                if st:
                    tdet[e] = {k: (round(v, 7) if isinstance(v, float) else v) for k, v in st.items()}
                    tdet[e]["raw_blocks"] = eng_times[e][(regime, sf.head)].get("raw")
                    tdet[e]["null_blocks"] = eng_times[e][(regime, sf.head)].get("null")
                    if st["cache_suspect"]:
                        suspect.append(e)
            rj[regime] = {
                "time_ms": tj,               # NET per-call ms (loop floor subtracted)
                "timing": tdet,              # raw / null / canary detail per engine
                "cache_suspect": suspect,    # engines whose canary says "memoised"
                "agree_digits": a.get("agree_digits"),
                "disagree": a.get("disagree", []),
                "coverage": a.get("coverage", {}),
                "python_lib": libs.get((regime, sf.head)),
            }
        cases_json[sf.head] = {"family": sf.family, "note": sf.note, "regimes": rj}

    result = {
        "meta": {
            "date": dt.date.today().isoformat(),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "engines": versions,
            "ap_dps": ap_dps,
            "kval": KVAL, "ktime": KTIME, "rblk": RBLK,
            "timing": "net = min(raw) - min(null loop); machine literals; ClearSystemCache per block; canary",
        },
        "cases": cases_json,
    }

    # ---- throughput
    if not args.no_throughput and {"mathilda", "wolfram", "python"} & want_engines:
        print("[shootout] throughput (10^6 buffer) ...", flush=True)
        thru = {"mathilda": {}, "wolfram": {}, "python": {}}
        wl = gen_wl_throughput(cases)
        if "mathilda" in want_engines:
            so, _, _ = run_proc([MATHILDA, "-file"], wl, ".m", args.timeout)
            for ln in so.splitlines():
                if ln.startswith("THRU\t"):
                    _, h, ms = ln.split("\t")[:3]
                    try:
                        thru["mathilda"][h] = float(ms)
                    except ValueError:
                        pass
        if "wolfram" in want_engines and os.path.exists(WOLFRAM):
            so, _, _ = run_proc([WOLFRAM, "-file"], wl, ".m", args.timeout)
            for ln in so.splitlines():
                if ln.startswith("THRU\t"):
                    _, h, ms = ln.split("\t")[:3]
                    try:
                        thru["wolfram"][h] = float(ms)
                    except ValueError:
                        pass
        if "python" in want_engines:
            thru["python"] = run_python_throughput(cases)
        result["throughput"] = thru

    # ---- exact
    if not args.no_exact:
        print("[shootout] exact arithmetic ...", flush=True)
        exact = {"mathilda": {"t": {}, "v": {}}, "wolfram": {"t": {}, "v": {}},
                 "maxima": {"t": {}, "v": {}}, "python": {"t": {}, "v": {}}}
        if {"mathilda", "wolfram"} & want_engines:
            wl = gen_wl_exact(exacts)
            for e, binpath in (("mathilda", MATHILDA), ("wolfram", WOLFRAM)):
                if e not in want_engines or (e == "wolfram" and not os.path.exists(WOLFRAM)):
                    continue
                so, _, _ = run_proc([binpath, "-file"], wl, ".m", args.timeout)
                _parse_exact(so, exact[e])
        if "maxima" in want_engines:
            so, _, _ = run_proc([MAXIMA, "--very-quiet", "-b"], gen_mx_exact(exacts), ".mac", args.timeout)
            _parse_exact(so, exact["maxima"])
        if "python" in want_engines:
            v, t = run_python_exact(exacts)
            exact["python"]["v"] = {f"{h}:{n}": val for (h, n), val in v.items()}
            exact["python"]["t"] = {f"{h}:{n}": tm for (h, n), tm in t.items()}
        result["exact"] = {"sizes": {ex.head: list(ex.sizes) for ex in exacts},
                           "notes": {ex.head: ex.note for ex in exacts}, "engines": exact}

    # ---- write
    out = args.out or os.path.join(HERE, "results", f"{dt.date.today().isoformat()}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(result, f, indent=1, default=str)
    print(f"[shootout] wrote {out}", flush=True)
    _print_summary(result, cases, want_regimes)


def _parse_exact(stdout, store):
    for ln in stdout.splitlines():
        if "\t" not in ln:
            continue
        p = ln.split("\t")
        if p[0] == "EXACTVAL" and len(p) >= 4:
            store["v"][f"{p[1]}:{p[2]}"] = p[3].strip()
        elif p[0] == "EXACTT" and len(p) >= 4:
            try:
                store["t"][f"{p[1]}:{p[2]}"] = float(p[3])
            except ValueError:
                pass


def _print_summary(result, cases, regimes):
    print("\n=== summary ===")
    for sf in cases:
        cj = result["cases"].get(sf.head, {})
        for regime, rj in cj.get("regimes", {}).items():
            cov = rj["coverage"]
            covs = "".join(e[0].upper() if cov.get(e) else "." for e in ENGINES)
            ad = rj["agree_digits"]
            dis = ",".join(d["engine"] for d in rj["disagree"])
            flag = f"  DISAGREE[{dis}]" if dis else ""
            print(f"  {sf.head:20s} {regime:3s} cov={covs} agree={ad}{flag}")


if __name__ == "__main__":
    main()
