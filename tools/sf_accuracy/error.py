"""Accuracy gate for the error-function family's machine kernels.

Heads: Erf, Erfc, Erfi, FresnelC, FresnelS, Sinc (machine real and machine
complex), InverseErf, InverseErfc (machine real).

Drives the Mathilda binary (`-file`, one evaluation per point, printed with
ToString[InputForm[SetPrecision[r, 17]]] so the double is recovered exactly)
on machine literals (trailing backtick; complex as Complex[a`, b`]) and
compares against mpmath at 30 digits.

Two error figures per point, both in ulps of the result:
  rel     |got - ref| / ulp(|ref|)
  scaled  |got - ref| / ulp(max(|ref|, |z f'(z)|))
`scaled` is the brief's "absolute error scaled by the local magnitude" near a
zero: |z f'(z)| is the size of the change a one-ulp perturbation of the
argument causes, i.e. the magnitude of the quantities the value is computed
from. Away from zeros the two coincide.

Usage (from the worktree root):
    /Users/user/Mathilda_dev/Mathilda/benchmarks/.venv/bin/python -I \
        tools/sf_accuracy/error.py [--bin ./Mathilda] [--n 2000] [--heads Erf,...]
"""
import argparse, os, re, subprocess, sys, tempfile
import numpy as np
import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sf_fit"))
import fitlib  # erfinv_mp / erfcinv_mp (Newton-polished, fast)


def ulp(a):
    a = mp.mpf(abs(a))
    if a == 0:
        return mp.mpf(2) ** -1074
    return mp.mpf(2) ** (int(mp.floor(mp.log(a, 2))) - 52)


def ref_and_scale(head, z):
    if head == "Erf":
        return mp.erf(z), abs(z * 2 / mp.sqrt(mp.pi) * mp.exp(-z * z))
    if head == "Erfc":
        return mp.erfc(z), abs(z * 2 / mp.sqrt(mp.pi) * mp.exp(-z * z))
    if head == "Erfi":
        return mp.erfi(z), abs(z * 2 / mp.sqrt(mp.pi) * mp.exp(z * z))
    if head == "FresnelC":
        return mp.fresnelc(z), abs(z * mp.cos(mp.pi * z * z / 2))
    if head == "FresnelS":
        return mp.fresnels(z), abs(z * mp.sin(mp.pi * z * z / 2))
    if head == "Sinc":
        if z == 0:
            return mp.mpf(1), mp.mpf(0)
        return mp.sin(z) / z, abs(mp.cos(z) - mp.sin(z) / z)
    if head == "InverseErf":
        with mp.workdps(40):
            v = fitlib.erfinv_mp(mp.re(z))
        return v, abs(z * mp.sqrt(mp.pi) / 2 * mp.exp(v * v))
    if head == "InverseErfc":
        with mp.workdps(40):
            v = fitlib.erfcinv_mp(mp.re(z))
        return v, abs(z * mp.sqrt(mp.pi) / 2 * mp.exp(v * v))
    raise ValueError(head)


# (head, regime, band label, sampler)
def bands(rng):
    U = lambda a, b: (lambda: float(rng.uniform(a, b)))
    C = lambda a, b, c, d: (lambda: complex(float(rng.uniform(a, b)), float(rng.uniform(c, d))))
    L = lambda a, b: (lambda: float(10 ** rng.uniform(a, b)))
    out = []
    for h in ("Erf", "Erfc", "Erfi"):
        out += [(h, "MR", "bench [-3,3]", U(-3, 3)), (h, "MR", "wide [-6,6]", U(-6, 6)),
                (h, "C", "bench [-2,2]^2", C(-2, 2, -2, 2)), (h, "C", "wide [-5,5]^2", C(-5, 5, -5, 5))]
    for h in ("FresnelC", "FresnelS"):
        out += [(h, "MR", "bench [-3,3]", U(-3, 3)), (h, "MR", "wide [-50,50]", U(-50, 50)),
                (h, "C", "bench [-2,2]^2", C(-2, 2, -2, 2)), (h, "C", "wide [-4,4]^2", C(-4, 4, -4, 4))]
    out += [("Sinc", "MR", "bench [0.2,10]", U(0.2, 10)), ("Sinc", "MR", "wide [-1e3,1e3]", U(-1e3, 1e3)),
            ("Sinc", "C", "bench [0.5,5]x[-3,3]", C(0.5, 5, -3, 3)), ("Sinc", "C", "wide [-20,20]^2", C(-20, 20, -20, 20))]
    out += [("InverseErf", "MR", "bench [-0.95,0.95]", U(-0.95, 0.95)),
            ("InverseErf", "MR", "wide (-1,1)", U(-0.9999999999, 0.9999999999)),
            ("InverseErfc", "MR", "bench [0.05,1.95]", U(0.05, 1.95)),
            ("InverseErfc", "MR", "tail 10^[-300,-1]", L(-300, -1))]
    return out


def lit(v):
    return repr(float(v)) + "`"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=os.path.join(HERE, "..", "..", "Mathilda"))
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--heads", default="")
    ap.add_argument("--seed", type=int, default=20261017)
    a = ap.parse_args()
    want = set(a.heads.split(",")) if a.heads else None
    rng = np.random.default_rng(a.seed)
    rows = []
    for head, reg, lab, samp in bands(rng):
        if want and head not in want:
            continue
        pts = [samp() for _ in range(a.n)]
        src = []
        for p in pts:
            arg = ("Complex[%s, %s]" % (lit(p.real), lit(p.imag))) if isinstance(p, complex) else lit(p)
            src.append("r = %s[%s]; Print[ToString[InputForm[SetPrecision[r, 17]]]];" % (head, arg))
        with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
            f.write("\n".join(src) + "\n")
            fn = f.name
        out = subprocess.run([a.bin, "-file", fn], capture_output=True, text=True).stdout.strip().split("\n")
        os.unlink(fn)
        rel, sc, far, bad = [], [], [], 0
        worst = (0.0, None)
        for p, line in zip(pts, out):
            line = line.strip().replace("*^", "e")
            m = re.fullmatch(r"(-?[0-9.e+\-]+)(?:\s*([+-])\s*([0-9.e+\-]+)\*I)?", line)
            if not m:
                bad += 1
                continue
            got = mp.mpf(m.group(1))
            if m.group(2):
                got = mp.mpc(got, mp.mpf(m.group(3)) * (1 if m.group(2) == "+" else -1))
            z = mp.mpc(p) if isinstance(p, complex) else mp.mpf(p)
            f_, loc = ref_and_scale(head, z)
            e = abs(got - f_)
            r_ = float(e / ulp(abs(f_)))
            s_ = float(e / ulp(max(abs(f_), loc)))
            rel.append(r_); sc.append(s_)
            if abs(f_) >= loc / 4:          # not near a zero of f
                far.append(r_)
            if s_ > worst[0]:
                worst = (s_, p)
        rel = np.array(rel); sc = np.array(sc); far = np.array(far) if far else np.zeros(1)
        rows.append((head, reg, lab, len(rel), bad, np.median(rel), far.max(), rel.max(),
                     np.median(sc), np.percentile(sc, 99), sc.max(), worst[1]))
        print("%-12s %-3s %-22s n=%-5d bad=%-2d | rel med %.2f max(away from zeros) %.2f max %.2f"
              " | scaled med %.2f p99 %.2f max %.2f  (worst scaled at %r)" % rows[-1], flush=True)


if __name__ == "__main__":
    main()
