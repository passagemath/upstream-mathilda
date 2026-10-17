#!/usr/bin/env python3.11
"""FactorInteger shootout harness.

For every (engine, group) pair, emit the engine's script, run it with a
timeout, parse the factorisations it prints and verify each one independently
(product of prime powers == n, every base passes sympy.isprime).  A *batch*
group runs in one process and is timed as a whole inside the engine; a
*single* group runs one process per case, so a timeout or crash is recorded
against that case alone.  Timing is always taken inside the engine, so process
start-up is never charged.

    python3.11 run.py [--engines mathilda,pari,...] [--groups qs,ecm]
                      [--timeout 120] [--tag baseline-v0450] [--method M]

--method passes a Method to the engines that can select one (Mathilda,
Mathematica, PARI flags, FLINT); others are run with their default.  Results
go to results/<tag>.json; report.py turns them into tables.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from math import prod

from sympy import isprime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MATHILDA = os.environ.get("MATHILDA", os.path.join(ROOT, "Mathilda"))
MATHILDA450 = os.environ.get("MATHILDA450", "")   # a v0.450 build, for the "before" column
WOLFRAM = "/Applications/Mathematica.app/Contents/MacOS/wolframscript"
PY = "/usr/local/opt/python@3.11/bin/python3.11"
FLINT = os.path.join(HERE, "drivers", "flint_factor")

# --------------------------------------------------------------------------
# Script emitters.  Each returns (filename-suffix, script text).  The output
# protocol every engine prints is:
#   batch : "#TOTAL\t<seconds>" then "<id>\t\t<p e p e ...>" per case
#   single: "<id>\t<seconds>\t<p e p e ...>"
# --------------------------------------------------------------------------


METHOD_WL = {   # Mathilda method specs needing suboptions for a fair comparison
    "TrialDivision": '{"TrialDivision", "Bound" -> 100000}',
}


def wl_script(cases, batch, method=None, engine="mathilda"):
    """Mathilda and Mathematica share one Wolfram-language script.  Mathematica
    has no Method option for FactorInteger, so it always runs its own driver."""
    if engine == "mathematica" or not method:
        opt = ""
    else:
        opt = f', Method -> {METHOD_WL.get(method, chr(34) + method + chr(34))}'

    lst = ",\n".join(f'{{"{c["id"]}", {c["n"]}}}' for c in cases)
    fmt = 'StringRiffle[ToString /@ Flatten[#], " "]&'
    if batch:
        body = (f"cases = {{\n{lst}\n}};\n"
                f"ns = cases[[All, 2]];\n"
                f"t = First[AbsoluteTiming[rs = Map[FactorInteger[#{opt}]&, ns];]];\n"
                f'Print["#TOTAL\\t", ToString[t, InputForm]];\n'
                f'Do[Print[cases[[i, 1]], "\\t\\t", ({fmt})[rs[[i]]]], {{i, Length[cases]}}];\n')
    else:
        body = (f"cases = {{\n{lst}\n}};\n"
                f"Do[{{t, r}} = AbsoluteTiming[FactorInteger[c[[2]]{opt}]];\n"
                f'  Print[c[[1]], "\\t", ToString[t, InputForm], "\\t", ({fmt})[r]], {{c, cases}}];\n')
    return ".m", body


# factorint flags: 1 no MPQS, 2 skip first-stage ECM, 4 no rho/SQUFOF, 8 no final ECM.
PARI_FLAGS = {None: None, "QuadraticSieve": 4 | 2 | 8, "ECM": 1, "PollardRho": 1 | 2 | 8,
              "ShanksSquareForms": 1 | 2 | 8}


def gp_script(cases, batch, method=None, engine="pari"):
    flag = PARI_FLAGS.get(method)
    call = f"factorint(n, {flag})" if flag is not None else "factor(n)"
    if method == "TrialDivision":
        call = "factor(n, 100000)"
    # gp reads its input line by line, so the case vector must be one line.
    lst = ", ".join(f'["{c["id"]}", {c["n"]}]' for c in cases)
    s = ["default(parisizemax, 2^32);", f"cases = [{lst}];",
         "fl(f) = my(v = List()); for(i=1, #f~, listput(v, Str(f[i,1])); listput(v, Str(f[i,2]))); strjoin(Vec(v), \" \");"]
    if batch:
        s += ["my(t0 = getabstime()); rs = vector(#cases, i, my(n = cases[i][2]); " + call + "); "
              "print(\"#TOTAL\\t\", (getabstime() - t0) / 1000.);",
              "for(i=1, #cases, print(cases[i][1], \"\\t\\t\", fl(rs[i])));"]
    else:
        s += ["for(i=1, #cases, my(n = cases[i][2], t0 = getabstime(), f = " + call + "); "
              "print(cases[i][1], \"\\t\", (getabstime() - t0) / 1000., \"\\t\", fl(f)));"]
    s.append("quit;")
    return ".gp", "\n".join(s) + "\n"


def maxima_script(cases, batch, method=None, engine="maxima"):
    lst = ",\n".join(f'["{c["id"]}", {c["n"]}]' for c in cases)
    s = [f"cases : [\n{lst}\n]$",
         "fl(f) := simplode(map(string, flatten(f)), \" \")$"]
    if batch:
        s += ["t0 : elapsed_real_time()$",
              "rs : map(lambda([c], ifactors(second(c))), cases)$",
              "printf(true, \"#TOTAL~a~a~%\", ascii(9), elapsed_real_time() - t0)$",
              "for i : 1 thru length(cases) do printf(true, \"~a~a~a~a~%\", first(cases[i]), ascii(9), ascii(9), fl(rs[i]))$"]
    else:
        s += ["for c in cases do block([t0 : elapsed_real_time(), f], f : ifactors(second(c)), "
              "printf(true, \"~a~a~a~a~a~%\", first(c), ascii(9), elapsed_real_time() - t0, ascii(9), fl(f)))$"]
    return ".mac", "\n".join(s) + "\n"


SYMPY_CALL = {
    None: "factorint(n)",
    "PollardRho": "split(n, pollard_rho(n))",
    "PollardPMinus1": "split(n, pollard_pm1(n, B=10**6))",
    "ECM": "dict.fromkeys(sorted(ecm(n)), 1)",
    "TrialDivision": "factorint(n, limit=10**5, use_rho=False, use_pm1=False, use_ecm=False)",
}


def sympy_script(cases, batch, method=None, engine="sympy"):
    data = json.dumps([[c["id"], c["n"]] for c in cases])
    call = SYMPY_CALL.get(method, "factorint(n)")
    return ".py", f"""
import json, time
from sympy import factorint
from sympy.ntheory import pollard_rho, pollard_pm1, ecm
def split(n, d):
    return {{n: 1}} if not d else {{d: 1, n // d: 1}}
def factorint_m(n):
    return {call}
cases = [(i, int(n)) for i, n in json.loads('''{data}''')]
fl = lambda f: " ".join(f"{{p}} {{e}}" for p, e in sorted(f.items()))
if {batch!r}:
    t0 = time.perf_counter(); rs = [factorint_m(n) for _, n in cases]
    print("#TOTAL\\t%.9f" % (time.perf_counter() - t0))
    for (i, _), f in zip(cases, rs): print(i + "\\t\\t" + fl(f))
else:
    for i, n in cases:
        t0 = time.perf_counter(); f = factorint_m(n); t = time.perf_counter() - t0
        print(f"{{i}}\\t{{t:.9f}}\\t{{fl(f)}}", flush=True)
"""


def flint_input(cases, batch, method=None, engine="flint"):
    return ".txt", "".join(f'{c["id"]} {c["n"]}\n' for c in cases)


FLINT_METHODS = {None: "auto", "QuadraticSieve": "qs", "ShanksSquareForms": "squfof",
                 "Lehman": "lehman", "WilliamsPPlus1": "pp1"}

# GMP-ECM command line, one process per case (its natural drivers: ECM with
# B1 climbing via -I, p-1 / p+1 at B1 = 10^6).
GMPECM_ARGS = {"ECM": ["-one", "-I", "1", "-c", "0", "2000"],
               "PollardPMinus1": ["-pm1", "1e6"],
               "WilliamsPPlus1": ["-pp1", "-c", "3", "1e6"]}


def gmpecm_input(cases, batch, method=None, engine="gmpecm"):
    return ".txt", "".join(c["n"] + "\n" for c in cases)


ENGINES = {
    "mathilda":    (wl_script,     lambda p, b, m: [MATHILDA, "-file", p]),
    "mathilda450": (wl_script,     lambda p, b, m: [MATHILDA450, "-file", p]),
    "mathematica": (wl_script,     lambda p, b, m: [WOLFRAM, "-file", p]),
    "pari":        (gp_script,     lambda p, b, m: ["gp", "-q", "-f", p]),
    "maxima":      (maxima_script, lambda p, b, m: ["maxima", "--very-quiet", "-b", p]),
    "sympy":       (sympy_script,  lambda p, b, m: [PY, p]),
    "flint":       (flint_input,   None),
    "gmpecm":      (gmpecm_input,  None),
}


def run_engine(engine, cases, batch, timeout, method):
    emit, cmd = ENGINES[engine]
    suffix, text = emit(cases, batch, method, engine)
    with tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False, dir=os.environ.get("TMPDIR")) as f:
        f.write(text)
        path = f.name
    try:
        if engine == "gmpecm":
            return run_gmpecm(cases, timeout, method)
        if engine == "flint":
            argv = [FLINT, "batch" if batch else "single", FLINT_METHODS.get(method, "auto")]
            stdin = open(path)
        else:
            argv = cmd(path, batch, method)
            stdin = subprocess.DEVNULL
        t0 = time.time()
        try:
            p = subprocess.run(argv, stdin=stdin, capture_output=True, text=True, timeout=timeout)
            out, status = p.stdout, ("ok" if p.returncode == 0 else f"exit{p.returncode}")
        except subprocess.TimeoutExpired as e:
            out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
            status = "timeout"
        wall = time.time() - t0
    finally:
        os.unlink(path)
    return parse(out), status, wall


def run_gmpecm(cases, timeout, method):
    """Factor each case with the ecm CLI, re-running on the cofactor until it
    is prime (the CLI reports one split per run).  Timed from outside."""
    from sympy import isprime as _isp
    rows, total, status = {}, 0.0, "ok"
    for c in cases:
        n = int(c["n"])
        stack, primes = [n], []
        t0 = time.time()
        try:
            while stack:
                m = stack.pop()
                if _isp(m):
                    primes.append(m)
                    continue
                p = subprocess.run(["ecm", "-q"] + GMPECM_ARGS.get(method, GMPECM_ARGS["ECM"]),
                                   input=f"{m}\n", capture_output=True, text=True,
                                   timeout=max(1, timeout - (time.time() - t0)))
                toks = [int(t) for t in p.stdout.split() if t.isdigit()]
                fs = [t for t in toks if 1 < t < m and m % t == 0]
                if not fs:
                    primes.append(m)        # no split: reported as-is (fails verification)
                    break
                d = fs[0]
                stack += [d, m // d]
        except subprocess.TimeoutExpired:
            status = "timeout"
        t = time.time() - t0
        total += t
        txt = " ".join(f"{q} 1" for q in sorted(primes)) if status == "ok" else ""
        rows[c["id"]] = (t, txt)
    return (total, rows), status, total


LINE = re.compile(r"^(\S+)\t([^\t]*)\t(.*)$")


def parse(out):
    total, rows = None, {}
    for line in out.splitlines():
        line = line.rstrip()
        if line.startswith("#TOTAL"):
            total = float(line.split("\t")[1].replace("`", "").split("*^")[0]) * (
                10 ** int(line.split("*^")[1]) if "*^" in line else 1)
            continue
        m = LINE.match(line)
        if m:
            t = m.group(2).strip()
            if "*^" in t:
                a, b = t.split("*^")
                t = float(a.rstrip("`")) * 10 ** int(b)
            rows[m.group(1)] = (float(str(t).rstrip("`")) if t else None, m.group(3).strip())
    return total, rows


def verify(n, text):
    toks = text.split()
    if not toks or len(toks) % 2:
        return False
    pairs = [(int(toks[i]), int(toks[i + 1])) for i in range(0, len(toks), 2)]
    if prod(p ** e for p, e in pairs) != n:
        return False
    return all(isprime(p) for p, _ in pairs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engines", default="mathilda,mathematica,pari,maxima,sympy,flint")
    ap.add_argument("--groups", default="")
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--tag", default="run")
    ap.add_argument("--method", default=None)
    a = ap.parse_args()
    corpus = json.load(open(os.path.join(HERE, "corpus.json")))
    want = set(a.groups.split(",")) if a.groups else None
    outpath = os.path.join(HERE, "results", a.tag + ".json")
    results = json.load(open(outpath)) if os.path.exists(outpath) else {}
    for g in corpus["groups"]:
        if want and g["name"] not in want:
            continue
        for eng in a.engines.split(","):
            key = f"{eng}|{g['name']}"
            if g["kind"] == "batch":
                (total, rows), status, wall = run_engine(eng, g["cases"], True, a.timeout, a.method)
                bad = [c["id"] for c in g["cases"]
                       if c["id"] not in rows or not verify(int(c["n"]), rows[c["id"]][1])]
                results[key] = {"kind": "batch", "status": status, "total": total,
                                "n": len(g["cases"]), "bad": bad, "wall": wall}
                print(f"{eng:12s} {g['name']:8s} {status:8s} total={total} bad={len(bad)}", flush=True)
            else:
                per = {}
                for c in g["cases"]:
                    (_, rows), status, wall = run_engine(eng, [c], False, a.timeout, a.method)
                    t, txt = rows.get(c["id"], (None, ""))
                    ok = bool(txt) and verify(int(c["n"]), txt)
                    per[c["id"]] = {"status": status, "t": t, "ok": ok, "wall": wall}
                    print(f"{eng:12s} {c['id']:18s} {status:8s} t={t} ok={ok}", flush=True)
                results[key] = {"kind": "single", "cases": per}
            json.dump(results, open(outpath, "w"), indent=1)


if __name__ == "__main__":
    main()
