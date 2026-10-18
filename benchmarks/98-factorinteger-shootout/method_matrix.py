#!/usr/bin/env python3
"""Method x number-shape sweep for the FactorInteger tutorial's performance page.

Each cell runs in its own Mathilda process (a hung method cannot poison the
next cell): AbsoluteTiming[FactorInteger[n, Method -> m]] under a WALL-clock
process timeout T (TimeConstrained counts CPU time summed over threads).
A cell is OK only if the product of the returned factors is n and every base
passes PrimeQ; PARTIAL if a composite is left; TIMEOUT if the budget ran out.
Writes JSON to argv[1].  Usage (from the repo root, after `make`):
    python3 benchmarks/98-factorinteger-shootout/method_matrix.py \
        benchmarks/98-factorinteger-shootout/results/method-matrix.json
The result is tabulated on site/docs/tutorials/factorinteger-methods/performance.md.
"""
import json, os, subprocess, sys, concurrent.futures as cf

BIN = os.environ.get("MATHILDA_BIN", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Mathilda"))
T = float(os.environ.get("FI_SWEEP_T", "20"))
JOBS = int(os.environ.get("FI_SWEEP_JOBS", "1"))

# Shapes are built inside Mathilda so the page can reproduce them verbatim.
SM = "Product[Prime[i], {i, 2, 18}]"
SHAPES = {
    "u64 semiprime (2×32-bit)": "NextPrime[2^31 + 12345] NextPrime[2^32 - 999]",
    "bal20": "NextPrime[3 10^9] NextPrime[7 10^9]",
    "bal30": "NextPrime[3 10^14] NextPrime[7 10^14]",
    "bal40": "NextPrime[3 10^19] NextPrime[7 10^19]",
    "bal60": "NextPrime[3 10^29] NextPrime[7 10^29]",
    "p10 x p50": "NextPrime[10^9] NextPrime[10^50]",
    "p20 x p80": "NextPrime[10^19] NextPrime[10^79]",
    "close (60d, |p-q|~1e6)": "NextPrime[10^30] NextPrime[10^30 + 10^6]",
    "smooth p-1 (p25 x p40)":
        f"(2 NestWhile[# + 1 &, 1, ! PrimeQ[2 # {SM} + 1] &] {SM} + 1) NextPrime[10^40]",
    "smooth p+1 (p25 x p40)":
        f"(2 NestWhile[# + 1 &, 1, ! PrimeQ[2 # {SM} - 1] &] {SM} - 1) NextPrime[10^40]",
    "RBD-structured (88d)":
        "742245361114170698080188160556986382238609541182118959142360194921632978579427957105431",
    "F7 = 2^128+1": "2^128 + 1",
}
METHODS = ["Automatic", "TrialDivision", "PollardRho", "Fermat", "Lehman",
           "ShanksSquareForms", "CFRAC", "Dixon", "QuadraticSieve",
           "PollardPMinus1", "WilliamsPPlus1", "ECM", "BlakeRationalBaseDescent"]


def cell(shape, method):
    opt = "" if method == "Automatic" else f', Method -> "{method}"'
    lines = [
        f"n = {SHAPES[shape]};",
        f"r = AbsoluteTiming[Quiet[FactorInteger[n{opt}]]];",
        'If[r === $Aborted, "TIMEOUT", '
        '{r[[1]], If[Times @@ (Power @@@ r[[2]]) == n && AllTrue[r[[2, All, 1]], PrimeQ], "OK", "PARTIAL"], '
        'IntegerLength[n]}]',
    ]
    reqs = [json.dumps({"id": i, "expr": l}) for i, l in enumerate(lines, 1)]
    reqs.append(json.dumps({"type": "quit"}))
    try:
        p = subprocess.run([BIN], input="\n".join(reqs) + "\n", capture_output=True,
                           text=True, errors="replace", timeout=T)
    except subprocess.TimeoutExpired:
        return shape, method, "TIMEOUT", None
    out = None
    for raw in p.stdout.splitlines():
        try:
            o = json.loads(raw)
        except ValueError:
            continue
        if o.get("id") == 3 and o.get("type") == "expr":
            out = o.get("payload")
    if out is None:
        return shape, method, "CRASH", None
    if "TIMEOUT" in out:
        return shape, method, "TIMEOUT", None
    body = out.strip("{} ").split(",")
    try:
        t = float(body[0].replace("*^", "e"))
    except ValueError:
        return shape, method, "ERR:" + out[:80], None
    return shape, method, body[1].strip().strip('"'), t


def main():
    cells = [(s, m) for s in SHAPES for m in METHODS]
    res = []
    with cf.ThreadPoolExecutor(JOBS) as ex:
        for shape, method, status, t in ex.map(lambda c: cell(*c), cells):
            print(f"{shape:28s} {method:26s} {status:8s} {t}", flush=True)
            res.append({"shape": shape, "method": method, "status": status, "seconds": t})
    json.dump({"budget_s": T, "shapes": SHAPES, "results": res}, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
