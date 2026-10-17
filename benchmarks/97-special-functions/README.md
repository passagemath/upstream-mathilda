# Special-functions shootout

A self-contained benchmark comparing **every numerically-evaluable special function in
Mathilda** against three external systems, over three numeric regimes, for both **speed**
and **accuracy**. It produces a standalone LaTeX study,
[`special_function_performance.tex`](special_function_performance.tex), with a rendered
PDF.

This is **not** part of the weekly `benchmarks/run_all.py` gap gate (that is a three-engine,
boolean-check, CI-adjacent job). This folder has no two-digit numeric prefix, so
`run_all.py`'s experiment discovery ignores it; it runs on its own driver.

## Engines

| Column | Backend | Regimes |
|---|---|---|
| Mathilda | `./Mathilda -file` | MR, C, AP |
| Mathematica | `wolframscript -file` (same `.m` as Mathilda) | MR, C, AP |
| Maxima | `maxima --very-quiet -b` (`float` / `rectform` / `bfloat`) | MR, C, AP (where it has a kernel) |
| Python | `scipy.special` (machine), `mpmath` (arbitrary precision) | MR, C, AP |

Regimes: **MR** machine real (IEEE double), **C** complex machine, **AP** arbitrary
precision (default 50 digits).

## What it measures

- **Speed, scalar latency** — per-call time in every regime, timed *inside* each engine
  over blocks of *distinct* arguments (so no engine serves a memoised result), min over
  blocks. One process per engine, so start-up is never charged per call.
- **Speed, vectorized throughput** — machine-real only: Mathilda's packed buffer vs
  SciPy's NumPy array over `10^6` elements.
- **Accuracy, agreement-only** — no oracle. For each point, the number of leading
  significant digits the engines mutually agree on; any engine that deviates from the
  consensus is flagged. Doubles as the cross-engine mapping validator.

## Files

| File | Role |
|---|---|
| `registry.py` | the function↔engine spine: per-head WL/Maxima templates, scipy/mpmath callables, domains, regime support |
| `run_shootout.py` | orchestrator: generates `.m`/`.mac`, runs the four engines, parses, computes agreement + timings → `results/<date>.json` |
| `report.py` | `results/*.json` → LaTeX fragments (`tab_*.tex`) + `macros.tex` |
| `special_function_performance.tex` | the standalone article (`\input`s the fragments) |
| `results/<date>.json` | raw per-case data + host/version block |

## Run

```bash
cd <mathilda-repo>
make -j                                         # build ./Mathilda
benchmarks/.venv/bin/python benchmarks/97-special-functions/run_shootout.py --ap-dps 50
benchmarks/.venv/bin/python benchmarks/97-special-functions/report.py
cd benchmarks/97-special-functions && latexmk -pdf special_function_performance.tex
```

Useful flags on `run_shootout.py`:

- `--only Gamma,BesselJ` — restrict to a subset of heads (smoke testing).
- `--regimes MR,C,AP` and `--engines mathilda,wolfram,maxima,python` — restrict axes.
- `--no-throughput` / `--no-exact` — skip those sections.
- `--ap-dps 100` — change the arbitrary-precision digit count.

Engine binaries are located via `MATHILDA_BIN`, `WOLFRAMSCRIPT`, `MAXIMA_BIN`, and the
benchmark venv; an engine not found is simply reported as absent (`—`).

## Adding a function

Add one `SF(...)` record to `NUMERIC` in `registry.py`: the Mathilda/WL call template
(shared with `wolframscript`), the Maxima template (`None` if Maxima lacks it), an mpmath
callable (used for AP and as the SciPy fallback), an optional scipy callable for MR/C, the
regimes it supports, and the per-slot sampling domains. The agreement check will tell you
immediately if a cross-engine convention was mapped wrong (the new head shows up in the
disagreements table).
