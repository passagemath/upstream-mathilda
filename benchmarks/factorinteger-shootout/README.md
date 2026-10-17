# FactorInteger shootout

Measures Mathilda's `FactorInteger` against Mathematica (`wolframscript`), PARI/GP
(`factor` / `factorint`), Maxima (`ifactors`), SymPy (`factorint`) and raw FLINT
(`fmpz_factor`, `qsieve_factor`). The corpus is generated with fixed seeds, and every
engine's answer is verified independently (the product equals n, and every base passes
`sympy.isprime`).

    python3.11 gen_corpus.py                       # writes corpus.json (seeded)
    cc -O2 -o drivers/flint_factor drivers/flint_factor.c $(pkg-config --cflags --libs flint)
    python3.11 run.py --tag <name> [--engines mathilda,pari,...] [--groups qs,ecm] [--timeout 120]
    python3.11 verify_methods.py --count 100       # randomised correctness sweep, every method

- **Batch groups** (`u32`, `u63`, `semi63`, `semi96`, `rand40`, `smooth`) run in one process
  per engine and are timed as a whole inside the engine.
- **Single groups** (`qs`, `ecm`, `pm1`, `pp1`, `fermat`, `power`, `special`) run one process
  per case, so a timeout or crash is recorded for that case alone.
- Process start-up is never included in the timings.

`results/baseline-v0450-*.json` is the state before the rewrite (Mathilda v0.450).
The full report is `REPORT.md`.
