#!/usr/bin/env python3.11
"""Turn results/<tag>-<engine>.json files into Markdown comparison tables.

    python3.11 report.py <tag> [engine ...]   -> prints Markdown to stdout

Batch groups show total seconds; single groups show seconds per case.  A
cell is "timeout", "wrong" (answer failed verification) or "fail" (crash /
no answer).  The fastest correct engine in each row is in bold.
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(tag):
    res = {}
    for path in sorted(glob.glob(os.path.join(HERE, "results", f"{tag}-*.json"))):
        res.update(json.load(open(path)))
    return res


def fmt(t):
    if t is None:
        return "—"
    if t < 0.001:
        return f"{t * 1e6:.0f} µs"
    if t < 1:
        return f"{t * 1e3:.1f} ms"
    return f"{t:.2f} s"


def main():
    tag = sys.argv[1]
    res = load(tag)
    engines = sys.argv[2:] or sorted({k.split("|")[0] for k in res})
    corpus = json.load(open(os.path.join(HERE, "corpus.json")))
    print(f"| case | " + " | ".join(engines) + " |")
    print("|---|" + "---|" * len(engines))
    for g in corpus["groups"]:
        rows = []
        if g["kind"] == "batch":
            cells = {}
            for e in engines:
                r = res.get(f"{e}|{g['name']}")
                if not r:
                    cells[e] = (None, "")
                elif r["status"] != "ok" or r["total"] is None:
                    cells[e] = (None, r["status"])
                elif r["bad"]:
                    cells[e] = (None, f"{len(r['bad'])}/{r['n']} wrong")
                else:
                    cells[e] = (r["total"], "")
            rows.append((f"{g['name']} ({len(g['cases'])}×, total)", cells))
        else:
            for c in g["cases"]:
                cells = {}
                for e in engines:
                    r = res.get(f"{e}|{g['name']}")
                    x = r["cases"].get(c["id"]) if r else None
                    if not x:
                        cells[e] = (None, "")
                    elif x["status"] == "timeout":
                        cells[e] = (None, "timeout")
                    elif not x["ok"]:
                        cells[e] = (None, "wrong" if x["t"] is not None else "fail")
                    else:
                        cells[e] = (x["t"], "")
                rows.append((c["id"], cells))
        for label, cells in rows:
            ok = [t for t, _ in cells.values() if t is not None]
            best = min(ok) if ok else None
            out = []
            for e in engines:
                t, note = cells[e]
                if t is None:
                    out.append(note or "—")
                else:
                    s = fmt(t)
                    out.append(f"**{s}**" if t == best else s)
            print(f"| {label} | " + " | ".join(out) + " |")


if __name__ == "__main__":
    main()
