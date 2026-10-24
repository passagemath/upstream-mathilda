#!/usr/bin/env python3
"""Run FeynCalc's upstream unit tests under Mathilda and ratchet the pass set.

    python3 tools/check_feyncalc.py                 # all unit-test directories
    python3 tools/check_feyncalc.py --dir Lorentz --dir Dirac
    python3 tools/check_feyncalc.py --file Contract.test
    python3 tools/check_feyncalc.py --integration   # the IntegrationTests drivers too
    python3 tools/check_feyncalc.py --update-baseline
    python3 tools/check_feyncalc.py --triage        # message / verdict histograms

The vendored package lives in packages/feyncalc/ (local, untracked); without it
this check prints a skip and exits 0.  Each upstream driver
(Tests/<Dir>/<Dir>.mt) is run UNMODIFIED by tools/feyncalc/run_fc_tests.m, a
minimal MUnit`Test.  The driver loads every .test file of its directory and
keeps the groups that match $OnlySubTest (the upstream onlySubTest hook), so a
run is sharded one process per .test file: shards run in parallel, a crash
loses one file's tests rather than a directory's, and FeynCalc loads in about a
second, which is what a shard costs on top of its tests.

Verdicts, per test id: PASS, FAIL (wrong value), MSG (right value, unexpected
messages), ABORT, TIMEOUT, SKIP (typesetting tests: box output is out of
scope), NOTRUN (an id in a .test file that no line reported: its shard crashed
or timed out).  Results are written to packages/feyncalc/harness/results/
<Dir>.tsv.  The ratchet: every id in harness/BASELINE_PASS.txt that was run
must still PASS (exit 1 otherwise); newly passing ids are listed so the
baseline can be raised with --update-baseline."""
import argparse
import collections
import concurrent.futures as cf
import os
import re
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(REPO, "packages", "feyncalc")
TESTS = os.path.join(PKG, "Tests")
HARNESS = os.path.join(PKG, "harness")
BASELINE = os.path.join(HARNESS, "BASELINE_PASS.txt")
RUNNER = os.path.join(REPO, "tools", "feyncalc", "run_fc_tests.m")
UNIT_DIRS = ["Shared", "NonCommAlgebra", "Lorentz", "Dirac", "SUN", "Pauli",
             "Feynman", "Tables", "LoopIntegrals", "ExportImport"]
VERDICTS = ["PASS", "FAIL", "MSG", "ABORT", "TIMEOUT", "SKIP", "NOTRUN"]

GROUP_RE = re.compile(r"^(Tests`[A-Za-z]+`[A-Za-z0-9]+)\s*=", re.M)
ID_RE = re.compile(r'\{\s*"([A-Za-z][A-Za-z0-9]*-ID[0-9]+)"\s*,')


def drivers(dirs, integration):
    """(label, .mt path, directory of its .test files) for each driver."""
    out = []
    for d in dirs:
        mt = os.path.join(TESTS, d, d + ".mt")
        if os.path.isfile(mt):
            out.append((d, mt, os.path.join(TESTS, d)))
        if integration:
            idir = os.path.join(TESTS, d, "IntegrationTests")
            for f in sorted(os.listdir(idir)) if os.path.isdir(idir) else []:
                if f.endswith(".mt"):
                    out.append((d + "/IntegrationTests", os.path.join(idir, f), idir))
    return out


def shards(label, mt, tdir, only_files):
    """One shard per .test file: its groups (for $OnlySubTest) and test ids."""
    out = []
    for f in sorted(os.listdir(tdir)):
        if not f.endswith(".test") or (only_files and f not in only_files):
            continue
        text = open(os.path.join(tdir, f), encoding="utf-8", errors="replace").read()
        groups = GROUP_RE.findall(text)
        if groups:
            out.append({"label": label, "mt": mt, "file": f, "groups": groups,
                        "ids": ID_RE.findall(text)})
    return out


def run_shard(sh, binary, timeout, shard_timeout):
    # `X$ after a backtick: the runner wraps each as .*re.*, and a bare group
    # name would also select every group it is a prefix of.
    only = "|".join("`" + g.split("`")[-1] + "$" for g in sh["groups"])
    env = dict(os.environ, FC_ROOT=PKG, FC_MT=sh["mt"], FC_ONLYSUB=only,
               FC_TIMEOUT=str(timeout))
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-file", RUNNER], cwd=REPO, env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=shard_timeout)
        out, rc = p.stdout.decode("utf-8", "replace"), p.returncode
    except subprocess.TimeoutExpired as e:
        out, rc = (e.stdout or b"").decode("utf-8", "replace"), "timeout"
    results, details, ended = {}, {}, False
    for line in out.splitlines():
        f = line.split("\t")
        if f[0] == "FC" and len(f) >= 3:
            key = f[2]
            n = 2
            while key in results:          # an id repeated across groups
                key, n = f"{f[2]}#{n}", n + 1
            results[key] = (f[1], f[3] if len(f) > 3 else "")
        elif f[0] == "FCD" and len(f) >= 2:
            details.setdefault(f[1], f[2:])
        elif line == "FCEND":
            ended = True
    for i in sh["ids"]:
        if i not in results:
            results[i] = ("NOTRUN", "")
    msgs = collections.Counter(m.group(1) for m in
                               re.finditer(r"^([A-Za-z$`]+::[A-Za-z0-9]+)", out, re.M))
    return {"shard": sh, "results": results, "details": details, "ended": ended,
            "rc": rc, "secs": time.time() - t0, "msgs": msgs,
            "tail": "\n".join(out.splitlines()[-5:])}


def load_baseline():
    base = set()
    if os.path.isfile(BASELINE):
        for line in open(BASELINE, encoding="utf-8"):
            line = line.rstrip("\n")
            if line and not line.startswith("#"):
                base.add(tuple(line.split("\t", 1)))
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", action="append", help="test directory (repeatable)")
    ap.add_argument("--file", action="append", help="only this .test file (repeatable)")
    ap.add_argument("--integration", action="store_true", help="include IntegrationTests drivers")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--timeout", type=int, default=60, help="per-test budget, seconds")
    ap.add_argument("--shard-timeout", type=int, default=3600, help="per-process budget, seconds")
    ap.add_argument("--binary", default=os.path.join(REPO, "Mathilda"))
    ap.add_argument("--update-baseline", action="store_true")
    ap.add_argument("--triage", action="store_true", help="print verdict and message histograms")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    if not os.path.isdir(TESTS):
        print("check-feyncalc: packages/feyncalc not present -- skipped")
        return 0
    dirs = a.dir or UNIT_DIRS
    work = [s for (lab, mt, td) in drivers(dirs, a.integration)
            for s in shards(lab, mt, td, set(a.file or []))]
    if not work:
        print("check-feyncalc: nothing selected")
        return 2
    # Longest first: the slowest shards should not start last.
    work.sort(key=lambda s: -len(s["ids"]))
    t0 = time.time()
    runs = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_shard, s, a.binary, a.timeout, a.shard_timeout) for s in work]
        for fu in cf.as_completed(futs):
            r = fu.result()
            runs.append(r)
            if not a.quiet:
                c = collections.Counter(v for v, _ in r["results"].values())
                bad = "" if r["ended"] else f"  [no FCEND, rc={r['rc']}]"
                print(f"  {r['shard']['label']}/{r['shard']['file']}: "
                      f"{c['PASS']}/{len(r['results'])} in {r['secs']:.0f}s{bad}", flush=True)

    # Results per label.
    os.makedirs(os.path.join(HARNESS, "results"), exist_ok=True)
    per = collections.defaultdict(dict)
    for r in runs:
        for i, (v, secs) in r["results"].items():
            d = r["details"].get(i, ["", "", ""])
            per[r["shard"]["label"]][i] = (v, secs, r["shard"]["file"], *d[:3])
    for lab, res in per.items():
        path = os.path.join(HARNESS, "results", lab.replace("/", "_") + ".tsv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("id\tverdict\tseconds\tfile\tgot\texpected\tmessages\n")
            for i in sorted(res):
                fh.write(i + "\t" + "\t".join(str(x) for x in res[i]) + "\n")

    # Summary.
    print(f"\n{'directory':<28}" + "".join(f"{v:>8}" for v in VERDICTS) + f"{'total':>8}{'pass%':>8}")
    tot = collections.Counter()
    for lab in sorted(per):
        c = collections.Counter(v[0] for v in per[lab].values())
        tot.update(c)
        n = sum(c.values())
        print(f"{lab:<28}" + "".join(f"{c[v]:>8}" for v in VERDICTS) +
              f"{n:>8}{100.0 * c['PASS'] / max(1, n - c['SKIP']):>7.1f}%")
    n = sum(tot.values())
    print(f"{'ALL':<28}" + "".join(f"{tot[v]:>8}" for v in VERDICTS) +
          f"{n:>8}{100.0 * tot['PASS'] / max(1, n - tot['SKIP']):>7.1f}%")
    print(f"wall {time.time() - t0:.0f}s, {len(runs)} shards, {a.jobs} jobs")

    if a.triage:
        msgs = collections.Counter()
        for r in runs:
            msgs.update(r["msgs"])
        print("\nmessages printed (all shards):")
        for m, k in msgs.most_common(40):
            print(f"  {k:>6}  {m}")
        crashed = [r for r in runs if not r["ended"]]
        if crashed:
            print("\nshards without FCEND:")
            for r in crashed:
                print(f"  {r['shard']['label']}/{r['shard']['file']} rc={r['rc']}: {r['tail'][-300:]}")

    # Ratchet.
    now_pass = {(lab, i) for lab, res in per.items() for i, v in res.items() if v[0] == "PASS"}
    ran = {(lab, i) for lab, res in per.items() for i in res}
    base = load_baseline()
    regressed = sorted(k for k in base if k in ran and k not in now_pass)
    gained = sorted(now_pass - base)
    if gained:
        print(f"\n{len(gained)} newly passing (not in the baseline)" +
              ("" if a.update_baseline else "; raise it with --update-baseline"))
    if a.update_baseline:
        keep = {k for k in base if k not in ran}
        os.makedirs(HARNESS, exist_ok=True)
        with open(BASELINE, "w", encoding="utf-8") as fh:
            fh.write("# FeynCalc unit tests that pass under Mathilda (tools/check_feyncalc.py)\n")
            for lab, i in sorted(keep | now_pass):
                fh.write(f"{lab}\t{i}\n")
        print(f"baseline: {len(keep | now_pass)} passing ids written to {os.path.relpath(BASELINE, REPO)}")
        return 0
    if regressed:
        print(f"\nREGRESSED: {len(regressed)} baseline ids no longer pass:")
        for lab, i in regressed[:50]:
            v = per[lab].get(i, ("?",))[0]
            print(f"  {lab}\t{i}\t{v}")
        return 1
    print("\nratchet OK" + ("" if base else " (no baseline yet)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
