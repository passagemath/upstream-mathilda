#!/usr/bin/env python3
"""
check_override.py -- run the numeric oracle on a file of hand-repaired entries.

Usage:
    python3 -I tools/mathar/check_override.py <entries.json> <shard.json> <out.json> [nkernels]

<shard.json> has the overrides.json shape:
    { "M1234": {"lhs": "Integrate[...]", "rhs": ["..."], "conds": ["a > 0"],
                "note": "typo: paper has x^2, intent x^3 (verified)"},
      "M1300": {"skip": "recurrence: integral on both sides", "status": "SKIP-RECURRENCE"} }

The entries in the shard (non-skip) are merged over <entries.json>, and
verify.wls is run on exactly those ids; its JSON verdicts land in <out.json>
and a one-line summary per id is printed.  Overrides are held to the same
oracle as automatic translations: a reading is only admitted when the paper's
value agrees with NIntegrate.
"""
import json
import os
import subprocess
import sys
import tempfile

WOLFRAM = '/Applications/Mathematica.app/Contents/MacOS/wolframscript'


def main():
    entries = json.load(open(sys.argv[1]))
    shard = json.load(open(sys.argv[2]))
    out = sys.argv[3]
    nk = sys.argv[4] if len(sys.argv) > 4 else '2'
    by = {e['id']: e for e in entries}
    ids = []
    for oid, o in shard.items():
        if 'skip' in o:
            continue
        e = dict(by[oid])
        e.update({k: o[k] for k in ('lhs', 'rhs', 'conds') if k in o})
        e.setdefault('conds', [])
        e['status'] = 'ok'
        by[oid] = e
        ids.append(oid)
    if not ids:
        print('nothing to verify')
        return
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as tf:
        json.dump([by[i] for i in ids], tf)
        tmp = tf.name
    here = os.path.dirname(os.path.abspath(__file__))
    subprocess.run([WOLFRAM, '-file', os.path.join(here, 'verify.wls'), tmp, out, nk],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.unlink(tmp)
    for r in json.load(open(out)):
        print(r['id'], r['status'], (r.get('why') or r.get('bad') or '')[:200])


if __name__ == '__main__':
    main()
