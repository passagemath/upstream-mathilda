#!/usr/bin/env python3
"""
triage.py -- summarise a Mathar corpus run for the triage report.

Usage:
    python3 -I tools/mathar/triage.py <results.jsonl> <corpus.m>... [--mma mma.json]

<results.jsonl> is what tests/test_mathar_corpus.c writes under
MATHAR_CORPUS_JSON.  The corpus files supply each record's integrand, limits
and section; --mma optionally adds whether Mathematica's Integrate returns a
closed form (tools/mathar/mma_integrate.wls), the "algorithmically solvable"
signal used to rank failures.

Prints: outcome counts per mode and per paper section, a feature table of the
failures (integrand head families x interval shape), and the WRONG / CRASH /
NONNUMERIC list in full.
"""
import json
import re
import sys
from collections import Counter, defaultdict

FAMILIES = [
    ('bessel', r'Bessel[JYIK]|Hankel|SphericalBessel'),
    ('elliptic', r'Elliptic[KEFP]|EllipticPi|JacobiSN'),
    ('gamma/zeta', r'Gamma|PolyGamma|Zeta|HurwitzZeta|BarnesG|Beta\['),
    ('polylog', r'PolyLog|LerchPhi'),
    ('hypergeo', r'Hypergeometric|MeijerG|AppellF1'),
    ('orthopoly', r'LegendreP|LegendreQ|HermiteH|LaguerreL|Chebyshev|Gegenbauer'),
    ('expint/erf', r'ExpIntegral|SinIntegral|CosIntegral|LogIntegral|Erf|Fresnel'),
    ('invtrig', r'Arc(Sin|Cos|Tan|Cot|Sec|Csc)\b|ArcTan\['),
    ('invhyp', r'Arc(Sinh|Cosh|Tanh|Coth|Sech|Csch)'),
    ('log', r'Log\['),
    ('trig', r'\b(Sin|Cos|Tan|Cot|Sec|Csc)\['),
    ('hyp', r'\b(Sinh|Cosh|Tanh|Coth|Sech|Csch)\['),
    ('exp', r'E\^|Exp\['),
    ('algebraic', r'Sqrt\[|\^\(-?\d+/\d+\)|\^\(1/'),
]


def families(f):
    fam = [name for name, rx in FAMILIES if re.search(rx, f)]
    return fam or ['rational']


def interval(lim):
    lim = lim.replace(' ', '')
    if lim == '':
        return 'indef'
    table = {'0,1': '[0,1]', '0,Infinity': '[0,inf)', '-Infinity,Infinity': '(-inf,inf)',
             '0,Pi/2': '[0,pi/2]', '0,Pi': '[0,pi]', '0,2*Pi': '[0,2pi]', '1,Infinity': '[1,inf)',
             '0,Pi/4': '[0,pi/4]', '-1,1': '[-1,1]'}
    return table.get(lim, 'other')


REC_RE = re.compile(r'\{"(M\d+i?)", "(definite|indefinite)",\n  Hold\[Integrate\[(.*?)\]\],\n', re.S)
SEC_RE = re.compile(r'\(\* ---- (.*?) ---- \*\)')


def load_corpus(path):
    txt = open(path).read()
    recs = {}
    pos_secs = [(m.start(), m.group(1)) for m in SEC_RE.finditer(txt)]
    for m in REC_RE.finditer(txt):
        sec = ''
        for p, s in pos_secs:
            if p < m.start():
                sec = s
        body = m.group(3)
        # split integrand / spec at the last top-level ", {" or ", x"
        recs[m.group(1)] = {'kind': m.group(2), 'call': body, 'section': sec.split(' / ')[0]}
    return recs


def main():
    args = sys.argv[1:]
    mma = {}
    if '--mma' in args:
        i = args.index('--mma')
        for r in json.load(open(args[i + 1])):
            mma[r['id']] = r
        args = args[:i] + args[i + 2:]
    res_path, corpora = args[0], args[1:]
    recs = {}
    for c in corpora:
        recs.update(load_corpus(c))
    results = defaultdict(dict)
    for line in open(res_path):
        r = json.loads(line)
        results[r['id']][r['mode']] = r
    by_mode = defaultdict(Counter)
    by_sec = defaultdict(Counter)
    fam_fail = Counter()
    fam_tot = Counter()
    bad = []
    for rid, modes in results.items():
        rec = recs.get(rid, {'kind': '?', 'call': '', 'section': '?'})
        for m, r in modes.items():
            by_mode[(rec['kind'], m)][r['status']] += 1
        st = modes.get('S', {}).get('status')
        by_sec[rec['section']][st] += 1
        key = (rec['kind'], '+'.join(families(rec['call'])))
        fam_tot[key] += 1
        if st != 'CORRECT':
            fam_fail[key] += 1
        for m, r in modes.items():
            if r['status'] in ('WRONG', 'CRASH', 'NONNUMERIC'):
                bad.append((rid, m, r['status'], rec['call'][:150], r['result'][:150]))
    print('== outcomes by kind/mode')
    for k in sorted(by_mode):
        print(' ', k, dict(by_mode[k]))
    print('== symbolic outcome by section')
    for s in sorted(by_sec):
        print('  %-60s %s' % (s[:60], dict(by_sec[s])))
    print('== failures by integrand family (kind, family): fail/total')
    for k, v in fam_fail.most_common(60):
        print('  %-55s %d/%d' % ('%s %s' % k, v, fam_tot[k]))
    if mma:
        solv = Counter()
        for rid, modes in results.items():
            st = modes.get('S', {}).get('status')
            ms = mma.get(rid, {}).get('status', '?')
            solv[(ms, st == 'CORRECT')] += 1
        print('== Mathematica closed form vs Mathilda CORRECT (symbolic):', dict(solv))
    print('== WRONG / CRASH / NONNUMERIC')
    for b in sorted(bad):
        print('  %s %s %-10s %s\n      -> %s' % b)


if __name__ == '__main__':
    main()
