#!/usr/bin/env python3
"""
build_corpus.py -- emit the Mathar Integrate corpora from oracle-verified entries.

Usage:
    python3 -I tools/mathar/build_corpus.py <entries.json> <verified.json>[,<more.json>...] \
        <tex sha256> <out-dir> [<mathilda names file>]

Inputs:
    entries.json   from extract.py (every integral of the paper, ok or SKIP)
    verified*.json from verify.wls; later files override earlier ones by id
    names file     optional: output of Mathilda `Names["*"]`, one per line; an
                   entry whose integrand uses a head Mathilda does not define is
                   recorded as SKIP-UNSUPPORTED-HEAD instead of entering the corpus

Outputs (in <out-dir>):
    MatharDefiniteIntegralsCorpus.m
    MatharIndefiniteIntegralsCorpus.m
    MatharSkipped.md      -- every entry not in a corpus, with its reason

Record shape (see tests/test_mathar_corpus.c):
    {id, kind, Hold[symbolic], Hold[instance], {inst...}, {ref...}, Hold[f], var}
"""
import json
import re
import sys
from collections import Counter, defaultdict

INT_NAMES = {'n', 'm', 'k', 'l', 'j', 'N'}
# Heads introduced by the translator or by Mathematica's InputForm that are
# not functions: skip them in the unsupported-head scan.
NON_HEADS = {'List', 'Rule', 'Hold', 'Integrate', 'Sqrt', 'Times', 'Plus', 'Power'}


def plain_number(s):
    """Mathematica InputForm number (with `prec marks) -> plain decimal text."""
    s = re.sub(r'`[0-9.]*', '', s)
    s = re.sub(r'\*\^', '*10^', s)

    def trim(m):
        txt = m.group(0)
        if '.' not in txt:
            return txt
        intpart, frac = txt.split('.')
        digits = len(intpart.lstrip('-0')) + len(frac)
        if digits > 16:
            keep = max(1, 16 - len(intpart.lstrip('-0')))
            frac = frac[:keep]
        return intpart + '.' + frac
    return re.sub(r'-?\d+\.\d*', trim, s)


def heads_of(expr_text):
    return set(re.findall(r'([A-Za-z][A-Za-z0-9]*)\[', expr_text))


def parse_inst(s):
    """'{a -> 3/8, n -> 2}' -> list of (name, value-text)."""
    s = s.strip()
    if s in ('{}', ''):
        return []
    body = s[1:-1]
    out = []
    for part in re.split(r',\s*(?=[A-Za-z][A-Za-z0-9]*\s*->)', body):
        k, v = part.split('->', 1)
        out.append((k.strip(), v.strip()))
    return out


def is_int_text(v):
    return re.fullmatch(r'-?\d+', v) is not None


def assumptions_for(ent, ver):
    conds = list(ent.get('conds') or [])
    params = ver.get('params') or []
    mentioned = set()
    for c in conds:
        mentioned |= set(re.findall(r'[A-Za-z][A-Za-z0-9]*', c))
    insts = [dict(parse_inst(g['inst'])) for g in ver.get('good', [])]
    for p in params:
        vals = [d[p] for d in insts if p in d]
        if not vals:
            continue
        if p in mentioned:
            continue
        if all(is_int_text(v) for v in vals) and p in INT_NAMES:
            lo = min(int(v) for v in vals)
            conds.append('Element[%s, Integers] && %s >= %d' % (p, p, min(lo, 1) if lo >= 1 else 0))
        elif all(not v.startswith('-') for v in vals):
            conds.append('%s > 0' % p)
    return ' && '.join('(%s)' % c for c in conds) if conds else ''


def main():
    entries = json.load(open(sys.argv[1]))
    ver = {}
    for path in sys.argv[2].split(','):
        for r in json.load(open(path)):
            ver[r['id']] = r
    tex_sha = sys.argv[3]
    outdir = sys.argv[4]
    names = None
    if len(sys.argv) > 5:
        names = set(l.strip() for l in open(sys.argv[5]) if l.strip())

    recs = {'definite': [], 'indefinite': []}
    skipped = []
    for e in entries:
        st = e['status']
        if st != 'ok':
            skipped.append((e, st, e.get('reason', '')))
            continue
        v = ver.get(e['id'])
        if v is None:
            skipped.append((e, 'SKIP-UNVERIFIED', 'not run through the oracle'))
            continue
        if v['status'] not in ('VERIFIED', 'VERIFIED1'):
            why = v.get('why') or v.get('bad') or ''
            skipped.append((e, 'SKIP-ORACLE-' + v['status'], why[:200]))
            continue
        f = v['integrand']
        if names is not None:
            missing = sorted(h for h in heads_of(f) if h not in names and h not in NON_HEADS)
            if missing:
                skipped.append((e, 'SKIP-UNSUPPORTED-HEAD', ', '.join(missing)))
                continue
        kind = v['kind']
        var = v['var']
        good = v['good']
        insts = [g['inst'] for g in good]
        refs = [plain_number(g['ref']) for g in good] if kind == 'definite' else []
        assum = assumptions_for(e, v)
        spec = var if kind == 'indefinite' else '{%s, %s}' % (var, v['limits'].strip()[1:-1])
        opt = ', Assumptions -> %s' % assum if assum else ''
        sym = 'Integrate[%s, %s%s]' % (f, spec, opt)
        # the numeric instance, substituted by Mathematica (instance.wls)
        if v.get('inst_integrand'):
            ispec = var if kind == 'indefinite' else '{%s, %s}' % (var, v['inst_limits'].strip()[1:-1])
            num = 'Integrate[%s, %s]' % (v['inst_integrand'], ispec)
        else:
            num = 'Integrate[%s, %s]' % (f, spec)
        rec = '{"%s", "%s",\n  Hold[%s],\n  Hold[%s],\n  {%s},\n  {%s},\n  Hold[%s], %s}' % (
            e['id'], kind, sym, num, ', '.join(insts), ', '.join(refs), f, var)
        recs[kind].append((e, rec, v))

    head = {}
    for kind in recs:
        head[kind] = (
            '(* Mathar%sIntegralsCorpus.m -- generated by tools/mathar/build_corpus.py.\n'
            '   DO NOT EDIT BY HAND: fix tools/mathar/overrides.json and regenerate.\n\n'
            '   Source: R. J. Mathar, "Yet Another Table of Integrals", arXiv:1207.5845\n'
            '   (TeX source sha256 %s).  Each record is one display equation of the\n'
            '   paper, id = M<TeX line>.  Every reference value was produced by\n'
            '   tools/mathar/verify.wls (Mathematica): the paper\'s value had to agree\n'
            '   with NIntegrate (no convergence messages) at two parameter\n'
            '   instantiations before the entry was admitted.  Entries the oracle\n'
            '   rejected (paper typos, ambiguous notation) are listed with reasons in\n'
            '   tests/MatharSkipped.md; corrected readings are in overrides.json.\n\n'
            '   Record: {id, kind, Hold[symbolic], Hold[instance], {inst..}, {ref..},\n'
            '            Hold[integrand], var}   -- runner: tests/test_mathar_corpus.c\n'
            '*)\n' % (kind.capitalize(), tex_sha))
    for kind, lst in recs.items():
        lst.sort(key=lambda t: t[0]['line'])
        body = []
        cur = None
        for e, rec, v in lst:
            sec = '%s / %s' % (e['section'], e['subsection'])
            if sec != cur:
                body.append('(* ---- %s ---- *)' % sec.replace('*)', '* )'))
                cur = sec
            body.append(rec)
        txt = head[kind] + '{\n' + '\n'.join(_join(body)) + '\n}\n'
        open('%s/Mathar%sIntegralsCorpus.m' % (outdir, kind.capitalize()), 'w').write(txt)

    # skipped report
    by = defaultdict(list)
    for e, st, why in skipped:
        by[st].append((e, why))
    lines = ['# Mathar integrals not in the corpus', '',
             'Generated by `tools/mathar/build_corpus.py`. Every `\\int` of arXiv:1207.5845 is',
             'either a corpus record or listed here with the reason it is not one.', '',
             '| Status | Count |', '|---|---|']
    for st in sorted(by):
        lines.append('| %s | %d |' % (st, len(by[st])))
    lines.append('')
    for st in sorted(by):
        lines.append('## %s' % st)
        lines.append('')
        for e, why in sorted(by[st], key=lambda t: t[0]['line']):
            tex = ' '.join(e['tex'].split())[:160].replace('|', '\\|')
            lines.append('- `%s` (%s): %s — `%s`' % (e['id'], e.get('cite', ''), why.replace('|', '\\|')[:200], tex))
        lines.append('')
    open('%s/MatharSkipped.md' % outdir, 'w').write('\n'.join(lines))
    print('definite:', len(recs['definite']), 'indefinite:', len(recs['indefinite']),
          'skipped:', len(skipped), dict(Counter(st for _, st, _ in skipped)))


def _join(body):
    """Insert commas only between records (comments carry none)."""
    out = []
    recs_idx = [i for i, b in enumerate(body) if not b.startswith('(*')]
    last = recs_idx[-1] if recs_idx else -1
    for i, b in enumerate(body):
        out.append(b)
    # commas: every record except the last gets a trailing comma
    res = []
    for i, b in enumerate(out):
        if b.startswith('(*'):
            res.append(b)
        else:
            res.append(b + (',' if i != last else ''))
    return res


if __name__ == '__main__':
    main()
