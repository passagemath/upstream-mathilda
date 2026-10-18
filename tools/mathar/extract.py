#!/usr/bin/env python3
"""
extract.py -- harvest the integrals of R. J. Mathar, "Yet Another Table of
Integrals" (arXiv:1207.5845) from its TeX source into a JSON list of entries.

Usage:
    python3 -I tools/mathar/extract.py <mathar.tex> <entries.json>

Every display block (\\be..\\ee, equation, multline, gather, eqnarray, align)
and every inline $..$ fragment that contains `\\int` becomes exactly one entry,
so the output accounts for every integral in the paper: an entry is either
translated (`status: "ok"`) or carries a SKIP status and the reason.

Entry fields:
    id         "M<line>"  -- the TeX line of the block (stable across runs)
    line, section, subsection, cite
    tex        the raw block
    kind       "indefinite" | "definite"
    lhs        InputForm of the side holding the integral (may be c*Integrate+d)
    rhs        list of InputForm candidates for the value (other '=' pieces)
    params     free symbols of the integrand (excluding the integration var)
    conds      InputForm conditions translated from the side text (best effort)
    conds_tex  raw condition text
    status     ok | SKIP-PARSE | SKIP-MULTIPLE | SKIP-NO-RHS | SKIP-INLINE | ...
    reason     for SKIP statuses
"""
import json
import re
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tex2m  # noqa: E402
from tex2m import TexError  # noqa: E402

BLOCK_RE = re.compile(
    r'\\be\b(?P<a>.*?)\\ee\b|'
    r'\\begin\{(?P<env>equation\*?|multline\*?|gather\*?|eqnarray\*?|align\*?)\}'
    r'(?P<b>.*?)\\end\{(?P=env)\}', re.S)

SEC_RE = re.compile(r'\\(section|subsection)\*?\{((?:[^{}]|\{[^{}]*\})*)\}')
CITE_RE = re.compile(r'\\cite(?:\[([^\]]*)\])?\{([^}]*)\}')


def line_of(s, pos, _cache={}):
    return s.count('\n', 0, pos) + 1


def trim_tail(tex):
    """Strip trailing punctuation and return (expr, condition-fragments).

    Conditions written inside the display after a top-level comma
    (`..., \\, a>0, \\, n=1,2,\\ldots`) are split off."""
    t = tex.strip()
    t = re.sub(r'(\\[,;!]|\\quad|\\qquad|\s)+$', '', t)
    t = re.sub(r'[.,;]\s*$', '', t)
    t = re.sub(r'(\\[,;!]|\\quad|\\qquad|\s)+$', '', t)
    # split at top-level commas
    parts, cur, depth, i = [], [], 0, 0
    while i < len(t):
        c = t[i]
        if c == '\\':
            m = re.match(r'\\([A-Za-z]+|.|$)', t[i:], re.S)
            tok = m.group(0)
            if tok in ('\\left', '\\right'):
                cur.append(tok)
                i += len(tok)
                continue
            if tok in ('\\quad', '\\qquad') and depth == 0:
                parts.append(''.join(cur))
                cur = []
                i += len(tok)
                continue
            cur.append(tok)
            i += len(tok)
            continue
        if c in '{([':
            depth += 1
        elif c in '})]':
            depth -= 1
        if c == ',' and depth == 0:
            parts.append(''.join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    parts.append(''.join(cur))
    expr = parts[0]
    conds = [p for p in parts[1:] if p.strip() and p.strip() not in ('\\,', '\\;')]
    return expr, conds


def has_relation_op(t):
    return bool(re.search(r'(?<![\\a-zA-Z])[<>]|\\(le|leq|ge|geq|ne|neq|in|notin)\b', t))


def cond_to_m(tex):
    """Translate one condition fragment to InputForm, or None."""
    t = tex.strip().strip('$').strip()
    t = re.sub(r'\\(?:text|mbox|rm)\{[^}]*\}', '', t)
    if not t:
        return None
    t = re.sub(r'\\Re\s*\(([^()]*)\)', r'REPART{\1}', t)
    t = re.sub(r'\\Re\s*([A-Za-z]|\\[a-z]+)', r'REPART{\1}', t)
    t = re.sub(r'\\Im\s*\(([^()]*)\)', r'IMPART{\1}', t)
    t = re.sub(r'\\Im\s*([A-Za-z]|\\[a-z]+)', r'IMPART{\1}', t)
    m = re.match(r'^\s*([A-Za-z]|\\[a-z]+)\s*=\s*([0-9\-]+)\s*,\s*([0-9\-]+)\s*,?\s*(\\ldots|\\dots|\.\.\.)', t)
    if m:
        v = tex2m.translate(m.group(1))
        lo = int(m.group(2))
        return 'Element[%s, Integers] && %s >= %d' % (v, v, lo)
    m = re.match(r'^\s*([A-Za-z]|\\[a-z]+)\s*\\in\s*\\mathbb\s*\{?\s*([NZRC])\s*\}?\s*(\^\s*\{?\s*([+*])\s*\}?)?\s*$', t)
    if m:
        v = tex2m.translate(m.group(1))
        dom = m.group(2)
        sup = m.group(4)
        if dom == 'N':
            return 'Element[%s, Integers] && %s >= %d' % (v, v, 1)
        if dom == 'Z':
            return 'Element[%s, Integers]' % v + (' && %s > 0' % v if sup else '')
        if dom == 'R':
            return 'Element[%s, Reals]' % v + (' && %s > 0' % v if sup else '')
        return None
    if not has_relation_op(t) and 'REPART' not in t:
        return None
    # chain of relations
    toks = re.split(r'(<|>|\\leq?\b|\\geq?\b|\\neq?\b)', t)
    if len(toks) < 3:
        return None
    ops = {'<': '<', '>': '>', '\\le': '<=', '\\leq': '<=', '\\ge': '>=', '\\geq': '>=',
           '\\ne': '!=', '\\neq': '!='}
    out = []
    for i in range(0, len(toks), 2):
        piece = toks[i].strip()
        piece = re.sub(r'REPART\{([^}]*)\}', r'\\Re(\1)', piece)
        piece = re.sub(r'IMPART\{([^}]*)\}', r'\\Im(\1)', piece)
        if not piece:
            return None
        s = tex2m.translate(piece)
        out.append(s)
        if i + 1 < len(toks):
            out.append(ops[toks[i + 1].strip()])
    expr = ' '.join(out)
    return expr.replace('!=', '!=')


def collect_conds(frags):
    res = []
    for f in frags:
        try:
            c = cond_to_m(f)
        except (TexError, KeyError, AttributeError):
            c = None
        if c:
            res.append(c)
    return res


def prose_conds(text):
    """Math fragments in the prose after a block: `for $a>0$, $n\\ge 1$`."""
    frags = []
    for m in re.finditer(r'\$([^$]+)\$', text):
        for piece in re.split(r',\s*(?![^{]*\})', m.group(1)):
            frags.append(piece)
    return frags


def count_ints(t):
    return len(re.findall(r'\\i+nt\b|\\oint\b', t))


def entry_for_block(tex, line, section, subsection, cite, after):
    ent = {'id': 'M%d' % line, 'line': line, 'section': section,
           'subsection': subsection, 'cite': cite, 'tex': tex.strip()}
    body = re.sub(r'(?<!\\)%[^\n]*', '', tex)
    body = re.sub(r'\\label\{[^}]*\}', '', body)
    body = body.replace('\\nonumber', '').replace('\\notag', '')
    # eqnarray / align rows: join (an aligned chain is one relation)
    body = body.replace('&', ' ').replace('\\\\', ' ')
    if re.search(r'\\iint|\\iiint|\\int\s*(?:_\{[^}]*\}|_\S)?\s*(?:\^\{[^}]*\}|\^\S)?\s*\\int', body):
        ent.update(status='SKIP-MULTIPLE', reason='nested/multiple integral')
        return ent
    if re.search(r'\\oint', body):
        ent.update(status='SKIP-CONTOUR', reason='contour integral')
        return ent
    expr_tex, inline_conds = trim_tail(body)
    pieces = tex2m.split_relation(expr_tex)
    if len(pieces) < 2:
        ent.update(status='SKIP-NO-RHS', reason='no top-level =')
        return ent
    int_pieces = [i for i, p in enumerate(pieces) if '\\int' in p]
    if not int_pieces:
        ent.update(status='SKIP-NO-INT', reason='integral not at top level')
        return ent
    li = int_pieces[0]
    lhs_tex = pieces[li]
    if count_ints(lhs_tex) > 1:
        # a sum/product of integrals, or an integral inside an integrand
        ent.update(status='SKIP-MULTIPLE', reason='more than one integral in one side')
        return ent
    ent['lhs_tex'] = lhs_tex
    # one reading of the letter e for the whole entry (Euler's e vs parameter)
    e_param = tex2m.bare_e(tex2m.tokenize(tex2m.preprocess(expr_tex)))
    try:
        tree = tex2m.parse_tex(lhs_tex, e_param)
    except TexError as e:
        ent.update(status='SKIP-PARSE', reason='lhs: %s' % e)
        return ent
    ints = []

    def find_ints(n):
        if isinstance(n, tuple):
            if n[0] == 'int':
                ints.append(n)
            for c in n[1:]:
                if isinstance(c, tuple):
                    find_ints(c)
                elif isinstance(c, list):
                    for x in c:
                        find_ints(x)
    find_ints(tree)
    if len(ints) != 1:
        ent.update(status='SKIP-PARSE', reason='integral count %d' % len(ints))
        return ent
    it = ints[0]
    var = it[2]
    params = sorted(tex2m.free_syms(it[1]) - {var})
    ent['var'] = tex2m.mm_symbol(var)
    ent['kind'] = 'indefinite' if it[3] is None else 'definite'
    int_syms = tex2m.free_syms(tree)
    ctx = {'gamma_const': 'gamma' not in int_syms, 'catalan': 'G' not in int_syms}
    ent['params'] = [tex2m.mm_symbol(p) for p in params]
    try:
        ent['lhs'] = tex2m.emit(tree, {'gamma_const': False, 'catalan': False})
    except TexError as e:
        ent.update(status='SKIP-PARSE', reason='lhs emit: %s' % e)
        return ent
    rhs, errs = [], []
    for j, p in enumerate(pieces):
        if j == li or '\\int' in p:
            continue
        try:
            rhs.append(tex2m.translate(p, ctx, e_param))
        except TexError as e:
            errs.append(str(e))
    ent['rhs'] = rhs
    ent['rhs_errors'] = errs
    ent['conds_tex'] = inline_conds + prose_conds(after)
    ent['conds'] = collect_conds(ent['conds_tex'])
    if not rhs:
        ent.update(status='SKIP-PARSE', reason='rhs: ' + '; '.join(errs))
        return ent
    ent['status'] = 'ok'
    return ent


def main():
    src, out = sys.argv[1], sys.argv[2]
    s = open(src, encoding='utf-8').read()
    # section index
    secs = [(m.start(), m.group(1), m.group(2)) for m in SEC_RE.finditer(s)]

    def sec_at(pos):
        sec = sub = ''
        for p, kind, title in secs:
            if p > pos:
                break
            if kind == 'section':
                sec, sub = title, ''
            else:
                sub = title
        return sec, sub

    cites = [(m.start(), (m.group(1) or '') + '|' + m.group(2)) for m in CITE_RE.finditer(s)]

    def cite_before(pos):
        best = ''
        for p, c in cites:
            if p > pos:
                break
            best = c
        return best

    entries = []
    spans = []
    for m in BLOCK_RE.finditer(s):
        body = m.group('a') if m.group('a') is not None else m.group('b')
        spans.append((m.start(), m.end()))
        if '\\int' not in body:
            continue
        line = line_of(s, m.start())
        sec, sub = sec_at(m.start())
        # prose after the block up to the next block / cite / blank line pair
        tail = s[m.end(): m.end() + 400]
        cut = re.search(r'\\cite|\\be\b|\\begin\{|\\(sub)?section|\n\s*\n', tail)
        after = tail[:cut.start()] if cut else tail
        entries.append(entry_for_block(body, line, sec, sub, cite_before(m.start()), after))
    # inline integrals outside display blocks
    mask = list(s)
    for a, b in spans:
        for k in range(a, b):
            mask[k] = ' '
    rest = ''.join(mask)
    for m in re.finditer(r'\\int', rest):
        line = line_of(s, m.start())
        sec, sub = sec_at(m.start())
        entries.append({'id': 'M%di' % line, 'line': line, 'section': sec, 'subsection': sub,
                        'cite': cite_before(m.start()),
                        'tex': s[max(0, m.start() - 80): m.start() + 160],
                        'status': 'SKIP-INLINE', 'reason': 'integral in running text (definition/remark)'})
    # Hand repairs (tools/mathar/overrides.json): a corrected reading of an
    # entry the translator could not parse or misread, or a reasoned SKIP.
    # Overridden entries go through the same numeric oracle as the rest.
    ov_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'overrides.json')
    if os.path.exists(ov_path):
        ov = json.load(open(ov_path, encoding='utf-8'))
        by_id = {e['id']: e for e in entries}
        for oid, o in ov.items():
            e = by_id.get(oid)
            if e is None:
                print('WARNING: override for unknown id', oid)
                continue
            e['override'] = True
            if 'skip' in o:
                e['status'] = o.get('status', 'SKIP-MANUAL')
                e['reason'] = o['skip']
                continue
            for k in ('lhs', 'rhs', 'conds', 'note'):
                if k in o:
                    e[k] = o[k]
            e.setdefault('conds', [])
            e['status'] = 'ok'
            e.pop('reason', None)
    json.dump(entries, open(out, 'w'), indent=1, ensure_ascii=False)
    from collections import Counter
    c = Counter(e['status'] for e in entries)
    print('entries:', len(entries), dict(c))
    k = Counter(e.get('kind') for e in entries if e['status'] == 'ok')
    print('ok kinds:', dict(k))


if __name__ == '__main__':
    main()
