#!/usr/bin/env python3
"""
tex2m.py -- a small LaTeX-math -> Mathematica/Mathilda InputForm translator.

Built for mining R. J. Mathar, "Yet Another Table of Integrals" (arXiv:1207.5845),
whose display equations are hand-written TeX in the classical table style:
implicit multiplication, un-parenthesised function arguments (`\\sin 2x`,
`\\log\\frac{a}{b}`), powers on function names (`\\sin^2 x`, `\\Gamma^3(1+j)`),
subscripted special functions (`J_\\nu(x)`, `{}_2F_1(a,b;c;z)`), and the integral
differential written anywhere in the integrand (`\\frac{x\\,dx}{1+x^2}`,
`\\int dk\\, J_0(k)`, trailing `dx`).

Mathematica's own `ToExpression[..., TeXForm]` cannot read these (it fails on
every `\\int` and reads `\\log` as a symbol), so this is a dedicated tokenizer +
Pratt parser that produces an AST and prints it in InputForm, which both
Mathematica and Mathilda parse.

Translation is deliberately *conservative*: anything it does not understand
raises `TexError` with the offending fragment, and the caller records the entry
as SKIP-PARSE rather than guessing.  Correctness of every translated entry is
established downstream by a numeric oracle (verify.wls), never assumed.

API:
    translate(tex)         -> InputForm string of a single expression
    split_relation(tex)    -> list of top-level '='-separated TeX pieces
"""
import re

class TexError(Exception):
    pass


VAR_NAMES = {'x', 't', 'u', 'v', 'y', 'z', 's', 'theta', 'phi', 'k', 'r', 'w', 'q'}

# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------
GREEK = {
    'alpha': 'alpha', 'beta': 'beta', 'gamma': 'gamma', 'delta': 'delta',
    'epsilon': 'epsilon', 'varepsilon': 'epsilon', 'zeta': 'zeta', 'eta': 'eta',
    'theta': 'theta', 'vartheta': 'theta', 'iota': 'iota', 'kappa': 'kappa',
    'lambda': 'lambda', 'mu': 'mu', 'nu': 'nu', 'xi': 'xi', 'rho': 'rho',
    'varrho': 'rho', 'sigma': 'sigma', 'tau': 'tau', 'upsilon': 'upsilon',
    'phi': 'phi', 'varphi': 'phi', 'chi': 'chi', 'psi': 'psi', 'omega': 'omega',
    'Delta': 'CapDelta', 'Theta': 'CapTheta', 'Lambda': 'CapLambda',
    'Xi': 'CapXi', 'Sigma': 'CapSigma', 'Phi': 'CapPhi', 'Psi': 'CapPsi',
    'Omega': 'CapOmega', 'ell': 'ell',
}

# Elementary functions written as TeX operators.  Value: (head, kind) where kind
# 'f' = ordinary one-argument function.
ELEM = {
    'sin': 'Sin', 'cos': 'Cos', 'tan': 'Tan', 'cot': 'Cot', 'sec': 'Sec',
    'csc': 'Csc', 'cosec': 'Csc', 'sinh': 'Sinh', 'cosh': 'Cosh', 'tanh': 'Tanh',
    'coth': 'Coth', 'sech': 'Sech', 'csch': 'Csch', 'sh': 'Sinh', 'ch': 'Cosh',
    'arcsin': 'ArcSin', 'arccos': 'ArcCos', 'arctan': 'ArcTan',
    'arccot': 'ArcCot', 'arcsec': 'ArcSec', 'arccsc': 'ArcCsc',
    'arcsinh': 'ArcSinh', 'arccosh': 'ArcCosh', 'arctanh': 'ArcTanh',
    'arccoth': 'ArcCoth', 'Arsh': 'ArcSinh', 'arsinh': 'ArcSinh',
    'arcosh': 'ArcCosh', 'artanh': 'ArcTanh', 'Arch': 'ArcCosh', 'Arth': 'ArcTanh',
    'ln': 'Log', 'log': 'Log', 'exp': 'Exp', 'erf': 'Erf', 'erfc': 'Erfc',
    'erfi': 'Erfi', 'Ei': 'ExpIntegralEi', 'Si': 'SinIntegral',
    'Ci': 'CosIntegral', 'Shi': 'SinhIntegral', 'Chi': 'CoshIntegral',
    'li': 'LogIntegral', 'sgn': 'Sign', 'sign': 'Sign', 'Re': 'Re', 'Im': 'Im',
    'arg': 'Arg',
}

# Named functions written as plain letters/macros followed by a parenthesised
# argument list.  Value: Mathematica head (argument order unchanged).
PAREN_FUNCS = {
    'Gamma': 'Gamma', 'psi': 'PolyGamma0', 'zeta': 'ZetaFn', 'Beta': 'Beta',
    'B': 'Beta', 'Li': 'PolyLogFn', 'Cl': 'ClausenFn', 'Ls': 'LsFn',
    'W': 'ProductLog', 'K': 'EllipticKFn', 'E': 'EllipticEFn', 'Ai': 'AiryAi',
    'Bi': 'AiryBi', 'erf': 'Erf', 'Phi': 'LerchPhi', 'G': 'MeijerGFn',
    'beta': 'DirichletBetaFn', 'eta': 'DirichletEtaFn', 'operatorname': None,
    'F': 'MatharF2', 'Pi': 'MatharPi',
}

# Subscripted special functions  X_{nu}(z).
SUB_FUNCS = {
    'J': 'BesselJ', 'Y': 'BesselY', 'N': 'BesselY', 'I': 'BesselI',
    'K': 'BesselK', 'H': 'HermiteH', 'P': 'LegendreP', 'Q': 'LegendreQ',
    'L': 'LaguerreL', 'T': 'ChebyshevT', 'U': 'ChebyshevU', 'C': 'GegenbauerC',
    'Li': 'PolyLog', 'j': 'SphericalBesselJ', 'y': 'SphericalBesselY',
    'B': 'BernoulliB', 'E': 'EulerE', 'Cl': 'ClausenSub', 'psi': 'PolyGammaSub',
    'zeta': 'ZetaSub', 'D': 'ParabolicCylinderD', 'F': 'HypergeometricSub',
    'Ei': 'ExpIntegralSub', 'Ls': 'LsSub', 'Ai': 'AiryAi',
}

INVERSE = {'Sin': 'ArcSin', 'Cos': 'ArcCos', 'Tan': 'ArcTan', 'Cot': 'ArcCot',
           'Sec': 'ArcSec', 'Csc': 'ArcCsc', 'Sinh': 'ArcSinh', 'Cosh': 'ArcCosh',
           'Tanh': 'ArcTanh', 'Coth': 'ArcCoth', 'Sech': 'ArcSech', 'Csch': 'ArcCsch'}

IGNORE_CMDS = {
    ',', ';', '!', ':', ' ', 'quad', 'qquad', 'displaystyle', 'textstyle',
    'nonumber', 'notag', 'limits', 'nolimits', 'big', 'Big', 'bigg', 'Bigg',
    'bigl', 'bigr', 'Bigl', 'Bigr', 'left', 'right', 'biggl', 'biggr', 'Biggl', 'Biggr', 'middle',
    'mathop', 'label', 'hfill', 'allowbreak', 'phantom', 'hspace', 'mbox', 'enspace',
    'thinspace', 'negthinspace', '/', 'relax', 'strut', 'cdots@', 'medskip',
}

RELATIONS = {'=', '<', '>', 'le', 'leq', 'ge', 'geq', 'neq', 'ne', 'approx', 'sim',
             'equiv', 'simeq', 'cong', 'to', 'rightarrow', 'Rightarrow', 'in',
             'notin', 'subset', 'll', 'gg', 'propto'}

# --------------------------------------------------------------------------
# Tokenizer
# --------------------------------------------------------------------------
TOK_RE = re.compile(r"""
    (?P<cmd>\\(?:[A-Za-z]+|.))      |
    (?P<num>\d+(?:\.\d+)?|\.\d+)    |
    (?P<ws>\s+)                     |
    (?P<ch>.)
""", re.X | re.S)


def tokenize(tex):
    toks = []
    pending = None
    for m in TOK_RE.finditer(tex):
        k = m.lastgroup
        v = m.group(k)
        if k == 'ws':
            continue
        if k == 'cmd':
            name = v[1:]
            if name in IGNORE_CMDS:
                continue
            if name in ('LEFT', 'RIGHT'):
                pending = name
                continue
            toks.append(('cmd', name))
        elif k == 'num':
            toks.append(('num', v))
        else:
            if v == '%':
                continue
            if pending:
                toks.append(('cmd', pending + v))
                pending = None
                continue
            toks.append(('ch', v))
    return toks


def strip_comments(tex):
    return re.sub(r'(?<!\\)%[^\n]*', '', tex)


def preprocess(tex):
    tex = strip_comments(tex)
    tex = re.sub(r'\\label\{[^}]*\}', '', tex)
    tex = re.sub(r'\\(?:text|mbox|mathrm|textrm|rm|operatorname|mathit)\s*\{\s*([A-Za-z]+)\s*\}',
                 lambda m: '\\' + m.group(1) if m.group(1) in ELEM or m.group(1) in PAREN_FUNCS
                 or m.group(1) in SUB_FUNCS else '\\TXT{' + m.group(1) + '}', tex)
    tex = re.sub(r'\\(?:hspace|vspace|phantom)\*?\{[^}]*\}', '', tex)
    tex = tex.replace('\\\\', ' ').replace('&', ' ')
    tex = re.sub(r'\\(?:begin|end)\{(?:array|aligned|split)\}(\{[^}]*\})?', ' ', tex)
    tex = tex.replace('\\dfrac', '\\frac').replace('\\tfrac', '\\frac')
    tex = tex.replace('\\left.', '').replace('\\right.', '')
    tex = re.sub(r'\\left\\?([(\[|{])', lambda m: '\\LEFT' + m.group(1), tex)
    tex = re.sub(r'\\right\\?([)\]|}])', lambda m: '\\RIGHT' + m.group(1), tex)
    tex = tex.replace('\\lbrace', '{').replace('\\rbrace', '}')
    tex = tex.replace('\\lvert', '|').replace('\\rvert', '|')
    tex = tex.replace('\\vert', '|')
    tex = tex.replace("\\prime", "'")
    tex = re.sub(r'\(\s*-\s*\)', '(-1)', tex)
    tex = re.sub(r'\\mathbf\s*\{\s*([KEF])\s*\}', r'\1', tex)
    tex = tex.replace('\\{', '(').replace('\\}', ')')
    tex = tex.replace('\\mid', '|')
    tex = re.sub(r'(^|[^}A-Za-z0-9])\s*_\s*(\d)\s*F\s*_', r'\1{}_\2F_', tex)
    return tex


# --------------------------------------------------------------------------
# AST
# --------------------------------------------------------------------------
# Nodes are tuples:
#   ('num', '3') ('sym', 'x') ('const', 'Pi') ('call', head, [args])
#   ('add', [terms]) ('mul', [factors]) ('neg', e) ('pow', b, e) ('div', n, d)
#   ('int', body, var, lo, hi) ('sum', body, var, lo, hi, kind)
#   ('fact', e) ('abs', e) ('deriv', head, order, sub, args)  ('d', var)
#   ('rel', op, a, b)

class Parser:
    def __init__(self, toks):
        self.t = toks
        self.i = 0
        self.bound = set()  # bound summation / integration variables

    # --- token helpers -------------------------------------------------
    def peek(self, k=0):
        j = self.i + k
        return self.t[j] if j < len(self.t) else ('eof', None)

    def next(self):
        tok = self.peek()
        self.i += 1
        return tok

    def at(self, kind, val=None):
        tk = self.peek()
        return tk[0] == kind and (val is None or tk[1] == val)

    def expect(self, kind, val):
        tk = self.next()
        if tk[0] != kind or tk[1] != val:
            raise TexError('expected %s got %r near %s' % (val, tk, self.ctx()))
        return tk

    def ctx(self):
        return ' '.join(str(v) for _, v in self.t[max(0, self.i - 4): self.i + 4])

    # --- grammar -------------------------------------------------------
    def parse_all(self):
        e = self.expr()
        if self.peek()[0] != 'eof':
            raise TexError('trailing tokens: ' + self.ctx())
        return e

    def group(self):
        """A TeX argument: {...} or a single token."""
        if self.at('ch', '{'):
            self.next()
            e = self.expr(stop={'}'})
            self.expect('ch', '}')
            return e
        tk = self.peek()
        if tk[0] == 'num' and len(tk[1]) > 1 and '.' not in tk[1]:
            self.t[self.i] = ('num', tk[1][1:])
            return ('num', tk[1][0])
        return self.atom(single=True)

    def raw_group_tokens(self):
        if not self.at('ch', '{'):
            return [self.next()]
        depth = 0
        out = []
        while True:
            tk = self.next()
            if tk[0] == 'eof':
                raise TexError('unbalanced {')
            if tk == ('ch', '{'):
                depth += 1
                if depth == 1:
                    continue
            elif tk == ('ch', '}'):
                depth -= 1
                if depth == 0:
                    return out
            out.append(tk)

    def is_stop(self, stop):
        tk = self.peek()
        if tk[0] == 'eof':
            return True
        if tk[0] == 'ch' and tk[1] in stop:
            return True
        if tk[0] == 'cmd' and tk[1] in stop:
            return True
        if tk[0] == 'cmd' and tk[1] in ('RIGHT)', 'RIGHT]', 'RIGHT|', 'RIGHT}'):
            return True
        if tk[0] == 'ch' and tk[1] in ('=', '<', '>') and '=' in stop:
            return True
        if tk[0] == 'cmd' and tk[1] in RELATIONS and '=' in stop:
            return True
        return False

    def expr(self, stop=frozenset()):
        stop = set(stop) | {'}'}
        terms = []
        sign = 1
        if self.at('ch', '+'):
            self.next()
        elif self.at('ch', '-'):
            self.next()
            sign = -1
        elif self.at('cmd', 'pm') or self.at('cmd', 'mp'):
            raise TexError('\\pm')
        t = self.term(stop)
        terms.append(t if sign == 1 else ('neg', t))
        while not self.is_stop(stop):
            if self.at('ch', '+'):
                self.next()
                terms.append(self.term(stop))
            elif self.at('ch', '-'):
                self.next()
                terms.append(('neg', self.term(stop)))
            elif self.at('cmd', 'pm') or self.at('cmd', 'mp'):
                raise TexError('\\pm in expression')
            elif self.at('ch', ',') or self.at('ch', ';') or self.at('ch', '.'):
                break
            else:
                raise TexError('unexpected %r near %s' % (self.peek(), self.ctx()))
        return terms[0] if len(terms) == 1 else ('add', terms)

    TERM_END_CH = {'+', '-', '=', '<', '>', ',', ';', ')', ']', '}', '|', '.'}

    def term(self, stop):
        facs = []
        while True:
            tk = self.peek()
            if self.is_stop(stop):
                break
            if tk[0] == 'ch' and tk[1] in self.TERM_END_CH:
                if tk[1] == '|' and self.abs_can_open(facs):
                    pass
                elif tk[1] == '.' and not facs:
                    raise TexError('stray .')
                else:
                    break
            if tk[0] == 'cmd' and tk[1] in ('pm', 'mp'):
                break
            if tk[0] == 'cmd' and tk[1] in ('cdot', 'times'):
                self.next()
                continue
            if tk == ('ch', '/'):
                self.next()
                if not facs:
                    raise TexError('leading /')
                d = self.factor(stop)
                num = facs[0] if len(facs) == 1 else ('mul', facs)
                facs = [('div', num, d)]
                continue
            if tk == ('ch', '*'):
                self.next()
                continue
            facs.append(self.factor(stop))
        if not facs:
            raise TexError('empty term near ' + self.ctx())
        return facs[0] if len(facs) == 1 else ('mul', facs)

    def abs_can_open(self, facs):
        return getattr(self, '_absdepth', 0) == 0 or not facs

    def factor(self, stop):
        base = self.atom(stop=stop)
        return self.postfix(base)

    def postfix(self, base):
        while True:
            if self.at('ch', '^'):
                self.next()
                ex = self.group()
                base = ('pow', base, ex)
            elif self.at('ch', '!'):
                self.next()
                if self.at('ch', '!'):
                    self.next()
                    base = ('call', 'Factorial2', [base])
                else:
                    base = ('fact', base)
            elif self.at('ch', "'") and base[0] == 'call':
                raise TexError("prime on call")
            else:
                return base

    # --- atoms ---------------------------------------------------------
    def atom(self, stop=frozenset(), single=False):
        tk = self.peek()
        kind, v = tk
        if kind == 'eof':
            raise TexError('unexpected end')
        if kind == 'num':
            self.next()
            return ('num', v)
        if kind == 'ch':
            if v == '(':
                self.next()
                e = self.expr(stop={')'})
                self.expect('ch', ')')
                return e
            if v == '[':
                self.next()
                e = self.expr(stop={']'})
                self.expect('ch', ']')
                return e
            if v == '{':
                self.next()
                if self.at('ch', '}'):   # {}_2F_1
                    self.next()
                    return self.prescript_func()
                e = self.expr(stop={'}'})
                self.expect('ch', '}')
                return e
            if v == '|':
                self.next()
                self._absdepth = getattr(self, '_absdepth', 0) + 1
                e = self.expr(stop={'|'})
                self._absdepth -= 1
                self.expect('ch', '|')
                return ('abs', e)
            if v.isalpha():
                return self.letter()
            if v == '-' and single:
                self.next()
                return ('neg', self.atom(single=True))
            raise TexError('unexpected char %r near %s' % (v, self.ctx()))
        # commands
        if v in ('LEFT(', 'LEFT[', 'LEFT{', 'LEFT|'):
            self.next()
            close = 'RIGHT' + {'(': ')', '[': ']', '{': '}', '|': '|'}[v[-1]]
            e = self.expr(stop={close})
            if not self.at('cmd', close):
                raise TexError('unmatched ' + v)
            self.next()
            return ('abs', e) if v == 'LEFT|' else e
        if v == 'frac':
            self.next()
            n = self.group()
            d = self.group()
            return ('div', n, d)
        if v == 'binom' or v == 'choose':
            self.next()
            a = self.group()
            b = self.group()
            return ('call', 'Binomial', [a, b])
        if v == 'sqrt':
            self.next()
            if self.at('ch', '['):
                self.next()
                idx = self.expr(stop={']'})
                self.expect('ch', ']')
                a = self.group()
                return ('pow', a, ('div', ('num', '1'), idx))
            return ('call', 'Sqrt', [self.group()])
        if v == 'surd':
            self.next()
            a = self.postfix(self.atom(single=True))
            return ('call', 'Sqrt', [a])
        if v == 'infty':
            self.next()
            return ('const', 'Infinity')
        if v == 'pi':
            self.next()
            return ('const', 'Pi')
        if v in ('int', 'intop'):
            return self.integral(stop)
        if v in ('sum', 'prod'):
            return self.summation(v, stop)
        if v in ('lfloor', 'lceil'):
            self.next()
            close = 'rfloor' if v == 'lfloor' else 'rceil'
            e = self.expr(stop={close})
            self.expect('cmd', close)
            return ('call', 'Floor' if v == 'lfloor' else 'Ceiling', [e])
        if v == 'TXT':
            raise TexError('text in formula')
        if v in ELEM and not (v in PAREN_FUNCS and self.peek(1) == ('ch', '(')
                              and v in ('Ei', 'Si', 'Ci', 'Shi', 'Chi', 'li', 'erf')):
            return self.elem_func(ELEM[v])
        if v in GREEK or v in PAREN_FUNCS or v in SUB_FUNCS:
            return self.named(v, is_cmd=True)
        if v == 'mathcal' or v == 'cal':
            raise TexError('mathcal')
        if v == 'ldots' or v == 'cdots' or v == 'dots':
            raise TexError('ellipsis')
        raise TexError('unknown command \\' + v)

    def letter(self):
        _, v = self.next()
        # differential: 'd' followed by a variable letter / greek, in integrand
        if v == 'd' and self.in_integrand and self.diff_follows():
            var = self.read_var()
            return ('d', var)
        if v == 'e' and getattr(self, 'e_param', False):
            return ('sym', 'ee')
        if v == 'e' and self.at('ch', '^'):
            self.next()
            ex = self.group()
            return ('call', 'Exp', [ex])
        if v == 'e':
            return ('const', 'E')
        if v == 'i' and 'i' not in self.bound:
            return ('const', 'I')
        self.i -= 1
        return self.named(v, is_cmd=False)

    in_integrand = False

    def diff_follows(self):
        tk = self.peek()
        if tk[0] == 'ch' and tk[1].isalpha():
            nxt = self.peek(1)
            return True if self.int_var_hint is None else tk[1] == self.int_var_hint
        if tk[0] == 'cmd' and tk[1] in GREEK:
            return self.int_var_hint is None or GREEK[tk[1]] == self.int_var_hint
        return False

    int_var_hint = None

    def read_var(self):
        tk = self.next()
        name = tk[1] if tk[0] == 'ch' else GREEK[tk[1]]
        # x_1 style
        if self.at('ch', '_'):
            self.next()
            sub = self.raw_group_tokens()
            name += ''.join(str(t[1]) for t in sub)
        return name

    def named(self, v, is_cmd):
        """A symbol or a named function (letter or macro)."""
        self.next()
        name = GREEK[v] if (is_cmd and v in GREEK) else v
        # Subscript: special function or indexed symbol
        if self.at('ch', '_'):
            self.next()
            sub = self.group()
            sup = None
            if self.at('ch', '^') and v in ('P', 'Q'):   # P_\nu^\mu
                self.next()
                sup = self.group()
            primes = self.primes()
            if (v in SUB_FUNCS) and (self.at('ch', '(') or self.at('cmd', 'LEFT(')):
                args = self.paren_args()
                return self.mk_sub_func(v, sub, sup, args, primes)
            if v in SUB_FUNCS and v in ('B', 'E') and not primes:
                return ('call', SUB_FUNCS[v], [sub])     # Bernoulli / Euler numbers
            if primes:
                raise TexError('primes on indexed symbol')
            return ('sym', name + '_' + flat_name(sub))
        if v == 'psi' and self.at('ch', '^') and self.peek(1) == ('ch', '{') and self.peek(2) == ('ch', '('):
            self.next(); self.next(); self.next()
            order = self.expr(stop={')'})
            self.expect('ch', ')')
            self.expect('ch', '}')
            args = self.paren_args()
            return ('call', 'PolyGamma', [order] + args)
        powfn = None
        if self.at('ch', '^') and v in PAREN_FUNCS and v not in ('E', 'B', 'K', 'G', 'beta', 'eta', 'W') and self.paren_after_pow():
            # \Gamma^3(1+j)
            self.next()
            powfn = self.group()
        primes = self.primes()
        if v in PAREN_FUNCS and (self.at('ch', '(') or self.at('cmd', 'LEFT(')):
            if is_cmd or v in ('B', 'W', 'K', 'E', 'G', 'Ai', 'Bi', 'F'):
                if v in ('K', 'E', 'G', 'W', 'B') and self.in_integrand_symbol(v):
                    pass
                else:
                    args = self.paren_args()
                    r = self.mk_paren_func(v, args, primes)
                    return ('pow', r, powfn) if powfn is not None else r
        if primes:
            if v == 'psi' or name == 'psi':
                args = self.paren_args()
                return ('call', 'PolyGamma', [('num', str(primes))] + args)
            raise TexError('primes on symbol ' + v)
        if powfn is not None:
            raise TexError('power on non-function')
        if is_cmd and v in PAREN_FUNCS and v not in GREEK:
            raise TexError('function without args: ' + v)
        if name == 'gamma':
            return ('sym', 'gamma')   # resolved later (EulerGamma vs parameter)
        return ('sym', name)

    def in_integrand_symbol(self, v):
        return False

    def paren_after_pow(self):
        # lookahead: ^ group then '('
        j = self.i + 1
        if j < len(self.t) and self.t[j] == ('ch', '{'):
            depth = 0
            while j < len(self.t):
                if self.t[j] == ('ch', '{'):
                    depth += 1
                elif self.t[j] == ('ch', '}'):
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            j += 1
        else:
            j += 1
        return j < len(self.t) and self.t[j] in (('ch', '('), ('cmd', 'LEFT('))

    def primes(self):
        n = 0
        while self.at('ch', "'"):
            self.next()
            n += 1
        return n

    def paren_args(self):
        if self.at('cmd', 'LEFT('):
            self.next()
            close = ('cmd', 'RIGHT)')
        else:
            self.expect('ch', '(')
            close = ('ch', ')')
        args = []
        groups = [[]]
        while True:
            stopset = {')', ',', ';', '|'} if close[0] == 'ch' else {'RIGHT)', ',', ';', '|'}
            e = self.expr(stop=stopset)
            groups[-1].append(e)
            tk = self.next()
            if tk == close:
                break
            if tk == ('ch', ','):
                continue
            if tk in (('ch', ';'), ('ch', '|')):
                groups.append([])
                continue
            raise TexError('bad arg list near ' + self.ctx())
        if len(groups) == 1:
            return groups[0]
        return [('list', g) for g in groups]

    def mk_paren_func(self, v, args, primes):
        flat = [a for a in args if a[0] != 'list']
        if any(a[0] == 'list' for a in args):
            if v == 'G':
                raise TexError('MeijerG args')
            flat = [x for a in args for x in (a[1] if a[0] == 'list' else [a])]
        head = PAREN_FUNCS[v]
        if v == 'psi':
            if primes:
                return ('call', 'PolyGamma', [('num', str(primes))] + flat)
            return ('call', 'PolyGamma', flat)
        if primes and v in ('K', 'E') and primes == 1 and len(flat) == 1:
            return ('call', 'MatharK1' if v == 'K' else 'MatharE1', flat)
        if primes and v in ('zeta', 'Gamma'):
            return ('call', 'DerivN', [('num', str(primes)), ('str', 'Zeta' if v == 'zeta' else 'Gamma')] + flat)
        if primes:
            raise TexError('primes on ' + v)
        if v == 'zeta':
            return ('call', 'Zeta' if len(flat) == 1 else 'HurwitzZeta', flat)
        if v == 'Li':
            raise TexError('Li without index')
        if v == 'K':
            if len(flat) != 1:
                raise TexError('K arity')
            return ('call', 'MatharK', flat)
        if v == 'E':
            if len(flat) == 2:
                return ('call', 'MatharE2', flat)
            if len(flat) != 1:
                raise TexError('E arity')
            return ('call', 'MatharE', flat)
        if v == 'beta':
            return ('call', 'DirichletBeta', flat)
        if v == 'eta':
            return ('call', 'DirichletEta', flat)
        if v == 'Cl':
            raise TexError('Cl without index')
        if v == 'G' and len(flat) == 1:
            return ('call', 'BarnesG', flat)
        if v == 'G':
            raise TexError('G arity')
        if v == 'Phi':
            return ('call', 'LerchPhi', flat)
        if v == 'F':
            if len(flat) != 2:
                raise TexError('F arity')
            return ('call', 'MatharF2', flat)
        if v == 'Pi':
            return ('call', 'MatharPi', flat)
        if v in ('B', 'Beta'):
            if len(flat) == 2:
                return ('call', 'Beta', flat)
            raise TexError('B arity')
        return ('call', head, flat)

    def mk_sub_func(self, v, sub, sup, args, primes):
        flat = [x for a in args for x in (a[1] if a[0] == 'list' else [a])]
        if primes:
            if v in ('J', 'Y', 'I', 'K', 'N') and primes == 1:
                head = SUB_FUNCS[v]
                return ('call', 'Derivative1', [('str', head), sub] + flat)
            raise TexError('primes on subscripted function')
        if v == 'F':
            # {}_pF_q handled by prescript_func; a bare F_q(...) is unknown
            raise TexError('F_q')
        if v == 'Li':
            return ('call', 'PolyLog', [sub] + flat)
        if v == 'Cl':
            # Cl_s(t) = Im Li_s(e^{it}) for even s, Re Li_s(e^{it}) for odd s.
            if sub[0] != 'num' or len(flat) != 1:
                raise TexError('Clausen order')
            part = 'Im' if int(sub[1]) % 2 == 0 else 'Re'
            return ('call', part, [('call', 'PolyLog', [sub, ('call', 'Exp', [('mul', [('const', 'I'), flat[0]])])])])
        if v == 'psi':
            raise TexError('psi_n')
        if v == 'P' and sup is not None:
            return ('call', 'LegendreP', [sub, sup] + flat)
        if v == 'L' and sup is None and len(flat) == 1:
            return ('call', 'LaguerreL', [sub] + flat)
        if v in ('B', 'E') and len(flat) == 1:
            return ('call', 'BernoulliB' if v == 'B' else 'EulerE', [sub] + flat)
        if v == 'C' and sup is None:
            raise TexError('Gegenbauer without lambda')
        if len(flat) != 1 and v not in ('L', 'C'):
            raise TexError('arity of ' + v)
        return ('call', SUB_FUNCS[v], [sub] + flat)

    def prescript_func(self):
        # {}_pF_q(a1..;b1..;z)
        self.expect('ch', '_')
        p = self.group()
        tk = self.next()
        if tk != ('ch', 'F'):
            raise TexError('prescript not F')
        self.expect('ch', '_')
        q = self.group()
        if self.at('ch', '['):
            raise TexError('pFq bracket form')
        if self.at('ch', '\\'):
            pass
        args = self.paren_args()
        lists = [a for a in args if a[0] == 'list']
        if len(args) != 3 or len(lists) != 3:
            # also accept (a,b;c;z) where groups parse into lists
            raise TexError('pFq arg shape')
        a = lists[0][1]
        b = lists[1][1]
        z = lists[2][1]
        if len(z) != 1:
            raise TexError('pFq z')
        return ('call', 'HypergeometricPFQ', [('list', a), ('list', b), z[0]])

    def elem_func(self, head):
        self.next()
        powe = None
        if self.at('ch', '^'):
            self.next()
            powe = self.group()
            if powe == ('neg', ('num', '1')):
                inv = INVERSE.get(head)
                if inv is None:
                    raise TexError('f^{-1} inverse notation')
                head, powe = inv, None
        if self.at('ch', '_') and head == 'Log':
            raise TexError('log base')
        arg = self.func_arg(head)
        r = ('call', head, [arg])
        if head == 'Exp':
            r = ('call', 'Exp', [arg])
        return ('pow', r, powe) if powe is not None else r

    def func_arg(self, head):
        """Argument of an operator function: parenthesised, or the implicit
        product of simple factors that follows (`\\sin 2x`, `\\log\\frac{a}{b}`)."""
        if self.at('ch', '(') or self.at('cmd', 'LEFT(') or self.at('cmd', 'LEFT['):
            # (..) then possibly ^ applies to the function, not argument:
            # \sin(x)^2  -> Sin[x]^2 handled by postfix of caller
            save = self.i
            a = self.atom()
            # \log(1+x) y : stop here
            return a
        facs = []
        if self.at('ch', '|') or self.at('cmd', 'LEFT|'):
            return self.atom()
        tk = self.peek()
        if tk[0] == 'cmd' and tk[1] in ELEM:
            return self.factor(set())
        while True:
            tk = self.peek()
            if tk[0] == 'eof':
                break
            if tk[0] == 'ch' and (tk[1] in self.TERM_END_CH or tk[1] in '(['):
                if facs or tk[1] not in '([':
                    break
            if tk[0] == 'cmd':
                if tk[1] in ELEM or tk[1] in ('int', 'sum', 'prod', 'cdot', 'times', 'pm', 'mp'):
                    break
                if tk[1].startswith('RIGHT') or tk[1] in RELATIONS:
                    break
                if tk[1] == 'frac' and facs:
                    break
                if tk[1] in ('Gamma', 'psi', 'zeta') and facs:
                    break
            if tk[0] == 'ch' and tk[1] == 'd' and self.in_integrand:
                nxt = self.peek(1)
                if (nxt[0] == 'ch' and nxt[1].isalpha() and (self.int_var_hint is None or nxt[1] == self.int_var_hint)) or \
                   (nxt[0] == 'cmd' and nxt[1] in GREEK and (self.int_var_hint is None or GREEK[nxt[1]] == self.int_var_hint)):
                    break
            if tk[0] == 'ch' and tk[1] == '{' and False:
                break
            facs.append(self.factor(set()))
            # A function argument ends after a factor followed by another
            # letter that begins a new function-like object is ambiguous; we
            # keep collecting simple factors only.
        if not facs:
            raise TexError('missing function argument near ' + self.ctx())
        return facs[0] if len(facs) == 1 else ('mul', facs)

    def integral(self, stop):
        self.next()
        lo = hi = None
        for _ in range(2):
            if self.at('ch', '_'):
                self.next()
                lo = self.group()
            elif self.at('ch', '^'):
                self.next()
                hi = self.group()
        old = (self.in_integrand, self.int_var_hint)
        self.in_integrand = True
        self.int_var_hint = None
        body = self.expr(stop=set(stop) | {'='})
        self.in_integrand, self.int_var_hint = old
        body, var = extract_differential(body)
        if var is None:
            # Some entries omit the differential (`\\int_0^{\\pi/2}\\frac{x}{\\sin x}`).
            # Accept only when exactly one conventional integration variable
            # occurs; the numeric oracle confirms the reading downstream.
            cands = free_syms(body) & VAR_NAMES
            if len(cands) != 1:
                raise TexError('no differential found')
            var = cands.pop()
        return ('int', body, var, lo, hi)

    def summation(self, kind, stop):
        self.next()
        lo = hi = None
        var = None
        for _ in range(2):
            if self.at('ch', '_'):
                self.next()
                toks = self.raw_group_tokens()
                sub = Parser(toks)
                sub.bound = self.bound
                if ('ch', '=') in toks:
                    k = toks.index(('ch', '='))
                    vt = toks[:k]
                    if len(vt) != 1:
                        raise TexError('sum var')
                    var = vt[0][1] if vt[0][0] == 'ch' else GREEK.get(vt[0][1], vt[0][1])
                    sp = Parser(toks[k + 1:])
                    lo = sp.parse_all()
                elif len(toks) >= 3 and toks[1] == ('cmd', 'ge'):
                    var = toks[0][1]
                    lo = Parser(toks[2:]).parse_all()
                    hi = ('const', 'Infinity')
                else:
                    raise TexError('sum lower spec')
            elif self.at('ch', '^'):
                self.next()
                hi = self.group()
        if var is None or lo is None or hi is None:
            raise TexError('incomplete sum')
        self.bound = self.bound | {var}
        # summand: the rest of the term (product of factors)
        body = self.term(set(stop) | {'='})
        return ('sum', body, var, lo, hi, kind)


def flat_name(node):
    if node[0] == 'num':
        return node[1]
    if node[0] == 'sym':
        return node[1]
    if node[0] == 'mul':
        return ''.join(flat_name(f) for f in node[1])
    raise TexError('complex subscript')


def extract_differential(body):
    """Remove exactly one ('d', var) factor from the integrand AST."""
    found = []

    def walk(n, top):
        k = n[0]
        if k == 'd':
            found.append(n[1])
            return ('num', '1')
        if k == 'mul':
            return ('mul', [walk(f, top) for f in n[1]])
        if k == 'div':
            return ('div', walk(n[1], top), walk(n[2], False))
        if k == 'neg':
            return ('neg', walk(n[1], top))
        if k == 'add':
            return ('add', [walk(t, top) for t in n[1]])
        return n
    nb = walk(body, True)
    if len(found) != 1:
        return body, None
    return nb, found[0]


# --------------------------------------------------------------------------
# Printer
# --------------------------------------------------------------------------
PREC = {'add': 10, 'neg': 15, 'mul': 20, 'div': 20, 'pow': 30, 'atom': 100}


def prec(n):
    k = n[0]
    if k in ('add',):
        return 10
    if k == 'neg':
        return 15
    if k in ('mul', 'div'):
        return 20
    if k == 'pow':
        return 30
    if k == 'fact':
        return 40
    if k == 'num' and n[1].startswith('-'):
        return 15
    return 100


def par(s, inner, outer):
    return '(' + s + ')' if inner < outer else s


SYM_RENAME = {'gamma': 'gamma'}


def emit(n, ctx):
    k = n[0]
    if k == 'num':
        v = n[1]
        if '.' in v:
            # exact rational for decimals
            a, b = v.split('.')
            return '(%s/%s)' % (int(a + b) if a + b else 0, 10 ** len(b))
        return v
    if k == 'sym':
        name = n[1]
        if name == 'gamma' and ctx.get('gamma_const'):
            return 'EulerGamma'
        if name == 'G' and ctx.get('catalan'):
            return 'Catalan'
        return mm_symbol(name)
    if k == 'const':
        return n[1]
    if k == 'str':
        return n[1]
    if k == 'add':
        out = []
        for i, t in enumerate(n[1]):
            if t[0] == 'neg':
                s = emit(t[1], ctx)
                s = par(s, prec(t[1]), 20)
                out.append(('-' if i == 0 else ' - ') + s)
            else:
                s = emit(t, ctx)
                out.append(s if i == 0 else ' + ' + s)
        return ''.join(out)
    if k == 'neg':
        return '-' + par(emit(n[1], ctx), prec(n[1]), 20)
    if k == 'mul':
        fs = [f for f in n[1] if f != ('num', '1')] or [('num', '1')]
        return '*'.join(par(emit(f, ctx), prec(f), 21) for f in fs)
    if k == 'div':
        return par(emit(n[1], ctx), prec(n[1]), 20) + '/' + par(emit(n[2], ctx), prec(n[2]), 21)
    if k == 'pow':
        return par(emit(n[1], ctx), prec(n[1]), 31) + '^' + par(emit(n[2], ctx), prec(n[2]), 31)
    if k == 'fact':
        return par(emit(n[1], ctx), prec(n[1]), 41) + '!'
    if k == 'abs':
        return 'Abs[' + emit(n[1], ctx) + ']'
    if k == 'list':
        return '{' + ', '.join(emit(a, ctx) for a in n[1]) + '}'
    if k == 'call':
        head, args = n[1], n[2]
        if head == 'Derivative1':
            return 'Derivative[0, 1][%s][%s]' % (args[0][1], ', '.join(emit(a, ctx) for a in args[1:]))
        if head == 'EllipticKmod':
            return 'EllipticK[(%s)^2]' % emit(args[0], ctx)
        if head == 'EllipticEmod':
            return 'EllipticE[(%s)^2]' % emit(args[0], ctx)
        if head == 'DerivN':
            return 'Derivative[%s][%s][%s]' % (emit(args[0], ctx), args[1][1], ', '.join(emit(a, ctx) for a in args[2:]))
        if head == 'PolyGamma0':
            head = 'PolyGamma'
        return head + '[' + ', '.join(emit(a, ctx) for a in args) + ']'
    if k == 'int':
        _, body, var, lo, hi = n
        v = mm_symbol(var)
        b = emit(body, ctx)
        if lo is None and hi is None:
            return 'Integrate[%s, %s]' % (b, v)
        if lo is None or hi is None:
            raise TexError('one-sided limits')
        return 'Integrate[%s, {%s, %s, %s}]' % (b, v, emit(lo, ctx), emit(hi, ctx))
    if k == 'sum':
        _, body, var, lo, hi, kind = n
        return '%s[%s, {%s, %s, %s}]' % ('Sum' if kind == 'sum' else 'Product', emit(body, ctx),
                                         mm_symbol(var), emit(lo, ctx), emit(hi, ctx))
    if k == 'd':
        raise TexError('stray differential')
    raise TexError('cannot emit ' + k)


def mm_symbol(name):
    s = name.replace('_', '')
    if not re.match(r'^[A-Za-z][A-Za-z0-9]*$', s):
        raise TexError('bad symbol ' + name)
    # Avoid clashing with Mathematica builtins single-letter: C, D, E, I, K, N, O
    if s in ('C', 'D', 'E', 'I', 'K', 'N', 'O'):
        return 'c' + s
    return s


# --------------------------------------------------------------------------
# Relation splitting & public API
# --------------------------------------------------------------------------

def split_relation(tex):
    """Split at top-level '=' (brace depth 0, paren depth 0)."""
    out, cur, depth = [], [], 0
    i = 0
    while i < len(tex):
        c = tex[i]
        if c == '\\':
            m = re.match(r'\\([A-Za-z]+|.|$)', tex[i:], re.S)
            cur.append(m.group(0))
            i += len(m.group(0))
            continue
        if c in '{([':
            depth += 1
        elif c in '})]':
            depth -= 1
        if c == '=' and depth == 0:
            out.append(''.join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    out.append(''.join(cur))
    return out


def bare_e(toks):
    """True when the letter e occurs other than as e^{..}: then it is a
    parameter (eccentricity, modulus) throughout the entry, not Euler's e."""
    for j, tk in enumerate(toks):
        if tk == ('ch', 'e'):
            nxt = toks[j + 1] if j + 1 < len(toks) else ('eof', None)
            if nxt != ('ch', '^'):
                return True
            # e^2, e^{2}: a squared parameter, not exp(2)
            n2 = toks[j + 2] if j + 2 < len(toks) else ('eof', None)
            n3 = toks[j + 3] if j + 3 < len(toks) else ('eof', None)
            n4 = toks[j + 4] if j + 4 < len(toks) else ('eof', None)
            if n2[0] == 'num' or (n2 == ('ch', '{') and n3[0] == 'num' and n4 == ('ch', '}')):
                return True
    return False


E_PARAM = False


def parse_tex(tex, e_param=None):
    toks = tokenize(preprocess(tex))
    p = Parser(toks)
    p.e_param = bare_e(toks) if e_param is None else e_param
    return p.parse_all()


def free_syms(n, acc=None, bound=frozenset()):
    if acc is None:
        acc = set()
    k = n[0]
    if k == 'sym':
        if n[1] not in bound:
            acc.add(n[1])
    elif k == 'int':
        free_syms(n[1], acc, bound | {n[2]})
        for lim in (n[3], n[4]):
            if lim is not None:
                free_syms(lim, acc, bound)
    elif k == 'sum':
        free_syms(n[1], acc, bound | {n[2]})
        free_syms(n[3], acc, bound)
        free_syms(n[4], acc, bound)
    else:
        for c in n[1:]:
            if isinstance(c, tuple):
                free_syms(c, acc, bound)
            elif isinstance(c, list):
                for x in c:
                    if isinstance(x, tuple):
                        free_syms(x, acc, bound)
    return acc


def translate(tex, ctx=None, e_param=None):
    return emit(parse_tex(tex, e_param), ctx or {})


if __name__ == '__main__':
    import sys
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            print(translate(line))
        except TexError as e:
            print('ERR', e)
