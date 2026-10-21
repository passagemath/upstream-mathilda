# Bitwise operations

Bit-level operations on integers. Implemented in `src/bitwise/` (one builtin per
file, registered by `bitwise_init()`).

## BitLength

- `BitLength[n]`: The number of binary bits needed to represent the integer `n`.

**Features**:
- `Protected`, `Listable`. Threads element-wise over a list of integers, e.g.
  `BitLength[{0, 1, 2, 7, 8, 255, 256}]`.
- For `n > 0`, `BitLength[n]` is an efficient exact version of
  `Floor[Log[2, n]] + 1` -- it never converts through floating-point and is
  exact for arbitrarily large `n` (via GMP's `mpz_sizeinbase` in base 2).
- `BitLength[0]` is `0`.
- For `n < 0`, `BitLength[n]` is equivalent to `BitLength[BitNot[n]]`. In two's
  complement `BitNot[n] = -n - 1`, so `BitLength[-1]` is `0`, `BitLength[-2]` is
  `1`, and `BitLength[-2^k]` is `k`.
- Works for both machine integers and bignums. The packed / `NDArray` fast path
  (`src/ndinteger.c`) and the `Compile[]` lowering (`OP_BLEN_I`) both handle the
  full `int64` range, including `INT64_MIN` (`BitLength[-2^63]` is `63`), which
  the complement-based kernel computes without overflow.

```mathematica
In[1]:= BitLength[255]
Out[1]= 8

In[2]:= BitLength[256]
Out[2]= 9

In[3]:= BitLength[2^100]
Out[3]= 101

In[4]:= BitLength[100!]
Out[4]= 525

In[5]:= BitLength[{0, 1, 2, 7, 8, 255, 256}]
Out[5]= {0, 1, 2, 3, 4, 8, 9}

In[6]:= BitLength[{-1, -2, -8, -256}]
Out[6]= {0, 1, 3, 8}
```

**Diagnostics**. A non-integer argument emits `BitLength::int` and echoes the
call back unevaluated; a wrong-arity call emits `BitLength::argx`. A symbolic
argument is left unevaluated silently.

```mathematica
In[7]:= BitLength[3.5]
        BitLength::int: Integer expected at position 1 in BitLength[3.5].
Out[7]= BitLength[3.5]

In[8]:= BitLength[1, 2]
        BitLength::argx: BitLength called with 2 arguments; 1 argument is expected.
Out[8]= BitLength[1, 2]
```

## BitXor

- `BitXor[n1, n2, ...]`: The bitwise XOR of the integers `ni` -- the integer
  whose binary representation has ones at the positions where an odd number of
  the `ni` have ones.

**Features**:
- `Flat`, `Listable`, `OneIdentity`, `Orderless`, `Protected`. Takes any number
  of arguments; `BitXor[]` is `0` and `BitXor[n]` is `n`. Threads over lists
  (`BitXor[Range[120], 5]`, `BitXor[{1, 2}, 3, {4, 5}]`), including nested and
  empty ones (`BitXor[{}]` is `{}`).
- Negative integers are taken in two's complement (`BitXor[-2, -3]` is `3`).
  Integers of any size are exact: bignums go through GMP's `mpz_xor`, read in
  place and XORed straight into a presized result (no operand or result
  copies), and a result that fits in a machine integer is returned as one
  (`BitXor[2^70, 2^70 + 5]` is `5`).
- Automatic symbolic simplification, in this order:
  1. identical arguments cancel in pairs (`BitXor[x, x]` is `0`,
     `BitXor[x, y, x, y, x]` is `x`, `BitXor[x, 2.5, 2.5]` is `x`);
  2. an odd number of `-1` arguments is pulled out as `BitNot`
     (`BitXor[-1, x]` is `BitNot[x]`, `BitXor[1, -1, x]` is
     `BitNot[BitXor[1, x]]`);
  3. the remaining integers combine into one constant, dropped when it is `0`,
     pulled out as `BitNot` when it is `-1` (`BitXor[1, 2, x]` is
     `BitXor[3, x]`, `BitXor[-2, 1, x]` is `BitNot[x]`); any other negative
     constant stays inline (`BitXor[-5, x]`).

  `BitNot` appears only as an inert head.
- **Packed / NDArray**: an exact int64 kernel (`ndk_BitXor_ii`, plus
  vectorisable whole-buffer arms for array ⊕ scalar and array ⊕ array) in
  `src/ndinteger.c`, on the
  `AWARE` and `INT64_OK` lists. It handles array ⊕ scalar and array ⊕ array, and
  the n-ary form (`BitXor[v, w, 7]`, `BitXor[v, 1, 2]`) stays on the buffer, one
  pass per operand. A visible `NDArray` gives a visible `NDArray`. Operands that
  no int64 buffer can hold (a bignum, a symbol) take the exact List path.
- **Compile[]**: `OP_BXOR_I` (total, never bails), folded left for any arity.
  It lowers at scalar shape, array ⊕ scalar, array ⊕ array and the n-ary array
  form, and `CompileDiagnostics` reports `Compiled -> True` for each. Scalar,
  List, packed, `NDArray` and compiled results agree bit for bit, `INT64_MIN`
  included.

```
In[1]:= BitXor[61, 15]
Out[1]= 50

In[2]:= BitXor[3^100, 5^100]
Out[2]= 7888609052210118054116793203790437104101594041461795136781437785275168

In[3]:= BitXor[3333, 5555, 7777, 9999]
Out[3]= 8664

In[4]:= {BitXor[x, x, x], BitXor[1, 2, x], BitXor[-1, x], BitXor[y, x]}
Out[4]= {x, BitXor[3, x], BitNot[x], BitXor[x, y]}

In[5]:= Outer[BitXor, {1, 0}, {1, 0}]
Out[5]= {{0, 1}, {1, 0}}
```

**Diagnostics**. A concrete non-integer number that survives pair cancellation
(a Real, Rational or Complex) emits `BitXor::int`, naming its position, and the
call stays unevaluated. Symbolic arguments are silent. The message routes
through the message funnel, so `Quiet` and `Check` see it.

```
In[6]:= BitXor[2.5, 3]
        BitXor::int: Integer expected at position 1 in BitXor[2.5, 3].
Out[6]= BitXor[2.5, 3]
```

## BitAnd

- `BitAnd[n1, n2, ...]`: The bitwise AND of the integers `ni` -- the integer
  whose binary representation has ones at the positions where all of the `ni`
  have ones.

**Features**:
- `Flat`, `Listable`, `OneIdentity`, `Orderless`, `Protected`. Takes any number
  of arguments; `BitAnd[]` is `-1` (the identity) and `BitAnd[n]` is `n`.
  Threads over lists (`BitAnd[Range[70, 130], 17]`, `BitAnd[{1, 2}, 3, {4, 5}]`),
  including nested and empty ones (`BitAnd[{}]` and `BitAnd[{}, {}]` are `{}`).
- Negative integers are taken in two's complement (`BitAnd[-2, -3]` is `-4`).
  Integers of any size are exact: bignums go through GMP's `mpz_and`, read in
  place and ANDed straight into a presized result (no operand or result
  copies), and a result that fits in a machine integer is returned as one
  (`BitAnd[2^70 + 5, 7]` is `5`).
- Automatic symbolic simplification (AND is idempotent, `-1` is its identity
  and `0` absorbs):
  1. identical arguments collapse to one (`BitAnd[x, y, y, x]` is
     `BitAnd[x, y]`, `BitAnd[2.5, 2.5]` is `2.5`);
  2. the integers combine into one constant (`BitAnd[3, 5, x]` is
     `BitAnd[1, x]`). A constant of `0` absorbs everything, symbols and
     non-integers included (`BitAnd[1, 2, x]` and `BitAnd[0, 2.5]` are `0`); a
     constant of `-1` is dropped (`BitAnd[-1, x]` is `x`); any other constant
     stays inline (`BitAnd[-5, x]`).

  There are no complement or absorption rules: `BitAnd[x, BitNot[x]]` and
  `BitAnd[x, BitOr[x, y]]` stay as they are.
- **Packed / NDArray**: an exact int64 kernel (`ndk_BitAnd_ii`, plus
  vectorisable whole-buffer arms for array & scalar and array & array) in
  `src/ndinteger.c`, on the `AWARE` and `INT64_OK` lists. The n-ary form
  (`BitAnd[v, w, 7]`, `BitAnd[v, 7, 14]`) stays on the buffer, one pass per
  operand; a combined constant of `-1` costs no pass. A visible `NDArray` gives
  a visible `NDArray`. Operands that no int64 buffer can hold (a bignum, a
  symbol) take the exact List path. At 10^5 elements `BitAnd[v, 255]` and
  `BitAnd[v, w]` take about 0.04 and 0.05 ms, level with NumPy's `bitwise_and`.
- **Compile[]**: `OP_BAND_I` (total, never bails), folded left for any arity.
  It lowers at scalar shape, array & scalar, array & array and the n-ary array
  form, and `CompileDiagnostics` reports `Compiled -> True` for each. Scalar,
  List, packed, `NDArray` and compiled results agree bit for bit, `INT64_MIN`
  included.

```
In[1]:= BitAnd[61, 15]
Out[1]= 13

In[2]:= BitAnd[3^100, 5^100]
Out[2]= 503913279078601980753385379100266097026814443729

In[3]:= BitAnd[3333, 5555, 7777, 9999]
Out[3]= 1025

In[4]:= {BitAnd[x, y, y, x], BitAnd[3, 5, x], BitAnd[1, 2, x], BitAnd[-1, x]}
Out[4]= {BitAnd[x, y], BitAnd[1, x], 0, x}

In[5]:= Table[Sign[BitAnd[i, 2^2 + 2^3]], {i, 16}]
Out[5]= {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}
```

**Diagnostics**. None. A non-integer number (a Real, Rational or Complex) is
left in the call silently, exactly as a symbol is (`BitAnd[2.5, 3]` stays
unevaluated and `Check` takes the value branch).

## BitOr

- `BitOr[n1, n2, ...]`: The bitwise OR of the integers `ni` -- the integer
  whose binary representation has ones at the positions where any of the `ni`
  have ones.

**Features**:
- `Flat`, `Listable`, `OneIdentity`, `Orderless`, `Protected`. Takes any number
  of arguments; `BitOr[]` is `0` (the identity) and `BitOr[n]` is `n`. Threads
  over lists (`BitOr[Range[-40, 40], 17]`, `BitOr[{1, 2}, 8, {4, 5}]`),
  including nested and empty ones (`BitOr[{}]` and `BitOr[{}, {}]` are `{}`).
- Negative integers are taken in two's complement (`BitOr[-2, -3]` is `-1`).
  Integers of any size are exact: bignums go through GMP's `mpz_ior`, read in
  place and ORed straight into a presized result (no operand or result
  copies). A result that fits in a machine integer is returned as one, which
  happens whenever a negative operand sign-extends over a larger positive one
  (`BitOr[-5, 2^70]` is `-5`).
- Automatic symbolic simplification, the dual of BitAnd (OR is idempotent, `0`
  is its identity and `-1`, all bits on, absorbs):
  1. identical arguments collapse to one (`BitOr[x, y, y, x]` is
     `BitOr[x, y]`, `BitOr[2.5, 2.5]` is `2.5`);
  2. the integers combine into one constant (`BitOr[1, 2, x]` is
     `BitOr[3, x]`). A constant of `-1` absorbs everything, symbols and
     non-integers included (`BitOr[-1, 2, x]` and `BitOr[-1, 2.5]` are `-1`); a
     constant of `0` is dropped (`BitOr[0, x]` is `x`); any other constant
     stays inline (`BitOr[-5, x]`).

  There are no complement or absorption rules: `BitOr[x, BitNot[x]]` and
  `BitOr[x, BitAnd[x, y]]` stay as they are.
- **Packed / NDArray**: an exact int64 kernel (`ndk_BitOr_ii`, plus
  vectorisable whole-buffer arms for array | scalar and array | array) in
  `src/ndinteger.c`, on the `AWARE` and `INT64_OK` lists. The n-ary form
  (`BitOr[v, w, 7]`, `BitOr[v, 8, 16]`) stays on the buffer, one pass per
  operand; a combined constant of `0` costs no pass. A visible `NDArray` gives
  a visible `NDArray`. Operands that no int64 buffer can hold (a bignum, a
  symbol) take the exact List path. At 10^5 elements `BitOr[v, 255]`,
  `BitOr[v, w]` and `BitOr[u, v, w]` take about 0.035, 0.044 and 0.089 ms,
  level with NumPy's `bitwise_or` (0.034, 0.043, 0.085 ms).
- **Compile[]**: `OP_BOR_I` (total, never bails), folded left for any arity.
  It lowers at scalar shape, array | scalar, array | array and the n-ary array
  form, and `CompileDiagnostics` reports `Compiled -> True` for each, so
  auto-compilation picks it up. Scalar, List, packed, `NDArray` and compiled
  results agree bit for bit, `INT64_MIN` included.

```
In[1]:= BitOr[61, 15]
Out[1]= 63

In[2]:= BitOr[3^100, 5^100]
Out[2]= 7888609052210118054117297117069515706082347426840895402878464599718897

In[3]:= BitOr[3333, 5555, 7777, 9999]
Out[3]= 16383

In[4]:= {BitOr[x, y, y, x], BitOr[1, 2, x], BitOr[-1, 2, x], BitOr[0, x]}
Out[4]= {BitOr[x, y], BitOr[3, x], -1, x}

In[5]:= Outer[BitOr, {1, 0}, {1, 0}]
Out[5]= {{1, 1}, {1, 0}}
```

**Diagnostics**. None. A non-integer number (a Real, Rational or Complex) is
left in the call silently, exactly as a symbol is (`BitOr[2.5, 3]` stays
unevaluated and `Check` takes the value branch).

## BitNot

- `BitNot[n]`: The bitwise NOT of the integer `n` -- ones become zeros and
  vice versa.

**Features**:
- `Listable`, `Protected`. Integers are taken in two's complement with an
  unlimited number of digits, so `BitNot[n]` is exactly `-1 - n`: `BitNot[0]`
  is `-1` (all bits on) and `BitNot[7]` is `-8`. Threads over lists, including
  nested and empty ones (`BitNot[{}]` is `{}`).
- Exact at any size. Machine integers complement in int64, which is total
  (`BitNot[-2^63]` is `2^63 - 1`). Bignums go through GMP's `mpz_com` into a
  fresh node, never touching the operand; a result that fits in a machine
  integer is returned as one (`BitNot[2^63]` is `-2^63 - 1`, `BitNot[-2^63 - 1]`
  is `2^63`).
- Automatic symbolic simplification: BitNot is an involution, so
  `BitNot[BitNot[x]]` is `x`. This also completes BitXor's `-1` extraction:
  `BitXor[-1, BitNot[x]]` is `x`. There are no other rules: `BitNot[-1 - x]`
  and `BitAnd[x, BitNot[x]]` stay as they are.
- **Packed / NDArray**: an exact int64 unary kernel (`ndk_BitNot_ii`, with a
  vectorisable whole-buffer arm `ndk_BitNot_whole`) in `src/ndinteger.c`, on
  the `AWARE` and `INT64_OK` lists. A packed array stays packed; a visible
  `NDArray` gives a visible `NDArray`. A real buffer has no integer answer, so
  its elements stay unevaluated, as the scalar does. At 10^5 elements
  `BitNot[v]` takes about 0.038 ms, level with NumPy's `invert` (0.034 ms).
  At 10^7 elements it takes 22 ms against 18.5 ms; both are bound by memory
  bandwidth.
- **Compile[]**: `OP_BNOT_I` (total, never bails). It lowers at scalar and
  array shape, and `CompileDiagnostics` reports `Compiled -> True` for both.
  Auto-compilation picks it up: `Table[x + BitNot[Floor[x]], {x, 0., 99999.}]`
  takes 1.9 ms compiled and 84 ms interpreted, with the same answer. Scalar,
  List, packed, `NDArray` and compiled results agree bit for bit, `INT64_MIN`
  included.

```
In[1]:= BitNot[7]
Out[1]= -8

In[2]:= BitNot[Range[-10, 10]]
Out[2]= {9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11}

In[3]:= BitNot[2^70]
Out[3]= -1180591620717411303425

In[4]:= {BitNot[BitNot[x]], BitXor[-1, BitNot[x]]}
Out[4]= {x, x}

In[5]:= BitNot[]
BitNot::argx: BitNot called with 0 arguments; 1 argument is expected.
Out[5]= BitNot[]
```

**Diagnostics**. `BitNot::argx` for a call that does not have exactly one
argument (`BitNot[]`, `BitNot[1, 2]`). It goes through the message funnel, so
`Quiet` silences it and `Check` sees it. A non-integer number (a Real,
Rational or Complex) is left in the call silently, exactly as a symbol is
(`BitNot[2.5]` stays unevaluated and `Check` takes the value branch).
