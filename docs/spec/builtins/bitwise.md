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
