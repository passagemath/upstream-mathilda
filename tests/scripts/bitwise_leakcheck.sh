#!/bin/sh
# bitwise_leakcheck.sh -- leak gate for the bit-level integer builtins.
#
# The C test binary cannot detect a leak on its own, so this is the real gate.
# Drives BitLength, BitXor, BitAnd, BitOr and BitNot (machine-int, bignum, negative,
# INT64_MIN, listable, packed, NDArray, n-ary buffer folds, symbolic
# cancellation / idempotence / 0- and -1-absorption, BitNot extraction, Compile, and
# error paths) under a leak checker and FAILS on any
# leaked byte.
#
#   macOS : leaks --atExit
#   Linux : valgrind --error-exitcode
#   neither present: SKIPS, loudly, with exit 0 -- an unavailable checker is not
#   a pass, so it says so rather than printing nothing (docs/adr/0004).
#
# Usage: bash tests/scripts/bitwise_leakcheck.sh [path-to-Mathilda]
set -e
BIN="${1:-./Mathilda}"
[ -x "$BIN" ] || { echo "bitwise_leakcheck: no Mathilda binary at $BIN — build first"; exit 2; }

WORK="$(mktemp -t bitleak.XXXXXX)"
cat > "$WORK" <<'MSCRIPT'
Do[BitLength[255], {200}];
Do[BitLength[-8], {200}];
Do[BitLength[0], {200}];
Do[BitLength[-9223372036854775808], {200}];
Do[BitLength[2^100], {200}];
Do[BitLength[2^100 - 1], {200}];
Do[BitLength[100!], {200}];
Do[BitLength[{0, 1, 2, 7, 8, 255, 256}], {200}];
Do[BitLength[{-1, -2, -8, -256}], {200}];
Do[BitLength[Range[64]], {200}];
Do[BitLength[NDArray[{-9223372036854775808, -1, 0, 1, 255}, DataType -> "int64"]], {200}];
Do[BitLength[3.5], {200}];
Do[BitLength[], {200}];
Do[BitLength[1, 2], {200}];
Do[BitLength[x], {200}];
Do[BitXor[61, 15], {200}];
Do[BitXor[3333, 5555, 7777, 9999], {200}];
Do[BitXor[-2, -3], {200}];
Do[BitXor[3^100, 5^100], {200}];
Do[BitXor[2^63, -1], {200}];
Do[BitXor[2^70, 2^70 + 5], {200}];
Do[BitXor[-3^50, 7, -5^30], {200}];
Do[BitXor[x, x, y], {200}];
Do[BitXor[x, y, x, y, x], {200}];
Do[BitXor[1, 2, x], {200}];
Do[BitXor[1, -1, x], {200}];
Do[BitXor[-2, 1, x], {200}];
Do[BitXor[-1, x, x], {200}];
Do[BitXor[-1, 2^70, x], {200}];
Do[BitXor[x, 2.5, 2.5], {200}];
Do[BitXor[], {200}];
Do[BitXor[{}, {}], {200}];
Do[BitXor[{1, 2}, 3, {4, 5}], {200}];
Do[BitXor[Range[64], 5], {200}];
Do[BitXor[Range[64], Reverse[Range[64]]], {200}];
Do[BitXor[Range[64], 1, 2, 3], {200}];
Do[BitXor[Range[64], Range[64], 9], {200}];
Do[BitXor[Range[8], 2^70], {200}];
Do[BitXor[Range[8], x, x], {200}];
Do[BitXor[NDArray[{-9223372036854775808, -1, 0, 1}, DataType -> "int64"], 7], {200}];
Do[BitXor[NDArray[{1, 2, 3}, DataType -> "int64"], NDArray[{3, 2, 1}, DataType -> "int64"], 1], {200}];
Do[Quiet[BitXor[{1, 2}, {1, 2, 3}]], {200}];
Do[Quiet[BitXor[2.5, 3]], {200}];
Do[Quiet[BitXor[1/2, x]], {200}];
Do[Quiet[BitXor[{1., 2.}, 3, 4]], {200}];
cfx = Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitXor[v, w, 6]];
Do[cfx[Range[32], Reverse[Range[32]]], {200}];
Do[Quiet[cfx[{1, 2}, {1, 2, 3}]], {200}];
Do[BitAnd[61, 15], {200}];
Do[BitAnd[3333, 5555, 7777, 9999], {200}];
Do[BitAnd[-2, -3], {200}];
Do[BitAnd[3^100, 5^100], {200}];
Do[BitAnd[3^100, 5^100, -7], {200}];
Do[BitAnd[2^63, -1], {200}];
Do[BitAnd[-2^63, -2^63 + 1], {200}];
Do[BitAnd[2^100 - 1, 2, 4, -8], {200}];
Do[BitAnd[-3^50, -5^40, x, x], {200}];
Do[BitAnd[x, y, y, x], {200}];
Do[BitAnd[3, 5, x], {200}];
Do[BitAnd[1, 2, x], {200}];
Do[BitAnd[0, x, 1/2], {200}];
Do[BitAnd[-1, x], {200}];
Do[BitAnd[-1, x, y], {200}];
Do[BitAnd[2^70, x, -1], {200}];
Do[BitAnd[x, 2^70 + 3, 2^70 + 1], {200}];
Do[BitAnd[x, 2.5, 2.5], {200}];
Do[BitAnd[2.5, 3], {200}];
Do[BitAnd[I, I, x], {200}];
Do[BitAnd[], {200}];
Do[BitAnd[{}, {}], {200}];
Do[BitAnd[{1, 2}, 3, {4, 5}], {200}];
Do[BitAnd[Range[64], 17], {200}];
Do[BitAnd[Range[64], Reverse[Range[64]]], {200}];
Do[BitAnd[Range[64], 7, 14], {200}];
Do[BitAnd[Range[64], Range[64], 9], {200}];
Do[BitAnd[Range[64], -1, -1], {200}];
Do[BitAnd[Range[8], 2^70 + 2], {200}];
Do[BitAnd[Range[8], x, x], {200}];
Do[BitAnd[NDArray[{-9223372036854775808, -1, 0, 1}, DataType -> "int64"], 7], {200}];
Do[BitAnd[NDArray[{1, 2, 3}, DataType -> "int64"], NDArray[{3, 2, 1}, DataType -> "int64"], 1], {200}];
Do[Quiet[BitAnd[{1, 2}, {1, 2, 3}]], {200}];
Do[Quiet[BitAnd[Range[2], Range[3], 1]], {200}];
cfa = Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitAnd[v, w, 6]];
Do[cfa[Range[32], Reverse[Range[32]]], {200}];
cfs = Compile[{{a, _Integer}, {b, _Integer}, {c, _Integer}}, BitAnd[a, b, c]];
Do[cfs[3333, 5555, -7], {200}];
Do[BitOr[61, 15], {200}];
Do[BitOr[3333, 5555, 7777, 9999], {200}];
Do[BitOr[-2, -3], {200}];
Do[BitOr[3^100, 5^100], {200}];
Do[BitOr[3^100, 5^100, 6], {200}];
Do[BitOr[2^63, -1], {200}];
Do[BitOr[-5, 2^70], {200}];
Do[BitOr[-2^70, 1], {200}];
Do[BitOr[-3^50, -5^40, x, x], {200}];
Do[BitOr[x, y, y, x], {200}];
Do[BitOr[1, 2, x], {200}];
Do[BitOr[-1, 2, x], {200}];
Do[BitOr[-1, x, 1/2], {200}];
Do[BitOr[0, x], {200}];
Do[BitOr[0, x, y], {200}];
Do[BitOr[2^70, x, 0], {200}];
Do[BitOr[x, 2^70 + 3, 2^70 + 1], {200}];
Do[BitOr[x, 2.5, 2.5], {200}];
Do[BitOr[2.5, 3], {200}];
Do[BitOr[I, I, x], {200}];
Do[BitOr[], {200}];
Do[BitOr[{}, {}], {200}];
Do[BitOr[{1, 2}, 8, {4, 5}], {200}];
Do[BitOr[Range[-40, 40], 17], {200}];
Do[BitOr[Range[64], Reverse[Range[64]]], {200}];
Do[BitOr[Range[64], 8, 16], {200}];
Do[BitOr[Range[64], Range[64], 9], {200}];
Do[BitOr[Range[64], 0, 0], {200}];
Do[BitOr[Range[8], 2^70], {200}];
Do[BitOr[Range[8], x, x], {200}];
Do[BitOr[NDArray[{-9223372036854775808, -1, 0, 1}, DataType -> "int64"], 7], {200}];
Do[BitOr[NDArray[{1, 2, 3}, DataType -> "int64"], NDArray[{3, 2, 1}, DataType -> "int64"], 16], {200}];
Do[Quiet[BitOr[{1, 2}, {1, 2, 3}]], {200}];
Do[Quiet[BitOr[Range[2], Range[3], 4]], {200}];
cfo = Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitOr[v, w, 6]];
Do[cfo[Range[32], Reverse[Range[32]]], {200}];
cfos = Compile[{{a, _Integer}, {b, _Integer}, {c, _Integer}}, BitOr[a, b, c]];
Do[cfos[3333, 5555, -7], {200}];
Do[BitNot[7], {200}];
Do[BitNot[-9223372036854775808], {200}];
Do[BitNot[2^63], {200}];
Do[BitNot[-2^63 - 1], {200}];
Do[BitNot[2^70], {200}];
Do[BitNot[100!], {200}];
Do[BitNot[BitNot[x]], {200}];
Do[BitNot[BitNot[BitNot[x]]], {200}];
bnz = BitNot[z]; Do[BitNot[bnz], {200}];
Do[BitXor[-1, BitNot[x]], {200}];
Do[BitNot[{1, {2, x}, 2.5, 2^70}], {200}];
Do[BitNot[Range[-10, 10]], {200}];
Do[BitNot[Partition[Range[12], 3]], {200}];
Do[BitNot[NDArray[{-9223372036854775808, -1, 0, 9223372036854775807}, DataType -> "int64"]], {200}];
Do[BitNot[NDArray[{1.5, 2.}]], {200}];
Do[BitNot[2.5], {200}];
Do[Quiet[BitNot[]], {200}];
Do[Quiet[BitNot[1, 2]], {200}];
cfn = Compile[{{v, _Integer, 1}}, BitNot[v]];
Do[cfn[Range[-50, 50]], {200}];
cfns = Compile[{{n, _Integer}}, BitNot[n]];
Do[cfns[12345], {200}];
MSCRIPT

if command -v leaks >/dev/null 2>&1; then
    OUT="$(leaks --atExit -- "$BIN" -file "$WORK" 2>/dev/null | grep 'total leaked bytes' | tail -1)"
    rm -f "$WORK"
    echo "bitwise_leakcheck: $OUT"
    case "$OUT" in
        *"0 leaks for 0 total leaked bytes"*) echo "bitwise_leakcheck: PASS"; exit 0 ;;
        "") echo "bitwise_leakcheck: FAIL — no verdict line from leaks"; exit 1 ;;
        *) echo "bitwise_leakcheck: FAIL — leaked memory above"; exit 1 ;;
    esac
elif command -v valgrind >/dev/null 2>&1; then
    valgrind --leak-check=full --errors-for-leak-kinds=definite \
             --error-exitcode=1 "$BIN" -file "$WORK" >/dev/null 2>/tmp/bitleak.valgrind
    rc=$?
    rm -f "$WORK"
    [ $rc -eq 0 ] && { echo "bitwise_leakcheck: PASS (valgrind)"; exit 0; }
    echo "bitwise_leakcheck: FAIL (valgrind) — see /tmp/bitleak.valgrind"; exit 1
else
    rm -f "$WORK"
    echo "bitwise_leakcheck: SKIPPED — neither 'leaks' nor 'valgrind' found."
    echo "bitwise_leakcheck: this is NOT a pass; the leak gate did not run."
    exit 0
fi
