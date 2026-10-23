/*
 * test_parallel.c -- the parallel builtins with the serial semantics of a
 * kernel that has no subkernels (src/parallel.c).
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>

static void test_parallel_map_table(void) {
    assert_eval_eq("ParallelMap[f, {1, 2}]", "{f[1], f[2]}", 0);
    assert_eval_eq("ParallelMap[f, {1, 2}, DistributedContexts -> None]", "{f[1], f[2]}", 0);
    assert_eval_eq("ParallelMap[f, {{1}, {2}}, {2}]", "{{f[1]}, {f[2]}}", 0);
    assert_eval_eq("ParallelMap[#^2 &, Range[4], Method -> \"ItemsPerEvaluation\" -> 2]", "{1, 4, 9, 16}", 0);
    assert_eval_eq("ParallelTable[i^2, {i, 3}]", "{1, 4, 9}", 0);
    assert_eval_eq("ParallelTable[i j, {i, 2}, {j, 2}, Method -> \"CoarsestGrained\"]", "{{1, 2}, {2, 4}}", 0);
    /* ParallelTable holds its arguments like Table. */
    assert_eval_eq("i = 7; ParallelTable[i, {i, 2}]", "{1, 2}", 0);
}

static void test_no_subkernels(void) {
    /* ParallelEvaluate runs on every subkernel -- none -- so expr is not evaluated. */
    assert_eval_eq("pe = 0; {ParallelEvaluate[pe = 1; pe, DistributedContexts -> None], pe}", "{{}, 0}", 0);
    assert_eval_eq("{Kernels[], DistributeDefinitions[f, g], $KernelCount, $KernelID}", "{{}, {}, 0, 0}", 0);
    assert_eval_eq("Quiet[LaunchKernels[4]]", "{}", 0);
    assert_eval_eq("Check[LaunchKernels[4], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("{Context[ParallelKernels], Context[DistributedContexts]}", "{\"System`\", \"System`\"}", 0);
    assert_eval_eq("{Attributes[ParallelMap], Attributes[ParallelTable], Attributes[ParallelEvaluate], "
                   "Attributes[DistributeDefinitions]}",
                   "{{Protected}, {HoldAll, Protected}, {HoldFirst, Protected}, {HoldAll, Protected}}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_parallel_map_table);
    TEST(test_no_subkernels);
    printf("All parallel tests passed.\n");
    return 0;
}
