# Ch.7 -- timing and pausing execution: the two clocks, Pause, $TimeUnit
$TimeUnit
AbsoluteTiming[Pause[0.3]]
Timing[Pause[0.3]]
AbsoluteTiming[Pause[1/4]]
SessionTime[]
TimeUsed[]
# TimeConstrained bounds wall-clock time (v0.463): a pause is charged against the budget
AbsoluteTiming[TimeConstrained[Pause[5], 0.5]]
Timing[TimeConstrained[Pause[5], 0.5]]
TimeConstrained[Pause[5], 0.5, "too slow"]
