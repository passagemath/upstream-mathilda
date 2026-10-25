/*
 * tc_unwind_safe.h -- may TimeConstrained's SIGPROF handler unwind the code
 * it interrupted?
 *
 * A siglongjmp out of the C library's allocator (or out of a stdio call)
 * abandons the lock it holds, and the next allocation traps or deadlocks.
 * These two functions let the handler ask, from the interrupted signal
 * context, whether any frame of the interrupted stack lies inside the
 * allocator's or the C library's code, so it can defer the jump instead.
 */
#ifndef MATHILDA_TC_UNWIND_SAFE_H
#define MATHILDA_TC_UNWIND_SAFE_H

/* Record the code ranges of the allocator and C library images.  Call once,
 * outside any signal handler, before the first SIGPROF can arrive.  Cheap to
 * call again (later calls return at once). */
void tc_unwind_safe_init(void);

/* From a SA_SIGINFO handler: nonzero when the interrupted PC, or a return
 * address on its frame-pointer chain, lies in a recorded range -- the jump
 * must be deferred.  Zero when the stack is clear, and also on a platform
 * where the context cannot be read (the handler then behaves as before).
 * Async-signal-safe: reads registers and stack words only. */
int tc_interrupted_in_libc(const void* ucontext);

#endif /* MATHILDA_TC_UNWIND_SAFE_H */
