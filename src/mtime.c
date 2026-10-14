/*
 * Portable timing shim.  See mtime.h for the contract.
 *
 * clock_gettime and nanosleep are POSIX, not C99: glibc hides them under
 * -std=c99 unless a feature-test macro is defined BEFORE the first include
 * (below it the header is already parsed with the wrong namespace and the macro
 * does nothing).  macOS exposes them implicitly, so getting this wrong compiles
 * clean there and breaks only on Linux -- the exact class SPEC.md §10 is about.
 * The whole POSIX path is skipped on _WIN32, where none of it exists.
 */
#ifndef _WIN32
#  ifndef _POSIX_C_SOURCE
#  define _POSIX_C_SOURCE 200809L
#  endif
#endif

#include "mtime.h"
#include <time.h>          /* clock_gettime, nanosleep, clock, CLOCKS_PER_SEC */

#ifdef _WIN32
#  include <windows.h>     /* QueryPerformance*, GetSystemTimePreciseAsFileTime, Sleep */
#else
#  include <errno.h>       /* EINTR for the nanosleep retry */
#endif

double mth_mono_seconds(void) {
#ifdef _WIN32
    LARGE_INTEGER freq, now;
    if (QueryPerformanceFrequency(&freq) && freq.QuadPart != 0 &&
        QueryPerformanceCounter(&now))
        return (double)now.QuadPart / (double)freq.QuadPart;
    /* No high-resolution counter: CPU time is at least a duration. */
    return (double)clock() / (double)CLOCKS_PER_SEC;
#else
    struct timespec ts;
    /* CLOCK_MONOTONIC is mandatory in POSIX.1-2001 and present on every host we
     * target (Linux/glibc, macOS 10.12+, WSL 1+2); treat a failure as "fall
     * back to CPU time" rather than propagate an error callers cannot use. */
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
        return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
    return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

double mth_wall_seconds(void) {
#ifdef _WIN32
    FILETIME ft;
    GetSystemTimePreciseAsFileTime(&ft);
    unsigned long long ticks = ((unsigned long long)ft.dwHighDateTime << 32)
                             |  (unsigned long long)ft.dwLowDateTime;
    /* FILETIME counts 100-ns ticks since 1601-01-01; 11644473600 s separate
     * that epoch from the Unix epoch (1970-01-01). */
    return (double)ticks / 1e7 - 11644473600.0;
#else
    struct timespec ts;
    if (clock_gettime(CLOCK_REALTIME, &ts) == 0)
        return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
    return (double)time(NULL);   /* whole-second last resort */
#endif
}

void mth_sleep_seconds(double seconds) {
    if (!(seconds > 0.0)) return;
#ifdef _WIN32
    /* Sleep() takes whole milliseconds; round up so we never undershoot, and
     * clamp below DWORD max (~49.7 days) so the cast cannot wrap. */
    double ms = seconds * 1000.0 + 0.5;
    if (ms > 4294967294.0) ms = 4294967294.0;
    Sleep((DWORD)ms);
#else
    struct timespec req, rem;
    req.tv_sec  = (time_t)seconds;
    req.tv_nsec = (long)((seconds - (double)req.tv_sec) * 1e9);
    if (req.tv_nsec < 0) req.tv_nsec = 0;
    if (req.tv_nsec > 999999999L) req.tv_nsec = 999999999L;
    /* A signal can interrupt nanosleep; resume with the remainder so the whole
     * interval elapses (matches the former Pause[] loop). */
    while (nanosleep(&req, &rem) != 0 && errno == EINTR)
        req = rem;
#endif
}
