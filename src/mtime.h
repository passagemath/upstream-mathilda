#ifndef MATHILDA_MTIME_H
#define MATHILDA_MTIME_H

/*
 * Portable timing shim -- the single seam where Mathilda's clocks and sleep
 * cross the OS boundary.
 *
 * POSIX (Linux/glibc, macOS) uses clock_gettime / nanosleep; native Windows
 * (MinGW UCRT64, with no Cygwin POSIX layer) uses QueryPerformanceCounter,
 * GetSystemTimePreciseAsFileTime and Sleep.  Before this module the
 * clock_gettime(CLOCK_MONOTONIC)-with-clock()-fallback idiom was copy-pasted in
 * five places (core.c's TimeConstrained deadline, datetime.c, simp_util.c,
 * integrate.c, integrate_derivdivides.c), none of which compiled on MinGW --
 * centralising them here gives the Windows port a single file to carry and
 * removes the duplication.  See docs/building.md (§Windows) and SPEC.md §10.
 */

/*
 * Seconds from an unspecified monotonic origin.  Never decreases and never
 * fails (falls back to CPU clock() if the OS monotonic clock is unavailable),
 * so callers need no error path.  For measuring elapsed intervals -- NOT a
 * wall-clock reading.
 */
double mth_mono_seconds(void);

/*
 * Wall-clock seconds since the Unix epoch (1970-01-01 00:00 UTC), with
 * sub-second resolution where the OS provides it.  Backs AbsoluteTime and
 * friends.  Can move backwards across an NTP step -- use mth_mono_seconds for
 * durations.
 */
double mth_wall_seconds(void);

/*
 * Sleep for approximately `seconds` (values <= 0 return immediately).  On POSIX
 * a signal can cut nanosleep short; the loop retries with the remaining time so
 * the full interval elapses.  On Windows the request is rounded up to whole
 * milliseconds so the sleep never undershoots.
 */
void mth_sleep_seconds(double seconds);

#endif /* MATHILDA_MTIME_H */
