/*
 * fi_par.c — the parallel region shared by the quadratic sieve, ECM, p+1 and
 * the CFRAC / Dixon relation collectors, and the work-claim counter the
 * latter use to balance candidates of uneven cost.
 *
 * fi_run_parallel(fn, jobs, stride, n) runs fn on n independent jobs, one
 * thread each.  The threading discipline (the same as src/graph/gmet_core.c):
 *
 *   - workers are created with every signal blocked, so an asynchronous
 *     signal (TimeConstrained's SIGPROF, SIGINT) is only ever delivered on
 *     the main thread;
 *   - the region is bracketed by tc_async_defer_push/pop, and the main
 *     thread does no work of its own: it waits in pthread_join, so a timeout
 *     arriving mid-region is recorded and taken at the main thread's next
 *     safe point after every worker has finished (the guarded GMP allocators
 *     only jump on the main thread; see tc_alloc_safepoint in src/core.c);
 *   - a job whose thread could not be created runs inline after the joins.
 *
 * Callers keep a region short (one polynomial family per QS worker, one
 * curve per ECM worker, one seed per p+1 worker, a few milliseconds of
 * candidates for CFRAC / Dixon) so a timeout is honoured promptly.  Without
 * MATHILDA_THREADS everything runs serially on the calling thread.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L        /* sysconf */
#endif
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE 1             /* _SC_NPROCESSORS_ONLN under _POSIX_C_SOURCE */
#endif

#include "fi_internal.h"
#include "core.h"

#include <stdlib.h>

#ifdef MATHILDA_THREADS
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#endif

int fi_nthreads(void) {
#ifdef MATHILDA_THREADS
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpu < 1) ncpu = 1;
    if (ncpu > FI_MAXTHREADS) ncpu = FI_MAXTHREADS;
    unsigned nt = (unsigned)ncpu;
    FI_TUNE_U("FI_THREADS", nt);
    if (nt < 1) nt = 1;
    if (nt > FI_MAXTHREADS) nt = FI_MAXTHREADS;
    return (int)nt;
#else
    return 1;
#endif
}

#ifdef MATHILDA_THREADS
typedef struct { void (*fn)(void*); void* arg; } FiJob;

static void* fi_job_run(void* p) {
    FiJob* j = (FiJob*)p;
    j->fn(j->arg);
    return NULL;
}
#endif

void fi_run_parallel(void (*fn)(void*), void* jobs, size_t stride, int n) {
#ifdef MATHILDA_THREADS
    if (n > 1) {
        pthread_t th[FI_MAXTHREADS];
        FiJob job[FI_MAXTHREADS];
        int started[FI_MAXTHREADS];
        if (n > FI_MAXTHREADS) n = FI_MAXTHREADS;
        sigset_t all, old;
        sigfillset(&all);
        pthread_sigmask(SIG_SETMASK, &all, &old);
        tc_async_defer_push();
        for (int i = 0; i < n; i++) {
            job[i].fn = fn;
            job[i].arg = (char*)jobs + (size_t)i * stride;
            started[i] = pthread_create(&th[i], NULL, fi_job_run, &job[i]) == 0;
        }
        pthread_sigmask(SIG_SETMASK, &old, NULL);
        for (int i = 0; i < n; i++) if (started[i]) pthread_join(th[i], NULL);
        for (int i = 0; i < n; i++) if (!started[i]) fn(job[i].arg);
        tc_async_defer_pop();
        return;
    }
#endif
    for (int i = 0; i < n; i++) fn((char*)jobs + (size_t)i * stride);
}

/* Work claiming.  A claim is one lock round trip per chunk; callers size the
 * chunk so that is noise (tens of microseconds of work per chunk). */
void fi_claim_init(FiClaim* q) {
    q->next = q->end = 0;
    q->lock = NULL;
#ifdef MATHILDA_THREADS
    pthread_mutex_t* m = malloc(sizeof(pthread_mutex_t));
    if (m && pthread_mutex_init(m, NULL) == 0) q->lock = m;
    else free(m);
#endif
}

void fi_claim_reset(FiClaim* q, size_t end) {
    q->next = 0;
    q->end = end;
}

int fi_claim_next(FiClaim* q, size_t chunk, size_t* lo, size_t* hi) {
#ifdef MATHILDA_THREADS
    if (q->lock) pthread_mutex_lock((pthread_mutex_t*)q->lock);
#endif
    size_t a = q->next, b = a + chunk < q->end ? a + chunk : q->end;
    q->next = b;
#ifdef MATHILDA_THREADS
    if (q->lock) pthread_mutex_unlock((pthread_mutex_t*)q->lock);
#endif
    *lo = a;
    *hi = b;
    return a < b;
}

void fi_claim_destroy(FiClaim* q) {
#ifdef MATHILDA_THREADS
    if (q->lock) {
        pthread_mutex_destroy((pthread_mutex_t*)q->lock);
        free(q->lock);
    }
#endif
    q->lock = NULL;
}
