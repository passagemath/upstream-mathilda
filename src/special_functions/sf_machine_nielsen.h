#ifndef SF_MACHINE_NIELSEN_H
#define SF_MACHINE_NIELSEN_H

#include <stdbool.h>

/* Nielsen's generalized polylogarithm PolyLog[n, p, z] at machine precision,
 * complex n, p, z as (re, im) pairs. The one kernel behind the scalar builtin,
 * the NDArray / packed map and Compile[]. It evaluates exactly the calls the
 * interpreter evaluates and returns false for the rest -- and also where its
 * own quadrature or series cannot reach machine accuracy, so the caller's
 * arbitrary-precision path answers instead. See sf_machine_nielsen.c. */
bool sf_nielsen_c(double nr, double ni, double pr, double pi,
                  double zr, double zi, double* rr, double* ri);

#endif /* SF_MACHINE_NIELSEN_H */
