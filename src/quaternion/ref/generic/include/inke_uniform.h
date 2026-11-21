// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_UNIFORM_H
#define INKE_UNIFORM_H

#include <quaternion.h>

enum inke_sample_status {
    INKE_SAMPLE_OK = 1,
    INKE_SAMPLE_RETRY = 0,
    INKE_SAMPLE_INVALID = -1,
    INKE_SAMPLE_RNG_ERROR = -2
};

typedef struct {
    unsigned long prime_candidates;
    unsigned long norm_candidates;
    unsigned long unit_candidates;
} inke_sample_counts;

/* Sample uniformly from {q : 2^(a-2) <= q <= 3*2^(a-2), q = 5 (mod 6)};
 * a >= 4 is even. This is the current INKE degree distribution. */
int inke_sample_q(ibz_t *q, unsigned int a);
int inke_q_is_admissible(const ibz_t *q, unsigned int a);

/* Construct (O gamma + d O) unit + d O. The caller supplies a primitive
 * gamma in O whose norm is a multiple of d, and a unit modulo d O.
 * Invalid input leaves out unchanged. out must be initialized; order and
 * algebra must remain alive for the lifetime of out. */
int inke_ideal_from_unit(quat_left_ideal_t *out,
                         const quat_alg_elem_t *gamma,
                         const quat_alg_elem_t *unit,
                         const ibz_t *d,
                         const quat_order_t *order,
                         const quat_alg_t *algebra);

/* One bounded sampling attempt with q fixed. This implements only ideal
 * sampling, not ideal-to-isogeny conversion or Keygen. A RETRY must be
 * retried with the same q. On success the returned ideal must be retained
 * throughout all subsequent conversion attempts.
 *
 * order is the maximal order Z<1,i,(i+j)/2,(1+ij)/2> of the manuscript.
 * norm_trials bounds the norm search separately for each auxiliary m.
 * The return value distinguishes budget exhaustion from invalid input and
 * RNG failure. All non-success returns leave out unchanged.
 * Auxiliary primality tests use GMP; exact norm and lattice checks do not
 * rely on a primality-test result being correct. */
int inke_sample_ideal(quat_left_ideal_t *out,
                      const ibz_t *q, unsigned int a,
                      const quat_order_t *order,
                      const quat_alg_t *algebra,
                      unsigned long prime_trials,
                      unsigned long norm_trials,
                      unsigned long unit_trials,
                      inke_sample_counts *counts);

#endif
