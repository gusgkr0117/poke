// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_FAST_KEYGEN_H
#define INKE_FAST_KEYGEN_H
#include "root-extension.h"
#include <inke_conversion.h>

/* Default work per Keygen step. A miss resumes fast search at the SAME q;
 * this is not a global cap and does not invoke an ideal-based fallback. */
#define B_FAST_DEFAULT_TRIALS 20000UL
typedef struct {
    unsigned long norm_candidates,witnesses,conversion_failures,
                  guard_misses,interval_misses;
} b_fast_counts;

/* Absolute cumulative norm-candidate limit, including earlier calls interrupted
 * by RNG failure. The caller increases it for each new work batch. Interval
 * sampling is exact with expected <2 draws; the Cornacchia search is bounded.
 * The output is committed only after exact norm, O0 membership and primitivity
 * checks. No division by content is permitted. A non-prime residual cannot
 * cause an invalid output because the final checks are exact.
 * OK/RETRY/INVALID/RNG_ERROR use the INKE_SAMPLE_* values. */
int b_fast_sample_theta(quat_alg_elem_t *theta,const ibz_t *target,
    unsigned long candidate_limit,b_fast_counts *counts,root_random_bytes rng,void *context);

/* For verified primitive theta of norm d*3^b, compute
 * phi=eta o theta/[3^b], ker eta=theta(E0[3^b]). Retain the exact configured
 * full 2-power basis images and verify their degree-d Weil pairing. */
int b_fast_convert_theta(inke_conversion_result *out,const quat_alg_elem_t *theta,
    const ibz_t *d,int *stage);

/* Fixed original q; no q sampling, root swap, model selection, or masks here.
 * All non-success returns preserve out. RETRY is a work-batch miss;
 * RNG_ERROR/INVALID must be surfaced, not interpreted as permission to choose
 * another q. The caller caches a successful conversion across finish retries. */
int b_fast_receiver(inke_conversion_result *out,const ibz_t *q,unsigned int a,
    unsigned long candidate_limit,b_fast_counts *counts,root_random_bytes rng,void *context,
    int *stage);
#endif
