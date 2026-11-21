// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_RESEARCH_UNIFORM_KEYGEN_WRAPPER_H
#define INKE_RESEARCH_UNIFORM_KEYGEN_WRAPPER_H
#include "b-gaussian-convert.h"
#include <inke_keygen.h>
/* The fallback must convert precisely I, not an equivalent ideal or a fresh
 * sampled ideal. Its own resumable state may live in context. A successful
 * output obeys the exact fixed-basis inke_conversion_result contract. This
 * callback boundary is not an assertion that a complete fallback is linked. */
typedef int (*b_same_ideal_fallback)(inke_conversion_result *out,
    const quat_left_ideal_t *I,void *context);
typedef struct {
    /* mixed_trials is a cumulative candidate cap for this original ideal,
     * including every failed witness conversion and resumed RNG call. */
    unsigned long prime_trials,norm_trials,unit_trials,mixed_trials,point_trials;
    b_same_ideal_fallback fallback;void *fallback_context;
    inke_product_fallback product_fallback;void *product_context;
} b_uniform_keygen_options;
enum {B_KEYGEN_NEEDS_FALLBACK=-10,B_KEYGEN_CONVERSION_ERROR=-11};
typedef struct {
    ibz_t q;quat_left_ideal_t ideal;inke_mixed_state mixed;
    b_gaussian_certificate gaussian;inke_conversion_result conversion;
    inke_sample_counts sampling;b_gaussian_counts gaussian_counts;
    unsigned long mixed_witnesses,mixed_failures,gaussian_attempts,gaussian_failures,
                  fallback_attempts,finish_attempts,finish_failures;
    unsigned int a;
    int has_q,has_ideal,mixed_ready,mixed_done,gaussian_prepared,
        gaussian_failed,has_conversion,conversion_origin,complete,phase,last_stage;
    b_gaussian_convert_stage conversion_stage;
} b_uniform_keygen;
void b_uniform_keygen_init(b_uniform_keygen *state);
void b_uniform_keygen_clear(b_uniform_keygen *state);
/* Every call resumes the SAME owned q and original ideal. Sampling commits
 * the ideal BEFORE any conversion preparation. Bounded-prefix misses move
 * to the next algorithm; they never reset the sampled input. A successful
 * conversion is cached across finish/masking retries. sk/pk are atomic.
 * Supply the same options throughout a run. All statuses except OK preserve
 * sampled q/I. NEEDS_FALLBACK explicitly means no complete fallback is linked.
 * Public plans/converter remain caller-owned and outlive this invocation. */
int b_uniform_keygen_step(inke_sk_t *sk,inke_pk_t *pk,b_uniform_keygen *state,
    const b_uniform_keygen_options *options,const b_gaussian_plans *plans,
    const b_gaussian_converter *converter,a_ec_rng rng,void *rng_context,
    a_kani_progress progress,void *progress_context);
#endif
