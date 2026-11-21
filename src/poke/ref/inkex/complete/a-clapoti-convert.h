#ifndef INKE_RESEARCH_CLAPOTI_CONVERT_H
#define INKE_RESEARCH_CLAPOTI_CONVERT_H
#include "a-clapoti-preparation.h"
#include "a-kani-geometry.h"
#include "b-poke-basis.h"
enum {A_CLAPOTI_CONVERT_OK=1,A_CLAPOTI_CONVERT_RETRY=0,
 A_CLAPOTI_CONVERT_INVALID=-1,A_CLAPOTI_CONVERT_RNG_ERROR=-2,
 A_CLAPOTI_CONVERT_GEOMETRY_ERROR=-3,A_CLAPOTI_CONVERT_FACTOR_ERROR=-4,
 A_CLAPOTI_CONVERT_VERIFY_ERROR=-5,A_CLAPOTI_CONVERT_MEMORY_ERROR=-6};
typedef struct {unsigned phase,detail;size_t components_completed;} a_clapoti_convert_stage;
/* The initialized arrays contain three consecutive eight-point vectors.
 * For each exact fixed-basis point P, the first four coordinates are
 * X^t e_0 [n1^-1] beta1(P), and the last four are zero. three may be NULL.
 * The borrowed field has the linked public prime. Outputs are atomic. */
int a_clapoti_fixed_payloads(a_ec_point two[24],a_ec_point three[24],
 const quat_left_ideal_t *I,const a_clapoti_certificate *certificate,
 const root_extension *E);
/* Complete certificate -> labelled genus-eight quotient -> original-I images.
 * The caller retains exactly I, q, and this certificate on auxiliary failure.
 * budget0 gives unbounded auxiliary torsion sampling; a finite budget returns
 * RETRY on the same problem. base_model is a borrowed public E0 marking.
 * Output two is the image of the FULL fixed BASIS_EVEN. Optional three is the
 * image of BASIS_THREE under the SAME factor map. Failures preserve outputs.
 * This adapter is not itself evidence that a full geometric run has passed. */
int a_clapoti_convert(inke_conversion_result *out,ec_basis_t *three,
 const quat_left_ideal_t *I,const a_clapoti_certificate *certificate,
 const root_ec_theta4 *base_model,a_ec_rng rng,void *rng_context,
 unsigned long budget,a_kani_progress progress,void *progress_context,
 a_clapoti_convert_stage *stage);
/* One-key state for b_same_ideal_fallback (matching its callback signature).
 * Clear/reinitialize before a different selected ideal. Algebraic preparation
 * and its certificate are cached across auxiliary geometric retries. */
typedef struct {
 a_clapoti_certificate certificate;a_clapoti_counts counts;
 const root_ec_theta4 *base_model;a_ec_rng rng;void *rng_context;
 unsigned long point_budget;a_kani_progress progress;void *progress_context;
 a_clapoti_convert_stage stage;int ready,prepared,last_status;
} a_clapoti_fallback;
int a_clapoti_fallback_init(a_clapoti_fallback *state,const root_ec_theta4 *base,
 a_ec_rng rng,void *rng_context,unsigned long point_budget,
 a_kani_progress progress,void *progress_context);
void a_clapoti_fallback_clear(a_clapoti_fallback *state);
/* Returns INKE_SAMPLE_*; detailed geometry/error status remains in state. */
int a_clapoti_fallback_convert(inke_conversion_result *out,
 const quat_left_ideal_t *I,void *state);
#endif
