// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_RESEARCH_GAUSSIAN_CONVERT_H
#define INKE_RESEARCH_GAUSSIAN_CONVERT_H
#include "b-poke-basis.h"
#include "a-kani-geometry.h"
/* One fixed-parameter public context; zero-initialize before init. It owns
 * Fp2 and one fixed theta marking of E0. Its setup cost is public precomputation,
 * and must be treated consistently in performance comparisons. */
typedef struct {root_extension field;root_ec_theta4 model;int field_ready,ready;} b_gaussian_converter;
enum {B_GAUSSIAN_CONVERT_OK=1,B_GAUSSIAN_CONVERT_RETRY=0,
 B_GAUSSIAN_CONVERT_INVALID=-1,B_GAUSSIAN_CONVERT_RNG_ERROR=-2,
 B_GAUSSIAN_CONVERT_GEOMETRY_ERROR=-3,B_GAUSSIAN_CONVERT_FACTOR_ERROR=-4,
 B_GAUSSIAN_CONVERT_VERIFY_ERROR=-5,B_GAUSSIAN_CONVERT_MEMORY_ERROR=-6};
typedef struct {unsigned phase,detail,components_completed;} b_gaussian_convert_stage;
int b_gaussian_converter_init(b_gaussian_converter *context);
void b_gaussian_converter_clear(b_gaussian_converter *context);
/* The input ideal and verified same-I certificate are immutable. The caller
 * must preserve I on every failure: RETRY means retry auxiliary torsion on
 * this same I/certificate; prefix MISS is handled before this function by
 * the complete converter on the SAME I. No sampler/Keygen state is changed.
 * A finite budget caps only auxiliary point sampling, never selects another I.
 * On success out->two is the image of the exact full linked BASIS_EVEN.
 * Optional three receives the image of exact BASIS_THREE on out->curve under
 * the SAME signed factor map. If NULL those extra payloads are not evaluated.
 * All output values are preserved on failure. This function is synchronous;
 * progress is called after quotient stages, without secret-value logging. */
int b_gaussian_convert(inke_conversion_result *out,ec_basis_t *three,
 const quat_left_ideal_t *I,const b_gaussian_certificate *certificate,
 const b_gaussian_converter *context,a_ec_rng rng,void *rng_context,
 unsigned long budget,a_kani_progress progress,void *progress_context,
 b_gaussian_convert_stage *stage);
#endif
