// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_RESEARCH_POKE_BASIS_H
#define INKE_RESEARCH_POKE_BASIS_H
#include "b-gaussian-preparation.h"
#include "a-extension-ec.h"
#include <inke_conversion.h>
/* The field prime must equal the linked POKE parameter prime. All fq and
 * a_ec_point outputs are initialized by the caller. Failures preserve outputs.
 * Lift P,Q with the relative sign fixed by the supplied P-Q. Only a common
 * sign remains; no basis change, restriction, or random sampling occurs. */
int b_poke_basis_lift(fq_t A,a_ec_point out[3],const ec_curve_t *curve,
                     const ec_basis_t *basis,const root_extension *E);
/* Export ONE model and the signed P,Q,P-Q triple, requiring F_(p^2) descent.
 * The curve is normalized to C=1 and its entire cache/storage is reset. */
int b_poke_basis_export(ec_curve_t *curve,ec_basis_t *basis,const fq_t A,
                       const a_ec_point points[3],const root_extension *E);
/* Evaluate [scale]alpha on the exact linked BASIS_EVEN (three=0) or
 * BASIS_THREE (three=1), using their precomputed signed O0 action matrices.
 * alpha must be integral in O0. This handles half-integral quaternions on
 * full 2-power torsion without requesting a point half in a larger field. */
int b_poke_fixed_action(a_ec_point out[3],int three,
                       const quat_alg_elem_t *alpha,const ibz_t *scale,
                       const root_extension *E);
/* Preparation boundary for a labelled same-I Kani converter. Verifies the
 * certificate against I, then returns [n1^-1]transport(P,Q,P-Q) on BOTH exact
 * fixed bases, at full TORSION_PLUS_EVEN_POWER and POWER_OF_3 respectively.
 * Feed each point as (point,O) into the SAME labelled quotient/factor map.
 * The two initialized output arrays must be disjoint. */
int b_gaussian_fixed_payloads(a_ec_point two[3],a_ec_point three[3],
    const quat_left_ideal_t *I,const b_gaussian_certificate *certificate,
    const root_extension *E);
/* Final boundary after applying that one labelled quotient/factor map.
 * Checks full fixed-basis 2-power orders and degree-scaled Weil pairing.
 * This check does not substitute for the caller's same-I labelled-map proof.
 * steps records the number of product-isogeny stages executed. */
int b_poke_conversion_export(inke_conversion_result *out,const fq_t A,
    const a_ec_point two[3],const ibz_t *degree,int steps,const root_extension *E);
#endif
