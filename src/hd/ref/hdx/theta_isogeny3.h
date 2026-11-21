// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_THETA_ISOGENY3_H
#define INKE_THETA_ISOGENY3_H

#include "theta_structure.h"

/* Generic (3,3) step, from Corte-Real Santos--Costello--Smith,
 * Research in Number Theory 11 (2025), Sections 4.1 and 4.2,
 * https://doi.org/10.1007/s40993-024-00600-y.
 * Coefficient coordinates/sign are checked against the authors' companion
 * implementation (SageMath/isogeny33.py in mariascrs/KummerIsogenies).
 * This C module is implemented in POKE's own field and theta types.
 *
 * The caller supplies points on the stated surface spanning an isotropic
 * subgroup. This routine is not an untrusted-point validation interface.
 * Unsupported coordinate charts return 0; they must not trigger resampling
 * of the ideal whose conversion is being attempted. */
typedef struct {
    fp2_t c[5];
    theta_structure_t codomain;
    int ready;
} inke_theta_isogeny3;

int inke_theta_equal(const theta_point_t *P, const theta_point_t *Q);
/* [3]P via diagonal biquadratic forms. Returns 0 for unsupported zero
 * coordinates, rather than passing an invalid projective point onward. */
int inke_theta_triple(theta_point_t *out, theta_structure_t *domain,
    const theta_point_t *P);
int inke_theta_isogeny3_compute(inke_theta_isogeny3 *out,
    const theta_structure_t *domain, const theta_point_t *R,
    const theta_point_t *S);
int inke_theta_isogeny3_eval(theta_point_t *out,
    const inke_theta_isogeny3 *phi, const theta_point_t *P);

#endif
