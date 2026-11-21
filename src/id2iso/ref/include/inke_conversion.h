// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_CONVERSION_H
#define INKE_CONVERSION_H
#include <inke_mixed.h>
#include <ec.h>

/* A verified fast conversion result. The full configured 2-primary basis
 * is retained; callers can then restrict it to the manuscript's 2^a. */
typedef struct { ec_curve_t curve; ec_basis_t two; int product_steps; } inke_conversion_result;

int inke_endomorphism_basis(ec_basis_t *out,const ec_basis_t *basis,
    const quat_alg_elem_t *endo,const ibz_t *modulus);
int inke_pairing_degree(const ec_basis_t *source,const ec_curve_t *E,
    const ec_basis_t *image,const ec_curve_t *F,int exponent,const ibz_t *degree);
int inke_convert_witness(inke_conversion_result *out,
    const inke_mixed_state *state,const inke_mixed_witness *witness,int *stage);
#endif
