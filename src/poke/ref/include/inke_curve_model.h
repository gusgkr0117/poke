// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_CURVE_MODEL_H
#define INKE_CURVE_MODEL_H
#include <ec.h>

#define INKE_MODEL_MAX_MAPS 6
typedef struct {
    ec_curve_t curve;
    ec_isom_t xmap;
} inke_model_map;

/* All Fp2 isomorphisms from source to normalized Montgomery models, modulo
 * the two choices of y sign. At most six x-coordinate maps occur.
 * Returns 1 on success; invalid input returns 0 without changing outputs. */
int inke_model_enumerate(inke_model_map maps[INKE_MODEL_MAX_MAPS],
    unsigned int *count,const ec_curve_t *source);

/* Retain precisely the maps to the model with lexicographically smallest
 * canonical fp2_encode(A). The caller samples uniformly among these maps.
 * Each map has exactly two y-coordinate lifts, so uniform selection gives
 * the x-only image of a uniform Fp2 isomorphism, including j=0 and j=1728. */
int inke_model_candidates(inke_model_map maps[INKE_MODEL_MAX_MAPS],
    unsigned int *count,const ec_curve_t *source);

/* Apply one selected map to the complete basis and normalize finite points
 * to Z=1. Supports out==basis. Invalid (0:0) data leave out unchanged. */
int inke_model_apply_basis(ec_basis_t *out,const ec_basis_t *basis,
    const inke_model_map *map);
#endif
