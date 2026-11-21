#ifndef INKE_RESEARCH_THETA_BASECHANGE_H
#define INKE_RESEARCH_THETA_BASECHANGE_H
#include "root-field-embedding.h"
#include "b-theta4-two.h"

/* Owned extension of an already constructed quotient.  The good affine
 * kernel lifts and trace chart are embedded literally: no new square roots,
 * kernel choices, normalizations, or markings are introduced. */
typedef struct {
    b_theta4_ctx source;
    b_theta4_two_ctx quotient;
    int ready;
} a_theta_two_basechange;

/* out must initially be zero; the given quotient's field must be L->sub.
 * On failure out remains empty. */
int a_theta_two_basechange_init(a_theta_two_basechange *out,
    const b_theta4_two_ctx *quotient,const root_field_embedding *L);
void a_theta_two_basechange_clear(a_theta_two_basechange *state);
#endif
