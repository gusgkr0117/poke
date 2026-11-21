#ifndef INKE_RESEARCH_THETA_MULTIFIELD_H
#define INKE_RESEARCH_THETA_MULTIFIELD_H
#include "b-theta4-two.h"
#include "root-theta-symplectic.h"

typedef struct {
    const fq_struct *const *points;
    size_t count;
} a_theta_payload_input;
typedef struct {
    fq_struct *points;
    size_t count;
} a_theta_payload_group;
typedef struct {
    const root_extension *base;
    b_theta4_ctx theta;
    fq_t zeta4,zeta8;
    fq_struct *generators;
    a_theta_payload_group *groups;
    size_t group_count;
    unsigned remaining,completed;
    int ready;
} a_theta_multifield;

/* The input null, lifted generators and all transported points use the same
 * level-four marking over Fp2. The generators lift an isotropic kernel to
 * order 2^(exponent+2). The caller owns the field and zero-initializes out. */
int a_theta_multifield_init(a_theta_multifield *out,unsigned g,unsigned exponent,
    const root_extension *base,
    const fq_struct *null,const fq_struct *const *generators,
    const a_theta_payload_input *payloads,size_t group_count,
    const fq_t zeta4,const fq_t zeta8);
void a_theta_multifield_clear(a_theta_multifield *state);
/* Compute one quotient and transport every group of points through it.
 * The curve and point data are committed together when the stage succeeds. */
int a_theta_multifield_step(a_theta_multifield *state);
#endif
