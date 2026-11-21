#ifndef INKE_RESEARCH_THETA_MULTIFIELD_H
#define INKE_RESEARCH_THETA_MULTIFIELD_H
#include "a-theta-basechange.h"
#include "root-theta-symplectic.h"

typedef struct {
    const root_extension *field;
    const fq_struct *const *points;
    size_t count;
} a_theta_payload_input;
typedef struct {
    const root_extension *field;
    root_field_embedding embedding;
    fq_struct *points;
    size_t count;
} a_theta_payload_group;
typedef struct {
    const root_extension *base,*kernel_field;
    root_field_embedding kernel_embedding;
    b_theta4_ctx theta;
    fq_t zeta4,zeta8;
    fq_struct *generators;
    a_theta_payload_group *groups;
    size_t group_count;
    unsigned remaining,completed;
    int ready,last_error;
} a_theta_multifield;

/* The theta null and roots are over base=Fp².  Generator vectors are over
 * kernel_field and lift the desired kernel to an isotropic direct summand
 * of order 2^(exponent+2).  Payloads use the literal base-changed marking.
 * All fields/points must satisfy the scalar-Frobenius and adapted-kernel
 * promise documented in reviewer-theta-rationality.md.  In particular the
 * current projective order-eight data descend, and rational good lifts
 * exist after symmetric compatibility.  Invalid promises are diagnosed;
 * no new ideal, kernel, field or random root is silently substituted.
 * The caller-owned fields must outlive the state. out initially is zero.
 * last_error=-2 identifies RNG failure during field embedding; -1 identifies
 * other failures, and 0 success. Failed init retains this after cleanup. */
int a_theta_multifield_init(a_theta_multifield *out,unsigned g,unsigned exponent,
    const root_extension *base,const root_extension *kernel_field,
    const fq_struct *null,const fq_struct *const *generators,
    const a_theta_payload_input *payloads,size_t group_count,
    const fq_t zeta4,const fq_t zeta8);
void a_theta_multifield_clear(a_theta_multifield *state);
/* One base-field quotient, cloned literally to every registered field.
 * No field compositum is used.  1=success,0=invariant/allocation failure or
 * no remaining stage. Failure preserves all point data and the current
 * target (apart from benign chart caching). */
int a_theta_multifield_step(a_theta_multifield *state);
#endif
