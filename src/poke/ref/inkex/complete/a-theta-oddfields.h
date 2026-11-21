#ifndef INKE_RESEARCH_THETA_ODDFIELDS_H
#define INKE_RESEARCH_THETA_ODDFIELDS_H
#include "a-theta-oddbasechange.h"
#include "a-theta-multifield.h"
typedef struct {
 const root_extension *field;
 const root_extension *evaluation_field;
 fq_struct *points;size_t count;
 root_extension common;int common_ready;
 root_field_embedding kernel_up,point_up;
} a_theta_odd_payload_group;
typedef struct {
 const root_extension *base,*kernel_field;
 root_field_embedding kernel_embedding;
 b_theta4_ctx theta;
 /* generators keep the remaining prime-power order; complements are the
  * transported initial [ell^(exponent-1)]L, hence always order ell. */
 fq_struct *generators,*complements,*pairing4;
 a_theta_odd_payload_group *groups;size_t group_count;
 ulong ell;unsigned exponent,remaining,completed;int ready,last_error;
} a_theta_oddfields;
/* Full isotropic K,L frames have order ell^exponent. The input pairing4 is
 * e_{ell^exponent}(K_i,L_j)^(4 ell^(exponent-1)), including the polarization
 * multiplier of preceding DIFFERENT prime steps. Caller supplies a genuine
 * symplectic frame and an Fp2-stable kernel. The weighted odd quotient induces
 * the base-defined level-four marking (reviewer-theta-descent.md).
 * Payload fields stay separate: each group uses only its pairwise compositum
 * with the kernel field. All caller-owned fields must outlive the state.
 * last_error is -2 for RNG failure during embedding/compositum construction,
 * -1 for other failures, and 0 on success; failed init retains the status. */
int a_theta_oddfields_init(a_theta_oddfields *out,unsigned g,ulong ell,unsigned exponent,
 const root_extension *base,const root_extension *kernel_field,const fq_struct *null,
 const fq_struct *const *generators,const fq_struct *const *complements,const fq_struct *pairing4,
 const a_theta_payload_input *payloads,size_t groups);
void a_theta_oddfields_clear(a_theta_oddfields *S);
/* Atomic point/target update. Root-free: uses no new randomness. */
int a_theta_oddfields_step(a_theta_oddfields *S);
#endif
