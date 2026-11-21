#ifndef INKE_RESEARCH_THETA4_ODD_H
#define INKE_RESEARCH_THETA4_ODD_H
#include "b-theta4.h"
/* Root-free excellent lifts for odd ell, LR2022 Algorithm1/Corollary4.6.
 * All outputs are preserved on failure. P must be ell-torsion. */
int b_theta4_odd_kernel_scale(fq_t scale_ell,const fq_struct *P,ulong ell,b_theta4_ctx *T);
int b_theta4_odd_coset_scale(fq_t scale_ell,const fq_struct *xP,const fq_struct *P,
    const fq_t kernel_scale_ell,const fq_struct *x,ulong ell,b_theta4_ctx *T);
/* One trace summand for the seed theta_0^ell. The caller must sum over the
 * entire maximal isotropic kernel and reject a zero total null. Complete
 * alternate-section selection is a separate interface, not asserted here. */
int b_theta4_odd_term(fq_struct *out,const fq_struct *xP,const fq_struct *P,
    const fq_struct *x,ulong ell,b_theta4_ctx *T);

typedef struct {
 b_theta4_ctx *source;
 ulong ell;
 size_t kernel_size, selected_chart;
 fq_struct *basis,*translate,*weights,*target_null;
 /* Optional bounded cache. -1 disables it and retains the complete streamed
  * computation; 0 is unattempted, 1 ready. No search law depends on it. */
 fq_struct *cached_kernel,*cached_scales;
 int cache_status;
} b_theta4_odd_ctx;
/* Complete odd-prime quotient. kernel/complement each contain g vectors;
 * together they must be symplectic bases of full ell-torsion, with each
 * half isotropic. pairing4[i*g+j] is e_ell(kernel_i,complement_j)^4, in the
 * principal-polarization convention. All torsion/roots lie in source.field.
 * The theta marking must be the same for every input.
 * At most ell^g complement charts are tried, with ell^g streamed trace
 * terms per chart. A bounded optional cache reuses kernel points and their
 * lift scales; larger inputs use the complete stream. No genericity assumption.
 * On failure Q is empty; success selects one fixed chart for every image. */
int b_theta4_odd_init(b_theta4_odd_ctx *Q,b_theta4_ctx *source,ulong ell,
 const fq_struct *const *kernel,const fq_struct *const *complement,
 const fq_struct *pairing4);
void b_theta4_odd_clear(b_theta4_odd_ctx *Q);
/* Evaluates using only the current level-four coordinates of x. */
int b_theta4_odd_evaluate(fq_struct *out,const fq_struct *x,b_theta4_odd_ctx *Q);
#endif
