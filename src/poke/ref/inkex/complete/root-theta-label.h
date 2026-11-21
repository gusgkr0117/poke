#ifndef INKE_RESEARCH_THETA_LABEL_H
#define INKE_RESEARCH_THETA_LABEL_H
#include "b-theta4.h"
typedef struct {
 const b_theta4_ctx *theta;
 uint32_t pivot,index[8];
 unsigned char inverse[8][8];
 fq_t zeta;
 fq_struct *ratios;
 int ready;
} root_theta_label_ctx;
/* The theta context must outlive this label context. Requires an actual
 * level4 embedding; support differences then span (Z/4)^g. */
int root_theta_label_init(root_theta_label_ctx *L,const b_theta4_ctx *T,const fq_t zeta);
void root_theta_label_clear(root_theta_label_ctx *L);
/* Find the unique alpha,beta with P projectively H(alpha,beta)null.
 * Enumerates only4^g shifts. The beta label is solved from g fixed probes;
 * every accepted result is verified against all4^g coordinates. No roots.
 * Outputs remain unchanged on failure, including for non-4-torsion input. */
int root_theta_label(uint32_t *alpha,uint32_t *beta,const fq_struct *P,
                       const root_theta_label_ctx *L);
#endif
