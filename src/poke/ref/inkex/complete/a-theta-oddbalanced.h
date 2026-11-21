#ifndef INKE_RESEARCH_THETA_ODDBALANCED_H
#define INKE_RESEARCH_THETA_ODDBALANCED_H
#include "a-theta-oddfields.h"
typedef struct {
 a_theta_oddfields state;unsigned *height,depth,capacity,max_depth;size_t stack_group;
 unsigned long long frame_multiplications,frame_evaluations;int ready;
} a_theta_oddbalanced;
int a_theta_oddbalanced_init(a_theta_oddbalanced *S,unsigned g,ulong ell,unsigned exponent,
 const root_extension *base,const root_extension *kernel_field,const fq_struct *null,
 const fq_struct *const *generators,const fq_struct *const *complements,const fq_struct *pairing4,
 const a_theta_payload_input *payloads,size_t groups);
int a_theta_oddbalanced_step(a_theta_oddbalanced *S);
void a_theta_oddbalanced_clear(a_theta_oddbalanced *S);
#endif
