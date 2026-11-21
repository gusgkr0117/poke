#ifndef INKE_RESEARCH_THETA_BALANCED_H
#define INKE_RESEARCH_THETA_BALANCED_H
#include "a-theta-multifield.h"
/* Binary strategy with current-image frames; no source-preimage recursion.
 * All contracts of a_theta_multifield_init apply. The embedded state is
 * owned by this wrapper and must not be stepped/cleared independently.
 * Original payload groups retain their indices; the last group is private.
 * Counter values count group-level point operations, not field operations. */
typedef struct {
    a_theta_multifield state;
    unsigned *height,depth,capacity,max_depth;
    size_t stack_group;
    unsigned long long frame_doublings,frame_evaluations;
    int ready;
} b_theta_balanced;
int b_theta_balanced_init(b_theta_balanced *out,unsigned g,unsigned exponent,
    const root_extension *base,
    const fq_struct *null,const fq_struct *const *generators,
    const a_theta_payload_input *payloads,size_t group_count,
    const fq_t zeta4,const fq_t zeta8);
/* One atomic quotient stage. A failed attempt may retain newly computed
 * strategy frames, but does not change the curve or any original payload.
 * Such cached frames are valid on that unchanged source; retry is safe. */
int b_theta_balanced_step(b_theta_balanced *state);
void b_theta_balanced_clear(b_theta_balanced *state);
#endif
