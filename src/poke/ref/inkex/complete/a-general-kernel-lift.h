#ifndef INKE_RESEARCH_GENERAL_KERNEL_LIFT_H
#define INKE_RESEARCH_GENERAL_KERNEL_LIFT_H
#include "a-extension-pairing.h"
/* One signed half on y^2=x^3+A*x^2+x. The curve must have all its 2-torsion
 * rational in E. This is deterministic with the field's square-root cache,
 * checks doubling exactly, and preserves out on failure (aliasing allowed). */
int a_ec_half_montgomery(a_ec_point *out,const a_ec_point *P,const fq_t A,
 const root_extension *E);
/* Input K[generator][elliptic factor] is the given rank-two graph in
 * (E1 x E2)[2^exponent]. Its FIRST elliptic components must form a basis.
 * Both curves are maximal Fp2 models; caller embeds them in a field containing
 * their full 2^(exponent+2)-torsion (relative degree4 always suffices when
 * 2^exponent divides p+1). No isogeny or quaternion is needed.
 * Output generates an isotropic subgroup isomorphic to
 * (Z/2^(exponent+2)Z)^2 (each generator has exact order2^(exponent+2)),
 * with [4]out=K exactly as SIGNED points. It is suitable for a length-exponent
 * level-four binary chain. Supplied out points must already be initialized.
 * RETRY/RNG_ERROR/INVALID preserve all outputs; budget0 is Las Vegas. */
int a_general_kernel_lift(a_ec_point out[2][2],const a_ec_point K[2][2],
 unsigned exponent,const fq_t A1,const fq_t A2,a_ec_rng rng,void *context,
 unsigned long budget,const root_extension *E);
#endif
