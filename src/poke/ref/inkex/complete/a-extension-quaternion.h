#ifndef INKE_RESEARCH_EXTENSION_QUATERNION_H
#define INKE_RESEARCH_EXTENSION_QUATERNION_H
#include "a-extension-ec.h"
#include <flint/fmpz_mat.h>
/* E0:y^2=x^3+x, i(x,y)=(-x,i*y), j=absolute p-Frobenius.
 * coeff[0..3]/den represents c0+c1*i+c2*j+c3*i*j.
 * den is 1 or 2, and the element must lie in O0.
 * order annihilates P. If den=2 and order is even, lift must satisfy
 * [2]lift=P and [2*order]lift=O; no division of an output point is used.
 * With odd order, lift may be NULL and scalar inversion is exact.
 * On invalid input, out is unchanged. */
int a_ec_quaternion(a_ec_point *out,const a_ec_point *P,const fmpz *coeff,
                     unsigned den,const fmpz_t order,const a_ec_point *lift,
                     const root_extension *E);

typedef struct {
    /* generator[2*j+s] has first block [n1]B_s in coordinate j,
     * and second block matrix(:,j)*theta(B_s), j=0..3,s=0,1. */
    a_ec_point generator[8][8];
    /* Complement with e(generator[a],complement[b])=zeta^(delta_ab),
     * using the same ordered source basis in all four factors. */
    a_ec_point complement[8][8];
    a_ec_point basis[2];
    fmpz_t order;
} a_clapoti_kernel;
void a_clapoti_kernel_init(a_clapoti_kernel *K,const root_extension *E);
void a_clapoti_kernel_clear(a_clapoti_kernel *K,const root_extension *E);
/* Construct a prime-power component of the reviewed graph kernel.
 * Caller supplies matrix=Y*X and n1=u*d1 from a verified preparation.
 * Require gcd(n1,ell)=1. For ell=2 and den=2 one extra torsion step
 * is sampled to evaluate the half-integral endomorphism exactly.
 * No geometric quotient or theta structure is produced here. */
int a_clapoti_kernel_component(a_clapoti_kernel *out,unsigned long ell,unsigned e,
                                const fmpz_t n1,const fmpz *theta,unsigned den,
                                const fmpz_mat_t matrix,a_ec_rng rng,void *context,
                                unsigned long budget,unsigned long *trials,
                                const root_extension *E);
/* Isotropic lift to ell^(e+extra). Here ell^e divides n1+n2, n1,n2
 * are ell-adic units. Multiply lifted generators by ell^extra to obtain
 * a basis of the original component. Q-inputs of eta are multiplied by
 * t=-n1/n2 modulo the enlarged order; t=1 modulo the original order.
 * The field must contain enough torsion, including any input-half step. */
int a_clapoti_kernel_lift(a_clapoti_kernel *out,unsigned long ell,unsigned e,
                           unsigned extra,const fmpz_t n1,const fmpz_t n2,
                           const fmpz *theta,unsigned den,const fmpz_mat_t matrix,
                           a_ec_rng rng,void *context,unsigned long budget,
                           unsigned long *trials,const root_extension *E);
#endif
