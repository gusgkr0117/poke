#ifndef INKE_RESEARCH_KANI_TWO_KERNEL_H
#define INKE_RESEARCH_KANI_TWO_KERNEL_H
#include "a-extension-quaternion.h"
typedef struct {
    a_ec_point generator[2][2],complement[2][2],basis[2];
    fmpz_t order;
} a_kani_two_kernel;
void a_kani_two_kernel_init(a_kani_two_kernel *K,const root_extension *E);
void a_kani_two_kernel_clear(a_kani_two_kernel *K,const root_extension *E);
/* eta is an integral O0 endomorphism of exact norm n1*n2, checked here.
 * ell^e divides n1+n2 and ell is prime to n1*n2. At order ell^(e+extra),
 * the isotropic lift is (n1 P,eta P),(n1 Q,t eta Q), t=-n1/n2.
 * Complement numbering matches these two generators and has diagonal
 * pairing equal to e(P,Q). For den2 and ell2, one extra input-half step
 * is required in the field. Output is preserved on failure. */
int a_kani_two_kernel_lift(a_kani_two_kernel *out,ulong ell,unsigned e,
    unsigned extra,const fmpz_t n1,const fmpz_t n2,const fmpz *eta,unsigned den,
    a_ec_rng rng,void *context,unsigned long budget,unsigned long *trials,const root_extension *E);
#endif
