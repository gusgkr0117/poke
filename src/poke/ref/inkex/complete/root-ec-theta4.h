#ifndef INKE_RESEARCH_EC_THETA4_H
#define INKE_RESEARCH_EC_THETA4_H
#include "a-extension-ec.h"
#include "b-theta4.h"
/* Signed genus-one level-four coordinates on y^2=x^3+A*x^2+x.
 * All roots defining the theta structure are checked. init returns zero
 * if this structure is not rational in the supplied field or A is singular.
 * The elliptic curve/field must outlive this context. Always clear after init. */
typedef struct {
    const root_extension *E;
    fq_t A,t,a,b,c,k0,k2;
    fq_struct *null;
    int ready;
} root_ec_theta4;
int root_ec_theta4_init(root_ec_theta4 *T,const fq_t A,const root_extension *E);
/* Reconstruct from a smooth symmetric theta null (a,b,c,b), with
 * ac(a^2+c^2)=2b^4. Uses only rational field operations, no square roots.
 * This does not perform a higher-dimensional factor/theta-structure change. */
int root_ec_theta4_from_null(root_ec_theta4 *T,const fq_struct *null,const root_extension *E);
void root_ec_theta4_clear(root_ec_theta4 *T);
/* Polynomial map, including infinity and every finite point. Invalid input
 * returns zero without changing the four initialized output elements. */
int root_ec_theta4_map(fq_struct *out,const a_ec_point *P,const root_ec_theta4 *T);
/* The inverse identifies the signed point, not only its Kummer image. */
int root_ec_theta4_unmap(a_ec_point *out,const fq_struct *P,const root_ec_theta4 *T);
/* Identify a four-torsion point in the 16-element Heisenberg orbit.
 * zeta is the chosen primitive fourth root; output labels are in0..3. */
int root_ec_theta4_label(unsigned *alpha,unsigned *beta,const a_ec_point *P,
                         const fq_t zeta,const root_ec_theta4 *T);
#endif
