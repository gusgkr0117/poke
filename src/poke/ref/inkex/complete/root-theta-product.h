#ifndef INKE_RESEARCH_THETA_PRODUCT_H
#define INKE_RESEARCH_THETA_PRODUCT_H
#include "root-theta-factor.h"
#include "root-theta-symplectic.h"
typedef struct {
 const root_extension *E;
 root_sp_word common,local[2];
 unsigned alpha[2],beta[2];
 root_theta_factor factors[2];
 fq_t eighth;
 size_t decompositions_tried;
 int ready;
} root_theta_product;
/* Exactly160 unordered symplectic decompositions of (Z/4)^4. Columns are
 * (u0,u1,v0,v1); each (ui,vi) spans its indicated elliptic factor. */
int root_theta_product_frames(root_sp_matrix frames[160]);
/* Complete finite marking/chart search for a genus2 principally polarized
 * PRODUCT with a symmetric full level4 theta null over E, and rational
 * compatible elliptic charts (in particular, the current maximal factors
 * with full rational16 torsion satisfy this hypothesis). Both factors and
 * their group origins are recovered in one consistent coordinate change.
 * P must be zero initialized. Failure preserves it. */
int root_theta_product_init(root_theta_product *P,const fq_struct *null,
 const fq_t eighth,const root_extension *E);
void root_theta_product_clear(root_theta_product *P);
/* Both outputs are signed; all coordinates undergo the identical stored
 * marking changes used for null. Failure preserves both initialized outputs. */
int root_theta_product_point(a_ec_point out[2],const fq_struct *point,
 const root_theta_product *P);
#endif
