#ifndef INKE_RESEARCH_THETA_FACTOR_H
#define INKE_RESEARCH_THETA_FACTOR_H
#include "root-ec-theta4.h"
#include "b-theta4.h"
typedef struct {
 root_ec_theta4 elliptic;
 const root_extension *E;
 unsigned g;
 size_t count;
 int ready;
} root_theta_factor;
/* The caller has already applied the theta change adapting the FIRST
 * elliptic symplectic pair to the desired factor. This checks every null
 * tensor identity and recovers a signed Montgomery model without new roots.
 * It does not guess or search for a desired elliptic factor. */
int root_theta_factor_init(root_theta_factor *F,const fq_struct *adapted_null,
                             unsigned g,const root_extension *E);
void root_theta_factor_clear(root_theta_factor *F);
/* Uses a nonzero column independently for each point, checking the full
 * tensor rank1 identity and the signed elliptic inverse. Zero columns are
 * skipped; outputs unchanged on failure. */
int root_theta_factor_point(a_ec_point *out,const fq_struct *adapted_point,
                              const root_theta_factor *F);
#endif
