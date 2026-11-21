#ifndef INKE_RESEARCH_EC_TATE_H
#define INKE_RESEARCH_EC_TATE_H
#include "a-extension-quaternion.h"
/* Exact two-dimensional logarithm in a supplied E[2^bits] basis, using a
 * recursive low/high-bit split. Curve and basis validity, final point
 * reconstruction, and coefficient range are checked. No integer factoring
 * or generic discrete logarithm is used. Outputs unchanged on failure. */
int root_ec_dlog2(fmpz_t x,fmpz_t y,const a_ec_point *P,const a_ec_point basis[2],
                    unsigned bits,const fq_t A,const root_extension *E);
/* Matrices of 1,i,(i+j)/2,(1+ij)/2 in the supplied E0[2^bits] basis.
 * half[0:2] are the compatible input halves. All output matrices are
 * initialized2x2; outputs unchanged on failure. */
int root_ec_maxorder_actions(fmpz_mat_struct out[4],const a_ec_point basis[2],
 const a_ec_point half[2],unsigned bits,const root_extension *E);
#endif
