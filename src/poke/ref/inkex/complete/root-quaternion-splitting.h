#ifndef INKE_RESEARCH_QUATERNION_SPLITTING_H
#define INKE_RESEARCH_QUATERNION_SPLITTING_H
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
/* Four initialized 2x2 matrices for 1,i,(i+j)/2,(1+ij)/2 modulo2^bits.
 * p must be 7mod8; i²=-1,j²=-p,ij=-ji.  This is an algebraic splitting,
 * not a claim that it already matches an elliptic torsion basis. */
int root_quaternion_splitting(fmpz_mat_struct out[4],const fmpz_t p,unsigned bits);
/* Known actions of ALL FOUR maximal-order basis elements must be supplied
 * in ONE actual torsion basis. Align the high-precision splitting with them
 * modulo2^known_bits, checking every resulting matrix.  No extra extension
 * field torsion is required. known_bits>=1, high_bits>=known_bits. */
int root_quaternion_splitting_align(fmpz_mat_struct out[4],
 const fmpz_mat_struct known[4],const fmpz_t p,unsigned known_bits,unsigned high_bits);
/* Convert (c0+c1*i+c2*j+c3*ij)/den to the same coherent matrix.  Supports
 * den1/2 and rejects elements outside the stated maximal order. */
int root_quaternion_splitting_element(fmpz_mat_t out,const fmpz *coeff,unsigned den,
 const fmpz_mat_struct basis[4],unsigned bits);
#endif
