#ifndef INKE_RESEARCH_THETA_SYMPLECTIC_H
#define INKE_RESEARCH_THETA_SYMPLECTIC_H
#include "b-theta4.h"
#define ROOT_SP_MAX_G 8
#define ROOT_SP_MAX_OPS 1024
enum {ROOT_SP_FOURIER=1,ROOT_SP_SWAP,ROOT_SP_SCALE,ROOT_SP_ADD,ROOT_SP_LOWER,ROOT_SP_UPPER};
typedef struct {int type,i,j,c;} root_sp_op;
typedef struct {unsigned g,count;root_sp_op ops[ROOT_SP_MAX_OPS];} root_sp_word;
/* Matrices have fixed stride16; labels are (alpha[0:g],beta[0:g]). */
typedef unsigned char root_sp_matrix[16][16];
void root_sp_identity(root_sp_matrix M,unsigned g);
int root_sp_valid(const root_sp_matrix M,unsigned g);
void root_sp_word_matrix(root_sp_matrix M,const root_sp_word *W);
/* Complete primitive isotropic columns K[:,0:g] to a symplectic matrix.
 * No division by2 or finite-field-only Gaussian elimination is used. */
int root_sp_complete_lagrangian(root_sp_matrix out,const root_sp_matrix K,unsigned g);
int root_sp_inverse(root_sp_matrix out,const root_sp_matrix M,unsigned g);
/* Decomposes any symplectic matrix over Z/4Z, with exact final check. */
int root_sp_decompose(root_sp_word *out,const root_sp_matrix M,unsigned g);
/* A map U satisfying U H(v) U^-1 proportional to H(Mv), in the exact
 * Heisenberg convention of b_theta4. Eighth root must have exact order8;
 * its square MUST equal the fourth root used for b_theta4_action and the
 * elliptic torsion-label identification (order alone does not fix its sign).
 * The same linear map MUST be applied to null and every evaluation point.
 * inverse applies the literal inverse word, including its phases/scales. */
int root_sp_apply(fq_struct *out,const fq_struct *in,const root_sp_word *W,
                  const fq_t eighth_root,int inverse,const fq_ctx_t field);
void root_sp_label(uint32_t *alpha,uint32_t *beta,const root_sp_matrix M,
                    uint32_t old_alpha,uint32_t old_beta,unsigned g);
/* DL Algorithm3: for an order8 point with 2P labelled(alpha,0),
 * compute H(alpha,0)Inv(P)=sign*P with sign in{-1,1}. */
int root_theta4_compatibility_sign(int *sign,const fq_struct *P,uint32_t alpha,
                                  const b_theta4_ctx *T);
/* For each basis point with sign=-1 set its bit in mask. Apply this SAME
 * marking change to null, basis, pair lifts, and every evaluation point. */
int root_theta4_compatibility_apply(fq_struct *out,const fq_struct *in,
                                   unsigned mask,unsigned g,const fq_ctx_t field);
#endif
