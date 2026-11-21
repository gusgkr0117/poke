#include "root-kani-factor.h"
#include "root-theta-label.h"
#include "root-theta-symplectic.h"
#include <stdlib.h>
int root_kani_extract(fq_t A,a_ec_point *points,size_t count,
 unsigned g,b_theta4_ctx *target,const fq_struct *lattice_points,
 const root_field_embedding *L,const fq_struct *payload,
 const fmpz_mat_t Pnum,const fmpz_mat_t eta,const fmpz_t n1,const fmpz_t N,
 const fq_t zeta4,const fq_t zeta8,const root_extension *E,unsigned *stage){
 if(stage)*stage=1;if(g<2||g>8||(g&1)||!target||target->g!=g||!lattice_points||(count&&(!points||!payload))||(L&&(!L->ready||L->sub!=E))||count>SIZE_MAX/sizeof(a_ec_point))return 0;
 size_t C=target->count;fq_struct *labels=b_theta4_vector_new(2*g*C,E->field),*adapted=b_theta4_vector_new(C,E->field);a_ec_point *answer=count?calloc(count,sizeof(*answer)):NULL;int ok=0,have_labels=0,have_factor=0;size_t initialized=0;root_theta_label_ctx label_ctx;root_theta_factor factor;root_tate_matrix H={{0}},e,S,F,T,inverse;
 if(!labels||!adapted||(count&&!answer))goto done;for(;initialized<count;initialized++)a_ec_init(answer+initialized,E);
 if(stage)*stage=2;for(unsigned j=0;j<2*g;j++){if(L){if(!root_field_projective_down(labels+j*C,lattice_points+j*C,C,L))goto done;}else b_theta4_copy(labels+j*C,lattice_points+j*C,target);}
 if(!root_theta_label_init(&label_ctx,target,zeta4))goto done;have_labels=1;
 for(unsigned j=0;j<2*g;j++){uint32_t alpha,beta;if(!root_theta_label(&alpha,&beta,labels+j*C,&label_ctx))goto done;for(unsigned i=0;i<g;i++){H[i][j]=(alpha>>(2*i))&3;H[i+g][j]=(beta>>(2*i))&3;}}
 if(stage)*stage=3;if(!root_tate_projector_dim(e,S,F,Pnum,eta,n1,N,g)||!root_tate_theta_frame_dim(T,H,F,g)||!root_sp_inverse(inverse,T,g))goto done;root_sp_word word;if(!root_sp_decompose(&word,inverse,g))goto done;
 if(stage)*stage=4;if(!root_sp_apply(adapted,target->null,&word,zeta8,0,E->field)||!root_theta_factor_init(&factor,adapted,g,E))goto done;have_factor=1;
 if(stage)*stage=5;for(size_t j=0;j<count;j++)if(!root_sp_apply(adapted,payload+j*C,&word,zeta8,0,E->field)||!root_theta_factor_point(answer+j,adapted,&factor)||!a_ec_on_curve(answer+j,factor.elliptic.A,E))goto done;
 fq_set(A,factor.elliptic.A,E->field);for(size_t j=0;j<count;j++)a_ec_set(points+j,answer+j,E);if(stage)*stage=0;ok=1;
done:if(have_factor)root_theta_factor_clear(&factor);if(have_labels)root_theta_label_clear(&label_ctx);for(size_t j=0;j<initialized;j++)a_ec_clear(answer+j,E);free(answer);b_theta4_vector_free(adapted,C,E->field);b_theta4_vector_free(labels,2*g*C,E->field);return ok;
}
