#include "b-product-fallback.h"
#include "root-theta-product.h"
#include "root-poke-field-bridge.h"
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <string.h>

int b_product_context_init(b_product_context *C,const root_extension *base,
 a_ec_rng rng,void *rng_context,unsigned long budget){
 if(!C||C->ready||!base||base->relative_degree!=1||!rng)return INKE_SAMPLE_INVALID;
 C->base=base;C->rng=rng;C->rng_context=rng_context;C->point_budget=budget;C->ready=1;return INKE_SAMPLE_OK;
}
void b_product_context_clear(b_product_context *C){
 if(C)memset(C,0,sizeof *C);
}
static int product_point(fq_struct *out,const a_ec_point *P,const a_ec_point *Q,
 const root_ec_theta4 T[2]){
 const fq_ctx_struct *F=T[0].E->field;fq_struct *small=b_theta4_vector_new(8,F);if(!small)return 0;
 int ok=root_ec_theta4_map(small,P,T)&&root_ec_theta4_map(small+4,Q,T+1);
 if(ok)for(unsigned j=0;j<4;j++)for(unsigned i=0;i<4;i++)fq_mul(out+i+4*j,small+i,small+4+j,F);
 b_theta4_vector_free(small,8,F);return ok;
}
int b_product_fallback(ec_curve_t *target,ec_basis_t *left_images,
 ec_basis_t *right_images,const ec_curve_t *left_curve,
 const ec_curve_t *right_curve,const ec_basis_t *left_two,
 const ec_basis_t *right_two,unsigned int exponent,const ibz_t *q,
 const ec_basis_t *left_points,const ec_basis_t *right_points,
 void *context,int *stage){
 b_product_context *C=context;if(stage)*stage=1;
 if(!C||!C->ready||!target||!left_curve||!right_curve||!left_two||!right_two||
    exponent<2||exponent>TORSION_PLUS_EVEN_POWER-2||!q||mpz_sgn(*q)<=0||!mpz_odd_p(*q)||
    (!!left_images!=!!left_points)||(!!right_images!=!!right_points))return INKE_SAMPLE_INVALID;
 C->calls++;
 const root_extension *E=C->base;
 int status=INKE_SAMPLE_INVALID;unsigned phase=2;
 size_t count=1+3*(!!left_points)+3*(!!right_points);
 root_ec_theta4 T[2]={{0}};root_theta_product split={0};b_theta_balanced chain={0};
 fq_t A[2],zeta8,dummy;for(unsigned i=0;i<2;i++)fq_init(A[i],E->field);
 fq_init(zeta8,E->field);fq_init(dummy,E->field);
 a_ec_point source[3],right[3],payload[7],zero,pair[2],images[2][7],multiple;
 for(unsigned j=0;j<3;j++){a_ec_init(source+j,E);a_ec_init(right+j,E);}
 for(unsigned j=0;j<7;j++){a_ec_init(payload+j,E);for(unsigned i=0;i<2;i++)a_ec_init(images[i]+j,E);}
 a_ec_init(&zero,E);a_ec_infinity(&zero,E);a_ec_init(&multiple,E);
 for(unsigned i=0;i<2;i++)a_ec_init(pair+i,E);
 fq_struct *null=b_theta4_vector_new(16,E->field),*points=b_theta4_vector_new(count*16,E->field),*generators=b_theta4_vector_new(32,E->field);
 const fq_struct *point_ptrs[7],*generator_ptrs[2];
 fmpz_t scalar,order;fmpz_init(scalar);fmpz_init(order);
 ec_curve_t result_curve={0};ec_basis_t result_left,result_right;
 if(!null||!points||!generators)goto done;
 if(!b_poke_basis_lift(A[0],source,left_curve,left_two,E)||
    !b_poke_basis_lift(A[1],right,right_curve,right_two,E))goto done;
 /* P, Q and P-Q all have exact order 4N precisely when P,Q form
  * a basis. Check this before invoking the complete theta formulas. */
 fmpz_one(order);fmpz_mul_2exp(order,order,exponent);
 fmpz_mul_ui(scalar,order,2);
 for(unsigned side=0;side<2;side++){
  const a_ec_point *basis=side?right:source;
  for(unsigned i=0;i<3;i++){
   a_ec_mul(&multiple,basis+i,scalar,A[side],E);
   if(multiple.infinity)goto done;
   a_ec_add(&multiple,&multiple,&multiple,A[side],E);
   if(!multiple.infinity)goto done;
  }
 }
 for(unsigned i=0;i<2;i++)if(!root_ec_theta4_init(T+i,A[i],E))goto done;
 if(!root_extension_sqrt(zeta8,E->i,E))goto done;
 /* N=2^exponent. On E[4N] the two left coefficients -q and N-q
  * give exponent -q(N-q), cancelling the pairing exponent of the
  * two right images. Quadrupling gives the original graph kernel.
  * The first lifted generator maps to order four on the desired factor. */
 fmpz_set_mpz(scalar,*q);if(fmpz_cmp(scalar,order)>=0)goto done;
 fmpz_neg(scalar,scalar);
 if(!product_point(null,&zero,&zero,T))goto done;
 for(unsigned i=0;i<2;i++){
  if(i)fmpz_add(scalar,scalar,order);
  a_ec_mul(&multiple,source+i,scalar,A[0],E);
  generator_ptrs[i]=generators+16*i;
  if(!product_point(generators+16*i,&multiple,right+i,T))goto done;
 }
 point_ptrs[0]=generators;
 size_t used=1;if(left_points){if(!b_poke_basis_lift(dummy,payload+used,left_curve,left_points,E))goto done;used+=3;}
 if(right_points){if(!b_poke_basis_lift(dummy,payload+used,right_curve,right_points,E))goto done;}
 for(size_t i=1;i<count;i++){
  point_ptrs[i]=points+16*i;int side=right_points&&i>=count-3;
  if(!product_point(points+16*i,side?&zero:payload+i,side?payload+i:&zero,T))goto done;
 }
 phase=4;a_theta_payload_input group={E,point_ptrs,count};
 if(!b_theta_balanced_init(&chain,2,exponent,E,E,null,generator_ptrs,&group,1,E->i,zeta8))goto done;
 while(chain.state.remaining){if(!b_theta_balanced_step(&chain))goto done;if(C->progress)C->progress(C->progress_context,0,2,chain.state.completed,exponent);}
 phase=5;if(!root_theta_product_init(&split,chain.state.theta.null,zeta8,E))goto done;
 for(size_t i=0;i<count;i++){
  if(!root_theta_product_point(pair,chain.state.groups[0].points+16*i,&split))goto done;
  for(unsigned j=0;j<2;j++)a_ec_set(images[j]+i,pair+j,E);
 }
 if(images[0][0].infinity==images[1][0].infinity)goto done;
 unsigned selected=images[0][0].infinity?1:0;
 const fq_struct *target_A=split.factors[selected].elliptic.A;
 a_ec_add(&multiple,&images[selected][0],&images[selected][0],target_A,E);
 if(multiple.infinity)goto done;
 a_ec_add(&multiple,&multiple,&multiple,target_A,E);
 if(!multiple.infinity||!root_to_poke_fp2(&result_curve.A,target_A,E))goto done;
 fp2_set_one(&result_curve.C);
 phase=6;used=1;ec_curve_t check_curve;
 if(left_points){if(!b_poke_basis_export(&check_curve,&result_left,split.factors[selected].elliptic.A,images[selected]+used,E))goto done;used+=3;}
 if(right_points&&!b_poke_basis_export(&check_curve,&result_right,split.factors[selected].elliptic.A,images[selected]+used,E))goto done;
 *target=result_curve;if(left_images)*left_images=result_left;if(right_images)*right_images=result_right;status=INKE_SAMPLE_OK;phase=0;
done:
 if(status==INKE_SAMPLE_INVALID&&phase==4&&chain.state.last_error==-2)status=INKE_SAMPLE_RNG_ERROR;
 if(stage)*stage=(int)phase;root_theta_product_clear(&split);b_theta_balanced_clear(&chain);
 fmpz_clear(order);fmpz_clear(scalar);b_theta4_vector_free(generators,32,E->field);b_theta4_vector_free(points,count*16,E->field);b_theta4_vector_free(null,16,E->field);
 a_ec_clear(&multiple,E);
 for(unsigned i=0;i<2;i++)a_ec_clear(pair+i,E);a_ec_clear(&zero,E);for(unsigned j=0;j<7;j++){for(unsigned i=0;i<2;i++)a_ec_clear(images[i]+j,E);a_ec_clear(payload+j,E);}for(unsigned j=0;j<3;j++){a_ec_clear(right+j,E);a_ec_clear(source+j,E);}
 for(unsigned i=0;i<2;i++)root_ec_theta4_clear(T+i);fq_clear(dummy,E->field);fq_clear(zeta8,E->field);for(unsigned i=0;i<2;i++)fq_clear(A[i],E->field);return status;
}
