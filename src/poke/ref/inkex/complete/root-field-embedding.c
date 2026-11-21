#include "root-field-embedding.h"
#include <flint/fq_embed.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fq_poly.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
/* The source modulus splits into distinct linear factors in the target.
 * Fresh uniform coefficients modulo the current factor give independent
 * values at its roots. Zero roots and the two quadratic-character classes
 * therefore give a proper factor with constant probability for degree>=2. */
static int embedding_generators(fq_t gsub,fq_t gsup,fmpz_mod_poly_t minpoly,
 const root_extension *sub,const root_extension *sup){
 const fq_ctx_struct *f=sup->field;int ok=0;
 fq_poly_t current,h,part,power,one;fq_poly_init(current,f);fq_poly_init(h,f);fq_poly_init(part,f);fq_poly_init(power,f);fq_poly_init(one,f);
 fq_t coefficient,value;fq_init(coefficient,f);fq_init(value,f);fmpz_t exponent;fmpz_init(exponent);fq_ctx_order(exponent,f);fmpz_sub_ui(exponent,exponent,1);fmpz_tdiv_q_2exp(exponent,exponent,1);
 fq_poly_set_fmpz_mod_poly(current,fq_ctx_modulus(sub->field),f);fq_poly_one(one,f);
 while(fq_poly_degree(current,f)>1){
  slong degree=fq_poly_degree(current,f);fq_poly_zero(h,f);
  for(slong j=0;j<degree;j++){
   if(!root_random_field(coefficient,f,sup->rng,sup->rng_context)){ok=-2;goto done;}
   fq_poly_set_coeff(h,j,coefficient,f);
  }
  fq_poly_gcd(part,current,h,f);slong d=fq_poly_degree(part,f);
  if(d>0&&d<degree){fq_poly_make_monic(current,part,f);continue;}
  if(d==degree)continue;
  fq_poly_powmod_fmpz_binexp(power,h,exponent,current,f);fq_poly_sub(power,power,one,f);fq_poly_gcd(part,current,power,f);d=fq_poly_degree(part,f);
  if(d>0&&d<degree)fq_poly_make_monic(current,part,f);
 }
 if(fq_poly_degree(current,f)!=1)goto done;
 fq_poly_get_coeff(gsup,current,0,f);fq_neg(gsup,gsup,f);fq_poly_get_coeff(coefficient,current,1,f);fq_inv(coefficient,coefficient,f);fq_mul(gsup,gsup,coefficient,f);
 fq_poly_set_fmpz_mod_poly(current,fq_ctx_modulus(sub->field),f);fq_poly_evaluate_fq(value,current,gsup,f);if(!fq_is_zero(value,f))goto done;
 fq_gen(gsub,sub->field);fmpz_mod_poly_set(minpoly,fq_ctx_modulus(sub->field),sub->field->ctxp);ok=1;
done:fmpz_clear(exponent);fq_clear(value,f);fq_clear(coefficient,f);fq_poly_clear(one,f);fq_poly_clear(power,f);fq_poly_clear(part,f);fq_poly_clear(h,f);fq_poly_clear(current,f);return ok;
}
static void apply(fq_t out,const fq_t in,const fmpz_mod_mat_t matrix,const root_extension *from,const root_extension *to){
 fmpz_poly_t a,b;fmpz_poly_init(a);fmpz_poly_init(b);fq_get_fmpz_poly(a,in,from->field);fmpz_t sum,t;fmpz_init(sum);fmpz_init(t);
 slong rows=fmpz_mat_nrows(matrix),cols=fmpz_mat_ncols(matrix);
 for(slong i=0;i<rows;i++){fmpz_zero(sum);for(slong j=0;j<cols;j++){fmpz_poly_get_coeff_fmpz(t,a,j);fmpz_addmul(sum,fmpz_mat_entry(matrix,i,j),t);}fmpz_mod(sum,sum,fq_ctx_prime(to->field));fmpz_poly_set_coeff_fmpz(b,i,sum);}
 fq_set_fmpz_poly(out,b,to->field);fmpz_clear(t);fmpz_clear(sum);fmpz_poly_clear(b);fmpz_poly_clear(a);
}
void root_field_embedding_clear(root_field_embedding *L){if(L->ready){fmpz_mod_mat_clear(L->project,L->sub->field->ctxp);fmpz_mod_mat_clear(L->embed,L->sub->field->ctxp);}memset(L,0,sizeof(*L));}
void root_field_up(fq_t out,const fq_t in,const root_field_embedding *L){apply(out,in,L->embed,L->sub,L->sup);}
int root_field_down(fq_t out,const fq_t in,const root_field_embedding *L){
 if(!L->ready)return 0;fq_t a,check;fq_init(a,L->sub->field);fq_init(check,L->sup->field);apply(a,in,L->project,L->sup,L->sub);root_field_up(check,a,L);int ok=fq_equal(check,in,L->sup->field);if(ok)fq_set(out,a,L->sub->field);fq_clear(check,L->sup->field);fq_clear(a,L->sub->field);return ok;
}
int root_field_embedding_init(root_field_embedding *L,const root_extension *sub,const root_extension *sup){
 if(!L)return 0;memset(L,0,sizeof(*L));L->last_error=-1;if(!sub||!sup||sub->relative_degree<1||sup->relative_degree%sub->relative_degree||!fmpz_equal(fq_ctx_prime(sub->field),fq_ctx_prime(sup->field)))return 0;
 L->sub=sub;L->sup=sup;fmpz_mod_mat_init(L->embed,2*sup->relative_degree,2*sub->relative_degree,sub->field->ctxp);fmpz_mod_mat_init(L->project,2*sub->relative_degree,2*sup->relative_degree,sub->field->ctxp);L->ready=1;
 fq_t gsub,gsup,check,minus;fq_init(gsub,sub->field);fq_init(gsup,sup->field);fq_init(check,sup->field);fq_init(minus,sup->field);fmpz_mod_poly_t minpoly;fmpz_mod_poly_init(minpoly,sub->field->ctxp);
 int result=embedding_generators(gsub,gsup,minpoly,sub,sup),error=result==-2?-2:-1,ok=result==1;if(!ok)goto cleanup;
 fq_embed_matrices(L->embed,L->project,gsub,sub->field,gsup,sup->field,minpoly);root_field_up(check,sub->i,L);
 if(!fq_equal(check,sup->i,sup->field)){
  fq_neg(minus,sup->i,sup->field);if(!fq_equal(check,minus,sup->field))ok=0;
  else{fq_frobenius(gsup,gsup,1,sup->field);fq_embed_matrices(L->embed,L->project,gsub,sub->field,gsup,sup->field,minpoly);root_field_up(check,sub->i,L);ok=fq_equal(check,sup->i,sup->field);}
 }
 /* Check the entire claimed left inverse over Fp, not only one example. */
 fmpz_mat_t product;fmpz_mat_init(product,2*sub->relative_degree,2*sub->relative_degree);fmpz_mat_mul(product,L->project,L->embed);
 for(slong i=0;i<2*sub->relative_degree;i++)for(slong j=0;j<2*sub->relative_degree;j++){
  fmpz_mod(fmpz_mat_entry(product,i,j),fmpz_mat_entry(product,i,j),fq_ctx_prime(sub->field));if(!fmpz_equal_ui(fmpz_mat_entry(product,i,j),i==j))ok=0;
 }
 fmpz_mat_clear(product);
cleanup:fmpz_mod_poly_clear(minpoly,sub->field->ctxp);fq_clear(minus,sup->field);fq_clear(check,sup->field);fq_clear(gsup,sup->field);fq_clear(gsub,sub->field);if(!ok)root_field_embedding_clear(L);L->last_error=ok?0:error;return ok;
}
int root_field_projective_down(fq_struct *out,const fq_struct *in,size_t count,const root_field_embedding *L){
 if(!L->ready||!count||count>SIZE_MAX/sizeof(fq_struct))return 0;size_t pivot=0;while(pivot<count&&fq_is_zero(in+pivot,L->sup->field))pivot++;if(pivot==count)return 0;
 fq_struct *v=malloc(count*sizeof(fq_struct));if(!v)return 0;for(size_t i=0;i<count;i++)fq_init(v+i,L->sub->field);fq_t inverse,x;fq_init(inverse,L->sup->field);fq_init(x,L->sup->field);fq_inv(inverse,in+pivot,L->sup->field);int ok=1;
 for(size_t i=0;i<count;i++){fq_mul(x,in+i,inverse,L->sup->field);if(!root_field_down(v+i,x,L)){ok=0;break;}}
 if(ok)for(size_t i=0;i<count;i++)fq_set(out+i,v+i,L->sub->field);fq_clear(x,L->sup->field);fq_clear(inverse,L->sup->field);for(size_t i=0;i<count;i++)fq_clear(v+i,L->sub->field);free(v);return ok;
}
