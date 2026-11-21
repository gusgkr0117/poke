#include "a-kani-two-kernel.h"
#include <limits.h>
void a_kani_two_kernel_init(a_kani_two_kernel *K,const root_extension *E){
 for(unsigned j=0;j<2;j++){a_ec_init(K->basis+j,E);for(unsigned i=0;i<2;i++){a_ec_init(K->generator[j]+i,E);a_ec_init(K->complement[j]+i,E);}}fmpz_init(K->order);
}
void a_kani_two_kernel_clear(a_kani_two_kernel *K,const root_extension *E){
 fmpz_clear(K->order);for(unsigned j=0;j<2;j++){a_ec_clear(K->basis+j,E);for(unsigned i=0;i<2;i++){a_ec_clear(K->generator[j]+i,E);a_ec_clear(K->complement[j]+i,E);}}
}
int a_kani_two_kernel_lift(a_kani_two_kernel *out,ulong ell,unsigned e,
 unsigned extra,const fmpz_t n1,const fmpz_t n2,const fmpz *eta,unsigned den,
 a_ec_rng rng,void *context,unsigned long budget,unsigned long *trials,const root_extension *E)
{
 if(!e||extra>=UINT_MAX-e||ell<2||(den!=1&&den!=2)||fmpz_sgn(n1)<=0||fmpz_sgn(n2)<=0||!fmpz_fdiv_ui(n1,ell)||!fmpz_fdiv_ui(n2,ell))return A_EC_INVALID;
 if(den==2&&(fmpz_fdiv_ui(eta,2)!=fmpz_fdiv_ui(eta+3,2)||fmpz_fdiv_ui(eta+1,2)!=fmpz_fdiv_ui(eta+2,2)))return A_EC_INVALID;
 fmpz_t scalar,base,norm,t,inv_n1,inv_product,two;
 fmpz_init(scalar);fmpz_init(base);fmpz_init(norm);fmpz_init(t);fmpz_init(inv_n1);fmpz_init(inv_product);fmpz_init_set_ui(two,2);
 int status=A_EC_INVALID;a_kani_two_kernel K;a_kani_two_kernel_init(&K,E);
 a_ec_point lifts[2],images[2];for(unsigned j=0;j<2;j++){a_ec_init(lifts+j,E);a_ec_init(images+j,E);}
 fq_t zero;fq_init(zero,E->field);fq_zero(zero,E->field);
 fmpz_mul(norm,eta+2,eta+2);fmpz_addmul(norm,eta+3,eta+3);fmpz_mul(norm,norm,fq_ctx_prime(E->field));fmpz_addmul(norm,eta,eta);fmpz_addmul(norm,eta+1,eta+1);
 fmpz_mul(scalar,n1,n2);fmpz_mul_ui(scalar,scalar,den*den);if(!fmpz_equal(norm,scalar))goto done;
 fmpz_set_ui(base,ell);fmpz_pow_ui(base,base,e);fmpz_add(scalar,n1,n2);if(!fmpz_divisible(scalar,base))goto done;
 fmpz_set_ui(K.order,ell);fmpz_pow_ui(K.order,K.order,e+extra);
 if(!fmpz_invmod(t,n2,K.order)||!fmpz_invmod(inv_n1,n1,K.order))goto done;
 fmpz_mul(t,t,n1);fmpz_neg(t,t);fmpz_mod(t,t,K.order);fmpz_mul(scalar,n1,n2);if(!fmpz_invmod(inv_product,scalar,K.order))goto done;
 int half=ell==2&&den==2;
 status=a_ec_torsion_basis(lifts,lifts+1,ell,e+extra+half,zero,rng,context,budget,trials,E);if(status!=A_EC_OK)goto done;
 for(unsigned j=0;j<2;j++){
  if(half)a_ec_mul(K.basis+j,lifts+j,two,zero,E);else a_ec_set(K.basis+j,lifts+j,E);
  status=a_ec_quaternion(images+j,K.basis+j,eta,den,K.order,half?lifts+j:NULL,E);if(status!=A_EC_OK)goto done;
  fmpz_mod(scalar,n1,K.order);a_ec_mul(K.generator[j],K.basis+j,scalar,zero,E);
  if(j)a_ec_mul(K.generator[j]+1,images+j,t,zero,E);else a_ec_set(K.generator[j]+1,images+j,E);
 }
 a_ec_infinity(K.complement[0],E);a_ec_mul(K.complement[0]+1,images+1,inv_product,zero,E);
 fmpz_neg(scalar,inv_n1);fmpz_mod(scalar,scalar,K.order);a_ec_mul(K.complement[1],K.basis,scalar,zero,E);a_ec_infinity(K.complement[1]+1,E);
 for(unsigned j=0;j<2;j++){a_ec_set(out->basis+j,K.basis+j,E);for(unsigned i=0;i<2;i++){a_ec_set(out->generator[j]+i,K.generator[j]+i,E);a_ec_set(out->complement[j]+i,K.complement[j]+i,E);}}fmpz_set(out->order,K.order);status=A_EC_OK;
done:
 fq_clear(zero,E->field);for(unsigned j=0;j<2;j++){a_ec_clear(images+j,E);a_ec_clear(lifts+j,E);}a_kani_two_kernel_clear(&K,E);
 fmpz_clear(two);fmpz_clear(inv_product);fmpz_clear(inv_n1);fmpz_clear(t);fmpz_clear(norm);fmpz_clear(base);fmpz_clear(scalar);return status;
}
