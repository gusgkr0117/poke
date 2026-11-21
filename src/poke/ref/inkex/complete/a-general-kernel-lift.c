#include "a-general-kernel-lift.h"
#include <limits.h>

/* Bekker--Zarhin, Section2, formula(6): choose ri^2=x-ei and
 * r0*r1*r2=-y, then Q=(x+s2,-y-s1*s2), [2]Q=P. */
int a_ec_half_montgomery(a_ec_point *out,const a_ec_point *P,const fq_t A,const root_extension *E){
 if(!out||!P||!E||!a_ec_curve_valid(A,E)||!a_ec_on_curve(P,A,E))return 0;
 fq_t roots[3],r[3],s1,s2,t,product;for(unsigned i=0;i<3;i++){fq_init(roots[i],E->field);fq_init(r[i],E->field);}
 fq_init(s1,E->field);fq_init(s2,E->field);fq_init(t,E->field);fq_init(product,E->field);
 a_ec_point Q,twice;a_ec_init(&Q,E);a_ec_init(&twice,E);int ok=0;
 if(P->infinity){a_ec_infinity(&Q,E);ok=1;goto success;}
 fq_zero(roots[0],E->field);fq_sqr(t,A,E->field);fq_set_ui(product,4,E->field);fq_sub(t,t,product,E->field);
 if(!root_extension_sqrt(t,t,E))goto done;
 fq_neg(roots[1],A,E->field);fq_add(roots[1],roots[1],t,E->field);fq_mul_fmpz(roots[1],roots[1],E->half,E->field);
 fq_neg(roots[2],A,E->field);fq_sub(roots[2],roots[2],t,E->field);fq_mul_fmpz(roots[2],roots[2],E->half,E->field);
 for(unsigned i=0;i<3;i++){fq_sub(t,P->x,roots[i],E->field);if(!root_extension_sqrt(r[i],t,E))goto done;}
 fq_mul(product,r[0],r[1],E->field);fq_mul(product,product,r[2],E->field);fq_neg(t,P->y,E->field);
 if(!fq_equal(product,t,E->field)){fq_neg(r[0],r[0],E->field);fq_neg(product,product,E->field);}if(!fq_equal(product,t,E->field))goto done;
 fq_add(s1,r[0],r[1],E->field);fq_add(s1,s1,r[2],E->field);fq_mul(s2,r[0],r[1],E->field);fq_mul(t,r[0],r[2],E->field);fq_add(s2,s2,t,E->field);fq_mul(t,r[1],r[2],E->field);fq_add(s2,s2,t,E->field);
 fq_add(Q.x,P->x,s2,E->field);fq_mul(Q.y,s1,s2,E->field);fq_add(Q.y,Q.y,P->y,E->field);fq_neg(Q.y,Q.y,E->field);Q.infinity=0;
 if(!a_ec_on_curve(&Q,A,E))goto done;a_ec_add(&twice,&Q,&Q,A,E);if(!a_ec_equal(&twice,P,E))goto done;ok=1;
success:if(ok)a_ec_set(out,&Q,E);
done:a_ec_clear(&twice,E);a_ec_clear(&Q,E);fq_clear(product,E->field);fq_clear(t,E->field);fq_clear(s2,E->field);fq_clear(s1,E->field);for(unsigned i=0;i<3;i++){fq_clear(r[i],E->field);fq_clear(roots[i],E->field);}return ok;
}

int a_general_kernel_lift(a_ec_point out[2][2],const a_ec_point K[2][2],unsigned exponent,
 const fq_t A1,const fq_t A2,a_ec_rng rng,void *context,unsigned long budget,const root_extension *E){
 if(!out||!K||!E||!rng||!exponent||exponent>UINT_MAX-2||!a_ec_curve_valid(A1,E)||!a_ec_curve_valid(A2,E))return A_EC_INVALID;
 const fq_struct *A[2]={A1,A2};int status=A_EC_INVALID;a_ec_point V[2][2],T;for(unsigned i=0;i<2;i++)for(unsigned j=0;j<2;j++)a_ec_init(&V[i][j],E);a_ec_init(&T,E);
 fmpz_t n,large,m,t;fmpz_init(n);fmpz_init(large);fmpz_init(m);fmpz_init(t);fmpz_one(n);fmpz_mul_2exp(n,n,exponent);fmpz_mul_ui(large,n,4);fmpz_fdiv_q_2exp(m,n,1);
 fq_t pair[2],rho,xi,test;fq_init(pair[0],E->field);fq_init(pair[1],E->field);fq_init(rho,E->field);fq_init(xi,E->field);fq_init(test,E->field);
 for(unsigned i=0;i<2;i++)for(unsigned j=0;j<2;j++){
  if(!a_ec_on_curve(&K[i][j],A[j],E))goto done;a_ec_mul(&T,&K[i][j],n,A[j],E);if(!T.infinity)goto done;
 }
 for(unsigned j=0;j<2;j++){status=a_ec_weil(pair[j],&K[0][j],&K[1][j],n,A[j],rng,context,budget,E);if(status!=A_EC_OK)goto done;}
 status=A_EC_INVALID;fq_pow(test,pair[0],m,E->field);if(fq_is_one(test,E->field))goto done;fq_mul(test,pair[0],pair[1],E->field);if(!fq_is_one(test,E->field))goto done;
 for(unsigned i=0;i<2;i++)for(unsigned j=0;j<2;j++){
  if(!a_ec_half_montgomery(&V[i][j],&K[i][j],A[j],E)||!a_ec_half_montgomery(&V[i][j],&V[i][j],A[j],E))goto done;
  a_ec_add(&T,&V[i][j],&V[i][j],A[j],E);a_ec_add(&T,&T,&T,A[j],E);if(!a_ec_equal(&T,&K[i][j],E))goto done;
 }
 for(unsigned j=0;j<2;j++){status=a_ec_weil(pair[j],&V[0][j],&V[1][j],large,A[j],rng,context,budget,E);if(status!=A_EC_OK)goto done;}
 status=A_EC_INVALID;fq_mul(rho,pair[0],pair[1],E->field);fq_pow_ui(test,rho,4,E->field);if(!fq_is_one(test,E->field))goto done;
 fq_pow(xi,pair[0],n,E->field);fq_sqr(test,xi,E->field);fq_neg(test,test,E->field);if(!fq_is_one(test,E->field))goto done;
 unsigned correction=0;fq_set(test,rho,E->field);while(correction<4&&!fq_is_one(test,E->field)){fq_mul(test,test,xi,E->field);correction++;}if(correction==4)goto done;
 fmpz_mul_ui(t,n,correction);fmpz_add_ui(t,t,1);a_ec_mul(&V[1][0],&V[1][0],t,A1,E);
 for(unsigned i=0;i<2;i++)for(unsigned j=0;j<2;j++){a_ec_add(&T,&V[i][j],&V[i][j],A[j],E);a_ec_add(&T,&T,&T,A[j],E);if(!a_ec_equal(&T,&K[i][j],E))goto done;}
 for(unsigned i=0;i<2;i++)for(unsigned j=0;j<2;j++)a_ec_set(&out[i][j],&V[i][j],E);status=A_EC_OK;
done:fq_clear(test,E->field);fq_clear(xi,E->field);fq_clear(rho,E->field);fq_clear(pair[1],E->field);fq_clear(pair[0],E->field);fmpz_clear(t);fmpz_clear(m);fmpz_clear(large);fmpz_clear(n);a_ec_clear(&T,E);for(unsigned i=0;i<2;i++)for(unsigned j=0;j<2;j++)a_ec_clear(&V[i][j],E);return status;
}
