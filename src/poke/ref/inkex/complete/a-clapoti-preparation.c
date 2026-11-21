#include "a-clapoti-preparation.h"
#include "../../../../quaternion/ref/generic/internal.h"
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#define A_FOUR_SQUARES_LIBRARY
#include "a-four-squares.c"
static int integral(const quat_lattice_t *L,const quat_order_t *O){quat_alg_elem_t x;quat_alg_elem_init(&x);mpz_set(x.denom,L->denom);int ok=mpz_sgn(x.denom)>0;for(unsigned j=0;ok&&j<4;j++){for(unsigned i=0;i<4;i++)mpz_set(x.coord[i],L->basis[i][j]);ok=quat_lattice_contains(NULL,O,&x,&QUATALG_PINFTY);}quat_alg_elem_finalize(&x);return ok;}
static int input_valid(const quat_left_ideal_t *I){if(!I||I->parent_order!=&MAXORD_O0||mpz_sgn(I->norm)<=0||!integral(&I->lattice,I->parent_order))return 0;ibz_t n,t;ibz_init(&n);ibz_init(&t);quat_lattice_index(&n,&I->lattice,I->parent_order);mpz_mul(t,I->norm,I->norm);int ok=!mpz_cmp(n,t);mpz_mul_ui(t,QUATALG_PINFTY.p,6);mpz_gcd(n,t,I->norm);ok=ok&&!mpz_cmp_ui(n,1);ibz_finalize(&t);ibz_finalize(&n);return ok;}
static int norm_int(ibz_t *n,const quat_alg_elem_t *x){ibq_t q;ibq_init(&q);quat_alg_norm(&q,x,&QUATALG_PINFTY);int ok=ibq_to_ibz(n,&q);ibq_finalize(&q);return ok;}
static void multiply(quat_lattice_t *out,const quat_lattice_t *in,const quat_alg_elem_t *x){ibz_mat_4x4_t M;ibz_mat_4x4_init(&M);quat_alg_rightmul_mat(&M,x,&QUATALG_PINFTY);ibz_mat_4x4_mul(&out->basis,&M,&in->basis);mpz_mul(out->denom,in->denom,x->denom);quat_lattice_reduce_denom(out,out);quat_lattice_hnf(out);ibz_mat_4x4_finalize(&M);}
static void element(quat_alg_elem_t *x,const ibz_mat_4x4_t *B,const ibz_t *den,const ibz_vec_4_t *v){mpz_set(x->denom,*den);for(unsigned i=0;i<4;i++){mpz_set_ui(x->coord[i],0);for(unsigned j=0;j<4;j++)mpz_addmul(x->coord[i],(*B)[i][j],(*v)[j]);}quat_alg_normalize(x);}
static int unit_norm(ibz_t *d,const quat_alg_elem_t *x,const quat_left_ideal_t *I,const ibz_t *mod){if(!norm_int(d,x)||!mpz_divisible_p(*d,I->norm))return -1;mpz_divexact(*d,*d,I->norm);if(mpz_sgn(*d)<=0)return 0;ibz_t g;ibz_init(&g);mpz_gcd(g,*d,*mod);int ok=!mpz_cmp_ui(g,1);ibz_finalize(&g);return ok;}
void a_clapoti_certificate_init(a_clapoti_certificate *C){memset(C,0,sizeof(*C));quat_alg_elem_init(&C->beta1);quat_alg_elem_init(&C->beta2);quat_alg_elem_init(&C->theta);ibz_init(&C->d1);ibz_init(&C->d2);ibz_init(&C->u);ibz_init(&C->v);ibz_init(&C->N);ibz_init(&C->n1);ibz_init(&C->n2);ibz_mat_4x4_init(&C->X);ibz_mat_4x4_init(&C->Y);}
void a_clapoti_certificate_clear(a_clapoti_certificate *C){free(C->extensions);free(C->exponents);free(C->primes);ibz_mat_4x4_finalize(&C->Y);ibz_mat_4x4_finalize(&C->X);ibz_finalize(&C->n2);ibz_finalize(&C->n1);ibz_finalize(&C->N);ibz_finalize(&C->v);ibz_finalize(&C->u);ibz_finalize(&C->d2);ibz_finalize(&C->d1);quat_alg_elem_finalize(&C->theta);quat_alg_elem_finalize(&C->beta2);quat_alg_elem_finalize(&C->beta1);memset(C,0,sizeof(*C));}
static int first(quat_alg_elem_t *x,ibz_t *d,const quat_left_ideal_t *I,const ibz_mat_4x4_t *B){
 ibz_vec_4_t v,base;ibz_vec_4_init(&v);ibz_vec_4_init(&base);ibz_t six,mod;ibz_init(&six);ibz_init(&mod);mpz_set_ui(six,6);mpz_mul_ui(mod,QUATALG_PINFTY.p,6);int found=0;
 for(unsigned k=0;k<1296&&!found;k++){unsigned z=k;for(unsigned j=0;j<4;j++){mpz_set_ui(v[j],z%6);z/=6;}element(x,B,&I->lattice.denom,&v);found=unit_norm(d,x,I,&six);if(found<0)goto done;}
 if(!found)goto done;for(unsigned j=0;j<4;j++)mpz_set(base[j],v[j]);found=unit_norm(d,x,I,&mod);if(found<0)goto done;
 for(unsigned i=0;i<4&&!found;i++)for(unsigned t=1;t<=2&&!found;t++){for(unsigned j=0;j<4;j++)mpz_set(v[j],base[j]);mpz_add_ui(v[i],v[i],6*t);element(x,B,&I->lattice.denom,&v);found=unit_norm(d,x,I,&mod);if(found<0)goto done;}
 for(unsigned i=0;i<4&&!found;i++)for(unsigned j=i+1;j<4&&!found;j++){for(unsigned k=0;k<4;k++)mpz_set(v[k],base[k]);mpz_add_ui(v[i],v[i],6);mpz_add_ui(v[j],v[j],6);element(x,B,&I->lattice.denom,&v);found=unit_norm(d,x,I,&mod);if(found<0)goto done;}
done:ibz_finalize(&mod);ibz_finalize(&six);ibz_vec_4_finalize(&base);ibz_vec_4_finalize(&v);return found==1?A_CLAPOTI_OK:A_CLAPOTI_INVALID;
}
static int second(quat_alg_elem_t *x,ibz_t *d,const quat_left_ideal_t *I,const ibz_mat_4x4_t *B,const ibz_t *d1,a_clapoti_counts *counts){
 ibz_vec_4_t v;ibz_vec_4_init(&v);ibz_t mod,lo,hi;ibz_init(&mod);ibz_init(&lo);ibz_init(&hi);mpz_mul(mod,*d1,QUATALG_PINFTY.p);mpz_mul_ui(mod,mod,6);mpz_set_ui(lo,0);mpz_set_ui(hi,31);int status=A_CLAPOTI_INVALID;unsigned prefix=0;
 for(;;){
  if(prefix==32){mpz_sub_ui(hi,mod,1);counts->used_complete_residue_box=1;prefix++;}
  for(unsigned j=0;j<4;j++)if(!ibz_rand_interval(&v[j],&lo,&hi)){status=A_CLAPOTI_RNG_ERROR;goto done;}
  if(prefix<32)prefix++;if(counts->second_candidates<ULONG_MAX)counts->second_candidates++;
  element(x,B,&I->lattice.denom,&v);int found=unit_norm(d,x,I,&mod);if(found<0)goto done;if(found){status=A_CLAPOTI_OK;break;}
 }
done:ibz_finalize(&hi);ibz_finalize(&lo);ibz_finalize(&mod);ibz_vec_4_finalize(&v);return status;
}
static unsigned degree_prime(unsigned long ell){unsigned long z=(ell-mpz_fdiv_ui(QUATALG_PINFTY.p,ell))%ell;uint64_t a=1;for(unsigned long k=1;k<ell;k++){a=a*z%ell;if(a==1)return (unsigned)k;}return 0;}
static int plan(a_clapoti_certificate *C,a_clapoti_counts *counts){
 ibz_t D,threshold,g,t;ibz_init(&D);ibz_init(&threshold);ibz_init(&g);ibz_init(&t);mpz_mul(D,C->d1,C->d2);mpz_mul_ui(threshold,D,6);size_t cap=mpz_sizeinbase(threshold,2)+3;int status=A_CLAPOTI_INVALID;
 if(cap>SIZE_MAX/sizeof(*C->primes))goto done;C->primes=calloc(cap,sizeof(*C->primes));C->exponents=calloc(cap,sizeof(*C->exponents));C->extensions=calloc(cap,sizeof(*C->extensions));if(!C->primes||!C->exponents||!C->extensions){status=A_CLAPOTI_NO_MEMORY;goto done;}
 C->prime_count=2;C->primes[0]=2;C->exponents[0]=TORSION_PLUS_EVEN_POWER-2;C->extensions[0]=2;C->primes[1]=3;C->extensions[1]=1;mpz_set(t,TORSION_PLUS_3POWER);while(mpz_divisible_ui_p(t,3)){mpz_divexact_ui(t,t,3);C->exponents[1]++;}if(mpz_cmp_ui(t,1)||!C->exponents[1])goto done;
 mpz_set_ui(C->N,1);mpz_mul_2exp(C->N,C->N,C->exponents[0]);mpz_mul(C->N,C->N,TORSION_PLUS_3POWER);
 for(unsigned long ell=5;mpz_cmp(C->N,threshold)<=0;ell+=2){
  if(ell>UINT_MAX-2||C->prime_count==cap)goto done;
  if(!small_prime(ell)||!mpz_fdiv_ui(D,ell)||!mpz_fdiv_ui(QUATALG_PINFTY.p,ell))continue;
  unsigned k=degree_prime(ell);if(!k)goto done;size_t i=C->prime_count++;C->primes[i]=ell;C->exponents[i]=1;C->extensions[i]=k;mpz_mul_ui(C->N,C->N,ell);
 }
 if(mpz_cmp_ui(C->d2,1)==0)mpz_set_ui(C->u,1);else{if(!mpz_invert(C->u,C->d1,C->d2))goto done;mpz_mul(C->u,C->u,C->N);mpz_mod(C->u,C->u,C->d2);}
 while(mpz_even_p(C->u)||mpz_divisible_ui_p(C->u,3))mpz_add(C->u,C->u,C->d2);
 mpz_mul(C->v,C->u,C->d1);mpz_sub(C->v,C->N,C->v);if(!mpz_divisible_p(C->v,C->d2))goto done;mpz_divexact(C->v,C->v,C->d2);if(mpz_sgn(C->v)<=0)goto done;
 mpz_gcd(g,C->u,C->v);mpz_divexact(C->u,C->u,g);mpz_divexact(C->v,C->v,g);mpz_divexact(C->N,C->N,g);
 size_t kept=0;for(size_t j=0;j<C->prime_count;j++){unsigned e=C->exponents[j];while(e&&mpz_divisible_ui_p(g,C->primes[j])){mpz_divexact_ui(g,g,C->primes[j]);e--;}if(e){C->primes[kept]=C->primes[j];C->exponents[kept]=e;C->extensions[kept]=C->extensions[j];kept++;}}C->prime_count=kept;if(mpz_cmp_ui(g,1))goto done;
 mpz_mul(C->n1,C->u,C->d1);mpz_mul(C->n2,C->v,C->d2);
 a_hurwitz x,y;h_init(&x);h_init(&y);int fs=a_four_squares(&x,&C->u,0,&counts->four_square_trials);if(fs==1)fs=a_four_squares(&y,&C->v,0,&counts->four_square_trials);
 if(fs==1){a_hamilton_matrix(C->X,&x,0);a_hamilton_matrix(C->Y,&y,1);status=A_CLAPOTI_OK;}else status=fs<0?A_CLAPOTI_RNG_ERROR:A_CLAPOTI_INVALID;h_clear(&y);h_clear(&x);
done:ibz_finalize(&t);ibz_finalize(&g);ibz_finalize(&threshold);ibz_finalize(&D);return status;
}
static int auxiliary(const quat_left_ideal_t *I,const quat_alg_elem_t *beta,const ibz_t *d){
 int ok=0;quat_alg_elem_t bar;quat_alg_elem_init(&bar);quat_lattice_t J,back,target;quat_lattice_init(&J);quat_lattice_init(&back);quat_lattice_init(&target);ibz_t t,n;ibz_init(&t);ibz_init(&n);
 if(!quat_lattice_contains(NULL,&I->lattice,beta,&QUATALG_PINFTY)||!norm_int(&t,beta))goto done;mpz_mul(n,*d,I->norm);if(mpz_cmp(n,t))goto done;
 quat_alg_conj(&bar,beta);mpz_mul(bar.denom,bar.denom,I->norm);multiply(&J,&I->lattice,&bar);if(!integral(&J,I->parent_order))goto done;quat_lattice_index(&t,&J,I->parent_order);mpz_mul(n,*d,*d);if(mpz_cmp(n,t))goto done;
 multiply(&back,&J,beta);ibz_mat_4x4_copy(&target.basis,&I->lattice.basis);mpz_set(target.denom,I->lattice.denom);for(unsigned i=0;i<4;i++)for(unsigned j=0;j<4;j++)mpz_mul(target.basis[i][j],target.basis[i][j],*d);quat_lattice_reduce_denom(&target,&target);quat_lattice_hnf(&target);ok=quat_lattice_equal(&back,&target);
done:ibz_finalize(&n);ibz_finalize(&t);quat_lattice_finalize(&target);quat_lattice_finalize(&back);quat_lattice_finalize(&J);quat_alg_elem_finalize(&bar);return ok;
}
int a_clapoti_certificate_verify(const a_clapoti_certificate *C,const quat_left_ideal_t *I){
 if(!C||!input_valid(I)||mpz_sgn(C->d1)<=0||mpz_sgn(C->d2)<=0||mpz_sgn(C->u)<=0||mpz_sgn(C->v)<=0||mpz_sgn(C->N)<=0||!C->prime_count||!C->primes||!C->exponents||!C->extensions||C->primes[0]!=2)return 0;
 const quat_alg_elem_t *q[3]={&C->beta1,&C->beta2,&C->theta};for(unsigned i=0;i<3;i++)if(mpz_sgn(q[i]->denom)<=0)return 0;
 int ok=0;ibz_t t,z,w,mod,product;ibz_init(&t);ibz_init(&z);ibz_init(&w);ibz_init(&mod);ibz_init(&product);quat_alg_elem_t bar,theta;quat_alg_elem_init(&bar);quat_alg_elem_init(&theta);
#define V(c) do{if(!(c))goto done;}while(0)
 V(auxiliary(I,&C->beta1,&C->d1)&&auxiliary(I,&C->beta2,&C->d2));mpz_gcd(t,C->d1,C->d2);V(!mpz_cmp_ui(t,1));mpz_mul(t,C->d1,C->d2);mpz_mul_ui(z,QUATALG_PINFTY.p,6);mpz_gcd(z,z,t);V(!mpz_cmp_ui(z,1));
 quat_alg_conj(&bar,&C->beta1);quat_alg_mul(&theta,&C->beta2,&bar,&QUATALG_PINFTY);mpz_mul(theta.denom,theta.denom,I->norm);quat_alg_normalize(&theta);for(unsigned i=0;i<4;i++){mpz_mul(z,theta.coord[i],C->theta.denom);mpz_mul(w,C->theta.coord[i],theta.denom);V(!mpz_cmp(z,w));}V(quat_lattice_contains(NULL,&MAXORD_O0,&C->theta,&QUATALG_PINFTY)&&norm_int(&z,&C->theta)&&!mpz_cmp(z,t));
 mpz_mul(t,C->u,C->d1);V(!mpz_cmp(t,C->n1));mpz_mul(t,C->v,C->d2);V(!mpz_cmp(t,C->n2));mpz_add(t,C->n1,C->n2);V(!mpz_cmp(t,C->N));mpz_gcd(t,C->n1,C->n2);V(!mpz_cmp_ui(t,1));mpz_mul(t,C->n1,C->n2);V(mpz_odd_p(t)&&!mpz_divisible_ui_p(t,3));
 for(unsigned i=0;i<4;i++)for(unsigned j=0;j<4;j++){
  mpz_set_ui(t,0);mpz_set_ui(z,0);for(unsigned k=0;k<4;k++){mpz_addmul(t,C->X[k][i],C->X[k][j]);mpz_addmul(z,C->Y[k][i],C->Y[k][j]);}V(i==j?(!mpz_cmp(t,C->u)&&!mpz_cmp(z,C->v)):(!mpz_sgn(t)&&!mpz_sgn(z)));
  mpz_set_ui(t,0);mpz_set_ui(z,0);mpz_set_ui(w,0);mpz_set_ui(mod,0);for(unsigned k=0;k<4;k++){mpz_addmul(t,C->X[i][k],C->Y[k][j]);mpz_addmul(z,C->Y[i][k],C->X[k][j]);mpz_addmul(w,C->X[i][k],C->Y[j][k]);mpz_addmul(mod,C->Y[k][i],C->X[k][j]);}V(!mpz_cmp(t,z)&&!mpz_cmp(w,mod));
 }
 mpz_set_ui(product,1);for(size_t j=0;j<C->prime_count;j++){
  unsigned long ell=C->primes[j];unsigned e=C->exponents[j];V(small_prime(ell)&&e&&C->extensions[j]&&e<=UINT_MAX-3);for(size_t k=0;k<j;k++)V(C->primes[k]!=ell);
  mpz_ui_pow_ui(mod,ell,e);mpz_mul(product,product,mod);mpz_ui_pow_ui(mod,ell,e+(j?0:3));mpz_neg(t,QUATALG_PINFTY.p);mpz_powm_ui(z,t,C->extensions[j],mod);V(!mpz_cmp_ui(z,1));
 }V(!mpz_cmp(product,C->N));ok=1;
done:
#undef V
 quat_alg_elem_finalize(&theta);quat_alg_elem_finalize(&bar);ibz_finalize(&product);ibz_finalize(&mod);ibz_finalize(&w);ibz_finalize(&z);ibz_finalize(&t);return ok;
}
int a_clapoti_complete_profile_verify(const a_clapoti_certificate *C,const quat_left_ideal_t *I){
 if(!a_clapoti_certificate_verify(C,I)||C->exponents[0]!=TORSION_PLUS_EVEN_POWER-2||!mpz_divisible_p(C->N,TORSION_PLUS_3POWER))return 0;
 for(size_t j=2;j<C->prime_count;j++)if(C->exponents[j]!=1||C->extensions[j]>=C->primes[j])return 0;
 return 1;
}
int a_clapoti_prepare(a_clapoti_certificate *out,const quat_left_ideal_t *I,a_clapoti_counts *counts){
 if(counts)memset(counts,0,sizeof(*counts));if(!out||!input_valid(I))return A_CLAPOTI_INVALID;
 a_clapoti_counts stats={0};a_clapoti_certificate C;a_clapoti_certificate_init(&C);ibz_mat_4x4_t B;ibz_mat_4x4_init(&B);quat_alg_elem_t bar;quat_alg_elem_init(&bar);int status=A_CLAPOTI_INVALID;
 if(quat_lattice_lll(&B,&I->lattice,&QUATALG_PINFTY.p))goto done;status=first(&C.beta1,&C.d1,I,&B);if(status!=A_CLAPOTI_OK)goto done;status=second(&C.beta2,&C.d2,I,&B,&C.d1,&stats);if(status!=A_CLAPOTI_OK)goto done;
 quat_alg_conj(&bar,&C.beta1);quat_alg_mul(&C.theta,&C.beta2,&bar,&QUATALG_PINFTY);mpz_mul(C.theta.denom,C.theta.denom,I->norm);quat_alg_normalize(&C.theta);
 status=plan(&C,&stats);if(status!=A_CLAPOTI_OK)goto done;if(!a_clapoti_complete_profile_verify(&C,I)){status=A_CLAPOTI_INVALID;goto done;}C.verified=1;
 a_clapoti_certificate old=*out;*out=C;C=old;status=A_CLAPOTI_OK;
done:if(counts)*counts=stats;quat_alg_elem_finalize(&bar);ibz_mat_4x4_finalize(&B);a_clapoti_certificate_clear(&C);return status;
}
