#include "a-kani-geometry.h"
#include "a-extension-pairing.h"
#include "root-ec-halving.h"
#include "root-quaternion-splitting.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
int a_kani_geometry_product(fq_struct *out,const a_ec_point *P,unsigned g,const root_ec_theta4 *T){
 if(!out||!P||!T||!T->ready||!g||g>8)return 0;size_t n=(size_t)1<<(2*g);const fq_ctx_struct *F=T->E->field;
 fq_struct *small=b_theta4_vector_new(4*g,F),*v=b_theta4_vector_new(n,F);if(!small||!v){b_theta4_vector_free(small,4*g,F);b_theta4_vector_free(v,n,F);return 0;}
 int ok=0;for(unsigned j=0;j<g;j++)if(!root_ec_theta4_map(small+4*j,P+j,T))goto done;
 fq_one(v,F);size_t length=1;for(unsigned j=0;j<g;j++){for(size_t k=length*4;k-->0;)fq_mul(v+k,v+k%length,small+4*j+k/length,F);length*=4;}
 for(size_t i=0;i<n;i++)fq_set(out+i,v+i,F);ok=1;
done:b_theta4_vector_free(v,n,F);b_theta4_vector_free(small,4*g,F);return ok;
}
static int prime(ulong n){if(n<2)return 0;for(ulong d=2;d<=n/d;d++)if(n%d==0)return 0;return 1;}
static int check_spec(const a_kani_geometry_spec *S,const fmpz_t p){
 if(!S||(S->dimension!=2&&S->dimension!=8)||!S->n1||!S->n2||!S->N||!S->endomorphism||(S->denominator!=1&&S->denominator!=2)||!S->prime_count||!S->primes||!S->exponents||!S->extensions||S->primes[0]!=2)return 0;
 if(fmpz_sgn(S->n1)<=0||fmpz_sgn(S->n2)<=0)return 0;const fmpz *c=S->endomorphism;if(S->denominator==2&&(fmpz_fdiv_ui(c,2)!=fmpz_fdiv_ui(c+3,2)||fmpz_fdiv_ui(c+1,2)!=fmpz_fdiv_ui(c+2,2)))return 0;
 fmpz_t t,n,mod,power,u,v;fmpz_init(t);fmpz_init(n);fmpz_init(mod);fmpz_init(power);fmpz_init(u);fmpz_init(v);int ok=0;
 fmpz_add(t,S->n1,S->n2);if(!fmpz_equal(t,S->N))goto done;fmpz_gcd(t,S->n1,S->n2);if(!fmpz_is_one(t))goto done;
 fmpz_one(n);for(size_t j=0;j<S->prime_count;j++){
  ulong ell=S->primes[j];unsigned e=S->exponents[j],guard=j?0:3;if(!e||e>UINT_MAX-guard||ell>ULONG_MAX||!prime(ell)||!S->extensions[j])goto done;for(size_t k=0;k<j;k++)if(S->primes[k]==ell)goto done;
  fmpz_set_ui(mod,ell);fmpz_pow_ui(mod,mod,e);fmpz_mul(n,n,mod);fmpz_set_ui(mod,ell);fmpz_pow_ui(mod,mod,e+guard);fmpz_neg(t,p);fmpz_powm_ui(power,t,S->extensions[j],mod);if(!fmpz_is_one(power))goto done;
 }
 if(!fmpz_equal(n,S->N))goto done;
 fmpz_mul(n,c+2,c+2);fmpz_addmul(n,c+3,c+3);fmpz_mul(n,n,p);fmpz_addmul(n,c,c);fmpz_addmul(n,c+1,c+1);
 if(S->dimension==2){fmpz_mul(t,S->n1,S->n2);fmpz_mul_ui(t,t,S->denominator*S->denominator);ok=fmpz_equal(n,t);goto done;}
 if(!S->X||!S->Y||!S->d1||!S->d2||fmpz_sgn(S->d1)<=0||fmpz_sgn(S->d2)<=0||fmpz_mat_nrows(S->X)!=4||fmpz_mat_ncols(S->X)!=4||fmpz_mat_nrows(S->Y)!=4||fmpz_mat_ncols(S->Y)!=4||!fmpz_divisible(S->n1,S->d1)||!fmpz_divisible(S->n2,S->d2))goto done;
 fmpz_mul(t,S->d1,S->d2);fmpz_mul_ui(t,t,S->denominator*S->denominator);if(!fmpz_equal(n,t))goto done;fmpz_divexact(u,S->n1,S->d1);fmpz_divexact(v,S->n2,S->d2);
 fmpz_mat_t A,B,C;fmpz_mat_init(A,4,4);fmpz_mat_init(B,4,4);fmpz_mat_init(C,4,4);ok=1;
 fmpz_mat_transpose(A,S->X);fmpz_mat_mul(B,A,S->X);for(unsigned i=0;i<4;i++)for(unsigned j=0;j<4;j++)if(i==j?!fmpz_equal(fmpz_mat_entry(B,i,j),u):!fmpz_is_zero(fmpz_mat_entry(B,i,j)))ok=0;
 fmpz_mat_transpose(A,S->Y);fmpz_mat_mul(B,A,S->Y);for(unsigned i=0;i<4;i++)for(unsigned j=0;j<4;j++)if(i==j?!fmpz_equal(fmpz_mat_entry(B,i,j),v):!fmpz_is_zero(fmpz_mat_entry(B,i,j)))ok=0;
 fmpz_mat_mul(B,S->X,S->Y);fmpz_mat_mul(C,S->Y,S->X);if(!fmpz_mat_equal(B,C))ok=0;fmpz_mat_mul(B,S->X,A);fmpz_mat_mul(C,A,S->X);if(!fmpz_mat_equal(B,C))ok=0;
 fmpz_mat_clear(C);fmpz_mat_clear(B);fmpz_mat_clear(A);
done:fmpz_clear(v);fmpz_clear(u);fmpz_clear(power);fmpz_clear(mod);fmpz_clear(n);fmpz_clear(t);return ok;
}
static int lattice_data(a_kani_geometry *G,const a_kani_geometry_spec *S,const a_ec_point basis[2]){
 a_kani_geometry_component *C=G->components;const root_extension *E=&C->field;unsigned g=G->dimension,e=S->exponents[0],bits=e+2;if(e>(UINT_MAX-4)/2)return 0;unsigned high=2*e+4;int ok=0;
 a_ec_point half[2],points[8],temp;for(unsigned i=0;i<2;i++)a_ec_init(half+i,E);for(unsigned i=0;i<g;i++)a_ec_init(points+i,E);a_ec_init(&temp,E);
 fmpz_mat_struct known[4],aligned[4];for(unsigned i=0;i<4;i++){fmpz_mat_init(known+i,2,2);fmpz_mat_init(aligned+i,2,2);}
 fmpz_mat_t theta,H;fmpz_mat_init(theta,2,2);fmpz_mat_init(H,g,g);fmpz_t h,inv,t;fmpz_init(h);fmpz_init(inv);fmpz_init(t);fmpz_one(h);fmpz_mul_2exp(h,h,e);
 for(unsigned i=0;i<2;i++)if(!root_ec_half_e0(half+i,basis+i,E))goto done;
 if(!root_ec_maxorder_actions(known,basis,half,bits,E)||!root_quaternion_splitting_align(aligned,known,fq_ctx_prime(E->field),bits,high)||!root_quaternion_splitting_element(theta,S->endomorphism,S->denominator,aligned,high))goto done;
 if(g==2){fmpz_mat_set(G->eta,theta);if(!root_tate_projector_input_g2(G->Pnum,G->eta,S->n1,S->n2))goto done;}
 else if(!root_tate_projector_input(G->Pnum,G->eta,S->X,S->Y,theta,S->d1,S->d2,0))goto done;
 if(!fmpz_invmod(inv,S->n1,h))goto done;for(unsigned i=0;i<g;i++)for(unsigned j=0;j<g;j++){fmpz_mul(t,inv,fmpz_mat_entry(G->eta,i,j));fmpz_mod(fmpz_mat_entry(H,i,j),t,h);}
 G->lattice=b_theta4_vector_new(2*g*G->count,E->field);if(!G->lattice)goto done;
 for(unsigned col=0;col<2*g;col++){
  for(unsigned i=0;i<g;i++)a_ec_infinity(points+i,E);
  if(col<g){a_ec_set(points+col/2,basis+(col&1),E);for(unsigned row=0;row<g;row++){a_ec_mul(&temp,basis+(row&1),fmpz_mat_entry(H,row,col),C->model.A,E);a_ec_add(points+g/2+row/2,points+g/2+row/2,&temp,C->model.A,E);}}
  else a_ec_mul(points+g/2+(col-g)/2,basis+((col-g)&1),h,C->model.A,E);
  if(!a_kani_geometry_product(G->lattice+col*G->count,points,g,&C->model))goto done;G->lattice_ptrs[col]=G->lattice+col*G->count;
 }
 ok=1;
done:fmpz_clear(t);fmpz_clear(inv);fmpz_clear(h);fmpz_mat_clear(H);fmpz_mat_clear(theta);for(unsigned i=0;i<4;i++){fmpz_mat_clear(aligned+i);fmpz_mat_clear(known+i);}a_ec_clear(&temp,E);for(unsigned i=0;i<g;i++)a_ec_clear(points+i,E);for(unsigned i=0;i<2;i++)a_ec_clear(half+i,E);return ok;
}
void a_kani_geometry_clear(a_kani_geometry *G){
 if(!G||!G->ready)return;a_kani_fullchain_clear(&G->chain);size_t n=G->count;unsigned g=G->dimension;
 if(G->components&&G->components[0].field_ready)b_theta4_vector_free(G->lattice,2*g*n,G->components[0].field.field);
 for(size_t j=0;j<G->prime_count;j++){
  a_kani_geometry_component *C=G->components+j;if(C->field_ready){b_theta4_vector_free(C->generators,g*n,C->field.field);b_theta4_vector_free(C->complements,g*n,C->field.field);b_theta4_vector_free(C->pairing4,g*g,C->field.field);root_ec_theta4_clear(&C->model);root_field_embedding_clear(&C->embedding);root_extension_clear(&C->field);}
 }
 free(G->prime_inputs);free(G->components);b_theta4_vector_free(G->null,n,G->base_model->E->field);fq_clear(G->zeta8,G->base_model->E->field);fq_clear(G->zeta4,G->base_model->E->field);fmpz_clear(G->N);fmpz_clear(G->n1);fmpz_mat_clear(G->eta);fmpz_mat_clear(G->Pnum);memset(G,0,sizeof(*G));
}
int a_kani_geometry_init(a_kani_geometry *G,const root_ec_theta4 *B,const a_kani_geometry_spec *S,a_ec_rng rng,void *context,unsigned long budget){
 if(G&&!G->ready)G->last_error=-1;
 if(!G||G->ready||!B||!B->ready||B->E->relative_degree!=1||!fq_is_zero(B->A,B->E->field)||!rng||!check_spec(S,fq_ctx_prime(B->E->field)))return A_EC_INVALID;
 unsigned g=S->dimension;G->base_model=B;G->dimension=g;G->count=(size_t)1<<(2*g);G->ready=1;fmpz_init_set(G->n1,S->n1);fmpz_init_set(G->N,S->N);fmpz_mat_init(G->Pnum,2*g,2*g);fmpz_mat_init(G->eta,g,g);fq_init(G->zeta4,B->E->field);fq_init(G->zeta8,B->E->field);fq_set(G->zeta4,B->E->i,B->E->field);int status=A_EC_INVALID;
 if(!root_extension_sqrt(G->zeta8,G->zeta4,B->E))goto fail;G->components=calloc(S->prime_count,sizeof(*G->components));G->prime_inputs=calloc(S->prime_count,sizeof(*G->prime_inputs));if(!G->components||!G->prime_inputs)goto fail;G->prime_count=S->prime_count;
 G->null=b_theta4_vector_new(G->count,B->E->field);if(!G->null)goto fail;a_ec_point zero[8];for(unsigned i=0;i<g;i++){a_ec_init(zero+i,B->E);a_ec_infinity(zero+i,B->E);}int mapped=a_kani_geometry_product(G->null,zero,g,B);for(unsigned i=0;i<g;i++)a_ec_clear(zero+i,B->E);if(!mapped)goto fail;
 fmpz_mat_t matrix;fmpz_mat_init(matrix,4,4);if(g==8)fmpz_mat_mul(matrix,S->Y,S->X);
 for(size_t k=0;k<S->prime_count;k++){
  a_kani_geometry_component *C=G->components+k;const root_extension *E=&C->field;
  if(!root_extension_init_rng(&C->field,fq_ctx_prime(B->E->field),S->extensions[k],rng,context)){status=C->field.last_error==-2?A_EC_RNG_ERROR:A_EC_INVALID;goto loop_fail;}C->field_ready=1;
  if(!root_field_embedding_init(&C->embedding,B->E,E)){status=C->embedding.last_error==-2?A_EC_RNG_ERROR:A_EC_INVALID;goto loop_fail;}
  if(!root_ec_theta4_base_change(&C->model,B,E))goto loop_fail;
  C->generators=b_theta4_vector_new(g*G->count,E->field);if(k){C->complements=b_theta4_vector_new(g*G->count,E->field);C->pairing4=b_theta4_vector_new(g*g,E->field);}if(!C->generators||(k&&(!C->complements||!C->pairing4)))goto loop_fail;
  a_kani_two_kernel K2;a_clapoti_kernel K8;if(g==2)a_kani_two_kernel_init(&K2,E);else a_clapoti_kernel_init(&K8,E);
  status=g==2?a_kani_two_kernel_lift(&K2,S->primes[k],S->exponents[k],k?0:2,S->n1,S->n2,S->endomorphism,S->denominator,rng,context,budget,NULL,E):a_clapoti_kernel_lift(&K8,(unsigned long)S->primes[k],S->exponents[k],k?0:2,S->n1,S->n2,S->endomorphism,S->denominator,matrix,rng,context,budget,NULL,E);
  if(status==A_EC_OK){
   const a_ec_point *basis=g==2?K2.basis:K8.basis;const fmpz *order=g==2?K2.order:K8.order;int good=1;
   for(unsigned j=0;good&&j<g;j++){good=a_kani_geometry_product(C->generators+j*G->count,g==2?K2.generator[j]:K8.generator[j],g,&C->model);if(good&&k)good=a_kani_geometry_product(C->complements+j*G->count,g==2?K2.complement[j]:K8.complement[j],g,&C->model);C->generator_ptrs[j]=C->generators+j*G->count;if(k)C->complement_ptrs[j]=C->complements+j*G->count;}
   if(good&&!k)good=lattice_data(G,S,basis);
   if(good&&k){fq_t z;fq_init(z,E->field);fmpz_t power;fmpz_init(power);status=a_ec_weil(z,basis,basis+1,order,C->model.A,rng,context,budget,E);good=status==A_EC_OK;if(good){fmpz_divexact_ui(power,order,S->primes[k]);fmpz_mul_ui(power,power,4);fq_pow(z,z,power,E->field);for(unsigned i=0;i<g;i++)for(unsigned j=0;j<g;j++){if(i==j)fq_set(C->pairing4+i*g+j,z,E->field);else fq_one(C->pairing4+i*g+j,E->field);}}fmpz_clear(power);fq_clear(z,E->field);}
   if(!good&&status==A_EC_OK)status=A_EC_INVALID;
  }
  if(g==2)a_kani_two_kernel_clear(&K2,E);else a_clapoti_kernel_clear(&K8,E);if(status!=A_EC_OK)goto loop_fail;
  G->prime_inputs[k]=(a_kani_prime_input){S->primes[k],S->exponents[k],E,C->generator_ptrs,C->complement_ptrs,C->pairing4};
 }
 fmpz_mat_clear(matrix);G->last_error=0;return A_EC_OK;
loop_fail:fmpz_mat_clear(matrix);
fail:a_kani_geometry_clear(G);G->last_error=status==A_EC_RNG_ERROR?-2:-1;return status==A_EC_OK?A_EC_INVALID:status;
}
int a_kani_geometry_start(a_kani_geometry *G,const a_theta_payload_input *payloads,size_t groups){
 if(G)G->last_error=-1;
 if(!G||!G->ready||G->started||(groups&&!payloads)||groups==SIZE_MAX||groups>SIZE_MAX/sizeof(a_theta_payload_input)-1)return 0;
 a_theta_payload_input *all=calloc(groups+1,sizeof(*all));if(!all)return 0;all[0]=(a_theta_payload_input){&G->components[0].field,G->lattice_ptrs,2*G->dimension};if(groups)memcpy(all+1,payloads,groups*sizeof(*all));int good=a_kani_fullchain_init(&G->chain,G->dimension,G->base_model->E,G->null,G->zeta4,G->zeta8,G->prime_inputs,G->prime_count,all,groups+1);free(all);if(good)G->started=1;G->last_error=good?0:G->chain.last_error;return good;
}
int a_kani_geometry_next(a_kani_geometry *G,a_kani_progress progress,void *context){if(!G)return 0;G->last_error=-1;if(!G->ready||!G->started)return 0;int ok=a_kani_fullchain_next(&G->chain,progress,context);G->last_error=ok?0:G->chain.last_error;return ok;}
