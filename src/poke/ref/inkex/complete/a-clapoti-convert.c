#include "a-clapoti-convert.h"
#include "root-kani-factor.h"
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <inke_uniform.h>

int a_clapoti_fixed_payloads(a_ec_point two[24],a_ec_point three[24],
 const quat_left_ideal_t *I,const a_clapoti_certificate *C,const root_extension *E){
 if(!two||!E||E->relative_degree!=1||!a_clapoti_certificate_verify(C,I))return 0;
 fmpz_t p;fmpz_init(p);fmpz_set_mpz(p,QUATALG_PINFTY.p);int ok=fmpz_equal(p,fq_ctx_prime(E->field));fmpz_clear(p);if(!ok)return 0;
 const unsigned batches=three?2:1;a_ec_point temp[48],fixed[3];
 for(unsigned i=0;i<48;i++)a_ec_init(temp+i,E);for(unsigned i=0;i<3;i++)a_ec_init(fixed+i,E);
 fq_t A;fq_init(A,E->field);fq_zero(A,E->field);ibz_t inverse;ibz_init(&inverse);fmpz_t coefficient;fmpz_init(coefficient);
 for(unsigned b=0;b<batches&&ok;b++){
  ok=mpz_invert(inverse,C->n1,b?TORSION_PLUS_3POWER:TORSION_PLUS_2POWER)&&b_poke_fixed_action(fixed,(int)b,&C->beta1,&inverse,E);
  for(unsigned i=0;i<3&&ok;i++)for(unsigned j=0;j<8;j++){
   a_ec_point *target=temp+24*b+8*i+j;
   if(j<4){fmpz_set_mpz(coefficient,C->X[0][j]);a_ec_mul(target,fixed+i,coefficient,A,E);}else a_ec_infinity(target,E);
  }
 }
 if(ok){for(unsigned i=0;i<24;i++)a_ec_set(two+i,temp+i,E);if(three)for(unsigned i=0;i<24;i++)a_ec_set(three+i,temp+24+i,E);}
 fmpz_clear(coefficient);ibz_finalize(&inverse);fq_clear(A,E->field);for(unsigned i=0;i<3;i++)a_ec_clear(fixed+i,E);for(unsigned i=0;i<48;i++)a_ec_clear(temp+i,E);return ok;
}

int a_clapoti_convert(inke_conversion_result *out,ec_basis_t *three,
 const quat_left_ideal_t *I,const a_clapoti_certificate *C,const root_ec_theta4 *B,
 a_ec_rng rng,void *rng_context,unsigned long budget,a_kani_progress progress,
 void *progress_context,a_clapoti_convert_stage *stage){
 a_clapoti_convert_stage state={1,0,0};if(stage)*stage=state;
 if(!out||!B||!B->ready||!rng||B->E->relative_degree!=1||!fq_is_zero(B->A,B->E->field)||!a_clapoti_certificate_verify(C,I))return A_CLAPOTI_CONVERT_INVALID;
 const root_extension *E=B->E;const size_t count=three?6:3,dimension=8,coordinates=65536;
 int status=A_CLAPOTI_CONVERT_INVALID;unsigned steps=0;fmpz_t n1,n2,N,d1,d2;
 fmpz_init(n1);fmpz_init(n2);fmpz_init(N);fmpz_init(d1);fmpz_init(d2);
 fmpz_set_mpz(n1,C->n1);fmpz_set_mpz(n2,C->n2);fmpz_set_mpz(N,C->N);fmpz_set_mpz(d1,C->d1);fmpz_set_mpz(d2,C->d2);
 fmpz theta[4];for(unsigned i=0;i<4;i++){fmpz_init(theta+i);fmpz_set_mpz(theta+i,C->theta.coord[i]);}
 fmpz_mat_t X,Y;fmpz_mat_init(X,4,4);fmpz_mat_init(Y,4,4);
 for(unsigned i=0;i<4;i++)for(unsigned j=0;j<4;j++){fmpz_set_mpz(fmpz_mat_entry(X,i,j),C->X[i][j]);fmpz_set_mpz(fmpz_mat_entry(Y,i,j),C->Y[i][j]);}
 ulong *primes=calloc(C->prime_count,sizeof(*primes));
 a_kani_geometry G={0};a_ec_point inputs[48],answers[6];for(unsigned i=0;i<48;i++)a_ec_init(inputs+i,E);for(unsigned i=0;i<6;i++)a_ec_init(answers+i,E);
 fq_t A;fq_init(A,E->field);fq_struct *payload=NULL;const fq_struct *payload_ptrs[6];
 if(!primes){status=A_CLAPOTI_CONVERT_MEMORY_ERROR;goto done;}
 for(size_t i=0;i<C->prime_count;i++){primes[i]=(ulong)C->primes[i];if(C->exponents[i]>(unsigned)INT_MAX-steps)goto done;steps+=C->exponents[i];}
 a_kani_geometry_spec spec={0};spec.dimension=8;spec.n1=n1;spec.n2=n2;spec.N=N;spec.d1=d1;spec.d2=d2;spec.X=X;spec.Y=Y;
 spec.endomorphism=theta;spec.denominator=(unsigned)mpz_get_ui(C->theta.denom);spec.prime_count=C->prime_count;spec.primes=primes;spec.exponents=C->exponents;spec.extensions=C->extensions;
 state.phase=2;if(!a_clapoti_fixed_payloads(inputs,three?inputs+24:NULL,I,C,E))goto done;
 payload=b_theta4_vector_new(count*coordinates,E->field);if(!payload){status=A_CLAPOTI_CONVERT_MEMORY_ERROR;goto done;}
 for(size_t i=0;i<count;i++){if(!a_kani_geometry_product(payload+i*coordinates,inputs+i*dimension,8,B))goto done;payload_ptrs[i]=payload+i*coordinates;}
 state.phase=3;status=a_kani_geometry_init(&G,B,&spec,rng,rng_context,budget);
 if(status!=A_EC_OK){status=status==A_EC_RETRY?A_CLAPOTI_CONVERT_RETRY:status==A_EC_RNG_ERROR?A_CLAPOTI_CONVERT_RNG_ERROR:A_CLAPOTI_CONVERT_GEOMETRY_ERROR;goto done;}
 a_theta_payload_input group={E,payload_ptrs,count};if(!a_kani_geometry_start(&G,&group,1)){status=G.last_error==-2?A_CLAPOTI_CONVERT_RNG_ERROR:A_CLAPOTI_CONVERT_GEOMETRY_ERROR;goto done;}
 state.phase=4;for(size_t i=0;i<C->prime_count;i++){if(!a_kani_geometry_next(&G,progress,progress_context)){status=G.last_error==-2?A_CLAPOTI_CONVERT_RNG_ERROR:A_CLAPOTI_CONVERT_GEOMETRY_ERROR;goto done;}state.components_completed++;}
 state.phase=5;if(!root_kani_extract(A,answers,count,8,&G.chain.theta,G.chain.groups[0].points,&G.components[0].embedding,
  G.chain.groups[1].points,G.Pnum,G.eta,G.n1,G.N,G.zeta4,G.zeta8,E,&state.detail)){status=A_CLAPOTI_CONVERT_FACTOR_ERROR;goto done;}
 state.phase=6;inke_conversion_result result;ec_basis_t result_three;ec_curve_t three_curve;
 if(!b_poke_conversion_export(&result,A,answers,&I->norm,(int)steps,E)){status=A_CLAPOTI_CONVERT_VERIFY_ERROR;goto done;}
 if(three&&(!b_poke_basis_export(&three_curve,&result_three,A,answers+3,E)||
  !test_point_order_threef(&result_three.P,&three_curve,POWER_OF_3)||!test_point_order_threef(&result_three.Q,&three_curve,POWER_OF_3)||!test_point_order_threef(&result_three.PmQ,&three_curve,POWER_OF_3))){status=A_CLAPOTI_CONVERT_VERIFY_ERROR;goto done;}
 *out=result;if(three)*three=result_three;status=A_CLAPOTI_CONVERT_OK;state.phase=0;state.detail=0;
done:
 if(stage)*stage=state;a_kani_geometry_clear(&G);b_theta4_vector_free(payload,count*coordinates,E->field);fq_clear(A,E->field);
 for(unsigned i=0;i<6;i++)a_ec_clear(answers+i,E);for(unsigned i=0;i<48;i++)a_ec_clear(inputs+i,E);
 free(primes);fmpz_mat_clear(Y);fmpz_mat_clear(X);for(unsigned i=0;i<4;i++)fmpz_clear(theta+i);
 fmpz_clear(d2);fmpz_clear(d1);fmpz_clear(N);fmpz_clear(n2);fmpz_clear(n1);return status;
}

int a_clapoti_fallback_init(a_clapoti_fallback *S,const root_ec_theta4 *B,
 a_ec_rng rng,void *context,unsigned long budget,a_kani_progress progress,void *progress_context){
 if(!S||S->ready||!B||!B->ready||!rng)return 0;
 memset(S,0,sizeof(*S));a_clapoti_certificate_init(&S->certificate);S->base_model=B;S->rng=rng;S->rng_context=context;S->point_budget=budget;S->progress=progress;S->progress_context=progress_context;S->ready=1;return 1;
}
void a_clapoti_fallback_clear(a_clapoti_fallback *S){if(!S||!S->ready)return;a_clapoti_certificate_clear(&S->certificate);memset(S,0,sizeof(*S));}
int a_clapoti_fallback_convert(inke_conversion_result *out,const quat_left_ideal_t *I,void *context){
 a_clapoti_fallback *S=context;if(!S||!S->ready||!out)return INKE_SAMPLE_INVALID;
 if(!S->prepared){
  S->last_status=a_clapoti_prepare(&S->certificate,I,&S->counts);
  if(S->last_status!=A_CLAPOTI_OK)return S->last_status==A_CLAPOTI_RNG_ERROR?INKE_SAMPLE_RNG_ERROR:INKE_SAMPLE_INVALID;
  S->prepared=1;
 }else if(!a_clapoti_certificate_verify(&S->certificate,I)){S->last_status=A_CLAPOTI_CONVERT_INVALID;return INKE_SAMPLE_INVALID;}
 S->last_status=a_clapoti_convert(out,NULL,I,&S->certificate,S->base_model,S->rng,S->rng_context,S->point_budget,S->progress,S->progress_context,&S->stage);
 return S->last_status==A_CLAPOTI_CONVERT_OK?INKE_SAMPLE_OK:
  S->last_status==A_CLAPOTI_CONVERT_RETRY?INKE_SAMPLE_RETRY:
  S->last_status==A_CLAPOTI_CONVERT_RNG_ERROR?INKE_SAMPLE_RNG_ERROR:INKE_SAMPLE_INVALID;
}
