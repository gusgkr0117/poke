// SPDX-License-Identifier: Apache-2.0
#include "b-gaussian-convert.h"
#include "root-kani-factor.h"
#include "root-poke-field-bridge.h"
#include <endomorphism_action.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <limits.h>
#include <string.h>

int b_gaussian_converter_init(b_gaussian_converter *C){
    if(C->ready||C->field_ready)return 0;
    fmpz_t p;fmpz_init(p);fmpz_set_mpz(p,QUATALG_PINFTY.p);
    int ok=root_extension_init(&C->field,p,1);fmpz_clear(p);if(!ok)return 0;C->field_ready=1;
    fq_t A,t;fq_init(A,C->field.field);fq_init(t,C->field.field);root_from_poke_fp2(A,&CURVE_E0.A,&C->field);root_from_poke_fp2(t,&CURVE_E0.C,&C->field);
    if(fq_is_zero(t,C->field.field))ok=0;else{fq_inv(t,t,C->field.field);fq_mul(A,A,t,C->field.field);ok=root_ec_theta4_init(&C->model,A,&C->field);}
    fq_clear(t,C->field.field);fq_clear(A,C->field.field);
    if(!ok){b_gaussian_converter_clear(C);return 0;}C->ready=1;return 1;
}
void b_gaussian_converter_clear(b_gaussian_converter *C){
    root_ec_theta4_clear(&C->model);if(C->field_ready)root_extension_clear(&C->field);memset(C,0,sizeof *C);
}
int b_gaussian_convert(inke_conversion_result *out,ec_basis_t *three,const quat_left_ideal_t *I,
 const b_gaussian_certificate *C,const b_gaussian_converter *context,a_ec_rng rng,void *rng_context,
 unsigned long budget,a_kani_progress progress,void *progress_context,b_gaussian_convert_stage *stage){
    b_gaussian_convert_stage state={1,0,0};if(stage)*stage=state;
    if(!context->ready||!rng||!b_gaussian_certificate_verify(C,I))return B_GAUSSIAN_CONVERT_INVALID;
    const root_extension *E=&context->field;const size_t count=three?6:3;
    int status=B_GAUSSIAN_CONVERT_INVALID;unsigned steps=0;
    fmpz_t n1,n2,N;fmpz_init(n1);fmpz_init(n2);fmpz_init(N);fmpz_set_mpz(n1,C->n1);fmpz_set_mpz(n2,C->n2);fmpz_set_mpz(N,C->N);
    fmpz eta[4];for(int i=0;i<4;i++){fmpz_init(eta+i);fmpz_set_mpz(eta+i,C->eta.coord[i]);}
    ulong primes[11];unsigned exponents[11],extensions[11];size_t prime_count=0;
    for(int i=0;i<11;i++)if(C->exponents[i]){primes[prime_count]=C->primes[i];exponents[prime_count]=C->exponents[i];extensions[prime_count]=C->extensions[i];steps+=C->exponents[i];prime_count++;}
    a_kani_geometry_spec spec={0};spec.dimension=2;spec.n1=n1;spec.n2=n2;spec.N=N;spec.endomorphism=eta;spec.denominator=(unsigned)mpz_get_ui(C->eta.denom);spec.prime_count=prime_count;spec.primes=primes;spec.exponents=exponents;spec.extensions=extensions;
    a_kani_geometry G={0};a_ec_point inputs[6],answers[6],pair[2];for(int i=0;i<6;i++){a_ec_init(inputs+i,E);a_ec_init(answers+i,E);}for(int i=0;i<2;i++)a_ec_init(pair+i,E);
    fq_t A;fq_init(A,E->field);fq_struct *payload=b_theta4_vector_new(count*16,E->field);const fq_struct *payload_ptrs[6];
    if(!payload){status=B_GAUSSIAN_CONVERT_MEMORY_ERROR;goto done;}
    state.phase=2;
    if(three){if(!b_gaussian_fixed_payloads(inputs,inputs+3,I,C,E))goto done;}
    else{ibz_t inverse;ibz_init(&inverse);int ok=mpz_invert(inverse,C->n1,TORSION_PLUS_2POWER)&&b_poke_fixed_action(inputs,0,&C->transport,&inverse,E);ibz_finalize(&inverse);if(!ok)goto done;}
    a_ec_infinity(pair+1,E);
    for(size_t i=0;i<count;i++){a_ec_set(pair,inputs+i,E);if(!a_kani_geometry_product(payload+i*16,pair,2,&context->model))goto done;payload_ptrs[i]=payload+i*16;}
    state.phase=3;status=a_kani_geometry_init(&G,&context->model,&spec,rng,rng_context,budget);
    if(status!=A_EC_OK){status=status==A_EC_RETRY?B_GAUSSIAN_CONVERT_RETRY:status==A_EC_RNG_ERROR?B_GAUSSIAN_CONVERT_RNG_ERROR:B_GAUSSIAN_CONVERT_GEOMETRY_ERROR;goto done;}
    a_theta_payload_input group={E,payload_ptrs,count};
    if(!a_kani_geometry_start(&G,&group,1)){status=G.last_error==-2?B_GAUSSIAN_CONVERT_RNG_ERROR:B_GAUSSIAN_CONVERT_GEOMETRY_ERROR;goto done;}
    state.phase=4;
    for(size_t i=0;i<prime_count;i++){if(!a_kani_geometry_next(&G,progress,progress_context)){status=G.last_error==-2?B_GAUSSIAN_CONVERT_RNG_ERROR:B_GAUSSIAN_CONVERT_GEOMETRY_ERROR;goto done;}state.components_completed++;}
    state.phase=5;
    if(!root_kani_extract(A,answers,count,2,&G.chain.theta,G.chain.groups[0].points,&G.components[0].embedding,
          G.chain.groups[1].points,G.Pnum,G.eta,G.n1,G.N,G.zeta4,G.zeta8,E,&state.detail)){status=B_GAUSSIAN_CONVERT_FACTOR_ERROR;goto done;}
    state.phase=6;inke_conversion_result result;ec_basis_t result_three;ec_curve_t three_curve;
    if(steps>INT_MAX||!b_poke_conversion_export(&result,A,answers,&I->norm,(int)steps,E)){status=B_GAUSSIAN_CONVERT_VERIFY_ERROR;goto done;}
    if(three){
        if(!b_poke_basis_export(&three_curve,&result_three,A,answers+3,E)||
           !test_point_order_threef(&result_three.P,&three_curve,POWER_OF_3)||
           !test_point_order_threef(&result_three.Q,&three_curve,POWER_OF_3)||
           !test_point_order_threef(&result_three.PmQ,&three_curve,POWER_OF_3)){
            status=B_GAUSSIAN_CONVERT_VERIFY_ERROR;goto done;
        }
    }
    *out=result;if(three)*three=result_three;status=B_GAUSSIAN_CONVERT_OK;state.phase=0;state.detail=0;
done:
    if(stage)*stage=state;a_kani_geometry_clear(&G);b_theta4_vector_free(payload,count*16,E->field);fq_clear(A,E->field);
    for(int i=0;i<2;i++)a_ec_clear(pair+i,E);for(int i=0;i<6;i++){a_ec_clear(answers+i,E);a_ec_clear(inputs+i,E);}
    for(int i=0;i<4;i++)fmpz_clear(eta+i);fmpz_clear(N);fmpz_clear(n2);fmpz_clear(n1);return status;
}
