// SPDX-License-Identifier: Apache-2.0
#include "b-uniform-keygen.h"
#include <endomorphism_action.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <string.h>
#ifdef INKE_MANUSCRIPT_PARAMETERS
#include <inke_field_constants.h>
#endif
static unsigned public_exponent(void){
#ifdef INKE_PUBLIC_TWO_EXPONENT
    return INKE_PUBLIC_TWO_EXPONENT;
#else
    return (unsigned)TORSION_PLUS_EVEN_POWER-2;
#endif
}
void b_uniform_keygen_init(b_uniform_keygen *S){
    memset(S,0,sizeof *S);ibz_init(&S->q);quat_left_ideal_init(&S->ideal);
    inke_mixed_state_init(&S->mixed);b_gaussian_certificate_init(&S->gaussian);S->a=public_exponent();
}
void b_uniform_keygen_clear(b_uniform_keygen *S){
    b_gaussian_certificate_clear(&S->gaussian);inke_mixed_state_clear(&S->mixed);
    quat_left_ideal_finalize(&S->ideal);ibz_finalize(&S->q);memset(S,0,sizeof *S);
}
int b_uniform_keygen_step(inke_sk_t *sk,inke_pk_t *pk,b_uniform_keygen *S,
    const b_uniform_keygen_options *O,const b_gaussian_plans *plans,
    const b_gaussian_converter *converter,a_ec_rng rng,void *rng_context,
    a_kani_progress progress,void *progress_context){
    if(S->complete||S->a!=public_exponent()||!O||!plans||!converter||!rng)return INKE_SAMPLE_INVALID;
    int status;
    if(!S->has_q){S->phase=1;status=inke_sample_q(&S->q,S->a);if(status!=INKE_SAMPLE_OK)return status;S->has_q=1;}
    if(!S->has_ideal&&!S->has_conversion){
        S->phase=2;status=inke_sample_ideal(&S->ideal,&S->q,S->a,&MAXORD_O0,&QUATALG_PINFTY,
              O->prime_trials,O->norm_trials,O->unit_trials,&S->sampling);
        if(status!=INKE_SAMPLE_OK)return status;
        /* Commit immediately: conversion/preparation has no authority to
         * change this original receiver ideal, even if a later stage fails. */
        S->has_ideal=1;
    }
    if(!S->has_conversion&&!S->mixed_done){
        S->phase=3;
        if(!O->mixed_trials)S->mixed_done=1;
        else{
            if(!S->mixed_ready){
                ibz_t target;ibz_init(&target);mpz_mul_2exp(target,TORSION_PLUS_3POWER,TORSION_PLUS_EVEN_POWER-2);
                status=inke_mixed_prepare(&S->mixed,&S->ideal,&target,&QUATALG_PINFTY);ibz_finalize(&target);
                if(status!=INKE_SAMPLE_OK)return status;S->mixed_ready=1;
            }
            inke_mixed_witness witness;inke_mixed_witness_init(&witness);
            status=INKE_SAMPLE_RETRY;
            while(!S->has_conversion&&S->mixed.attempts<O->mixed_trials){
                unsigned long before=S->mixed.attempts;
                status=inke_mixed_search(&witness,&S->mixed,O->mixed_trials-before);
                if(status!=INKE_SAMPLE_OK&&status!=INKE_SAMPLE_RETRY)break;
                if(S->mixed.attempts<=before||S->mixed.attempts>O->mixed_trials){
                    status=INKE_SAMPLE_INVALID;break;
                }
                if(status==INKE_SAMPLE_OK){
                    S->mixed_witnesses++;inke_conversion_result result;
                    if(inke_convert_witness(&result,&S->mixed,&witness,&S->last_stage)){
                        S->conversion=result;S->has_conversion=1;S->conversion_origin=1;
                    }else S->mixed_failures++;
                }
            }
            inke_mixed_witness_clear(&witness);
            if(status!=INKE_SAMPLE_OK&&status!=INKE_SAMPLE_RETRY)return status;
            /* The public budget counts ALL sampled candidates, including
             * failed geometric witnesses and earlier calls interrupted by
             * RNG failure. A witness is not a successful conversion. */
            S->mixed_done=1;
        }
    }
    if(!S->has_conversion&&!S->gaussian_prepared){
        S->phase=4;status=b_gaussian_prepare(&S->gaussian,&S->ideal,plans,&S->gaussian_counts);
        if(status==B_GAUSSIAN_OK)S->gaussian_prepared=1;
        else if(status==B_GAUSSIAN_MISS)S->gaussian_prepared=-1;
        else return B_KEYGEN_CONVERSION_ERROR;
    }
    if(!S->has_conversion&&S->gaussian_prepared==1&&!S->gaussian_failed){
        S->phase=5;S->gaussian_attempts++;inke_conversion_result result;
        status=b_gaussian_convert(&result,NULL,&S->ideal,&S->gaussian,converter,rng,rng_context,
                    O->point_trials,progress,progress_context,&S->conversion_stage);
        if(status==B_GAUSSIAN_CONVERT_OK){S->conversion=result;S->has_conversion=1;S->conversion_origin=2;}
        else if(status==B_GAUSSIAN_CONVERT_RETRY)return INKE_SAMPLE_RETRY;
        else if(status==B_GAUSSIAN_CONVERT_RNG_ERROR)return INKE_SAMPLE_RNG_ERROR;
        else{
            /* A fallback can still convert this SAME I. Keep a visible
             * diagnostic count; never call this an ordinary sampler retry. */
            S->gaussian_failures++;S->gaussian_failed=1;
        }
    }
    if(!S->has_conversion){
        S->phase=6;if(!O->fallback)return B_KEYGEN_NEEDS_FALLBACK;
        inke_conversion_result result;S->fallback_attempts++;
        status=O->fallback(&result,&S->ideal,O->fallback_context);
        if(status!=INKE_SAMPLE_OK)return status;
        S->conversion=result;S->has_conversion=1;S->conversion_origin=3;
    }
    S->phase=7;S->finish_attempts++;
    status=inke_keygen_from_conversion_with_fallback(sk,pk,&S->q,&S->conversion,
                O->product_fallback,O->product_context,&S->last_stage);
    if(status==INKE_SAMPLE_OK){S->complete=1;S->phase=0;}
    else S->finish_failures++;
    return status;
}
