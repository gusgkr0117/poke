// SPDX-License-Identifier: Apache-2.0
#include <inke_conversion.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c); exit(1); } } while(0)

int main(void)
{
    unsigned int a=(unsigned int)TORSION_PLUS_EVEN_POWER-2;
    ibz_t q,target; ibz_init(&q); ibz_init(&target);
    mpz_set_ui(q,5); mpz_mul_2exp(target,TORSION_PLUS_3POWER,a);
    quat_left_ideal_t I; quat_left_ideal_init(&I);
    inke_sample_counts counts={0};
    CHECK(inke_sample_ideal(&I,&q,a,&MAXORD_O0,&QUATALG_PINFTY,
        200000,20000,1024,&counts)==INKE_SAMPLE_OK);
    inke_mixed_state state; inke_mixed_state_init(&state);
    inke_mixed_witness w; inke_mixed_witness_init(&w);
    CHECK(inke_mixed_prepare(&state,&I,&target,&QUATALG_PINFTY)==INKE_SAMPLE_OK);
    int success=0;
    for (int attempt=0;attempt<20 && !success;++attempt) {
        int status=inke_mixed_search(&w,&state,100000);
        CHECK(status==INKE_SAMPLE_OK || status==INKE_SAMPLE_RETRY);
        if (status==INKE_SAMPLE_RETRY) { puts("norm budget exhausted; retaining I"); fflush(stdout); continue; }
        int stage=0; inke_conversion_result result;
        printf("witness %d: evaluating mixed chain\n",attempt); fflush(stdout);
        success=inke_convert_witness(&result,&state,&w,&stage);
        printf("witness %d: conversion=%d stage=%d; same ideal=%d\n",attempt,success,stage,
            quat_lideal_equals(&I,&state.input,&QUATALG_PINFTY)); fflush(stdout);
        CHECK(quat_lideal_equals(&I,&state.input,&QUATALG_PINFTY));
    }
    inke_mixed_witness_clear(&w); inke_mixed_state_clear(&state);
    quat_left_ideal_finalize(&I); ibz_finalize(&target); ibz_finalize(&q);
    CHECK(success);
    puts("Fixed-ideal conversion and full 2-primary pairing verified (upstream test field).");
    return 0;
}
