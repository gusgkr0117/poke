// SPDX-License-Identifier: Apache-2.0
#include <inke_mixed.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c); exit(1); \
} } while (0)

static void init_order(quat_order_t *order)
{
    quat_lattice_init(order); mpz_set_ui(order->denom,2);
    mpz_set_ui(order->basis[0][0],2); mpz_set_ui(order->basis[1][1],2);
    mpz_set_ui(order->basis[1][2],1); mpz_set_ui(order->basis[2][2],1);
    mpz_set_ui(order->basis[0][3],1); mpz_set_ui(order->basis[3][3],1);
}

int main(int argc, char **argv)
{
    int full = argc==2 && strcmp(argv[1],"--full")==0;
    CHECK(argc==1 || full);
    const unsigned int a[]={256,384,512}, b[]={162,243,324};
    const unsigned int f2[]={257,386,512}, f3[]={163,244,324}, c[]={11,23,517};
    ibz_t p, power, q, target, saved_q;
    ibz_init(&p); ibz_init(&power); ibz_init(&q); ibz_init(&target); ibz_init(&saved_q);
    unsigned int limited=0, passed=0;
    for (int level=0;level<3;++level) {
        mpz_ui_pow_ui(p,3,f3[level]); mpz_mul_2exp(p,p,f2[level]);
        mpz_mul_ui(p,p,c[level]); mpz_sub_ui(p,p,1);
        quat_alg_t alg; quat_alg_init_set(&alg,&p);
        quat_order_t order; init_order(&order);
        mpz_set_ui(power,1); mpz_mul_2exp(power,power,a[level]);
        mpz_ui_pow_ui(target,3,b[level]); mpz_mul_2exp(target,target,f2[level]-2);
        for (int sample=0;sample<(full?3:1);++sample) {
            if (full && sample==0) mpz_set_ui(q,5);
            else if (full && sample==1) mpz_sub_ui(q,power,5);
            else CHECK(inke_sample_q(&q,a[level])==INKE_SAMPLE_OK);
            mpz_set(saved_q,q);
            quat_left_ideal_t input; quat_left_ideal_init(&input);
            inke_sample_counts counts={0};
            CHECK(inke_sample_ideal(&input,&q,a[level],&order,&alg,
                                    200000,20000,1024,&counts)==INKE_SAMPLE_OK);
            inke_mixed_state state; inke_mixed_state_init(&state);
            inke_mixed_witness witness; inke_mixed_witness_init(&witness);
            CHECK(inke_mixed_prepare(&state,&input,&target,&alg)==INKE_SAMPLE_OK);
            CHECK(inke_mixed_search(&witness,&state,0)==INKE_SAMPLE_RETRY);
            CHECK(state.attempts==0 && quat_lideal_equals(&input,&state.input,&alg));
            int status=inke_mixed_search(&witness,&state,1);
            CHECK(status==INKE_SAMPLE_OK || status==INKE_SAMPLE_RETRY);
            CHECK(state.attempts==1 && quat_lideal_equals(&input,&state.input,&alg));
            /* Force a small nonzero budget before resuming, exercising a
             * retry after actual work as well as the zero-budget path. */
            /* Explicit bounded windows; every retry retains the original I. */
            for (int window=0;window<10 && status==INKE_SAMPLE_RETRY;++window) {
                status=inke_mixed_search(&witness,&state,100000);
                CHECK(status==INKE_SAMPLE_OK || status==INKE_SAMPLE_RETRY);
                CHECK(quat_lideal_equals(&input,&state.input,&alg));
                CHECK(mpz_cmp(q,saved_q)==0);
                printf("a=%u sample=%d window=%d status=%s attempts=%lu generators=%lu positive=%lu squares=%lu\n",
                    a[level],sample,window,status==INKE_SAMPLE_OK?"verified":"limit",
                    state.attempts,state.generators,state.positive_targets,state.square_solutions);
                fflush(stdout);
            }
            if (status==INKE_SAMPLE_OK) {
                CHECK(inke_mixed_verify(&state,&witness));
                /* Negative tests: accepting a norm equation alone is insufficient. */
                mpz_add_ui(witness.transport.coord[0],witness.transport.coord[0],1);
                CHECK(!inke_mixed_verify(&state,&witness));
                mpz_sub_ui(witness.transport.coord[0],witness.transport.coord[0],1);
                mpz_add_ui(witness.d1,witness.d1,1);
                CHECK(!inke_mixed_verify(&state,&witness));
                mpz_sub_ui(witness.d1,witness.d1,1);
                mpz_add_ui(witness.mu1.coord[2],witness.mu1.coord[2],1);
                CHECK(!inke_mixed_verify(&state,&witness));
                mpz_sub_ui(witness.mu1.coord[2],witness.mu1.coord[2],1);
                CHECK(inke_mixed_verify(&state,&witness));
                ++passed;
            } else {
                /* A search limit is a reported limitation, never a success. */
                ++limited;
            }
            inke_mixed_witness_clear(&witness); inke_mixed_state_clear(&state);
            quat_left_ideal_finalize(&input);
        }
        quat_lattice_finalize(&order); quat_alg_finalize(&alg);
    }
    ibz_finalize(&saved_q); ibz_finalize(&target); ibz_finalize(&q);
    ibz_finalize(&power); ibz_finalize(&p);
    printf("C mixed norm stage: %u verified, %u search limits; no ideal resampling.\n",passed,limited);
    puts("This is not an isogeny evaluation or a revised Keygen benchmark.");
    return limited?2:0;
}
