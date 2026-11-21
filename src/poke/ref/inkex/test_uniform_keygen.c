// SPDX-License-Identifier: Apache-2.0
#include <inke_keygen.h>
#include <inke_conversion.h>
#include <quaternion_data.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inke_curve_model.h"
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c); exit(1); } } while(0)

static void check_model(const ec_curve_t *E,const ec_basis_t *B)
{
    inke_model_map maps[6]; unsigned int count=0;
    CHECK(inke_model_candidates(maps,&count,E) && count>0);
    CHECK(fp2_is_equal(&maps[0].curve.A,&E->A) && fp2_is_one(&E->C));
    CHECK(!E->is_A24_computed_and_normalized);
    CHECK(fp2_is_one(&B->P.z) && fp2_is_one(&B->Q.z) && fp2_is_one(&B->PmQ.z));
}

int main(int argc,char **argv)
{
    int count=argc>1?atoi(argv[1]):3;
    CHECK(count>0);
    for (int iteration=0;iteration<count;++iteration) {
        inke_keygen_context ctx; inke_keygen_init(&ctx);
        inke_sk_t sk; inke_pk_t pk; memset(&sk,0xa5,sizeof(sk)); memset(&pk,0xa5,sizeof(pk));
        inke_sk_t original_sk=sk; inke_pk_t original_pk=pk;
        CHECK(inke_keygen_uniform(&sk,&pk,&ctx,0)==INKE_SAMPLE_RETRY);
        CHECK(!ctx.has_q && !ctx.has_ideal);
        CHECK(!memcmp(&sk,&original_sk,sizeof(sk)) && !memcmp(&pk,&original_pk,sizeof(pk)));
        /* Test the actual sampler first, then the two admissible boundaries. */
        if (iteration%3==1) { mpz_set_ui(ctx.q,5); ctx.has_q=1; }
        if (iteration%3==2) {
            mpz_set_ui(ctx.q,1); mpz_mul_2exp(ctx.q,ctx.q,ctx.a); mpz_sub_ui(ctx.q,ctx.q,5); ctx.has_q=1;
        }
        int success=0,have_ideal=0;
        quat_left_ideal_t fixed; quat_left_ideal_init(&fixed); ibz_t fixed_q; ibz_init(&fixed_q);
        for (int attempt=0;attempt<100 && !success;++attempt) {
            int status=inke_keygen_uniform(&sk,&pk,&ctx,100000);
            CHECK(status==INKE_SAMPLE_OK || status==INKE_SAMPLE_RETRY);
            if (have_ideal) {
                CHECK(mpz_cmp(fixed_q,ctx.q)==0);
                CHECK(quat_lideal_equals(&fixed,&ctx.mixed.input,&QUATALG_PINFTY));
            } else if (ctx.has_ideal) {
                ibz_mat_4x4_copy(&fixed.lattice.basis,&ctx.mixed.input.lattice.basis);
                mpz_set(fixed.lattice.denom,ctx.mixed.input.lattice.denom);
                mpz_set(fixed.norm,ctx.mixed.input.norm); fixed.parent_order=ctx.mixed.input.parent_order;
                mpz_set(fixed_q,ctx.q); have_ideal=1;
            }
            success=status==INKE_SAMPLE_OK;
            if (!success) CHECK(!memcmp(&sk,&original_sk,sizeof(sk)) && !memcmp(&pk,&original_pk,sizeof(pk)));
            printf("key %d attempt %d: success=%d stage=%d witnesses=%lu same_ideal=%d\n",
                iteration,attempt,success,ctx.last_stage,ctx.witnesses,have_ideal); fflush(stdout);
        }
        CHECK(success);
        check_model(&pk.EA,&pk.PQ2);
        check_model(&pk.EA,&pk.PQ3);
        check_model(&pk.EA1,&pk.PQA13);
        printf("key %d: product-prefix steps=%d\n",iteration,ctx.product_steps);
        CHECK(inke_keygen_uniform(&sk,&pk,&ctx,100)==INKE_SAMPLE_INVALID);
        CHECK(test_point_order_twof(&pk.PQ2.P,&pk.EA,ctx.a+2));
        CHECK(test_point_order_twof(&pk.PQ2.Q,&pk.EA,ctx.a+2));
        CHECK(test_point_order_threef(&pk.PQ3.P,&pk.EA,POWER_OF_3));
        CHECK(test_point_order_threef(&pk.PQA13.P,&pk.EA1,POWER_OF_3));
        ec_basis_t original=BASIS_EVEN; ec_curve_t E0=CURVE_E0;
        int extra=TORSION_PLUS_EVEN_POWER-(ctx.a+2);
        ec_dbl_iter(&original.P,extra,&E0,&original.P); ec_dbl_iter(&original.Q,extra,&E0,&original.Q);
        ec_dbl_iter(&original.PmQ,extra,&E0,&original.PmQ);
        ibz_t degree,mask; ibz_init(&degree); ibz_init(&mask);
        ibz_copy_digits(&degree,sk.alpha,NWORDS_ORDER); ibz_copy_digits(&mask,sk.beta,NWORDS_ORDER);
        mpz_mul(degree,degree,mask); mpz_mul(degree,degree,ctx.mixed.input.norm);
        CHECK(inke_pairing_degree(&original,&E0,&pk.PQ2,&pk.EA,ctx.a+2,&degree));
        ibz_finalize(&mask); ibz_finalize(&degree);
        for (int trial=0;trial<8;++trial) {
            unsigned char m[32],recovered[32],seed[32];
            for (int k=0;k<32;++k) { m[k]=(unsigned char)(17*iteration+trial+k); seed[k]=(unsigned char)(trial*29+k); }
            inke_ct_t ct,ct_again;
            CHECK(inke_encrypt_uniform(&ct,&pk,m,seed,sizeof(seed)));
            CHECK(inke_encrypt_uniform(&ct_again,&pk,m,seed,sizeof(seed)));
            CHECK(!memcmp(&ct,&ct_again,sizeof(ct)));
            CHECK(inke_decrypt_uniform(recovered,&ct,&sk));
            CHECK(!memcmp(m,recovered,32));
            printf("key %d: encryption/decryption %d verified\n",iteration,trial); fflush(stdout);
        }
        quat_left_ideal_finalize(&fixed); ibz_finalize(&fixed_q); inke_keygen_clear(&ctx);
    }
#ifdef INKE_MANUSCRIPT_PARAMETERS
    puts("Uniform-ideal C Keygen and PKE round trips passed (MANUSCRIPT PARAMETERS; EXPERIMENTAL MIXED CONVERSION).");
#else
    puts("Uniform-ideal C Keygen and PKE round trips passed (UPSTREAM TEST FIELDS; NOT MANUSCRIPT PARAMETERS).");
#endif
    return 0;
}
