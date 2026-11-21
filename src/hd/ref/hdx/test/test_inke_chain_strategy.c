// SPDX-License-Identifier: Apache-2.0
#include "inke_chain.h"
#include <inke_conversion.h>
#include <quaternion_data.h>
#include <endomorphism_action.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef INKE_MANUSCRIPT_PARAMETERS
#include <inke_field_constants.h>
#endif

#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c); exit(1); } } while (0)

static int make_kernel(theta_couple_jac_point_t K[2],const ec_basis_t *B,
    const ibz_t *modulus,const inke_mixed_witness *w)
{
    ec_basis_t image,base=*B;
    if (!inke_endomorphism_basis(&image,B,&w->theta,modulus)) return 0;
    ec_curve_t E=CURVE_E0;
    jac_point_t left[2],right[2];
    lift_basis(&left[0],&left[1],&base,&E);
    lift_basis(&right[0],&right[1],&image,&E);
    ibz_t negative; ibz_init(&negative); mpz_neg(negative,w->d1);
    mpz_mod(negative,negative,*modulus);
    digit_t scalar[NWORDS_ORDER]={0}; ibz_to_digits(scalar,&negative);
    for (int j=0;j<2;++j) {
        MUL_generic(&K[j].P1,&left[j],scalar,NWORDS_ORDER,&E);
        K[j].P2=right[j];
    }
    ibz_finalize(&negative); return 1;
}

static void compare_chains(const inke_product_chain *balanced,
    const inke_product_chain *reference,const ec_basis_t *right_two,
    const ec_basis_t *right_three)
{
    CHECK(balanced->prefix_length==reference->prefix_length);
    CHECK(inke_theta_equal(&balanced->splitting.B.null_point,&reference->splitting.B.null_point));
    CHECK(balanced->len2==reference->len2 && balanced->len3==reference->len3);
    const ec_point_t *points[2][6]={
        {&BASIS_EVEN.P,&BASIS_EVEN.Q,&BASIS_EVEN.PmQ,
         &BASIS_THREE.P,&BASIS_THREE.Q,&BASIS_THREE.PmQ},
        {&right_two->P,&right_two->Q,&right_two->PmQ,
         &right_three->P,&right_three->Q,&right_three->PmQ}};
    for (int side=0;side<2;++side) for (int k=0;k<6;++k) {
        theta_couple_point_t b,r;
        CHECK(inke_chain_eval_one(&b,balanced,points[side][k],side));
        CHECK(inke_chain_eval_one(&r,reference,points[side][k],side));
        CHECK(ec_is_equal(&b.P1,&r.P1)); CHECK(ec_is_equal(&b.P2,&r.P2));
    }
    CHECK(balanced->kernel_doublings<=reference->kernel_doublings);
    CHECK(balanced->kernel_triplings<=reference->kernel_triplings);
    printf("strategy equality: 2 length %d, 3 length %d, product prefix %d; kernel doublings %lu/%lu, triplings %lu/%lu\n",
        balanced->len2,balanced->len3,balanced->prefix_length,
        balanced->kernel_doublings,reference->kernel_doublings,
        balanced->kernel_triplings,reference->kernel_triplings); fflush(stdout);
}

static void check_two_power_chains(const inke_conversion_result *conversion,unsigned int a,
    const ibz_t *q)
{
    theta_couple_curve_t E={CURVE_E0,conversion->curve};
    ec_basis_t bases[2]={BASIS_EVEN,conversion->two},right_three;
    ec_curve_to_basis_3(&right_three,&E.E2);
    for (int side=0;side<2;++side) {
        ec_curve_t *curve=side?&E.E2:&E.E1;
        int doublings=TORSION_PLUS_EVEN_POWER-(a+2);
        ec_dbl_iter(&bases[side].P,doublings,curve,&bases[side].P);
        ec_dbl_iter(&bases[side].Q,doublings,curve,&bases[side].Q);
        ec_dbl_iter(&bases[side].PmQ,doublings,curve,&bases[side].PmQ);
    }
    jac_point_t left[2],right[2];
    lift_basis(&left[0],&left[1],&bases[0],&E.E1);
    lift_basis(&right[0],&right[1],&bases[1],&E.E2);
    ibz_t modulus,negative; ibz_init(&modulus); ibz_init(&negative);
    mpz_set_ui(modulus,1); mpz_mul_2exp(modulus,modulus,a+2);
    mpz_neg(negative,*q); mpz_mod(negative,negative,modulus);
    digit_t scalar[NWORDS_ORDER]={0}; ibz_to_digits(scalar,&negative);
    theta_couple_jac_point_t K[2];
    for (int k=0;k<2;++k) {
        MUL_generic(&K[k].P1,&left[k],scalar,NWORDS_ORDER,&E.E1); K[k].P2=right[k];
    }
    inke_product_chain balanced,reference;
    inke_chain_init(&balanced); inke_chain_init(&reference);
    theta_couple_point_t marker;
    CHECK(inke_chain_compute_with_lift(&balanced,&E,K,a,&marker));
    /* Identify each factor independently using the restriction degree
     * modulo four. This pairing is only a test oracle for the marker. */
    ec_basis_t four=bases[0];
    ec_dbl_iter(&four.P,a,&E.E1,&four.P);
    ec_dbl_iter(&four.Q,a,&E.E1,&four.Q);
    ec_dbl_iter(&four.PmQ,a,&E.E1,&four.PmQ);
    theta_couple_point_t im[3];
    CHECK(inke_chain_eval_one(&im[0],&balanced,&bases[0].P,0));
    CHECK(inke_chain_eval_one(&im[1],&balanced,&bases[0].Q,0));
    CHECK(inke_chain_eval_one(&im[2],&balanced,&bases[0].PmQ,0));
    int degree_matches=0;
    for(int side=0;side<2;side++) {
        ec_curve_t *curve=side?&balanced.codomain.E2:&balanced.codomain.E1;
        ec_point_t *point=side?&marker.P2:&marker.P1;
        ec_basis_t image;
        ec_dbl_iter(&image.P,a,curve,side?&im[0].P2:&im[0].P1);
        ec_dbl_iter(&image.Q,a,curve,side?&im[1].P2:&im[1].P1);
        ec_dbl_iter(&image.PmQ,a,curve,side?&im[2].P2:&im[2].P1);
        if(inke_pairing_degree(&four,&E.E1,&image,curve,2,q)) {
            CHECK(test_point_order_twof(point,curve,2));
            degree_matches++;
        } else {
            CHECK(ec_is_zero(point));
        }
    }
    CHECK(degree_matches==1);
    for (int k=0;k<2;++k) double_couple_jac_point_iter(&K[k],2,&E,&K[k]);
    CHECK(inke_chain_compute_exact(&balanced,&E,K,a));
    CHECK(inke_chain_compute_exact_reference(&reference,&E,K,a));
    compare_chains(&balanced,&reference,&conversion->two,&right_three);
    inke_chain_clear(&balanced); inke_chain_clear(&reference);
    ibz_finalize(&negative); ibz_finalize(&modulus);
}

int main(void)
{
    unsigned int n2=TORSION_PLUS_EVEN_POWER-2;
#ifdef INKE_PUBLIC_TWO_EXPONENT
    unsigned int a=INKE_PUBLIC_TWO_EXPONENT;
#else
    unsigned int a=n2;
#endif
    ibz_t q,target; ibz_init(&q); ibz_init(&target);
    mpz_set_ui(q,5); mpz_mul_2exp(target,TORSION_PLUS_3POWER,n2);
    quat_left_ideal_t I; quat_left_ideal_init(&I);
    inke_sample_counts counts={0};
    CHECK(inke_sample_ideal(&I,&q,a,&MAXORD_O0,&QUATALG_PINFTY,
        200000,20000,1024,&counts)==INKE_SAMPLE_OK);
    inke_mixed_state state; inke_mixed_state_init(&state);
    inke_mixed_witness w; inke_mixed_witness_init(&w);
    CHECK(inke_mixed_prepare(&state,&I,&target,&QUATALG_PINFTY)==INKE_SAMPLE_OK);
    inke_product_chain balanced,reference;
    inke_chain_init(&balanced); inke_chain_init(&reference);
    theta_couple_curve_t E={CURVE_E0,CURVE_E0};
    int successes=0,exact_checked=0,prefix_checked=0;
    /* Prefix occurrence is determined by the input and the witness. It is
     * coverage information, not a condition for a fixed ideal to pass. */
    for (int attempt=0;attempt<100 && (successes<2 || !exact_checked);++attempt) {
        int status=inke_mixed_search(&w,&state,100000);
        CHECK(status==INKE_SAMPLE_OK || status==INKE_SAMPLE_RETRY);
        CHECK(quat_lideal_equals(&I,&state.input,&QUATALG_PINFTY));
        if (status==INKE_SAMPLE_RETRY) continue;
        theta_couple_jac_point_t K2[2],K3[2];
        if (!make_kernel(K2,&BASIS_EVEN,&TORSION_PLUS_2POWER,&w) ||
            !make_kernel(K3,&BASIS_THREE,&TORSION_PLUS_3POWER,&w)) continue;
        int b=inke_chain_compute(&balanced,&E,K2,n2,K3,POWER_OF_3);
        int r=inke_chain_compute_reference(&reference,&E,K2,n2,K3,POWER_OF_3);
        CHECK(b==r);
        if (!b) continue;
        compare_chains(&balanced,&reference,&BASIS_EVEN,&BASIS_THREE); ++successes;
        prefix_checked|=balanced.prefix_length>0;
        if (!exact_checked) {
            inke_conversion_result conversion; int stage;
            if (inke_convert_witness(&conversion,&state,&w,&stage)) {
                check_two_power_chains(&conversion,a,&q); exact_checked=1;
            }
        }
    }
    CHECK(successes>=2 && exact_checked);
    inke_chain_clear(&balanced); inke_chain_clear(&reference);
    inke_mixed_witness_clear(&w); inke_mixed_state_clear(&state);
    quat_left_ideal_finalize(&I); ibz_finalize(&target); ibz_finalize(&q);
    printf("Balanced and sequential chains agree, including exact-order tails; product prefix observed: %s.\n",
        prefix_checked?"yes":"no (not required for each input ideal)");
    return 0;
}
