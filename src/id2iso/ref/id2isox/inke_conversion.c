// SPDX-License-Identifier: Apache-2.0
#include <inke_conversion.h>
#include <id2iso.h>
#include <endomorphism_action.h>
#include <quaternion_data.h>
#include <biextension.h>
#include "inke_chain.h"

static int action_matrix(ibz_mat_2x2_t *M,const quat_alg_elem_t *alpha,const ibz_t *modulus)
{
    ibz_vec_4_t c; ibz_vec_4_init(&c);
    ibz_t t,norm,det; ibz_init(&t); ibz_init(&norm); ibz_init(&det);
    ibq_t nq; ibq_init(&nq);
    int ok=quat_lattice_contains(&c,&MAXORD_O0,alpha,&QUATALG_PINFTY);
    if (!ok) goto done;
    const ibz_mat_2x2_t *actions[3]={&ACTION_GEN2,&ACTION_GEN3,&ACTION_GEN4};
    for (int i=0;i<2;++i) for (int j=0;j<2;++j) {
        mpz_set_ui((*M)[i][j],0);
        if (i==j) mpz_set((*M)[i][j],c[0]);
        for (int k=0;k<3;++k) mpz_addmul((*M)[i][j],c[k+1],(*actions[k])[i][j]);
        mpz_mod((*M)[i][j],(*M)[i][j],*modulus);
    }
    mpz_mul(det,(*M)[0][0],(*M)[1][1]); mpz_submul(det,(*M)[0][1],(*M)[1][0]);
    quat_alg_norm(&nq,alpha,&QUATALG_PINFTY);
    ok=ibq_to_ibz(&norm,&nq);
    mpz_sub(t,det,norm); ok=ok && mpz_divisible_p(t,*modulus);
done:
    ibq_finalize(&nq); ibz_finalize(&det); ibz_finalize(&norm); ibz_finalize(&t);
    ibz_vec_4_finalize(&c); return ok;
}
static void scalar_jac(jac_point_t *out,const jac_point_t *P,const ibz_t *k,
    const ibz_t *modulus,const ec_curve_t *E)
{
    ibz_t t; ibz_init(&t); mpz_mod(t,*k,*modulus);
    digit_t s[NWORDS_ORDER]={0}; ibz_to_digits(s,&t);
    MUL_generic(out,P,s,NWORDS_ORDER,E); ibz_finalize(&t);
}
static void matrix_jac(jac_point_t out[2],const jac_point_t base[2],
    const ibz_mat_2x2_t *M,const ec_curve_t *E)
{
    for (int j=0;j<2;++j) {
        digit_t s[2][NWORDS_ORDER]={{0}};
        ibz_to_digits(s[0],&(*M)[0][j]); ibz_to_digits(s[1],&(*M)[1][j]);
        DBLMUL_generic(&out[j],&base[0],s[0],&base[1],s[1],E,NWORDS_ORDER);
    }
}
int inke_endomorphism_basis(ec_basis_t *out,const ec_basis_t *B,
    const quat_alg_elem_t *alpha,const ibz_t *modulus)
{
    ibz_mat_2x2_t M; ibz_mat_2x2_init(&M);
    int ok=action_matrix(&M,alpha,modulus);
    if (ok) {
        ec_curve_t E=CURVE_E0; ec_basis_t basis=*B;
        jac_point_t base[2],image[2],negative,difference;
        lift_basis(&base[0],&base[1],&basis,&E); matrix_jac(image,base,&M,&E);
        negative=image[1]; fp2_neg(&negative.y,&negative.y);
        ADD(&difference,&image[0],&negative,&E);
        jac_to_xz(&out->P,&image[0]); jac_to_xz(&out->Q,&image[1]);
        jac_to_xz(&out->PmQ,&difference);
    }
    ibz_mat_2x2_finalize(&M); return ok;
}
static void pairing(fp2_t *r,ec_basis_t basis,const ec_curve_t *E,int f)
{
    ec_point_t ac={E->A,E->C},a24;
    A24_from_AC(&a24,&ac); weil(r,f,&basis.P,&basis.Q,&basis.PmQ,&a24);
}
int inke_pairing_degree(const ec_basis_t *B,const ec_curve_t *E,
    const ec_basis_t *image,const ec_curve_t *F,int f,const ibz_t *degree)
{
    const ec_point_t *pts[3]={&image->P,&image->Q,&image->PmQ};
    for (int i=0;i<3;++i) if (!test_point_order_twof(pts[i],F,f)) return 0;
    fp2_t source,target,power;
    pairing(&source,*B,E,f); pairing(&target,*image,F,f);
    ibz_t modulus,d; ibz_init(&modulus); ibz_init(&d);
    mpz_set_ui(modulus,1); mpz_mul_2exp(modulus,modulus,f); mpz_mod(d,*degree,modulus);
    digit_t exponent[NWORDS_ORDER]={0}; ibz_to_digits(exponent,&d);
    fp2_pow_vartime(&power,&source,exponent,NWORDS_ORDER);
    int ok=fp2_is_equal(&power,&target)!=0;
    /* A primitive pairing is also needed, not only the degree identity. */
    for (int i=1;i<f;++i) fp2_sqr(&source,&source);
    fp2_t one; fp2_set_one(&one);
    ok=ok && !fp2_is_equal(&source,&one);
    ibz_finalize(&d); ibz_finalize(&modulus); return ok;
}
static int kernel(theta_couple_jac_point_t K[2],const ec_basis_t *B,
    const ibz_t *modulus,const inke_mixed_witness *w)
{
    ibz_mat_2x2_t M; ibz_mat_2x2_init(&M);
    int ok=action_matrix(&M,&w->theta,modulus);
    if (ok) {
        jac_point_t base[2],right[2]; ec_basis_t basis=*B; ec_curve_t E=CURVE_E0;
        lift_basis(&base[0],&base[1],&basis,&E); matrix_jac(right,base,&M,&E);
        ibz_t negative; ibz_init(&negative); mpz_neg(negative,w->d1);
        for (int j=0;j<2;++j) {
            scalar_jac(&K[j].P1,&base[j],&negative,modulus,&E); K[j].P2=right[j];
        }
        ibz_finalize(&negative);
    }
    ibz_mat_2x2_finalize(&M); return ok;
}
int inke_convert_witness(inke_conversion_result *out,const inke_mixed_state *s,
    const inke_mixed_witness *w,int *stage)
{
    *stage=10;
    if (mpz_cmp(s->algebra->p,QUATALG_PINFTY.p)!=0 || !inke_mixed_verify(s,w)) return 0;
    int n2=(int)TORSION_PLUS_EVEN_POWER-2,n3=POWER_OF_3;
    ibz_t target; ibz_init(&target);
    mpz_mul_2exp(target,TORSION_PLUS_3POWER,n2);
    int ok=mpz_cmp(target,s->target)==0; ibz_finalize(&target);
    if (!ok) return 0;
    theta_couple_jac_point_t K2[2],K3[2];
    *stage=11;
    if (!kernel(K2,&BASIS_EVEN,&TORSION_PLUS_2POWER,w) ||
        !kernel(K3,&BASIS_THREE,&TORSION_PLUS_3POWER,w)) return 0;
    inke_product_chain chain; inke_chain_init(&chain);
    theta_couple_curve_t E={CURVE_E0,CURVE_E0};
    ok=inke_chain_compute(&chain,&E,K2,n2,K3,n3);
    if (!ok) { *stage=20+chain.failed_stage; goto done; }
    ec_basis_t four=BASIS_EVEN, images[2];
    const ec_point_t *full[3]={&BASIS_EVEN.P,&BASIS_EVEN.Q,&BASIS_EVEN.PmQ};
    ec_point_t *input[3]={&four.P,&four.Q,&four.PmQ};
    ec_point_t *im[2][3]={{&images[0].P,&images[0].Q,&images[0].PmQ},
                        {&images[1].P,&images[1].Q,&images[1].PmQ}};
    *stage=30; ok=0;
    for (int k=0;k<3;++k) {
        ec_dbl_iter(input[k],n2,&E.E1,input[k]);
        theta_couple_point_t pair;
        /* Evaluate high-order points before reducing to 4-torsion. This
         * avoids a gluing coordinate chart vanishing on a low-order point. */
        if (!inke_chain_eval_one(&pair,&chain,full[k],0)) goto done;
        *im[0][k]=pair.P1; *im[1][k]=pair.P2;
        ec_dbl_iter(im[0][k],n2,&chain.codomain.E1,im[0][k]);
        ec_dbl_iter(im[1][k],n2,&chain.codomain.E2,im[1][k]);
    }
    int first=inke_pairing_degree(&four,&E.E1,&images[0],&chain.codomain.E1,2,&w->d1);
    int second=inke_pairing_degree(&four,&E.E1,&images[1],&chain.codomain.E2,2,&w->d1);
    *stage=31;
    if (first==second) goto done;
    out->curve=first?chain.codomain.E1:chain.codomain.E2;
    ec_basis_t mapped;
    *stage=32;
    if (!inke_endomorphism_basis(&mapped,&BASIS_EVEN,&w->transport,&TORSION_PLUS_2POWER)) goto done;
    ibz_t inverse; ibz_init(&inverse);
    mpz_invert(inverse,w->d1,TORSION_PLUS_2POWER);
    ec_point_t *transported[3]={&mapped.P,&mapped.Q,&mapped.PmQ};
    ec_point_t *result[3]={&out->two.P,&out->two.Q,&out->two.PmQ};
    for (int k=0;k<3;++k) {
        theta_couple_point_t pair;
        /* Precompose tau/d1 before an even-degree product chain; this
         * retains the original full 2-primary action. */
        ec_mul_ibz(transported[k],&E.E1,&inverse,transported[k]);
        if (!inke_chain_eval_one(&pair,&chain,transported[k],0)) {
            ibz_finalize(&inverse); goto done;
        }
        *result[k]=first?pair.P1:pair.P2;
    }
    ibz_finalize(&inverse); *stage=33;
    ok=inke_pairing_degree(&BASIS_EVEN,&E.E1,&out->two,&out->curve,
                            TORSION_PLUS_EVEN_POWER,&s->input.norm);
    if (ok) { *stage=0; out->product_steps=chain.prefix_length; }
done:
    inke_chain_clear(&chain); return ok;
}
