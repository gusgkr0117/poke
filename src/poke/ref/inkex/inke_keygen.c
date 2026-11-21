// SPDX-License-Identifier: Apache-2.0
#include <inke_keygen.h>
#include <inke_conversion.h>
#include <id2iso.h>
#include <quaternion_data.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <biextension.h>
#include <fips202.h>
#include <encoded_sizes.h>
#include "inke_chain.h"
#include "inke_curve_model.h"
#include <string.h>
#ifdef INKE_MANUSCRIPT_PARAMETERS
#include <inke_field_constants.h>
#endif

static unsigned int public_exponent(void)
{
#ifdef INKE_PUBLIC_TWO_EXPONENT
    return INKE_PUBLIC_TWO_EXPONENT;
#else
    return TORSION_PLUS_EVEN_POWER-2;
#endif
}

static void modulus_two(ibz_t *m,unsigned int a)
{ mpz_set_ui(*m,1); mpz_mul_2exp(*m,*m,a); }
static void restrict_two(ec_basis_t *B,const ec_curve_t *E,unsigned int a)
{
    ec_curve_t curve=*E;
    ec_dbl_iter(&B->P,TORSION_PLUS_EVEN_POWER-a,&curve,&B->P);
    ec_dbl_iter(&B->Q,TORSION_PLUS_EVEN_POWER-a,&curve,&B->Q);
    ec_dbl_iter(&B->PmQ,TORSION_PLUS_EVEN_POWER-a,&curve,&B->PmQ);
}
static void pair_to_basis(ec_basis_t *B,const jac_point_t P[2],const ec_curve_t *E)
{
    jac_point_t negative=P[1],difference;
    fp2_neg(&negative.y,&negative.y); ADD(&difference,&P[0],&negative,E);
    jac_to_xz(&B->P,&P[0]); jac_to_xz(&B->Q,&P[1]); jac_to_xz(&B->PmQ,&difference);
}
static void scale_basis(ec_basis_t *out,const ec_basis_t *B,const ec_curve_t *E,
    const ibz_t *x,const ibz_t *y)
{
    ec_basis_t basis=*B; ec_curve_t curve=*E;
    jac_point_t P[2]; digit_t s[NWORDS_ORDER]={0};
    lift_basis(&P[0],&P[1],&basis,&curve);
    ibz_to_digits(s,x); MUL_generic(&P[0],&P[0],s,NWORDS_ORDER,&curve);
    memset(s,0,sizeof(s)); ibz_to_digits(s,y); MUL_generic(&P[1],&P[1],s,NWORDS_ORDER,&curve);
    pair_to_basis(out,P,&curve);
}
static void combine_basis(ec_basis_t *out,const ec_basis_t *B,const ec_curve_t *E,
    digit_t s[4][NWORDS_ORDER])
{
    ec_basis_t basis=*B; ec_curve_t curve=*E;
    jac_point_t P[2],R[2]; lift_basis(&P[0],&P[1],&basis,&curve);
    DBLMUL_generic(&R[0],&P[0],s[0],&P[1],s[1],&curve,NWORDS_ORDER);
    DBLMUL_generic(&R[1],&P[0],s[2],&P[1],s[3],&curve,NWORDS_ORDER);
    pair_to_basis(out,R,&curve);
}
static int equal_basis(const ec_basis_t *B,const ec_basis_t *C)
{ return ec_is_equal(&B->P,&C->P) && ec_is_equal(&B->Q,&C->Q) && ec_is_equal(&B->PmQ,&C->PmQ); }
static int eval_basis(ec_basis_t *out,const inke_product_chain *chain,
    const ec_basis_t *B,int side,int factor)
{
    const ec_point_t *p[3]={&B->P,&B->Q,&B->PmQ};
    ec_point_t *q[3]={&out->P,&out->Q,&out->PmQ};
    for (int k=0;k<3;++k) {
        theta_couple_point_t im;
        if (!inke_chain_eval_one(&im,chain,p[k],side)) return 0;
        *q[k]=factor?im.P2:im.P1;
    }
    return 1;
}
static void graph_kernel(theta_couple_jac_point_t K[2],const ec_basis_t *left,
    const ec_basis_t *right,const theta_couple_curve_t *E,const ibz_t *q,
    const ibz_t *modulus)
{
    ec_basis_t L=*left,R=*right; ec_curve_t E1=E->E1,E2=E->E2;
    jac_point_t P[2],Q[2]; lift_basis(&P[0],&P[1],&L,&E1); lift_basis(&Q[0],&Q[1],&R,&E2);
    ibz_t negative; ibz_init(&negative); mpz_neg(negative,*q); mpz_mod(negative,negative,*modulus);
    digit_t s[NWORDS_ORDER]={0}; ibz_to_digits(s,&negative);
    for (int k=0;k<2;++k) { MUL_generic(&K[k].P1,&P[k],s,NWORDS_ORDER,&E1); K[k].P2=Q[k]; }
    ibz_finalize(&negative);
}
/* The lift (-qP,phi(P)) maps to (-2^a phi_1(P),0) in labelled
 * coordinates. The splitter can return the two factors in either order. */
static int marked_factor(const theta_couple_point_t *image,const theta_couple_curve_t *curves)
{
    if (ec_is_zero(&image->P2) && test_point_order_twof(&image->P1,&curves->E1,2)) return 0;
    if (ec_is_zero(&image->P1) && test_point_order_twof(&image->P2,&curves->E2,2)) return 1;
    return -1;
}
static int random_residue(ibz_t *v,const ibz_t *modulus,shake256ctx *rng)
{
    ibz_t upper; ibz_init(&upper); mpz_sub_ui(upper,*modulus,1);
    int ok=rng?ibz_rand_interval_with_state(v,&ibz_const_zero,&upper,rng):
               ibz_rand_interval(v,&ibz_const_zero,&upper);
    ibz_finalize(&upper); return ok;
}
static int random_unit(ibz_t *v,const ibz_t *modulus,shake256ctx *rng)
{
    ibz_t g; ibz_init(&g); int ok;
    do { ok=random_residue(v,modulus,rng); if (!ok) break; mpz_gcd(g,*v,*modulus); }
    while (mpz_cmp_ui(g,1));
    ibz_finalize(&g); return ok;
}
/* Return 1 on success, -2 on RNG failure, -3 for an invalid model.
 * Each listed x-map has two full lifts differing by simultaneous negation,
 * which has identical (x(P),x(Q),x(P-Q)) data. Uniform selection therefore
 * realizes the x-only projection of the paper's uniform isomorphism. */
static int canonicalize_receiver(ec_curve_t *E,ec_basis_t **bases,size_t count)
{
    inke_model_map maps[6]; unsigned int choices=0;
    if (!inke_model_candidates(maps,&choices,E) || !choices) return -3;
    ibz_t modulus,index; ibz_init(&modulus); ibz_init(&index);
    mpz_set_ui(modulus,choices);
    int ok=random_residue(&index,&modulus,NULL);
    unsigned long selected=ok?mpz_get_ui(index):0;
    ibz_finalize(&index); ibz_finalize(&modulus);
    if (!ok) return -2;
    for (size_t k=0;k<count;++k)
        if (!inke_model_apply_basis(bases[k],bases[k],&maps[selected])) return -3;
    *E=maps[selected].curve;
    return 1;
}
static int full_three(const ec_basis_t *B,const ec_curve_t *E)
{
    if (!test_point_order_threef(&B->P,E,POWER_OF_3) ||
        !test_point_order_threef(&B->Q,E,POWER_OF_3) ||
        !test_point_order_threef(&B->PmQ,E,POWER_OF_3)) return 0;
    ec_point_t P=B->P,Q=B->Q;
    ibz_t cofactor; ibz_init(&cofactor); mpz_divexact_ui(cofactor,TORSION_PLUS_3POWER,3);
    ec_mul_ibz(&P,E,&cofactor,&P); ec_mul_ibz(&Q,E,&cofactor,&Q); ibz_finalize(&cofactor);
    return !ec_is_equal(&P,&Q);
}
static int finish_keygen(inke_sk_t *sk,inke_pk_t *pk,inke_keygen_context *ctx,
    const inke_conversion_result *conversion,inke_product_fallback fallback,
    void *fallback_context)
{
    int ok=0; unsigned int a=ctx->a; inke_pk_t public_key={0}; inke_sk_t secret_key={0};
    ibz_t modulus,other,alpha,beta,gamma,gamma1;
    ibz_init(&modulus); ibz_init(&other); ibz_init(&alpha); ibz_init(&beta); ibz_init(&gamma); ibz_init(&gamma1);
    modulus_two(&other,a); mpz_sub(other,other,ctx->q);
    modulus_two(&modulus,a+2);
    theta_couple_curve_t E={CURVE_E0,conversion->curve};
    ec_basis_t source=BASIS_EVEN,image=conversion->two;
    restrict_two(&source,&E.E1,a+2); restrict_two(&image,&E.E2,a+2);
    inke_product_chain chain; inke_chain_init(&chain); theta_couple_jac_point_t K[2];
    graph_kernel(K,&source,&image,&E,&ctx->q,&modulus);
    public_key.EA=conversion->curve;
    ec_basis_t R,dualR,phi1,reconstructed;
    ec_curve_to_basis_3(&R,&public_key.EA);
    digit_t s[4][NWORDS_ORDER]={{0}};
    int factor=-1,geometry_ok=0;
    ctx->last_stage=40;
#ifndef INKE_TEST_FORCE_PRODUCT_FALLBACK
    theta_couple_point_t marker;
    if (inke_chain_compute_with_lift(&chain,&E,K,a,&marker)) {
        ctx->last_stage=45; factor=marked_factor(&marker,&chain.codomain);
        if (factor>=0) {
            public_key.EA1=factor?chain.codomain.E2:chain.codomain.E1;
            ctx->last_stage=46;
            geometry_ok=eval_basis(&phi1,&chain,&BASIS_THREE,0,factor) &&
                        eval_basis(&dualR,&chain,&R,1,factor);
        }
    } else ctx->last_stage+=chain.failed_stage;
#endif
    /* A failed deterministic fast computation is not cured by resampling
     * q, the original ideal, or masks. The optional implementation receives
     * the very same graph and every point needed by the finishing step. */
    for (int attempt=0;attempt<2;++attempt) {
        if (attempt) {
            if (!fallback) goto done;
            int detail=0;
            int status=fallback(&public_key.EA1,&phi1,&dualR,&E.E1,&E.E2,
                &source,&image,a,&ctx->q,&BASIS_THREE,&R,fallback_context,&detail);
            ctx->last_stage=status==INKE_SAMPLE_RNG_ERROR?-2:
                            status==INKE_SAMPLE_INVALID?-3:500+detail;
            if (status!=INKE_SAMPLE_OK) goto done;
            geometry_ok=1;
        }
        if (!geometry_ok || !full_three(&phi1,&public_key.EA1) ||
            !full_three(&dualR,&public_key.EA1)) { geometry_ok=0; continue; }
        memset(s,0,sizeof s);
        ec_dlog_weil_3(s[0],s[1],s[2],s[3],&dualR,&phi1,&public_key.EA1);
        combine_basis(&reconstructed,&dualR,&public_key.EA1,s);
        ctx->last_stage=47;
        if (equal_basis(&reconstructed,&phi1)) { geometry_ok=1; break; }
        geometry_ok=0;
    }
    if (!geometry_ok) goto done;
    combine_basis(&public_key.PQ3,&R,&public_key.EA,s);
    /* phi1 = -dual(phi2) o u, so phiI = [-(2^a-q)] u on odd torsion. */
    mpz_neg(other,other); mpz_mod(other,other,TORSION_PLUS_3POWER);
    scale_basis(&public_key.PQ3,&public_key.PQ3,&public_key.EA,&other,&other);
    ctx->last_stage=48;
    if (!full_three(&public_key.PQ3,&public_key.EA)) goto done;
    if (!random_unit(&alpha,&modulus,NULL) || !random_unit(&beta,&modulus,NULL) ||
        !random_unit(&gamma,&TORSION_PLUS_3POWER,NULL) || !random_unit(&gamma1,&TORSION_PLUS_3POWER,NULL)) {
        ctx->last_stage=-2; goto done;
    }
    scale_basis(&public_key.PQ2,&image,&public_key.EA,&alpha,&beta);
    scale_basis(&public_key.PQ3,&public_key.PQ3,&public_key.EA,&gamma,&gamma);
    scale_basis(&public_key.PQA13,&phi1,&public_key.EA1,&gamma1,&gamma1);
    ec_basis_t *on_EA[2]={&public_key.PQ2,&public_key.PQ3};
    ec_basis_t *on_EA1[1]={&public_key.PQA13};
    int model_status=canonicalize_receiver(&public_key.EA,on_EA,2);
    if (model_status!=1) { ctx->last_stage=model_status; goto done; }
    model_status=canonicalize_receiver(&public_key.EA1,on_EA1,1);
    if (model_status!=1) { ctx->last_stage=model_status; goto done; }
    ibz_to_digits(secret_key.deg,&ctx->q); ibz_to_digits(secret_key.alpha,&alpha); ibz_to_digits(secret_key.beta,&beta);
    *sk=secret_key; *pk=public_key; ok=1; ctx->last_stage=0;
done:
    inke_chain_clear(&chain);
    ibz_finalize(&gamma1); ibz_finalize(&gamma); ibz_finalize(&beta); ibz_finalize(&alpha);
    ibz_finalize(&other); ibz_finalize(&modulus); return ok;
}
void inke_keygen_init(inke_keygen_context *ctx)
{
    memset(ctx,0,sizeof(*ctx)); ibz_init(&ctx->q); inke_mixed_state_init(&ctx->mixed);
    ctx->a=public_exponent();
}
void inke_keygen_clear(inke_keygen_context *ctx)
{ inke_mixed_state_clear(&ctx->mixed); ibz_finalize(&ctx->q); memset(ctx,0,sizeof(*ctx)); }
int inke_keygen_from_conversion_with_fallback(inke_sk_t *sk,inke_pk_t *pk,const ibz_t *q,
    const inke_conversion_result *conversion,inke_product_fallback fallback,
    void *fallback_context,int *stage)
{
    if(stage)*stage=-1;
    inke_keygen_context ctx;inke_keygen_init(&ctx);
    ibz_t limit,degree;ibz_init(&limit);ibz_init(&degree);modulus_two(&limit,ctx.a);
    int status=INKE_SAMPLE_INVALID;
    if(!inke_q_is_admissible(q,ctx.a)||
       fp2_is_zero(&conversion->curve.C))goto done;
    fp2_t discriminant,t;fp2_sqr(&discriminant,&conversion->curve.A);
    fp2_sqr(&t,&conversion->curve.C);fp2_add(&t,&t,&t);fp2_add(&t,&t,&t);
    fp2_sub(&discriminant,&discriminant,&t);if(fp2_is_zero(&discriminant))goto done;
    mpz_sub(degree,limit,*q);mpz_mul(degree,degree,*q);
    if(!inke_pairing_degree(&BASIS_EVEN,&CURVE_E0,&conversion->two,
                           &conversion->curve,TORSION_PLUS_EVEN_POWER,&degree))goto done;
    mpz_set(ctx.q,*q);
    if(finish_keygen(sk,pk,&ctx,conversion,fallback,fallback_context))status=INKE_SAMPLE_OK;
    else status=ctx.last_stage==-2?INKE_SAMPLE_RNG_ERROR:
                ctx.last_stage==-3?INKE_SAMPLE_INVALID:INKE_SAMPLE_RETRY;
    if(stage)*stage=ctx.last_stage;
done:
    ibz_finalize(&degree);ibz_finalize(&limit);inke_keygen_clear(&ctx);return status;
}
int inke_keygen_from_conversion(inke_sk_t *sk,inke_pk_t *pk,const ibz_t *q,
    const inke_conversion_result *conversion,int *stage)
{
    return inke_keygen_from_conversion_with_fallback(sk,pk,q,conversion,NULL,NULL,stage);
}
int inke_keygen_uniform(inke_sk_t *sk,inke_pk_t *pk,inke_keygen_context *ctx,unsigned long budget)
{
    if (ctx->complete || ctx->a!=public_exponent()) return INKE_SAMPLE_INVALID;
    if (!budget) return INKE_SAMPLE_RETRY;
    int status;
    if (!ctx->has_q) {
        status=inke_sample_q(&ctx->q,ctx->a);
        if (status!=INKE_SAMPLE_OK) return status;
        ctx->has_q=1;
    }
    if (!ctx->has_ideal) {
        quat_left_ideal_t I; quat_left_ideal_init(&I);
        status=inke_sample_ideal(&I,&ctx->q,ctx->a,&MAXORD_O0,&QUATALG_PINFTY,200000,20000,1024,&ctx->sampling);
        if (status==INKE_SAMPLE_OK) {
            ibz_t target; ibz_init(&target);
            mpz_mul_2exp(target,TORSION_PLUS_3POWER,TORSION_PLUS_EVEN_POWER-2);
            status=inke_mixed_prepare(&ctx->mixed,&I,&target,&QUATALG_PINFTY);
            ibz_finalize(&target); if (status==INKE_SAMPLE_OK) ctx->has_ideal=1;
        }
        quat_left_ideal_finalize(&I); if (status!=INKE_SAMPLE_OK) return status;
    }
    inke_mixed_witness w; inke_mixed_witness_init(&w);
    status=inke_mixed_search(&w,&ctx->mixed,budget);
    if (status==INKE_SAMPLE_OK) {
        ++ctx->witnesses; inke_conversion_result conversion;
        if (!inke_convert_witness(&conversion,&ctx->mixed,&w,&ctx->last_stage)) {
            ++ctx->conversion_failures; status=INKE_SAMPLE_RETRY;
        } else if (!finish_keygen(sk,pk,ctx,&conversion,NULL,NULL)) {
            ++ctx->kani_failures;
            status=ctx->last_stage==-2?INKE_SAMPLE_RNG_ERROR:
                   ctx->last_stage==-3?INKE_SAMPLE_INVALID:INKE_SAMPLE_RETRY;
        } else { ctx->complete=1; ctx->product_steps=conversion.product_steps; }
    }
    inke_mixed_witness_clear(&w); return status;
}

static void curve_key(unsigned char key[32],const ec_curve_t *E)
{
    fp2_t j; ec_j_inv(&j,E);
    unsigned char encoded[2*FP_ENCODED_BYTES];
    fp2_encode(encoded,&j);
    SHAKE256(key,32,encoded,sizeof(encoded));
}
static void odd_isogeny(ec_isog_odd_t *phi,const ec_curve_t *E,const ec_basis_t *B,
    digit_t r[NWORDS_ORDER],digit_t s[NWORDS_ORDER])
{
    memset(phi,0,sizeof(*phi)); phi->curve=*E; phi->degree[0]=POWER_OF_3;
    ec_set_zero(&phi->ker_minus);
    ec_biscalar_mul_bounded(&phi->ker_plus,E,r,s,B,TORSION_3POWER_BYTES*8);
}
int inke_encrypt_uniform(inke_ct_t *ct,const inke_pk_t *pk,const unsigned char m[32],
    const unsigned char *seed,size_t seed_len)
{
    int ok=0; unsigned int a=public_exponent();
    inke_ct_t ciphertext={0}; ibz_t modulus,r,s,omega,inverse;
    ibz_init(&modulus); ibz_init(&r); ibz_init(&s); ibz_init(&omega); ibz_init(&inverse);
    modulus_two(&modulus,a+2); shake256ctx rng; shake256ctx *prng=NULL;
    if (seed) { shake256_absorb(&rng,seed,seed_len); prng=&rng; }
    /* Primitive pairs give every cyclic subgroup exactly phi(3^b) times. */
    do { if (!random_residue(&r,&TORSION_PLUS_3POWER,prng) ||
             !random_residue(&s,&TORSION_PLUS_3POWER,prng)) goto done; }
    while (mpz_divisible_ui_p(r,3) && mpz_divisible_ui_p(s,3));
    if (!random_unit(&omega,&modulus,prng)) goto done;
    mpz_invert(inverse,omega,modulus);
    digit_t rd[NWORDS_ORDER]={0},sd[NWORDS_ORDER]={0}; ibz_to_digits(rd,&r); ibz_to_digits(sd,&s);
    ec_isog_odd_t phi; ec_basis_t basis=BASIS_EVEN; ec_curve_t shared;
    restrict_two(&basis,&CURVE_E0,a+2);
    odd_isogeny(&phi,&CURVE_E0,&BASIS_THREE,rd,sd);
    ec_eval_three(&ciphertext.EB,&phi,(ec_point_t *)&basis,3);
    scale_basis(&ciphertext.PQ2_B,&basis,&ciphertext.EB,&omega,&inverse);
    basis=pk->PQ2; odd_isogeny(&phi,&pk->EA,&pk->PQ3,rd,sd);
    ec_eval_three(&ciphertext.EAB,&phi,(ec_point_t *)&basis,3);
    scale_basis(&ciphertext.PQ2_AB,&basis,&ciphertext.EAB,&omega,&inverse);
    odd_isogeny(&phi,&pk->EA1,&pk->PQA13,rd,sd);
    ec_eval_three(&shared,&phi,(ec_point_t *)&basis,0);
    unsigned char key[32]; curve_key(key,&shared);
    for (int k=0;k<32;++k) ciphertext.ct[k]=m[k]^key[k];
    *ct=ciphertext; ok=1;
done:
    if (prng) shake256_ctx_release(prng);
    ibz_finalize(&inverse); ibz_finalize(&omega); ibz_finalize(&s); ibz_finalize(&r); ibz_finalize(&modulus);
    return ok;
}
int inke_decrypt_uniform_with_fallback(unsigned char m[32],const inke_ct_t *ct,const inke_sk_t *sk,
    inke_product_fallback fallback,void *fallback_context,int *stage)
{
    int ok=0; unsigned int a=public_exponent();
    if (stage) *stage=-1;
    ibz_t modulus,q,alpha,beta; ibz_init(&modulus); ibz_init(&q); ibz_init(&alpha); ibz_init(&beta);
    modulus_two(&modulus,a+2); ibz_copy_digits(&q,sk->deg,NWORDS_ORDER);
    ibz_copy_digits(&alpha,sk->alpha,NWORDS_ORDER); ibz_copy_digits(&beta,sk->beta,NWORDS_ORDER);
    inke_product_chain chain; inke_chain_init(&chain);
    if (!mpz_invert(alpha,alpha,modulus) || !mpz_invert(beta,beta,modulus)) goto done;
    ec_basis_t image; theta_couple_curve_t E={ct->EB,ct->EAB};
    scale_basis(&image,&ct->PQ2_AB,&E.E2,&alpha,&beta);
    theta_couple_jac_point_t K[2]; graph_kernel(K,&ct->PQ2_B,&image,&E,&q,&modulus);
    int factor=-1,detail=40; ec_curve_t shared;
#ifndef INKE_TEST_FORCE_PRODUCT_FALLBACK
    theta_couple_point_t marker;
    if (inke_chain_compute_with_lift(&chain,&E,K,a,&marker)) {
        detail=45; factor=marked_factor(&marker,&chain.codomain);
    } else detail+=chain.failed_stage;
#endif
    if (factor>=0) shared=factor?chain.codomain.E2:chain.codomain.E1;
    else {
        if (!fallback) { if(stage)*stage=detail; goto done; }
        detail=0;
        int status=fallback(&shared,NULL,NULL,&E.E1,&E.E2,&ct->PQ2_B,&image,
            a,&q,NULL,NULL,fallback_context,&detail);
        if (status!=INKE_SAMPLE_OK) {
            if(stage)*stage=status==INKE_SAMPLE_RNG_ERROR?-2:
                           status==INKE_SAMPLE_INVALID?-3:500+detail;
            goto done;
        }
    }
    unsigned char key[32]; curve_key(key,&shared);
    for (int k=0;k<32;++k) m[k]=ct->ct[k]^key[k]; ok=1;
    if(stage)*stage=0;
done:
    inke_chain_clear(&chain); ibz_finalize(&beta); ibz_finalize(&alpha); ibz_finalize(&q); ibz_finalize(&modulus);
    return ok;
}
int inke_decrypt_uniform(unsigned char m[32],const inke_ct_t *ct,const inke_sk_t *sk)
{
    return inke_decrypt_uniform_with_fallback(m,ct,sk,NULL,NULL,NULL);
}
