// SPDX-License-Identifier: Apache-2.0
#include "b-poke-basis.h"
#include "root-poke-field-bridge.h"
#include <endomorphism_action.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <string.h>

static int same_prime(const root_extension *E){
    fmpz_t p;fmpz_init(p);fmpz_set_mpz(p,QUATALG_PINFTY.p);
    int ok=fmpz_equal(p,fq_ctx_prime(E->field));fmpz_clear(p);return ok;
}
static int xz_import(a_ec_point *out,const ec_point_t *P,const fq_t A,const root_extension *E){
    fq_t z,t;fq_init(z,E->field);fq_init(t,E->field);
    root_from_poke_fp2(out->x,&P->x,E);root_from_poke_fp2(z,&P->z,E);
    int ok=1;
    if(fq_is_zero(z,E->field)){
        ok=!fq_is_zero(out->x,E->field);if(ok)a_ec_infinity(out,E);
    }else{
        fq_inv(z,z,E->field);fq_mul(out->x,out->x,z,E->field);
        fq_add(t,out->x,A,E->field);fq_mul(t,t,out->x,E->field);
        fq_one(z,E->field);fq_add(t,t,z,E->field);fq_mul(t,t,out->x,E->field);
        ok=root_extension_sqrt(out->y,t,E);out->infinity=0;
    }
    fq_clear(t,E->field);fq_clear(z,E->field);return ok;
}
static int same_x(const a_ec_point *P,const a_ec_point *Q,const root_extension *E){
    return P->infinity||Q->infinity?P->infinity&&Q->infinity:fq_equal(P->x,Q->x,E->field);
}
int b_poke_basis_lift(fq_t A,a_ec_point out[3],const ec_curve_t *C,const ec_basis_t *B,const root_extension *E){
    if(!same_prime(E))return 0;
    fq_t a,c;fq_init(a,E->field);fq_init(c,E->field);root_from_poke_fp2(a,&C->A,E);root_from_poke_fp2(c,&C->C,E);
    a_ec_point P[3],negative,difference;for(int i=0;i<3;i++)a_ec_init(P+i,E);a_ec_init(&negative,E);a_ec_init(&difference,E);
    int ok=!fq_is_zero(c,E->field);if(!ok)goto done;
    fq_inv(c,c,E->field);fq_mul(a,a,c,E->field);if(!a_ec_curve_valid(a,E)) {ok=0;goto done;}
    ok=xz_import(P,&B->P,a,E)&&xz_import(P+1,&B->Q,a,E)&&xz_import(P+2,&B->PmQ,a,E);if(!ok)goto done;
    a_ec_neg(&negative,P+1,E);a_ec_add(&difference,P,&negative,a,E);
    if(!same_x(&difference,P+2,E)){
        a_ec_neg(P+1,P+1,E);a_ec_neg(&negative,P+1,E);a_ec_add(&difference,P,&negative,a,E);
        if(!same_x(&difference,P+2,E)){ok=0;goto done;}
    }
    a_ec_set(P+2,&difference,E);fq_set(A,a,E->field);for(int i=0;i<3;i++)a_ec_set(out+i,P+i,E);
done:
    a_ec_clear(&difference,E);a_ec_clear(&negative,E);for(int i=0;i<3;i++)a_ec_clear(P+i,E);fq_clear(c,E->field);fq_clear(a,E->field);return ok;
}
static int xz_export(ec_point_t *P,const a_ec_point *Q,const root_extension *E){
    if(Q->infinity){fp2_set_one(&P->x);fp2_set_zero(&P->z);return 1;}
    if(!root_to_poke_fp2(&P->x,Q->x,E))return 0;
    /* Signed descent matters before discarding y: a nonrational signed point
     * on a rational x-coordinate must not silently pass this boundary. */
    fp2_t y;if(!root_to_poke_fp2(&y,Q->y,E))return 0;fp2_set_one(&P->z);return 1;
}
int b_poke_basis_export(ec_curve_t *curve,ec_basis_t *basis,const fq_t A,const a_ec_point P[3],const root_extension *E){
    if(!same_prime(E)||!a_ec_curve_valid(A,E))return 0;
    ec_curve_t C;ec_basis_t B;memset(&C,0,sizeof C);memset(&B,0,sizeof B);
    a_ec_point neg,difference;a_ec_init(&neg,E);a_ec_init(&difference,E);
    int ok=a_ec_on_curve(P,A,E)&&a_ec_on_curve(P+1,A,E)&&a_ec_on_curve(P+2,A,E);
    if(ok){a_ec_neg(&neg,P+1,E);a_ec_add(&difference,P,&neg,A,E);ok=a_ec_equal(&difference,P+2,E);}
    if(ok)ok=root_to_poke_fp2(&C.A,A,E)&&xz_export(&B.P,P,E)&&xz_export(&B.Q,P+1,E)&&xz_export(&B.PmQ,P+2,E);
    if(ok){fp2_set_one(&C.C);*curve=C;*basis=B;}
    a_ec_clear(&difference,E);a_ec_clear(&neg,E);return ok;
}
int b_poke_fixed_action(a_ec_point out[3],int three,const quat_alg_elem_t *alpha,const ibz_t *scale,const root_extension *E){
    if((three!=0&&three!=1)||!same_prime(E))return 0;
    const ibz_t *modulus=three?&TORSION_PLUS_3POWER:&TORSION_PLUS_2POWER;
    const ibz_mat_2x2_t *actions[3]={&ACTION_GEN2,&ACTION_GEN3,&ACTION_GEN4};
    ibz_vec_4_t coeff;ibz_vec_4_init(&coeff);ibz_t entry;ibz_init(&entry);fmpz_t n;fmpz_init(n);
    fq_t A;fq_init(A,E->field);a_ec_point source[3],image[3],term;for(int i=0;i<3;i++){a_ec_init(source+i,E);a_ec_init(image+i,E);}a_ec_init(&term,E);
    int ok=quat_lattice_contains(&coeff,&MAXORD_O0,alpha,&QUATALG_PINFTY)&&b_poke_basis_lift(A,source,&CURVE_E0,three?&BASIS_THREE:&BASIS_EVEN,E);
    if(!ok)goto done;
    for(int j=0;j<2;j++){
        a_ec_infinity(image+j,E);
        for(int i=0;i<2;i++){
            mpz_set_ui(entry,0);if(i==j)mpz_set(entry,coeff[0]);
            for(int k=0;k<3;k++)mpz_addmul(entry,coeff[k+1],(*actions[k])[i][j]);
            mpz_mul(entry,entry,*scale);mpz_mod(entry,entry,*modulus);fmpz_set_mpz(n,entry);
            a_ec_mul(&term,source+i,n,A,E);a_ec_add(image+j,image+j,&term,A,E);
        }
    }
    a_ec_neg(&term,image+1,E);a_ec_add(image+2,image,&term,A,E);for(int i=0;i<3;i++)a_ec_set(out+i,image+i,E);
done:
    a_ec_clear(&term,E);for(int i=0;i<3;i++){a_ec_clear(image+i,E);a_ec_clear(source+i,E);}fq_clear(A,E->field);fmpz_clear(n);ibz_finalize(&entry);ibz_vec_4_finalize(&coeff);return ok;
}
int b_gaussian_fixed_payloads(a_ec_point two[3],a_ec_point three[3],const quat_left_ideal_t *I,const b_gaussian_certificate *C,const root_extension *E){
    if(!b_gaussian_certificate_verify(C,I))return 0;
    a_ec_point images[6];for(int i=0;i<6;i++)a_ec_init(images+i,E);ibz_t inverse;ibz_init(&inverse);
    int ok=mpz_invert(inverse,C->n1,TORSION_PLUS_2POWER)&&b_poke_fixed_action(images,0,&C->transport,&inverse,E);
    if(ok)ok=mpz_invert(inverse,C->n1,TORSION_PLUS_3POWER)&&b_poke_fixed_action(images+3,1,&C->transport,&inverse,E);
    if(ok)for(int i=0;i<3;i++){a_ec_set(two+i,images+i,E);a_ec_set(three+i,images+3+i,E);}
    ibz_finalize(&inverse);for(int i=0;i<6;i++)a_ec_clear(images+i,E);return ok;
}
int b_poke_conversion_export(inke_conversion_result *out,const fq_t A,const a_ec_point two[3],const ibz_t *degree,int steps,const root_extension *E){
    if(steps<0||mpz_sgn(*degree)<=0)return 0;
    inke_conversion_result result;memset(&result,0,sizeof result);
    int ok=b_poke_basis_export(&result.curve,&result.two,A,two,E);
    if(ok)ok=inke_pairing_degree(&BASIS_EVEN,&CURVE_E0,&result.two,&result.curve,TORSION_PLUS_EVEN_POWER,degree);
    if(ok){result.product_steps=steps;*out=result;}return ok;
}
