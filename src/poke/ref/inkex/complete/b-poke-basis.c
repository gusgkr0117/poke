// SPDX-License-Identifier: Apache-2.0
#include "b-poke-basis.h"
#include "root-poke-field-bridge.h"
#include <quaternion_data.h>
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
