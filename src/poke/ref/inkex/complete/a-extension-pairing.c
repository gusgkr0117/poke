#include "a-extension-pairing.h"
/* g_{P,Q}(T)=line(P,Q)(T)/vertical(P+Q)(T). Returns zero when the
 * selected evaluation meets a divisor's support and must be retried. */
static int line_ratio(fq_t value,a_ec_point *sum,const a_ec_point *P,const a_ec_point *Q,const a_ec_point *T,const fq_t A,const root_extension *E){
    if(T->infinity)return 0;
    if(P->infinity){a_ec_set(sum,Q,E);fq_one(value,E->field);return 1;}
    if(Q->infinity){a_ec_set(sum,P,E);fq_one(value,E->field);return 1;}
    fq_t slope,den,num,tmp,one;fq_t *vars[]={&slope,&den,&num,&tmp,&one};for(int i=0;i<5;i++)fq_init(*vars[i],E->field);fq_one(one,E->field);int ok=0;
    if(fq_equal(P->x,Q->x,E->field)){
        if(!fq_equal(P->y,Q->y,E->field)||fq_is_zero(P->y,E->field)){
            fq_sub(value,T->x,P->x,E->field);a_ec_infinity(sum,E);ok=!fq_is_zero(value,E->field);goto done;
        }
        fq_sqr(num,P->x,E->field);fq_mul_ui(num,num,3,E->field);fq_mul(tmp,A,P->x,E->field);fq_mul_ui(tmp,tmp,2,E->field);fq_add(num,num,tmp,E->field);fq_add(num,num,one,E->field);fq_mul_ui(den,P->y,2,E->field);
    }else{fq_sub(num,Q->y,P->y,E->field);fq_sub(den,Q->x,P->x,E->field);}
    fq_inv(slope,den,E->field);fq_mul(slope,slope,num,E->field);a_ec_add(sum,P,Q,A,E);
    fq_sub(num,T->x,P->x,E->field);fq_mul(num,num,slope,E->field);fq_sub(tmp,T->y,P->y,E->field);fq_sub(num,tmp,num,E->field);
    if(sum->infinity)fq_one(den,E->field);else fq_sub(den,T->x,sum->x,E->field);
    if(fq_is_zero(num,E->field)||fq_is_zero(den,E->field))goto done;fq_inv(tmp,den,E->field);fq_mul(value,num,tmp,E->field);ok=1;
done:for(int i=4;i>=0;i--)fq_clear(*vars[i],E->field);return ok;
}
static int miller(fq_t out,const a_ec_point *P,const a_ec_point *T,const fmpz_t n,const fq_t A,const root_extension *E){
    if(T->infinity)return 0;if(P->infinity){fq_one(out,E->field);return 1;}
    a_ec_point R,S;a_ec_init(&R,E);a_ec_init(&S,E);a_ec_set(&R,P,E);fq_t f,g;fq_init(f,E->field);fq_init(g,E->field);fq_one(f,E->field);int ok=0;
    for(slong bit=fmpz_bits(n)-2;bit>=0;bit--){
        if(!line_ratio(g,&S,&R,&R,T,A,E))goto done;fq_sqr(f,f,E->field);fq_mul(f,f,g,E->field);a_ec_set(&R,&S,E);
        if(fmpz_tstbit(n,(ulong)bit)){if(!line_ratio(g,&S,&R,P,T,A,E))goto done;fq_mul(f,f,g,E->field);a_ec_set(&R,&S,E);}
    }fq_set(out,f,E->field);ok=1;
done:fq_clear(g,E->field);fq_clear(f,E->field);a_ec_clear(&S,E);a_ec_clear(&R,E);return ok;
}
int a_ec_weil(fq_t out,const a_ec_point *P,const a_ec_point *Q,const fmpz_t n,const fq_t A,a_ec_rng rng,void *context,unsigned long budget,const root_extension *E){
    if(fmpz_sgn(n)<=0||fmpz_divisible(n,fq_ctx_prime(E->field))||!a_ec_curve_valid(A,E)||!a_ec_on_curve(P,A,E)||!a_ec_on_curve(Q,A,E)||!rng)return A_EC_INVALID;
    a_ec_point R,negR,QR,PR,check;a_ec_point *pts[]={&R,&negR,&QR,&PR,&check};for(int i=0;i<5;i++)a_ec_init(pts[i],E);int status=A_EC_INVALID;
    fq_t f1,f2,f3,f4,value;fq_t *vars[]={&f1,&f2,&f3,&f4,&value};for(int i=0;i<5;i++)fq_init(*vars[i],E->field);
    a_ec_mul(&check,P,n,A,E);if(!check.infinity)goto done;a_ec_mul(&check,Q,n,A,E);if(!check.infinity)goto done;
    if(P->infinity||Q->infinity){fq_one(out,E->field);status=A_EC_OK;goto done;}
    status=A_EC_RETRY;
    for(unsigned long k=0;!budget||k<budget;k++){
        status=a_ec_random_point(&R,A,rng,context,1,NULL,E);if(status<0)break;if(status==A_EC_RETRY)continue;
        a_ec_neg(&negR,&R,E);a_ec_add(&QR,Q,&R,A,E);a_ec_add(&PR,P,&negR,A,E);
        status=A_EC_RETRY;
        if(!miller(f1,P,&QR,n,A,E)||!miller(f2,P,&R,n,A,E)||!miller(f3,Q,&negR,n,A,E)||!miller(f4,Q,&PR,n,A,E))continue;
        fq_mul(f1,f1,f3,E->field);fq_mul(f2,f2,f4,E->field);fq_inv(value,f2,E->field);fq_mul(value,value,f1,E->field);
        fq_pow(f1,value,n,E->field);if(!fq_is_one(f1,E->field)){status=A_EC_INVALID;break;}
        fq_set(out,value,E->field);status=A_EC_OK;break;
    }
done:for(int i=4;i>=0;i--)fq_clear(*vars[i],E->field);for(int i=4;i>=0;i--)a_ec_clear(pts[i],E);return status;
}
