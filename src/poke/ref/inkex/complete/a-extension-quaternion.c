#include "a-extension-quaternion.h"
#include <limits.h>
static void action_i(a_ec_point *R,const a_ec_point *P,const root_extension *E){
    fq_neg(R->x,P->x,E->field);fq_mul(R->y,P->y,E->i,E->field);R->infinity=P->infinity;
}
int a_ec_quaternion(a_ec_point *out,const a_ec_point *P,const fmpz *coeff,unsigned den,const fmpz_t order,const a_ec_point *lift,const root_extension *E){
    if((den!=1&&den!=2)||fmpz_sgn(order)<=0)return A_EC_INVALID;
    if(den==2&&(fmpz_fdiv_ui(coeff+0,2)!=fmpz_fdiv_ui(coeff+3,2)||fmpz_fdiv_ui(coeff+1,2)!=fmpz_fdiv_ui(coeff+2,2)))return A_EC_INVALID;
    fq_t zero;fq_init(zero,E->field);fq_zero(zero,E->field);int status=A_EC_INVALID;
    a_ec_point basis[4],acc,tmp,check;for(int j=0;j<4;j++)a_ec_init(&basis[j],E);a_ec_init(&acc,E);a_ec_init(&tmp,E);a_ec_init(&check,E);
    fmpz_t modulus,c,two;fmpz_init(modulus);fmpz_init(c);fmpz_init(two);fmpz_set_ui(two,2);
    if(!a_ec_on_curve(P,zero,E))goto done;a_ec_mul(&check,P,order,zero,E);if(!check.infinity)goto done;
    fmpz_set(modulus,order);
    if(den==2&&!fmpz_fdiv_ui(order,2)){
        if(!lift||!a_ec_on_curve(lift,zero,E))goto done;a_ec_mul(&check,lift,two,zero,E);if(!a_ec_equal(&check,P,E))goto done;
        fmpz_mul_ui(modulus,order,2);a_ec_mul(&check,lift,modulus,zero,E);if(!check.infinity)goto done;a_ec_set(&basis[0],lift,E);
    }else a_ec_set(&basis[0],P,E);
    action_i(&basis[1],&basis[0],E);a_ec_frobenius(&basis[2],&basis[0],1,E);action_i(&basis[3],&basis[2],E);a_ec_infinity(&acc,E);
    for(int j=0;j<4;j++){fmpz_mod(c,coeff+j,modulus);a_ec_mul(&tmp,&basis[j],c,zero,E);a_ec_add(&acc,&acc,&tmp,zero,E);}
    if(den==2&&fmpz_fdiv_ui(order,2)){if(!fmpz_invmod(c,two,order))goto done;a_ec_mul(&acc,&acc,c,zero,E);}
    if(!a_ec_on_curve(&acc,zero,E))goto done;a_ec_mul(&check,&acc,order,zero,E);if(!check.infinity)goto done;a_ec_set(out,&acc,E);status=A_EC_OK;
done:fmpz_clear(two);fmpz_clear(c);fmpz_clear(modulus);a_ec_clear(&check,E);a_ec_clear(&tmp,E);a_ec_clear(&acc,E);for(int j=3;j>=0;j--)a_ec_clear(&basis[j],E);fq_clear(zero,E->field);return status;
}
void a_clapoti_kernel_init(a_clapoti_kernel *K,const root_extension *E){
    for(int i=0;i<8;i++)for(int j=0;j<8;j++){a_ec_init(&K->generator[i][j],E);a_ec_init(&K->complement[i][j],E);}
    for(int i=0;i<2;i++)a_ec_init(&K->basis[i],E);fmpz_init(K->order);
}
void a_clapoti_kernel_clear(a_clapoti_kernel *K,const root_extension *E){
    fmpz_clear(K->order);for(int i=1;i>=0;i--)a_ec_clear(&K->basis[i],E);
    for(int i=7;i>=0;i--)for(int j=7;j>=0;j--){a_ec_clear(&K->complement[i][j],E);a_ec_clear(&K->generator[i][j],E);}
}
int a_clapoti_kernel_lift(a_clapoti_kernel *out,unsigned long ell,unsigned e,unsigned extra,const fmpz_t n1,const fmpz_t n2,const fmpz *theta,unsigned den,const fmpz_mat_t matrix,a_ec_rng rng,void *context,unsigned long budget,unsigned long *trials,const root_extension *E){
    if(!e||extra>=UINT_MAX-e||ell<2||fmpz_fdiv_ui(n1,ell)==0||fmpz_fdiv_ui(n2,ell)==0||fmpz_mat_nrows(matrix)!=4||fmpz_mat_ncols(matrix)!=4||(den!=1&&den!=2))return A_EC_INVALID;
    a_clapoti_kernel K;a_clapoti_kernel_init(&K,E);a_ec_point lifts[2],images[2],scaled[2];for(int i=0;i<2;i++){a_ec_init(&lifts[i],E);a_ec_init(&images[i],E);a_ec_init(&scaled[i],E);}
    fq_t zero;fq_init(zero,E->field);fq_zero(zero,E->field);
    fmpz_t scalar,two,base,t,inv_product,inv_n1;fmpz_init(scalar);fmpz_init(two);fmpz_init(base);fmpz_init(t);fmpz_init(inv_product);fmpz_init(inv_n1);
    fmpz_set_ui(two,2);fmpz_set_ui(base,ell);fmpz_pow_ui(base,base,e);fmpz_add(scalar,n1,n2);int status=A_EC_INVALID;if(!fmpz_divisible(scalar,base))goto done;
    fmpz_set_ui(K.order,ell);fmpz_pow_ui(K.order,K.order,e+extra);
    if(!fmpz_invmod(t,n2,K.order)||!fmpz_invmod(inv_n1,n1,K.order))goto done;
    fmpz_mul(t,t,n1);fmpz_neg(t,t);fmpz_mod(t,t,K.order);fmpz_mul(scalar,n1,n2);if(!fmpz_invmod(inv_product,scalar,K.order))goto done;
    int half=(ell==2&&den==2);status=a_ec_torsion_basis(&lifts[0],&lifts[1],ell,e+extra+half,zero,rng,context,budget,trials,E);if(status!=A_EC_OK)goto done;
    for(int i=0;i<2;i++){
        if(half)a_ec_mul(&K.basis[i],&lifts[i],two,zero,E);else a_ec_set(&K.basis[i],&lifts[i],E);
        status=a_ec_quaternion(&images[i],&K.basis[i],theta,den,K.order,half?&lifts[i]:NULL,E);if(status!=A_EC_OK)goto done;
        fmpz_mod(scalar,n1,K.order);a_ec_mul(&scaled[i],&K.basis[i],scalar,zero,E);
    }
    for(int j=0;j<4;j++)for(int s=0;s<2;s++){
        int idx=2*j+s;for(int row=0;row<8;row++){a_ec_infinity(&K.generator[idx][row],E);a_ec_infinity(&K.complement[idx][row],E);}
        a_ec_set(&K.generator[idx][j],&scaled[s],E);
        for(int row=0;row<4;row++){
            fmpz_set(scalar,fmpz_mat_entry(matrix,row,j));if(s)fmpz_mul(scalar,scalar,t);fmpz_mod(scalar,scalar,K.order);
            a_ec_mul(&K.generator[idx][row+4],&images[s],scalar,zero,E);
        }
        if(!s){
            for(int row=0;row<4;row++){fmpz_mul(scalar,fmpz_mat_entry(matrix,row,j),inv_product);fmpz_mod(scalar,scalar,K.order);a_ec_mul(&K.complement[idx][row+4],&images[1],scalar,zero,E);}
        }else{fmpz_neg(scalar,inv_n1);fmpz_mod(scalar,scalar,K.order);a_ec_mul(&K.complement[idx][j],&K.basis[0],scalar,zero,E);}
    }
    for(int i=0;i<8;i++)for(int j=0;j<8;j++){a_ec_set(&out->generator[i][j],&K.generator[i][j],E);a_ec_set(&out->complement[i][j],&K.complement[i][j],E);}
    for(int i=0;i<2;i++)a_ec_set(&out->basis[i],&K.basis[i],E);fmpz_set(out->order,K.order);status=A_EC_OK;
done:fmpz_clear(inv_n1);fmpz_clear(inv_product);fmpz_clear(t);fmpz_clear(base);fmpz_clear(two);fmpz_clear(scalar);fq_clear(zero,E->field);
    for(int i=1;i>=0;i--){a_ec_clear(&scaled[i],E);a_ec_clear(&images[i],E);a_ec_clear(&lifts[i],E);}a_clapoti_kernel_clear(&K,E);return status;
}
int a_clapoti_kernel_component(a_clapoti_kernel *out,unsigned long ell,unsigned e,const fmpz_t n1,const fmpz *theta,unsigned den,const fmpz_mat_t matrix,a_ec_rng rng,void *context,unsigned long budget,unsigned long *trials,const root_extension *E){
    fmpz_t n2;fmpz_init(n2);fmpz_neg(n2,n1);int status=a_clapoti_kernel_lift(out,ell,e,0,n1,n2,theta,den,matrix,rng,context,budget,trials,E);fmpz_clear(n2);return status;
}