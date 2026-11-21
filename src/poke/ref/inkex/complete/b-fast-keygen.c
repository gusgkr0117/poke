// SPDX-License-Identifier: Apache-2.0
#include "b-fast-keygen.h"
#include <id2iso.h>
#include <quaternion_data.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <stdlib.h>
#include <string.h>

static void erase_bytes(unsigned char *p,size_t n){
    volatile unsigned char *v=p;while(n--)*v++=0;
}
/* Exact interval rejection sampling; its expected number of trials is <2.
 * Random-source failure is reported without changing the sampled q. */
static int interval(ibz_t *out,const ibz_t *lo,const ibz_t *hi,
    b_fast_counts *counts,a_ec_rng rng,void *context){
    if(mpz_cmp(*lo,*hi)>0)return INKE_SAMPLE_INVALID;
    ibz_t width,x;ibz_init(&width);ibz_init(&x);mpz_sub(width,*hi,*lo);
    int status=INKE_SAMPLE_RETRY;
    if(!mpz_sgn(width)){mpz_set(*out,*lo);status=INKE_SAMPLE_OK;goto done;}
    size_t bits=mpz_sizeinbase(width,2),bytes=(bits+7)/8;
    unsigned char *buffer=malloc(bytes);if(!buffer){status=INKE_SAMPLE_INVALID;goto done;}
    for(;;){
        if(!rng(context,buffer,bytes)){status=INKE_SAMPLE_RNG_ERROR;break;}
        if(bits%8)buffer[bytes-1]&=(unsigned char)((1u<<(bits%8))-1u);
        mpz_import(x,bytes,-1,1,0,0,buffer);
        if(mpz_cmp(x,width)<=0){mpz_add(*out,x,*lo);status=INKE_SAMPLE_OK;break;}
    }
    (void)counts;
    erase_bytes(buffer,bytes);free(buffer);
done:ibz_finalize(&x);ibz_finalize(&width);return status;
}
static int two_squares(ibz_t *x,ibz_t *y,const ibz_t *r,b_fast_counts *counts,
    a_ec_rng rng,void *context){
    if(mpz_cmp_ui(*r,2)==0){mpz_set_ui(*x,1);mpz_set_ui(*y,1);return INKE_SAMPLE_OK;}
    if(mpz_cmp_ui(*r,2)<0||mpz_fdiv_ui(*r,4)!=1)return INKE_SAMPLE_RETRY;
    ibz_t hi,e,w,t,u,v,s;ibz_init(&hi);ibz_init(&e);ibz_init(&w);
    ibz_init(&t);ibz_init(&u);ibz_init(&v);ibz_init(&s);
    mpz_sub_ui(hi,*r,1);mpz_fdiv_q_2exp(e,hi,2);
    int status=INKE_SAMPLE_RETRY,found=0;
    size_t budget=mpz_sizeinbase(*r,2)+128;
    for(size_t j=0;j<budget;j++){
        status=interval(&w,&ibz_const_one,&hi,counts,rng,context);
        if(status!=INKE_SAMPLE_OK)goto done;
        mpz_powm(t,w,e,*r);mpz_mul(s,t,t);mpz_mod(s,s,*r);
        if(mpz_cmp(s,hi)==0){found=1;break;}
    }
    status=INKE_SAMPLE_RETRY;if(!found)goto done;
    mpz_set(u,*r);mpz_set(v,t);
    for(;;){
        mpz_mul(s,v,v);if(mpz_cmp(s,*r)<=0)break;
        mpz_mod(t,u,v);mpz_set(u,v);mpz_set(v,t);
    }
    mpz_sub(s,*r,s);
    if(mpz_perfect_square_p(s)){mpz_set(*x,v);mpz_sqrt(*y,s);status=INKE_SAMPLE_OK;}
done:ibz_finalize(&s);ibz_finalize(&v);ibz_finalize(&u);ibz_finalize(&t);
    ibz_finalize(&w);ibz_finalize(&e);ibz_finalize(&hi);return status;
}
static int exact_theta(const quat_alg_elem_t *theta,const ibz_t *target){
    if(mpz_sgn(theta->denom)<=0)return 0;
    ibq_t n;ibq_init(&n);ibz_t exact;ibz_init(&exact);
    quat_alg_norm(&n,theta,&QUATALG_PINFTY);
    int ok=ibq_to_ibz(&exact,&n)&&mpz_cmp(exact,*target)==0&&
        quat_lattice_contains(NULL,&MAXORD_O0,theta,&QUATALG_PINFTY)&&
        quat_alg_is_primitive(theta,&MAXORD_O0,&QUATALG_PINFTY);
    ibz_finalize(&exact);ibq_finalize(&n);return ok;
}
int b_fast_sample_theta(quat_alg_elem_t *out,const ibz_t *target,
    unsigned long limit,b_fast_counts *counts,a_ec_rng rng,void *context){
    if(!out||!counts||!rng||mpz_sgn(*target)<=0||mpz_even_p(*target)||
       mpz_fdiv_ui(QUATALG_PINFTY.p,4)!=3)return INKE_SAMPLE_INVALID;
    if(!limit||counts->norm_candidates>=limit||mpz_cmp(*target,QUATALG_PINFTY.p)<=0)
        return INKE_SAMPLE_RETRY;
    /* Bounded prime-residual norm search in the standard maximal order.
     * Same coordinate bounds as GeneralizedRepresentInteger at q_form=1:
     * z<=floor(sqrt(floor(4M/p)-1)), |t|<=floor(sqrt(floor(4M/p)-z^2)).
     * The fast law includes these choices; no output uniformity is assumed. */
    ibz_t four,sq,zmax,tmax,lo,tcount,residual,tmp;ibz_init(&four);ibz_init(&sq);
    ibz_init(&zmax);ibz_init(&tmax);ibz_init(&lo);ibz_init(&residual);ibz_init(&tmp);
    ibz_init(&tcount);
    quat_alg_elem_t theta;quat_alg_elem_init(&theta);
    mpz_mul_2exp(four,*target,2);mpz_fdiv_q(sq,four,QUATALG_PINFTY.p);
    mpz_sub_ui(zmax,sq,1);mpz_sqrt(zmax,zmax);
    int status=INKE_SAMPLE_RETRY;
    while(counts->norm_candidates<limit){
        counts->norm_candidates++;
        status=interval(&theta.coord[2],&ibz_const_one,&zmax,counts,rng,context);
        if(status!=INKE_SAMPLE_OK)break;
        mpz_mul(tmp,theta.coord[2],theta.coord[2]);
        mpz_sub(tmax,sq,tmp);mpz_sqrt(tmax,tmax);mpz_neg(lo,tmax);
        /* Uniformly select t of parity opposite to z in the same interval. */
        if(mpz_tstbit(lo,0)==mpz_tstbit(theta.coord[2],0))mpz_add_ui(lo,lo,1);
        mpz_sub(tcount,tmax,lo);mpz_fdiv_q_2exp(tcount,tcount,1);
        status=interval(&theta.coord[3],&ibz_const_zero,&tcount,counts,rng,context);
        if(status!=INKE_SAMPLE_OK)break;
        mpz_mul_2exp(theta.coord[3],theta.coord[3],1);
        mpz_add(theta.coord[3],theta.coord[3],lo);
        mpz_mul(residual,theta.coord[3],theta.coord[3]);mpz_add(residual,residual,tmp);
        mpz_mul(residual,residual,QUATALG_PINFTY.p);mpz_sub(residual,four,residual);
        status=INKE_SAMPLE_RETRY;
        if(mpz_cmp_ui(residual,2)<0||!mpz_probab_prime_p(residual,64))continue;
        status=two_squares(&theta.coord[0],&theta.coord[1],&residual,counts,rng,context);
        if(status==INKE_SAMPLE_RNG_ERROR||status==INKE_SAMPLE_INVALID)break;
        if(status!=INKE_SAMPLE_OK)continue;
        /* In O0, x=t and y=z modulo 2. Exchanging x,y preserves the norm. */
        if(mpz_tstbit(theta.coord[0],0)!=mpz_tstbit(theta.coord[3],0))
            mpz_swap(theta.coord[0],theta.coord[1]);
        mpz_set_ui(theta.denom,2);status=INKE_SAMPLE_RETRY;
        if(!exact_theta(&theta,target))continue;
        quat_alg_normalize(&theta);for(unsigned j=0;j<4;j++)ibz_copy(&out->coord[j],&theta.coord[j]);
        ibz_copy(&out->denom,&theta.denom);
        status=INKE_SAMPLE_OK;break;
    }
    quat_alg_elem_finalize(&theta);ibz_finalize(&tmp);ibz_finalize(&residual);
    ibz_finalize(&tcount);ibz_finalize(&lo);ibz_finalize(&tmax);ibz_finalize(&zmax);ibz_finalize(&sq);
    ibz_finalize(&four);return status;
}
int b_fast_convert_theta(inke_conversion_result *out,const quat_alg_elem_t *theta,
    const ibz_t *d,int *stage){
    if(stage)*stage=10;
    if(!out||!theta||mpz_sgn(*d)<=0||mpz_even_p(*d)||mpz_divisible_ui_p(*d,3))
        return INKE_SAMPLE_INVALID;
    ibz_t target,inverse;ibz_init(&target);ibz_init(&inverse);
    mpz_mul(target,*d,TORSION_PLUS_3POWER);
    int status=INKE_SAMPLE_INVALID;
    if(!exact_theta(theta,&target))goto done;
    ec_basis_t image2,image3;ec_curve_t E0=CURVE_E0;
    if(!inke_endomorphism_basis(&image2,&BASIS_EVEN,theta,&TORSION_PLUS_2POWER)||
       !inke_endomorphism_basis(&image3,&BASIS_THREE,theta,&TORSION_PLUS_3POWER))goto done;
    if(stage)*stage=11;
    ec_isog_odd_t eta={0};eta.curve=E0;eta.degree[0]=POWER_OF_3;
    fp2_set_one(&eta.ker_minus.x);fp2_set_zero(&eta.ker_minus.z);
    if(test_point_order_threef(&image3.P,&E0,POWER_OF_3))eta.ker_plus=image3.P;
    else if(test_point_order_threef(&image3.Q,&E0,POWER_OF_3))eta.ker_plus=image3.Q;
    else goto done; /* A verified primitive theta has a primitive image column. */
    inke_conversion_result result={0};ec_eval_three(&result.curve,&eta,(ec_point_t *)&image2,3);
    if(!mpz_invert(inverse,TORSION_PLUS_3POWER,TORSION_PLUS_2POWER))goto done;
    if(stage)*stage=12;
    status=INKE_SAMPLE_RETRY;
    if(fp2_is_zero(&result.curve.C))goto done;
    fp2_t disc,t;fp2_sqr(&disc,&result.curve.A);fp2_sqr(&t,&result.curve.C);
    fp2_add(&t,&t,&t);fp2_add(&t,&t,&t);fp2_sub(&disc,&disc,&t);
    if(fp2_is_zero(&disc))goto done;
    ec_mul_ibz(&image2.P,&result.curve,&inverse,&image2.P);
    ec_mul_ibz(&image2.Q,&result.curve,&inverse,&image2.Q);
    ec_mul_ibz(&image2.PmQ,&result.curve,&inverse,&image2.PmQ);
    if(!inke_pairing_degree(&BASIS_EVEN,&E0,&image2,&result.curve,TORSION_PLUS_EVEN_POWER,d))goto done;
    result.two=image2;result.product_steps=0;*out=result;status=INKE_SAMPLE_OK;
    if(stage)*stage=0;
done:ibz_finalize(&inverse);ibz_finalize(&target);return status;
}
int b_fast_receiver(inke_conversion_result *out,const ibz_t *q,unsigned int a,
    unsigned long limit,b_fast_counts *counts,a_ec_rng rng,void *context,int *stage){
    if(stage)*stage=1;
    if(!out||!counts||!rng||!inke_q_is_admissible(q,a))
        return INKE_SAMPLE_INVALID;
    ibz_t A,d,target;ibz_init(&A);ibz_init(&d);ibz_init(&target);
    mpz_set_ui(A,1);mpz_mul_2exp(A,A,a);int status=INKE_SAMPLE_INVALID;
    if(mpz_cmp(*q,A)>=0)goto done;
    mpz_sub(d,A,*q);mpz_mul(d,d,*q);mpz_mul(target,d,TORSION_PLUS_3POWER);
    status=INKE_SAMPLE_RETRY;
    if(!limit||counts->norm_candidates>=limit)goto done;
    if(mpz_cmp(target,QUATALG_PINFTY.p)<=0){counts->guard_misses++;goto done;}
    quat_alg_elem_t theta;quat_alg_elem_init(&theta);
    while(counts->norm_candidates<limit){
        if(stage)*stage=2;
        status=b_fast_sample_theta(&theta,&target,limit,counts,rng,context);
        if(status!=INKE_SAMPLE_OK)break;
        counts->witnesses++;status=b_fast_convert_theta(out,&theta,&d,stage);
        if(status!=INKE_SAMPLE_RETRY)break;
        counts->conversion_failures++;
    }
    quat_alg_elem_finalize(&theta);
done:ibz_finalize(&target);ibz_finalize(&d);ibz_finalize(&A);return status;
}
