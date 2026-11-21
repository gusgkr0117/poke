// SPDX-License-Identifier: Apache-2.0
#include "inke_chain.h"
#include <stdlib.h>
#include <string.h>

void theta_product_structure_to_elliptic_product(theta_couple_curve_t *, theta_structure_t *);
void apply_isomorphism(theta_point_t *,const theta_splitting_t *,const theta_point_t *);
void theta_isogeny_comput4(theta_isogeny_t *,const theta_structure_t *,
    const theta_point_t *,const theta_point_t *,int,int);
void theta_isogeny_comput2(theta_isogeny_t *,const theta_structure_t *,
    const theta_point_t *,const theta_point_t *,int,int);

static int valid_point(const theta_point_t *P)
{
    return !(fp2_is_zero(&P->x) && fp2_is_zero(&P->y) &&
             fp2_is_zero(&P->z) && fp2_is_zero(&P->t));
}
static int valid_surface(const theta_structure_t *A)
{
    theta_point_t h;
    to_squared_theta(&h,&A->null_point);
    return !fp2_is_zero(&A->null_point.x) && !fp2_is_zero(&A->null_point.y) &&
           !fp2_is_zero(&A->null_point.z) && !fp2_is_zero(&A->null_point.t) &&
           !fp2_is_zero(&h.x) && !fp2_is_zero(&h.y) &&
           !fp2_is_zero(&h.z) && !fp2_is_zero(&h.t);
}
static int normalize_point(theta_point_t *P)
{
    fp2_t *coordinates[4]={&P->x,&P->y,&P->z,&P->t};
    int first=0;
    while (first<4 && fp2_is_zero(coordinates[first])) ++first;
    if (first==4) return 0;
    fp2_t inverse=*coordinates[first]; fp2_inv(&inverse);
    for (int j=0;j<4;++j) fp2_mul(coordinates[j],coordinates[j],&inverse);
    return 1;
}
void inke_chain_init(inke_product_chain *c) { memset(c,0,sizeof(*c)); }
void inke_chain_clear(inke_product_chain *c)
{
    free(c->two); free(c->three); free(c->product_prefix); memset(c,0,sizeof(*c));
}

static void add_pair(theta_couple_jac_point_t *out,
    const theta_couple_jac_point_t *P,const theta_couple_jac_point_t *Q,
    const theta_couple_curve_t *E,int subtract)
{
    theta_couple_jac_point_t R=*Q;
    if (subtract) { fp2_neg(&R.P1.y,&R.P1.y); fp2_neg(&R.P2.y,&R.P2.y); }
    ADD(&out->P1,&P->P1,&R.P1,&E->E1); ADD(&out->P2,&P->P2,&R.P2,&E->E2);
}

/* On y^2=x^3+x, i(x,y)=(-x,sqrt(-1)y). The map
 * (P,Q) -> (P+uQ,P-uQ), u=1 or i, has multiplier 2 and kills
 * the graph of u on E[2]. Its codomain remains E x E. */
static void product_step(theta_couple_jac_point_t *out,
    const theta_couple_jac_point_t *P,const theta_couple_curve_t *E,int use_i)
{
    jac_point_t Q=P->P2,negative; theta_couple_jac_point_t R;
    if (use_i) {
        fp2_t imaginary; fp2_set_zero(&imaginary); fp_set_one(&imaginary.im);
        fp2_neg(&Q.x,&Q.x); fp2_mul(&Q.y,&Q.y,&imaginary);
    }
    negative=Q; fp2_neg(&negative.y,&negative.y);
    ADD(&R.P1,&P->P1,&Q,&E->E1); ADD(&R.P2,&P->P1,&negative,&E->E1); *out=R;
}
static int product_kernel(const theta_couple_jac_point_t K[2],
    const theta_couple_curve_t *E,int exponent)
{
    if (!fp2_is_zero(&E->E1.A) || !fp2_is_zero(&E->E2.A)) return -1;
    theta_couple_jac_point_t T[2]; ec_point_t P[2],Q[2];
    for (int k=0;k<2;++k) {
        double_couple_jac_point_iter(&T[k],exponent-1,E,&K[k]);
        jac_to_xz(&P[k],&T[k].P1); jac_to_xz(&Q[k],&T[k].P2);
        if (fp2_is_zero(&P[k].z) || fp2_is_zero(&Q[k].z)) return -1;
    }
    for (int use_i=0;use_i<2;++use_i) {
        if (ec_is_equal(&P[0],&Q[0]) && ec_is_equal(&P[1],&Q[1])) return use_i;
        fp2_neg(&Q[0].x,&Q[0].x); fp2_neg(&Q[1].x,&Q[1].x);
    }
    return -1;
}

/* The generators in this routine have order 2^(count+2). Splitting the
 * remaining chain preserves a high-order copy for its second part, rather
 * than deriving every order-eight kernel from the original generators. */
static int balanced_two(inke_product_chain *c,theta_structure_t *A,
    theta_point_t q[2],int count,int first)
{
    if (!count) return 1;
    if (count>1) {
        int second=count/2,initial=count-second;
        theta_point_t lower[2];
        theta_precomputation(A);
        for (int j=0;j<2;++j) {
            double_iter(&lower[j],A,&q[j],second);
            c->kernel_doublings+=(unsigned long)second;
        }
        if (!balanced_two(c,A,lower,initial,first)) return 0;
        for (int k=first;k<first+initial;++k) for (int j=0;j<2;++j) {
            theta_isogeny_eval(&q[j],&c->two[k],&q[j]);
            if (!valid_point(&q[j])) return 0;
        }
        return balanced_two(c,A,q,second,first+initial);
    }
    theta_precomputation(A);
    theta_isogeny_comput(&c->two[first],A,&q[0],&q[1],
        !c->len3 && first==c->len2-2,c->len3 || first<c->len2-3);
    theta_point_t origin;
    theta_isogeny_eval(&origin,&c->two[first],&A->null_point);
    if (!inke_theta_equal(&origin,&c->two[first].codomain.null_point)) return 0;
    for (int j=0;j<2;++j) {
        theta_isogeny_eval(&q[j],&c->two[first],&q[j]);
        if (!valid_point(&q[j])) return 0;
    }
    *A=c->two[first].codomain;
    return !(c->len3 || first<c->len2-3) || valid_surface(A);
}

/* All four nonzero classes of the rank-two 3-primary kernel are retained.
 * The split only changes when multiples and images are computed, not the
 * kernel or the choice between the six supported coefficient charts. */
static int balanced_three(inke_product_chain *c,theta_structure_t *A,
    theta_point_t q[4],int count,int first)
{
    if (!count) return 1;
    if (count>1) {
        int second=count/2,initial=count-second;
        theta_point_t lower[4];
        for (int j=0;j<4;++j) {
            lower[j]=q[j];
            for (int k=0;k<second;++k) {
                if (!inke_theta_triple(&lower[j],A,&lower[j])) return 0;
                ++c->kernel_triplings;
            }
        }
        if (!balanced_three(c,A,lower,initial,first)) return 0;
        for (int k=first;k<first+initial;++k) for (int j=0;j<4;++j)
            if (!inke_theta_isogeny3_eval(&q[j],&c->three[k],&q[j])) return 0;
        return balanced_three(c,A,q,second,first+initial);
    }
    int found=0;
    for (int j=0;j<4 && !found;++j) for (int k=j+1;k<4 && !found;++k)
        found=inke_theta_isogeny3_compute(&c->three[first],A,&q[j],&q[k]);
    if (!found) return 0;
    for (int j=0;j<4;++j) {
        theta_point_t killed;
        if (!inke_theta_isogeny3_eval(&killed,&c->three[first],&q[j]) ||
            !inke_theta_equal(&killed,&c->three[first].codomain.null_point)) return 0;
        q[j]=killed;
    }
    *A=c->three[first].codomain;
    return 1;
}

static int split_image(theta_couple_point_t *out,const inke_product_chain *c,theta_point_t P);

static int compute(inke_product_chain *c,const theta_couple_curve_t *E,
    const theta_couple_jac_point_t K2[2],int len2,
    const theta_couple_jac_point_t K3[2],int len3,int exact,int balanced,
    theta_couple_point_t *lift_image)
{
    inke_chain_clear(c);
    if (len2<1 || len3<0 || (exact && (len2<3 || len3))) return 0;
    c->domain=*E; c->len2=len2; c->len3=len3; c->failed_stage=1;
    /* Jacobian arithmetic below uses the affine Montgomery coefficient A/C. */
    ec_curve_t *curves[2]={&c->domain.E1,&c->domain.E2};
    for (int k=0;k<2;++k) {
        fp2_t inverse=curves[k]->C; if (fp2_is_zero(&inverse)) return 0;
        fp2_inv(&inverse); fp2_mul(&curves[k]->A,&curves[k]->A,&inverse);
        fp2_set_one(&curves[k]->C); curves[k]->is_A24_computed_and_normalized=false;
    }
    E=&c->domain;
    c->two=calloc(len2>1?len2-1:1,sizeof(*c->two));
    c->three=calloc(len3?len3:1,sizeof(*c->three));
    c->product_prefix=calloc(len2,sizeof(*c->product_prefix));
    if (!c->two || !c->three || !c->product_prefix) return 0;
    theta_couple_jac_point_t input2[2]={K2[0],K2[1]},input3[2];
    if (len3) { input3[0]=K3[0]; input3[1]=K3[1]; }
    int type;
    while ((type=product_kernel(input2,E,len2+2-2*exact))>=0) {
        if (len2<=3) return 0;
        c->product_prefix[c->prefix_length++]=type;
        for (int k=0;k<2;++k) {
            product_step(&input2[k],&input2[k],E,type);
            if (len3) product_step(&input3[k],&input3[k],E,type);
        }
        --len2;
    }
    c->len2=len2; K2=input2; if (len3) K3=input3;
    theta_couple_jac_point_t eight[2];
    for (int j=0;j<2;++j) double_couple_jac_point_iter(&eight[j],len2-1-2*exact,E,&K2[j]);
    if (!gluing_comput_checked(&c->glue,&c->domain,&eight[0],&eight[1])) return 0;
    theta_structure_t A={0}; A.null_point=c->glue.codomain;
    if (!valid_surface(&A)) return 0;
    /* A product here requires an explicit product transition, not the
     * generic gluing formulas used below. Reject it without changing I. */
    theta_splitting_t early;
    if (splitting_comput(&early,&A)) return 0;
    theta_point_t q2[2],q3[4];
    for (int j=0;j<2;++j) {
        gluing_eval_point_no_help(&q2[j],&K2[j],E,&c->glue);
        if (!valid_point(&q2[j])) return 0;
    }
    if (len3) {
        theta_couple_jac_point_t input[4]={K3[0],K3[1]};
        add_pair(&input[2],&K3[0],&K3[1],E,0);
        add_pair(&input[3],&K3[0],&K3[1],E,1);
        for (int j=0;j<4;++j) {
            gluing_eval_point_no_help(&q3[j],&input[j],E,&c->glue);
            if (!valid_point(&q3[j])) return 0;
        }
    }
    int normal=len2-1-2*exact;
    if (balanced) {
        c->failed_stage=2;
        if (!balanced_two(c,&A,q2,normal,0)) return 0;
        if (len3) for (int step=0;step<normal;++step) for (int j=0;j<4;++j) {
            theta_isogeny_eval(&q3[j],&c->two[step],&q3[j]);
            if (!valid_point(&q3[j])) return 0;
        }
    }
    for (int step=balanced?normal:0;step<len2-1;++step) {
        c->failed_stage=2;
        theta_point_t r[2];
        /* A square-root choice on a projective representative must not
         * change the final curve model when the traversal changes. */
        if (exact && step>=normal) {
            if (!normalize_point(&A.null_point)) return 0;
            A.precomputation=0;
        }
        theta_precomputation(&A);
        if (exact && step==len2-3)
            theta_isogeny_comput4(&c->two[step],&A,&q2[0],&q2[1],0,0);
        else if (exact && step==len2-2)
            theta_isogeny_comput2(&c->two[step],&A,&q2[0],&q2[1],1,0);
        else {
            for (int j=0;j<2;++j) {
                int doublings=len2-2-step-2*exact;
                double_iter(&r[j],&A,&q2[j],doublings);
                c->kernel_doublings+=(unsigned long)doublings;
            }
            theta_isogeny_comput(&c->two[step],&A,&r[0],&r[1],
                !len3 && step==len2-2,len3 || step<len2-3);
        }
        theta_point_t origin;
        theta_isogeny_eval(&origin,&c->two[step],&A.null_point);
        if (!inke_theta_equal(&origin,&c->two[step].codomain.null_point)) return 0;
        for (int j=0;j<2 && !(exact && step==len2-2);++j) {
            theta_isogeny_eval(&q2[j],&c->two[step],&q2[j]);
            if (!valid_point(&q2[j])) return 0;
        }
        if (len3) for (int j=0;j<4;++j) {
            theta_isogeny_eval(&q3[j],&c->two[step],&q3[j]);
            if (!valid_point(&q3[j])) return 0;
        }
        A=c->two[step].codomain;
        if ((len3 || step<len2-3) && !valid_surface(&A)) return 0;
    }
    if (balanced && len3) {
        c->failed_stage=3;
        if (!balanced_three(c,&A,q3,len3,0)) return 0;
    }
    for (int step=0;!balanced && step<len3;++step) {
        c->failed_stage=3;
        theta_point_t r[4];
        for (int j=0;j<4;++j) {
            r[j]=q3[j];
            for (int k=step+1;k<len3;++k) {
                if (!inke_theta_triple(&r[j],&A,&r[j])) return 0;
                ++c->kernel_triplings;
            }
        }
        int found=0;
        for (int j=0;j<4 && !found;++j) for (int k=j+1;k<4 && !found;++k)
            found=inke_theta_isogeny3_compute(&c->three[step],&A,&r[j],&r[k]);
        if (!found) return 0;
        for (int j=0;j<4;++j) {
            theta_point_t killed;
            if (!inke_theta_isogeny3_eval(&killed,&c->three[step],&r[j]) ||
                !inke_theta_equal(&killed,&c->three[step].codomain.null_point) ||
                !inke_theta_isogeny3_eval(&q3[j],&c->three[step],&q3[j])) return 0;
        }
        A=c->three[step].codomain;
    }
    c->failed_stage=4;
    if (!splitting_comput(&c->splitting,&A)) return 0;
    theta_product_structure_to_elliptic_product(&c->codomain,&c->splitting.B);
    if (fp2_is_zero(&c->codomain.E1.C) || fp2_is_zero(&c->codomain.E2.C)) return 0;
    if (lift_image && !split_image(lift_image,c,q2[0])) return 0;
    c->ready=1; c->failed_stage=0;
    return 1;
}
int inke_chain_compute(inke_product_chain *c,const theta_couple_curve_t *E,
    const theta_couple_jac_point_t K2[2],int len2,
    const theta_couple_jac_point_t K3[2],int len3)
{ return compute(c,E,K2,len2,K3,len3,0,1,NULL); }
int inke_chain_compute_with_lift(inke_product_chain *c,const theta_couple_curve_t *E,
    const theta_couple_jac_point_t lifted[2],int len2,theta_couple_point_t *lift_image)
{ return compute(c,E,lifted,len2,NULL,0,0,1,lift_image); }
int inke_chain_compute_exact(inke_product_chain *c,const theta_couple_curve_t *E,
    const theta_couple_jac_point_t K2[2],int len2)
{ return compute(c,E,K2,len2,NULL,0,1,1,NULL); }
int inke_chain_compute_reference(inke_product_chain *c,const theta_couple_curve_t *E,
    const theta_couple_jac_point_t K2[2],int len2,
    const theta_couple_jac_point_t K3[2],int len3)
{ return compute(c,E,K2,len2,K3,len3,0,0,NULL); }
int inke_chain_compute_exact_reference(inke_product_chain *c,const theta_couple_curve_t *E,
    const theta_couple_jac_point_t K2[2],int len2)
{ return compute(c,E,K2,len2,NULL,0,1,0,NULL); }

static int finish(theta_couple_point_t *out,const inke_product_chain *c,theta_point_t P)
{
    for (int k=0;k<c->len2-1;++k) {
        theta_isogeny_eval(&P,&c->two[k],&P);
        if (!valid_point(&P)) return 0;
    }
    for (int k=0;k<c->len3;++k)
        if (!inke_theta_isogeny3_eval(&P,&c->three[k],&P)) return 0;
    return split_image(out,c,P);
}
static int split_image(theta_couple_point_t *out,const inke_product_chain *c,theta_point_t P)
{
    apply_isomorphism(&P,&c->splitting,&P);
    const theta_point_t *O=&c->splitting.B.null_point;
    fp2_t u,v;
    /* Two Segre coordinate charts per factor, including the chart needed
     * when a 2-torsion point has a whole zero row/column. */
    if (!fp2_is_zero(&P.x) || !fp2_is_zero(&P.y)) {
        fp2_mul(&u,&O->y,&P.x); fp2_mul(&v,&O->x,&P.y);
    } else { fp2_mul(&u,&O->t,&P.z); fp2_mul(&v,&O->z,&P.t); }
    fp2_add(&out->P1.x,&u,&v); fp2_sub(&out->P1.z,&v,&u);
    if (!fp2_is_zero(&P.y) || !fp2_is_zero(&P.t)) {
        fp2_mul(&u,&O->t,&P.y); fp2_mul(&v,&O->y,&P.t);
    } else { fp2_mul(&u,&O->z,&P.x); fp2_mul(&v,&O->x,&P.z); }
    fp2_add(&out->P2.x,&u,&v); fp2_sub(&out->P2.z,&v,&u);
    return !(fp2_is_zero(&out->P1.x) && fp2_is_zero(&out->P1.z)) &&
           !(fp2_is_zero(&out->P2.x) && fp2_is_zero(&out->P2.z));
}
int inke_chain_eval(theta_couple_point_t *out,const inke_product_chain *c,
    const theta_couple_jac_point_t *P)
{
    if (!c->ready) return 0;
    theta_point_t image; theta_couple_jac_point_t transformed=*P;
    for (int k=0;k<c->prefix_length;++k)
        product_step(&transformed,&transformed,&c->domain,c->product_prefix[k]);
    gluing_eval_point_no_help(&image,&transformed,&c->domain,&c->glue);
    return valid_point(&image) && finish(out,c,image);
}
int inke_chain_eval_one(theta_couple_point_t *out,const inke_product_chain *c,
    const ec_point_t *P,int side)
{
    if (!c->ready || (side!=0 && side!=1)) return 0;
    if (c->prefix_length) {
        theta_couple_jac_point_t pair;
        ec_point_t point=*P; theta_couple_curve_t curves=c->domain;
        /* ADD uses the upstream marker (0:1:0) for the identity. */
        fp2_set_zero(&pair.P1.x); fp2_set_one(&pair.P1.y); fp2_set_zero(&pair.P1.z);
        pair.P2=pair.P1;
        if (side) lift_point(&pair.P2,&point,&curves.E2);
        else lift_point(&pair.P1,&point,&curves.E1);
        return inke_chain_eval(out,c,&pair);
    }
    theta_point_t image;
    theta_couple_point_t pair;
    ec_set_zero(&pair.P1); ec_set_zero(&pair.P2);
    if (side) pair.P2=*P; else pair.P1=*P;
    gluing_eval_point_special_case(&image,&pair,&c->glue);
    return valid_point(&image) && finish(out,c,image);
}
