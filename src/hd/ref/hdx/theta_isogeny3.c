// SPDX-License-Identifier: Apache-2.0
#include "theta_isogeny3.h"
#include <string.h>

static void unpack(fp2_t p[4], const theta_point_t *P)
{
    p[0]=P->x; p[1]=P->y; p[2]=P->z; p[3]=P->t;
}

static int nonzero(const fp2_t p[4])
{
    for (int i=0;i<4;++i) if (!fp2_is_zero(&p[i])) return 1;
    return 0;
}

int inke_theta_equal(const theta_point_t *P, const theta_point_t *Q)
{
    fp2_t p[4],q[4],l,r;
    unpack(p,P); unpack(q,Q);
    if (!nonzero(p) || !nonzero(q)) return 0;
    for (int i=0;i<4;++i) for (int j=i+1;j<4;++j) {
        fp2_mul(&l,&p[i],&q[j]); fp2_mul(&r,&p[j],&q[i]);
        if (!fp2_is_equal(&l,&r)) return 0;
    }
    return 1;
}

/* The paper's Hadamard rows 2 and 3 are interchanged relative to POKE.
 * Keep the paper's h_i ordering explicitly for its coefficient formulas. */
static void squared_hadamard(fp2_t h[4], const fp2_t p[4])
{
    fp2_t sq[4],u,v,w,z;
    for (int i=0;i<4;++i) fp2_sqr(&sq[i],&p[i]);
    fp2_add(&u,&sq[0],&sq[1]); fp2_sub(&v,&sq[0],&sq[1]);
    fp2_add(&w,&sq[2],&sq[3]); fp2_sub(&z,&sq[2],&sq[3]);
    fp2_add(&h[0],&u,&w); fp2_sub(&h[1],&u,&w);
    fp2_add(&h[2],&v,&z); fp2_sub(&h[3],&v,&z);
}

/* Diagonal biquadratic coordinates H(H(P^2)/H(O^2)). The accompanying
 * reference implementation uses these coordinates in the coefficients;
 * H(P^2) alone, as written in the article's h_i definition, is insufficient.
 * Use one common projective denominator, avoiding inversions. */
static void biquadratic_coordinates(fp2_t h[4], const fp2_t p[4],
    const fp2_t h0[4])
{
    fp2_t hscaled[4],product,u,v,w,z;
    squared_hadamard(hscaled,p);
    for (int i=0;i<4;++i) {
        fp2_set_one(&product);
        for (int j=0;j<4;++j) if (i!=j) fp2_mul(&product,&product,&h0[j]);
        fp2_mul(&hscaled[i],&hscaled[i],&product);
    }
    fp2_add(&u,&hscaled[0],&hscaled[1]); fp2_sub(&v,&hscaled[0],&hscaled[1]);
    fp2_add(&w,&hscaled[2],&hscaled[3]); fp2_sub(&z,&hscaled[2],&hscaled[3]);
    fp2_add(&h[0],&u,&w); fp2_sub(&h[1],&u,&w);
    fp2_add(&h[2],&v,&z); fp2_sub(&h[3],&v,&z);
}

static void gamma_and_f(fp2_t *gamma, fp2_t f[3], const fp2_t p[4],
    const fp2_t pair[6], const fp2_t dd[3])
{
    fp2_t u,v,l,r,z;
    fp2_set_zero(gamma);
    for (int k=0;k<3;++k) {
        /* Complementary pairs: 01|23, 02|13, 03|12. */
        int j=k+1, m=(j==1?2:1), n=6-j-m;
        fp2_mul(&u,&p[0],&p[j]); fp2_mul(&v,&p[m],&p[n]);
        fp2_mul(&l,&u,&pair[k]); fp2_mul(&r,&v,&pair[k+3]);
        fp2_sub(&z,&l,&r); fp2_mul(&z,&z,&dd[k]);
        fp2_add(gamma,gamma,&z);
        fp2_mul(&l,&v,&pair[k]); fp2_mul(&r,&u,&pair[k+3]);
        fp2_sub(&f[k],&l,&r);
    }
}

int inke_theta_triple(theta_point_t *out, theta_structure_t *domain,
    const theta_point_t *P)
{
    /* B_ii(2P,P) is proportional to (3P)_i P_i, with a common factor.
     * Projective products clear the denominators without inversion. */
    theta_point_t twice, hp,hq,h0,result;
    fp2_t p[4],dual[4],prod;
    unpack(p,P);
    to_squared_theta(&h0,&domain->null_point); unpack(dual,&h0);
    for (int i=0;i<4;++i)
        if (fp2_is_zero(&p[i]) || fp2_is_zero(&dual[i])) return 0;
    double_point(&twice,domain,P);
    to_squared_theta(&hp,P); to_squared_theta(&hq,&twice);
    fp2_t x[4],y[4]; unpack(x,&hp); unpack(y,&hq);
    for (int i=0;i<4;++i) {
        fp2_set_one(&prod);
        for (int j=0;j<4;++j) if (i!=j) fp2_mul(&prod,&prod,&dual[j]);
        fp2_mul(&x[i],&x[i],&y[i]); fp2_mul(&x[i],&x[i],&prod);
    }
    result.x=x[0]; result.y=x[1]; result.z=x[2]; result.t=x[3];
    hadamard(&result,&result); unpack(x,&result);
    for (int i=0;i<4;++i) {
        fp2_set_one(&prod);
        for (int j=0;j<4;++j) if (i!=j) fp2_mul(&prod,&prod,&p[j]);
        fp2_mul(&x[i],&x[i],&prod);
    }
    if (!nonzero(x)) return 0;
    out->x=x[0]; out->y=x[1]; out->z=x[2]; out->t=x[3];
    return 1;
}

int inke_theta_isogeny3_eval(theta_point_t *out,
    const inke_theta_isogeny3 *phi, const theta_point_t *P)
{
    if (!phi->ready) return 0;
    fp2_t p[4],sq[4],v[4],sum,t,product;
    unpack(p,P);
    if (!nonzero(p)) return 0;
    for (int j=0;j<4;++j) fp2_sqr(&sq[j],&p[j]);
    for (int i=0;i<4;++i) {
        fp2_set_zero(&sum); fp2_set_one(&product);
        for (int j=0;j<4;++j) {
            fp2_mul(&t,&phi->c[i^j],&sq[j]); fp2_add(&sum,&sum,&t);
            if (j!=i) fp2_mul(&product,&product,&p[j]);
        }
        fp2_mul(&sum,&sum,&p[i]); fp2_mul(&product,&product,&phi->c[4]);
        fp2_add(&v[i],&sum,&product);
    }
    if (!nonzero(v)) return 0;
    out->x=v[0]; out->y=v[1]; out->z=v[2]; out->t=v[3];
    return 1;
}

int inke_theta_isogeny3_compute(inke_theta_isogeny3 *out,
    const theta_structure_t *domain, const theta_point_t *R,
    const theta_point_t *S)
{
    memset(out,0,sizeof(*out));
    fp2_t o[4],r[4],s[4],pair[6],d[3],dd[3],hr[4],hs[4],h0[4];
    fp2_t gr,gs,fr[3],fs[3],beta1,beta2,u,v,w;
    unpack(o,&domain->null_point); unpack(r,R); unpack(s,S);
    squared_hadamard(h0,o);
    for (int i=0;i<4;++i)
        if (fp2_is_zero(&o[i]) || fp2_is_zero(&h0[i])) return 0;
    if (!nonzero(r) || !nonzero(s) || inke_theta_equal(R,S) ||
        inke_theta_equal(R,&domain->null_point) ||
        inke_theta_equal(S,&domain->null_point)) return 0;
    theta_structure_t A=*domain;
    theta_point_t doubled;
    double_point(&doubled,&A,R);
    if (!inke_theta_equal(&doubled,R)) return 0;
    double_point(&doubled,&A,S);
    if (!inke_theta_equal(&doubled,S)) return 0;
    for (int k=0;k<3;++k) {
        int j=k+1, m=(j==1?2:1), n=6-j-m;
        fp2_mul(&pair[k],&o[0],&o[j]); fp2_mul(&pair[k+3],&o[m],&o[n]);
        fp2_sqr(&u,&pair[k]); fp2_sqr(&v,&pair[k+3]);
        fp2_sub(&d[k],&u,&v);
        if (fp2_is_zero(&d[k])) return 0;
    }
    fp2_mul(&dd[0],&d[1],&d[2]); fp2_mul(&dd[1],&d[0],&d[2]);
    fp2_mul(&dd[2],&d[0],&d[1]);
    biquadratic_coordinates(hr,r,h0); biquadratic_coordinates(hs,s,h0);
    gamma_and_f(&gr,fr,r,pair,dd); gamma_and_f(&gs,fs,s,pair,dd);
    fp2_mul(&u,&gr,&fs[0]); fp2_mul(&v,&gs,&fr[0]);
    fp2_sub(&beta1,&u,&v); fp2_mul(&beta1,&beta1,&dd[0]);
    fp2_mul(&u,&hr[0],&hs[1]); fp2_mul(&v,&hr[1],&hs[0]);
    /* Sign agrees with the diagonal/off-diagonal biquadratic convention. */
    fp2_sub(&beta2,&v,&u);
    fp2_mul(&u,&hr[0],&hs[0]); fp2_mul(&u,&u,&beta1);
    fp2_add(&out->c[0],&u,&u);
    for (int k=0;k<3;++k) {
        fp2_mul(&u,&hr[0],&hs[k+1]); fp2_mul(&v,&hr[k+1],&hs[0]);
        fp2_add(&u,&u,&v); fp2_mul(&u,&u,&beta1);
        fp2_mul(&v,&gr,&fs[k]); fp2_mul(&w,&gs,&fr[k]);
        fp2_add(&v,&v,&w); fp2_mul(&v,&v,&beta2); fp2_mul(&v,&v,&dd[k]);
        fp2_add(&out->c[k+1],&u,&v);
    }
    fp2_mul(&u,&beta2,&gr); fp2_mul(&u,&u,&gs); fp2_add(&out->c[4],&u,&u);
    out->ready=1;
    if (!inke_theta_isogeny3_eval(&out->codomain.null_point,out,&domain->null_point))
        goto reject;
    /* Exact necessary checks also reject a vanishing coefficient chart. */
    if (!inke_theta_isogeny3_eval(&doubled,out,R) ||
        !inke_theta_equal(&doubled,&out->codomain.null_point)) goto reject;
    if (!inke_theta_isogeny3_eval(&doubled,out,S) ||
        !inke_theta_equal(&doubled,&out->codomain.null_point)) goto reject;
    return 1;
reject:
    out->ready=0;
    return 0;
}
