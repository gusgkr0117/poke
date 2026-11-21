/* Exact bounded prefix using the linked INKE E0 parameters. The static
 * lattice and two-square routines below are copied from the tested POKE C
 * implementation; no runtime script, Sage, Rust or Python is required. */
#include "b-gaussian-preparation.h"
#include <inke_mixed.h>
#include "../../../../quaternion/ref/generic/internal.h"
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
static void elem_copy(quat_alg_elem_t *out, const quat_alg_elem_t *in)
{
    ibz_copy(&out->denom, &in->denom);
    for (int i=0; i<4; ++i) ibz_copy(&out->coord[i], &in->coord[i]);
}

static int elem_equal(const quat_alg_elem_t *a, const quat_alg_elem_t *b)
{
    ibz_t left, right;
    ibz_init(&left); ibz_init(&right);
    int equal = 1;
    for (int i=0; i<4 && equal; ++i) {
        mpz_mul(left,a->coord[i],b->denom);
        mpz_mul(right,b->coord[i],a->denom);
        equal = mpz_cmp(left,right)==0;
    }
    ibz_finalize(&right); ibz_finalize(&left);
    return equal;
}

static void ideal_copy(quat_left_ideal_t *out, const quat_left_ideal_t *in)
{
    ibz_mat_4x4_copy(&out->lattice.basis, &in->lattice.basis);
    ibz_copy(&out->lattice.denom, &in->lattice.denom);
    ibz_copy(&out->norm, &in->norm);
    out->parent_order = in->parent_order;
}

static int norm_integer(ibz_t *norm, const quat_alg_elem_t *x, const quat_alg_t *alg)
{
    ibq_t rational;
    ibq_init(&rational);
    quat_alg_norm(&rational, x, alg);
    int ok = ibq_to_ibz(norm, &rational);
    ibq_finalize(&rational);
    return ok;
}

/* Direct multiplication of a lattice by a quaternion. Avoids the random
 * generator search in the older POKE quat_lideal_mul() helper. */
static void multiply_lattice(quat_lattice_t *out, const quat_lattice_t *in,
                              const quat_alg_elem_t *right, const quat_alg_t *alg)
{
    ibz_mat_4x4_t matrix;
    ibz_mat_4x4_init(&matrix);
    quat_alg_rightmul_mat(&matrix, right, alg);
    ibz_mat_4x4_mul(&out->basis, &matrix, &in->basis);
    ibz_mul(&out->denom, &in->denom, &right->denom);
    quat_lattice_reduce_denom(out, out);
    quat_lattice_hnf(out);
    ibz_mat_4x4_finalize(&matrix);
}

static int integral_lattice(const quat_lattice_t *lattice,
                            const quat_order_t *order, const quat_alg_t *alg)
{
    quat_alg_elem_t element;
    quat_alg_elem_init(&element);
    ibz_copy(&element.denom, &lattice->denom);
    int ok = 1;
    for (int j=0; j<4 && ok; ++j) {
        for (int i=0; i<4; ++i) ibz_copy(&element.coord[i], &lattice->basis[i][j]);
        ok = quat_lattice_contains(NULL, order, &element, alg);
    }
    quat_alg_elem_finalize(&element);
    return ok;
}


/* Euclidean Cornacchia with at most 198 square-root candidates. Even if a
 * probable-prime test accepts a composite, acceptance requires x^2+y^2=n. */
static int prime_squares(ibz_t *x, ibz_t *y, const ibz_t *n)
{
    if (mpz_cmp_ui(*n, 2) == 0) { mpz_set_ui(*x,1); mpz_set_ui(*y,1); return 1; }
    if (mpz_cmp_ui(*n, 2) < 0 || mpz_fdiv_ui(*n, 4) != 1 ||
        !mpz_probab_prime_p(*n, 32)) return 0;
    ibz_t exponent, base, u, v, t;
    ibz_init(&exponent); ibz_init(&base); ibz_init(&u); ibz_init(&v); ibz_init(&t);
    mpz_sub_ui(exponent, *n, 1); mpz_fdiv_q_2exp(exponent, exponent, 2);
    int ok = 0;
    for (unsigned long z=2; z<200; ++z) {
        mpz_set_ui(base,z); mpz_powm(v,base,exponent,*n);
        mpz_mul(t,v,v); mpz_add_ui(t,t,1); mpz_mod(t,t,*n);
        if (mpz_sgn(t) == 0) { ok = 1; break; }
    }
    if (!ok) goto done;
    mpz_set(u,*n); mpz_mul(t,v,v);
    while (mpz_cmp(t,*n) > 0) {
        mpz_mod(t,u,v); mpz_set(u,v); mpz_set(v,t); mpz_mul(t,v,v);
    }
    mpz_sub(t,*n,t); mpz_sqrt(*y,t); mpz_mul(u,*y,*y);
    ok = mpz_cmp(u,t) == 0;
    if (ok) mpz_set(*x,v);
done:
    ibz_finalize(&t); ibz_finalize(&v); ibz_finalize(&u);
    ibz_finalize(&base); ibz_finalize(&exponent);
    return ok;
}

static void complex_multiply(ibz_t *a, ibz_t *b, const ibz_t *x, const ibz_t *y)
{
    ibz_t re, im, tmp;
    ibz_init(&re); ibz_init(&im); ibz_init(&tmp);
    mpz_mul(re,*a,*x); mpz_mul(tmp,*b,*y); mpz_sub(re,re,tmp);
    mpz_mul(im,*a,*y); mpz_mul(tmp,*b,*x); mpz_add(im,im,tmp);
    mpz_set(*a,re); mpz_set(*b,im);
    ibz_finalize(&tmp); ibz_finalize(&im); ibz_finalize(&re);
}

static int two_squares(ibz_t *x, ibz_t *y, const ibz_t *value)
{
    static const unsigned short primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,
        47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,
        151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,
        251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,
        359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,
        461,463,467,479,487,491,499};
    if (mpz_sgn(*value) <= 0) return 0;
    ibz_t n, prime, scale, u, v, check;
    ibz_init(&n); ibz_init(&prime); ibz_init(&scale);
    ibz_init(&u); ibz_init(&v); ibz_init(&check);
    mpz_set(n,*value); mpz_set_ui(*x,1); mpz_set_ui(*y,0);
    int ok = 0;
    for (unsigned int i=0; i<sizeof(primes)/sizeof(primes[0]); ++i) {
        unsigned long ell=primes[i];
        if (!mpz_divisible_ui_p(n,ell)) continue;
        mpz_set_ui(prime,ell);
        mp_bitcnt_t exponent = mpz_remove(n,n,prime);
        if (ell%4 == 3 && (exponent&1)) goto done;
        mpz_pow_ui(scale,prime,(unsigned long)(exponent/2));
        mpz_mul(*x,*x,scale); mpz_mul(*y,*y,scale);
        if (exponent&1) {
            if (!prime_squares(&u,&v,&prime)) goto done;
            complex_multiply(x,y,&u,&v);
        }
    }
    if (mpz_cmp_ui(n,1) != 0) {
        if (!prime_squares(&u,&v,&n)) goto done;
        complex_multiply(x,y,&u,&v);
    }
    mpz_mul(check,*x,*x); mpz_mul(u,*y,*y); mpz_add(check,check,u);
    ok = mpz_cmp(check,*value) == 0;
    mpz_abs(*x,*x); mpz_abs(*y,*y);
done:
    ibz_finalize(&check); ibz_finalize(&v); ibz_finalize(&u);
    ibz_finalize(&scale); ibz_finalize(&prime); ibz_finalize(&n);
    return ok;
}


static void vector_element(quat_alg_elem_t *out,const ibz_mat_4x4_t *basis,const ibz_t *den,const ibz_vec_4_t *c){ibz_copy(&out->denom,den);for(unsigned i=0;i<4;i++){mpz_set_ui(out->coord[i],0);for(unsigned j=0;j<4;j++)mpz_addmul(out->coord[i],(*basis)[i][j],(*c)[j]);}quat_alg_normalize(out);}
static int unit_at(ibz_t *d,const quat_alg_elem_t *beta,const quat_left_ideal_t *I,const ibz_t *mod){if(!norm_integer(d,beta,&QUATALG_PINFTY)||!mpz_divisible_p(*d,I->norm))return 0;mpz_divexact(*d,*d,I->norm);if(mpz_sgn(*d)<=0)return 0;ibz_t g;ibz_init(&g);mpz_gcd(g,*d,*mod);int ok=mpz_cmp_ui(g,1)==0;ibz_finalize(&g);return ok;}
typedef unsigned __int128 a_plan_cost;
typedef struct {ibz_t d; int c[4];} b_candidate;
typedef struct {ibz_t N; unsigned e[11],k[11]; a_plan_cost cost;} b_plan;
static const unsigned b_primes[11]={2,3,5,7,11,13,17,19,23,29,31};

static unsigned order_bound(const mpz_t p,unsigned ell,unsigned e,unsigned cap){
    mpz_t mod,a,x;mpz_inits(mod,a,x,NULL);mpz_ui_pow_ui(mod,ell,e);mpz_neg(a,p);mpz_mod(a,a,mod);mpz_set_ui(x,1);
    unsigned result=0;for(unsigned k=1;k<=cap;k++){mpz_mul(x,x,a);mpz_mod(x,x,mod);if(mpz_cmp_ui(x,1)==0){result=k;break;}}
    mpz_clears(mod,a,x,NULL);return result;
}

static int b_by_norm(const void*a,const void*b){int c=mpz_cmp(((const b_candidate*)a)->d,((const b_candidate*)b)->d);return c?c:memcmp(((const b_candidate*)a)->c,((const b_candidate*)b)->c,4*sizeof(int));}
static int b_by_cost(const void*a,const void*b){const b_plan*x=a,*y=b;if(x->cost!=y->cost)return x->cost<y->cost?-1:1;return mpz_cmp(x->N,y->N);}
static unsigned b_candidates(b_candidate*out,const quat_left_ideal_t*I,const ibz_mat_4x4_t*basis){
    unsigned n=0;quat_alg_elem_t beta;quat_alg_elem_init(&beta);ibz_vec_4_t c;ibz_vec_4_init(&c);ibz_t d,mod;ibz_init(&d);ibz_init(&mod);mpz_mul_ui(mod,QUATALG_PINFTY.p,6);
    for(unsigned word=0;word<625;word++){unsigned t=word;int values[4];for(unsigned j=0;j<4;j++){values[j]=(int)(t%5)-2;t/=5;mpz_set_si(c[j],values[j]);}vector_element(&beta,basis,&I->lattice.denom,&c);if(!unit_at(&d,&beta,I,&mod))continue;int seen=0;for(unsigned j=0;j<n;j++)if(!mpz_cmp(d,out[j].d)){seen=1;break;}if(seen)continue;ibz_init(&out[n].d);mpz_set(out[n].d,d);memcpy(out[n].c,values,sizeof(values));++n;}
    qsort(out,n,sizeof(*out),b_by_norm);ibz_finalize(&mod);ibz_finalize(&d);ibz_vec_4_finalize(&c);quat_alg_elem_finalize(&beta);return n;
}
static unsigned b_plans(b_plan*out,unsigned base3){
    unsigned n=0;const unsigned maxima[11]={7,5,3,2,2,1,1,1,1,1,1};ibz_t temp;ibz_init(&temp);
    for(unsigned word=0;word<8*6*4*3*3*64;word++){unsigned t=word,e[11],k[11];int ok=1;a_plan_cost cost=0;for(unsigned j=0;j<11;j++){e[j]=t%(maxima[j]+1);t/=maxima[j]+1;if(j==0)e[j]+=TORSION_PLUS_EVEN_POWER-2;if(j==1)e[j]+=base3;k[j]=e[j]?order_bound(QUATALG_PINFTY.p,b_primes[j],e[j]+(j==0?3:0),64):1;if(!k[j]){ok=0;break;}a_plan_cost z=4*b_primes[j];z=z*z;cost+=(a_plan_cost)e[j]*k[j]*k[j]*z;}if(!ok)continue;unsigned common=1;for(unsigned j=0;j<11;j++){unsigned a=common,b=k[j];while(b){unsigned r=a%b;a=b;b=r;}common=(common/a)*k[j];if(common>64){ok=0;break;}}if(!ok)continue;cost=0;for(unsigned j=0;j<11;j++)if(e[j]){a_plan_cost weight=1;for(unsigned l=0;l<11;l++)if(e[l]){unsigned a=k[j],b=k[l];while(b){unsigned r=a%b;a=b;b=r;}unsigned pair=(k[j]/a)*k[l];weight+=(a_plan_cost)pair*pair;}a_plan_cost z=4*b_primes[j];cost+=(a_plan_cost)e[j]*z*z*weight;}
        ibz_init(&out[n].N);mpz_set_ui(out[n].N,1);for(unsigned j=0;j<11;j++){mpz_ui_pow_ui(temp,b_primes[j],e[j]);mpz_mul(out[n].N,out[n].N,temp);out[n].e[j]=e[j];out[n].k[j]=k[j];}out[n++].cost=cost;
    }qsort(out,n,sizeof(*out),b_by_cost);ibz_finalize(&temp);return n;
}
static void b_set_beta(quat_alg_elem_t*out,const b_candidate*c,const quat_left_ideal_t*I,const ibz_mat_4x4_t*basis){ibz_vec_4_t coeff;ibz_vec_4_init(&coeff);for(unsigned j=0;j<4;j++)mpz_set_si(coeff[j],c->c[j]);vector_element(out,basis,&I->lattice.denom,&coeff);ibz_vec_4_finalize(&coeff);}

int b_gaussian_plans_init(b_gaussian_plans *P){
    if(!P||P->ready)return B_GAUSSIAN_INVALID;
    b_plan *plans=calloc(110592,sizeof(*plans));if(!plans)return B_GAUSSIAN_NO_MEMORY;
    unsigned b=0;ibz_t n;ibz_init(&n);mpz_set(n,TORSION_PLUS_3POWER);
    while(mpz_divisible_ui_p(n,3)&&mpz_sgn(n)>0){mpz_divexact_ui(n,n,3);b++;}
    int valid=mpz_cmp_ui(n,1)==0;ibz_finalize(&n);if(!valid){free(plans);return B_GAUSSIAN_INVALID;}
    P->count=b_plans(plans,b);P->plans=plans;P->ready=1;return B_GAUSSIAN_OK;
}
void b_gaussian_plans_clear(b_gaussian_plans *P){if(!P)return;if(P->ready){b_plan*plans=P->plans;for(unsigned i=0;i<P->count;i++)ibz_finalize(&plans[i].N);free(plans);}memset(P,0,sizeof(*P));}
void b_gaussian_certificate_init(b_gaussian_certificate *C){
    memset(C,0,sizeof(*C));quat_alg_elem_init(&C->beta1);quat_alg_elem_init(&C->beta2);quat_alg_elem_init(&C->x);quat_alg_elem_init(&C->y);quat_alg_elem_init(&C->eta);quat_alg_elem_init(&C->transport);quat_left_ideal_init(&C->Ia);ibz_init(&C->d1);ibz_init(&C->d2);ibz_init(&C->u);ibz_init(&C->v);ibz_init(&C->N);ibz_init(&C->n1);ibz_init(&C->n2);
}
void b_gaussian_certificate_clear(b_gaussian_certificate *C){
    ibz_finalize(&C->n2);ibz_finalize(&C->n1);ibz_finalize(&C->N);ibz_finalize(&C->v);ibz_finalize(&C->u);ibz_finalize(&C->d2);ibz_finalize(&C->d1);quat_left_ideal_finalize(&C->Ia);quat_alg_elem_finalize(&C->transport);quat_alg_elem_finalize(&C->eta);quat_alg_elem_finalize(&C->y);quat_alg_elem_finalize(&C->x);quat_alg_elem_finalize(&C->beta2);quat_alg_elem_finalize(&C->beta1);memset(C,0,sizeof(*C));
}
static int b_input_valid(const quat_left_ideal_t *I){
    if(!I||I->parent_order!=&MAXORD_O0||mpz_sgn(I->norm)<=0||mpz_sgn(I->lattice.denom)<=0||!integral_lattice(&I->lattice,I->parent_order,&QUATALG_PINFTY))return 0;
    ibz_t index,norm;ibz_init(&index);ibz_init(&norm);quat_lattice_index(&index,&I->lattice,I->parent_order);mpz_mul(norm,I->norm,I->norm);int ok=mpz_cmp(index,norm)==0;mpz_gcd(index,I->norm,QUATALG_PINFTY.p);ok=ok&&mpz_cmp_ui(index,1)==0&&mpz_odd_p(I->norm)&&!mpz_divisible_ui_p(I->norm,3);ibz_finalize(&norm);ibz_finalize(&index);return ok;
}
static void b_assemble(b_gaussian_certificate *C,const quat_left_ideal_t *I){
    quat_alg_elem_t bar,theta,tmp;quat_alg_elem_init(&bar);quat_alg_elem_init(&theta);quat_alg_elem_init(&tmp);
    quat_alg_conj(&bar,&C->beta1);quat_alg_mul(&theta,&C->beta2,&bar,&QUATALG_PINFTY);mpz_mul(theta.denom,theta.denom,I->norm);quat_alg_normalize(&theta);
    quat_alg_mul(&tmp,&theta,&C->x,&QUATALG_PINFTY);quat_alg_mul(&C->eta,&C->y,&tmp,&QUATALG_PINFTY);quat_alg_normalize(&C->eta);
    quat_alg_conj(&bar,&C->x);quat_alg_mul(&C->transport,&bar,&C->beta1,&QUATALG_PINFTY);quat_alg_normalize(&C->transport);
    quat_lattice_t J;quat_lattice_init(&J);quat_alg_conj(&bar,&C->beta1);mpz_mul(bar.denom,bar.denom,I->norm);multiply_lattice(&J,&I->lattice,&bar,&QUATALG_PINFTY);multiply_lattice(&C->Ia.lattice,&J,&C->x,&QUATALG_PINFTY);mpz_set(C->Ia.norm,C->n1);C->Ia.parent_order=I->parent_order;quat_lattice_finalize(&J);
    quat_alg_elem_finalize(&tmp);quat_alg_elem_finalize(&theta);quat_alg_elem_finalize(&bar);
}
int b_gaussian_certificate_verify(const b_gaussian_certificate *C,const quat_left_ideal_t *I){
    if(!C||!b_input_valid(I)||C->Ia.parent_order!=I->parent_order||mpz_sgn(C->Ia.lattice.denom)<=0||mpz_sgn(C->d1)<=0||mpz_sgn(C->d2)<=0||mpz_sgn(C->u)<=0||mpz_sgn(C->v)<=0||mpz_sgn(C->N)<=0||mpz_sgn(C->n1)<=0||mpz_sgn(C->n2)<=0)return 0;
    const quat_alg_elem_t *elems[6]={&C->beta1,&C->beta2,&C->x,&C->y,&C->eta,&C->transport};for(unsigned i=0;i<6;i++)if(mpz_sgn(elems[i]->denom)<=0)return 0;
    int ok=0;ibz_t z,w,g;ibz_init(&z);ibz_init(&w);ibz_init(&g);b_gaussian_certificate D;b_gaussian_certificate_init(&D);quat_lattice_t back,target;quat_lattice_init(&back);quat_lattice_init(&target);
#define B_VERIFY(q) do{if(!(q))goto finish;}while(0)
    B_VERIFY(quat_lattice_contains(NULL,&I->lattice,&C->beta1,&QUATALG_PINFTY)&&quat_lattice_contains(NULL,&I->lattice,&C->beta2,&QUATALG_PINFTY));
    mpz_mul(w,C->d1,I->norm);B_VERIFY(norm_integer(&z,&C->beta1,&QUATALG_PINFTY)&&!mpz_cmp(z,w));mpz_mul(w,C->d2,I->norm);B_VERIFY(norm_integer(&z,&C->beta2,&QUATALG_PINFTY)&&!mpz_cmp(z,w));
    B_VERIFY(!mpz_cmp_ui(C->x.denom,1)&&!mpz_cmp_ui(C->y.denom,1)&&!mpz_sgn(C->x.coord[2])&&!mpz_sgn(C->x.coord[3])&&!mpz_sgn(C->y.coord[2])&&!mpz_sgn(C->y.coord[3]));
    B_VERIFY(norm_integer(&z,&C->x,&QUATALG_PINFTY)&&!mpz_cmp(z,C->u));B_VERIFY(norm_integer(&z,&C->y,&QUATALG_PINFTY)&&!mpz_cmp(z,C->v));
    mpz_mul(z,C->u,C->d1);B_VERIFY(!mpz_cmp(z,C->n1));mpz_mul(z,C->v,C->d2);B_VERIFY(!mpz_cmp(z,C->n2));mpz_add(z,C->n1,C->n2);B_VERIFY(!mpz_cmp(z,C->N));mpz_gcd(g,C->n1,C->n2);B_VERIFY(!mpz_cmp_ui(g,1));mpz_mul(z,C->n1,C->n2);B_VERIFY(mpz_odd_p(z)&&!mpz_divisible_ui_p(z,3));mpz_gcd(g,C->N,QUATALG_PINFTY.p);B_VERIFY(!mpz_cmp_ui(g,1));
    mpz_set(D.n1,C->n1);elem_copy(&D.beta1,&C->beta1);elem_copy(&D.beta2,&C->beta2);elem_copy(&D.x,&C->x);elem_copy(&D.y,&C->y);b_assemble(&D,I);B_VERIFY(elem_equal(&D.eta,&C->eta)&&elem_equal(&D.transport,&C->transport)&&quat_lattice_equal(&D.Ia.lattice,&C->Ia.lattice)&&!mpz_cmp(C->Ia.norm,C->n1));
    B_VERIFY(integral_lattice(&C->Ia.lattice,I->parent_order,&QUATALG_PINFTY)&&quat_lattice_contains(NULL,&C->Ia.lattice,&C->eta,&QUATALG_PINFTY));quat_lattice_index(&z,&C->Ia.lattice,I->parent_order);mpz_mul(w,C->n1,C->n1);B_VERIFY(!mpz_cmp(z,w));mpz_mul(w,C->n1,C->n2);B_VERIFY(norm_integer(&z,&C->eta,&QUATALG_PINFTY)&&!mpz_cmp(z,w));mpz_mul(w,C->n1,I->norm);B_VERIFY(norm_integer(&z,&C->transport,&QUATALG_PINFTY)&&!mpz_cmp(z,w));
    multiply_lattice(&back,&C->Ia.lattice,&C->transport,&QUATALG_PINFTY);ibz_mat_4x4_copy(&target.basis,&I->lattice.basis);mpz_set(target.denom,I->lattice.denom);for(unsigned i=0;i<4;i++)for(unsigned j=0;j<4;j++)mpz_mul(target.basis[i][j],target.basis[i][j],C->n1);quat_lattice_reduce_denom(&target,&target);quat_lattice_hnf(&target);B_VERIFY(quat_lattice_equal(&back,&target));
    mpz_set_ui(z,1);unsigned common=1;for(unsigned j=0;j<11;j++){B_VERIFY(C->primes[j]==b_primes[j]);unsigned e=C->exponents[j];B_VERIFY(e<4096);unsigned k=e?order_bound(QUATALG_PINFTY.p,b_primes[j],e+(j==0?3:0),64):1;B_VERIFY(k&&k==C->extensions[j]);mpz_ui_pow_ui(w,b_primes[j],e);mpz_mul(z,z,w);unsigned a=common,b=k;while(b){unsigned r=a%b;a=b;b=r;}common=common/a*k;B_VERIFY(common<=64);}B_VERIFY(!mpz_cmp(z,C->N)&&common==C->common_extension);B_VERIFY(C->exponents[0]>=TORSION_PLUS_EVEN_POWER-2&&mpz_divisible_p(C->N,TORSION_PLUS_3POWER));ok=1;
finish:
#undef B_VERIFY
    quat_lattice_finalize(&target);quat_lattice_finalize(&back);b_gaussian_certificate_clear(&D);ibz_finalize(&g);ibz_finalize(&w);ibz_finalize(&z);return ok;
}
static void b_certificate_copy(b_gaussian_certificate *out,const b_gaussian_certificate *C){
    elem_copy(&out->beta1,&C->beta1);elem_copy(&out->beta2,&C->beta2);elem_copy(&out->x,&C->x);elem_copy(&out->y,&C->y);elem_copy(&out->eta,&C->eta);elem_copy(&out->transport,&C->transport);ideal_copy(&out->Ia,&C->Ia);mpz_set(out->d1,C->d1);mpz_set(out->d2,C->d2);mpz_set(out->u,C->u);mpz_set(out->v,C->v);mpz_set(out->N,C->N);mpz_set(out->n1,C->n1);mpz_set(out->n2,C->n2);memcpy(out->primes,C->primes,sizeof(C->primes));memcpy(out->exponents,C->exponents,sizeof(C->exponents));memcpy(out->extensions,C->extensions,sizeof(C->extensions));out->common_extension=C->common_extension;out->verified=C->verified;
}
int b_gaussian_prepare(b_gaussian_certificate *out,const quat_left_ideal_t *I,const b_gaussian_plans *P,b_gaussian_counts *counts){
    if(counts)memset(counts,0,sizeof(*counts));if(!out||!P||!P->ready||!P->plans||!b_input_valid(I))return B_GAUSSIAN_INVALID;
    ibz_mat_4x4_t basis;ibz_mat_4x4_init(&basis);if(quat_lattice_lll(&basis,&I->lattice,&QUATALG_PINFTY.p)!=0){ibz_mat_4x4_finalize(&basis);return B_GAUSSIAN_INVALID;}
    b_candidate cand[625];unsigned nc=b_candidates(cand,I,&basis),limit=nc<32?nc:32;const b_plan *plans=P->plans;b_gaussian_certificate C;b_gaussian_certificate_init(&C);ibz_t D,tmp,g;ibz_init(&D);ibz_init(&tmp);ibz_init(&g);unsigned long tests=0,first=0;unsigned selected=0;int status=B_GAUSSIAN_MISS;
    for(unsigned h=0;h<P->count&&status==0&&tests<200000;h++)for(unsigned i=0;i<limit&&status==0&&tests<200000;i++)for(unsigned j=i+1;j<limit&&status==0&&tests<200000;j++){
        mpz_set(C.d1,cand[i].d);mpz_set(C.d2,cand[j].d);if((mpz_fdiv_ui(C.d1,4)+mpz_fdiv_ui(C.d2,4))%4)continue;mpz_mul(D,C.d1,C.d2);mpz_add(tmp,C.d1,C.d2);if(mpz_cmp(plans[h].N,tmp)<0)continue;mpz_gcd(g,D,plans[h].N);if(mpz_cmp_ui(g,1))continue;mpz_gcd(g,C.d1,C.d2);if(mpz_cmp_ui(g,1))continue;
        if(!mpz_invert(C.u,C.d1,C.d2)){status=B_GAUSSIAN_INVALID;break;}mpz_mul(C.u,C.u,plans[h].N);mpz_mod(C.u,C.u,C.d2);mpz_mul(C.v,C.u,C.d1);mpz_sub(C.v,plans[h].N,C.v);mpz_divexact(C.v,C.v,C.d2);
        for(unsigned shift=0;shift<64&&mpz_sgn(C.v)>0&&tests<200000;shift++,mpz_add(C.u,C.u,C.d2),mpz_sub(C.v,C.v,C.d1)){
            if(mpz_fdiv_ui(C.u,4)!=1||mpz_divisible_ui_p(C.u,3)||mpz_fdiv_ui(C.v,4)!=1||mpz_divisible_ui_p(C.v,3))continue;mpz_gcd(g,C.u,C.v);if(mpz_cmp_ui(g,1))continue;++tests;if(!two_squares(&C.x.coord[0],&C.x.coord[1],&C.u))continue;++first;if(!two_squares(&C.y.coord[0],&C.y.coord[1],&C.v))continue;
            mpz_set_ui(C.x.denom,1);mpz_set_ui(C.y.denom,1);mpz_set(C.N,plans[h].N);mpz_mul(C.n1,C.u,C.d1);mpz_mul(C.n2,C.v,C.d2);b_set_beta(&C.beta1,cand+i,I,&basis);b_set_beta(&C.beta2,cand+j,I,&basis);b_assemble(&C,I);C.common_extension=1;for(unsigned l=0;l<11;l++){C.primes[l]=b_primes[l];C.exponents[l]=plans[h].e[l];C.extensions[l]=plans[h].k[l];unsigned a=C.common_extension,b=C.extensions[l];while(b){unsigned r=a%b;a=b;b=r;}C.common_extension=C.common_extension/a*C.extensions[l];}selected=h;
            if(!b_gaussian_certificate_verify(&C,I)){status=B_GAUSSIAN_INVALID;break;}C.verified=1;b_certificate_copy(out,&C);status=B_GAUSSIAN_OK;break;
        }
    }
    if(counts){counts->candidate_norms=nc;counts->selected_plan=selected;counts->two_square_tests=tests;counts->first_square_successes=first;}for(unsigned i=0;i<nc;i++)ibz_finalize(&cand[i].d);ibz_finalize(&g);ibz_finalize(&tmp);ibz_finalize(&D);b_gaussian_certificate_clear(&C);ibz_mat_4x4_finalize(&basis);return status;
}
