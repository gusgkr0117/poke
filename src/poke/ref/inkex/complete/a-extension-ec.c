#include "a-extension-ec.h"
#include <stdlib.h>
#include <string.h>
#include <flint/fmpz_poly.h>

void a_ec_init(a_ec_point *P,const root_extension *E){fq_init(P->x,E->field);fq_init(P->y,E->field);P->infinity=1;}
void a_ec_clear(a_ec_point *P,const root_extension *E){fq_clear(P->y,E->field);fq_clear(P->x,E->field);}
void a_ec_infinity(a_ec_point *P,const root_extension *E){fq_zero(P->x,E->field);fq_one(P->y,E->field);P->infinity=1;}
void a_ec_set(a_ec_point *P,const a_ec_point *Q,const root_extension *E){fq_set(P->x,Q->x,E->field);fq_set(P->y,Q->y,E->field);P->infinity=Q->infinity;}
int a_ec_equal(const a_ec_point *P,const a_ec_point *Q,const root_extension *E){return P->infinity||Q->infinity?P->infinity&&Q->infinity:fq_equal(P->x,Q->x,E->field)&&fq_equal(P->y,Q->y,E->field);}
int a_ec_curve_valid(const fq_t A,const root_extension *E){fq_t t,four;fq_init(t,E->field);fq_init(four,E->field);fq_sqr(t,A,E->field);fq_set_ui(four,4,E->field);int ok=!fq_equal(t,four,E->field);fq_clear(four,E->field);fq_clear(t,E->field);return ok;}

static void rhs(fq_t y,const fq_t x,const fq_t A,const root_extension *E){
    fq_t t;fq_init(t,E->field);fq_add(t,x,A,E->field);fq_mul(t,t,x,E->field);
    fq_one(y,E->field);fq_add(t,t,y,E->field);fq_mul(y,t,x,E->field);fq_clear(t,E->field);
}
int a_ec_on_curve(const a_ec_point *P,const fq_t A,const root_extension *E){
    if(P->infinity)return 1;fq_t u,v;fq_init(u,E->field);fq_init(v,E->field);rhs(u,P->x,A,E);fq_sqr(v,P->y,E->field);
    int ok=fq_equal(u,v,E->field);fq_clear(v,E->field);fq_clear(u,E->field);return ok;
}
void a_ec_neg(a_ec_point *R,const a_ec_point *P,const root_extension *E){a_ec_set(R,P,E);fq_neg(R->y,R->y,E->field);}

/* Jacobian coordinates x=X/Z^2,y=Y/Z^3. Z=0 is infinity. */
typedef struct {fq_t X,Y,Z;} jac;
static void j_init(jac *P,const root_extension *E){fq_init(P->X,E->field);fq_init(P->Y,E->field);fq_init(P->Z,E->field);}
static void j_clear(jac *P,const root_extension *E){fq_clear(P->Z,E->field);fq_clear(P->Y,E->field);fq_clear(P->X,E->field);}
static void j_inf(jac *P,const root_extension *E){fq_zero(P->X,E->field);fq_one(P->Y,E->field);fq_zero(P->Z,E->field);}
static void j_set(jac *P,const jac *Q,const root_extension *E){fq_set(P->X,Q->X,E->field);fq_set(P->Y,Q->Y,E->field);fq_set(P->Z,Q->Z,E->field);}
static void j_from(jac *P,const a_ec_point *Q,const root_extension *E){if(Q->infinity){j_inf(P,E);return;}fq_set(P->X,Q->x,E->field);fq_set(P->Y,Q->y,E->field);fq_one(P->Z,E->field);}
static void j_to(a_ec_point *P,const jac *Q,const root_extension *E){
    if(fq_is_zero(Q->Z,E->field)){a_ec_infinity(P,E);return;}
    fq_t z,t;fq_init(z,E->field);fq_init(t,E->field);fq_inv(z,Q->Z,E->field);fq_sqr(t,z,E->field);fq_mul(P->x,Q->X,t,E->field);fq_mul(t,t,z,E->field);fq_mul(P->y,Q->Y,t,E->field);P->infinity=0;fq_clear(t,E->field);fq_clear(z,E->field);
}
static void j_double(jac *R,const jac *P,const fq_t A,const root_extension *E){
    if(fq_is_zero(P->Z,E->field)||fq_is_zero(P->Y,E->field)){j_inf(R,E);return;}
    jac T;j_init(&T,E);fq_t z2,y2,m,s,t,u;fq_init(z2,E->field);fq_init(y2,E->field);fq_init(m,E->field);fq_init(s,E->field);fq_init(t,E->field);fq_init(u,E->field);
    fq_sqr(z2,P->Z,E->field);fq_sqr(y2,P->Y,E->field);
    fq_sqr(m,P->X,E->field);fq_mul_ui(m,m,3,E->field);fq_mul(t,A,P->X,E->field);fq_mul(t,t,z2,E->field);fq_mul_ui(t,t,2,E->field);fq_add(m,m,t,E->field);fq_sqr(t,z2,E->field);fq_add(m,m,t,E->field);
    fq_mul(s,P->X,y2,E->field);fq_mul_ui(s,s,4,E->field);
    fq_mul(T.Z,P->Y,P->Z,E->field);fq_mul_ui(T.Z,T.Z,2,E->field);
    fq_sqr(T.X,m,E->field);fq_mul_ui(t,s,2,E->field);fq_sub(T.X,T.X,t,E->field);fq_sqr(t,T.Z,E->field);fq_mul(t,t,A,E->field);fq_sub(T.X,T.X,t,E->field);
    fq_sub(t,s,T.X,E->field);fq_mul(T.Y,m,t,E->field);fq_sqr(u,y2,E->field);fq_mul_ui(u,u,8,E->field);fq_sub(T.Y,T.Y,u,E->field);
    j_set(R,&T,E);fq_clear(u,E->field);fq_clear(t,E->field);fq_clear(s,E->field);fq_clear(m,E->field);fq_clear(y2,E->field);fq_clear(z2,E->field);j_clear(&T,E);
}
static void j_add(jac *R,const jac *P,const jac *Q,const fq_t A,const root_extension *E){
    if(fq_is_zero(P->Z,E->field)){j_set(R,Q,E);return;}if(fq_is_zero(Q->Z,E->field)){j_set(R,P,E);return;}
    fq_t z1,z2,u1,u2,s1,s2,h,r,h2,h3,t;fq_t *vars[]={&z1,&z2,&u1,&u2,&s1,&s2,&h,&r,&h2,&h3,&t};for(int i=0;i<11;i++)fq_init(*vars[i],E->field);
    jac T;j_init(&T,E);fq_sqr(z1,P->Z,E->field);fq_sqr(z2,Q->Z,E->field);fq_mul(u1,P->X,z2,E->field);fq_mul(u2,Q->X,z1,E->field);
    fq_mul(t,z2,Q->Z,E->field);fq_mul(s1,P->Y,t,E->field);fq_mul(t,z1,P->Z,E->field);fq_mul(s2,Q->Y,t,E->field);
    fq_sub(h,u2,u1,E->field);fq_sub(r,s2,s1,E->field);
    if(fq_is_zero(h,E->field)){if(fq_is_zero(r,E->field))j_double(R,P,A,E);else j_inf(R,E);goto done;}
    fq_sqr(h2,h,E->field);fq_mul(h3,h2,h,E->field);fq_mul(T.Z,P->Z,Q->Z,E->field);fq_mul(T.Z,T.Z,h,E->field);
    fq_sqr(T.X,r,E->field);fq_add(t,u1,u2,E->field);fq_mul(t,t,h2,E->field);fq_sub(T.X,T.X,t,E->field);fq_sqr(t,T.Z,E->field);fq_mul(t,t,A,E->field);fq_sub(T.X,T.X,t,E->field);
    fq_mul(t,u1,h2,E->field);fq_sub(t,t,T.X,E->field);fq_mul(T.Y,r,t,E->field);fq_mul(t,s1,h3,E->field);fq_sub(T.Y,T.Y,t,E->field);j_set(R,&T,E);
done:j_clear(&T,E);for(int i=10;i>=0;i--)fq_clear(*vars[i],E->field);
}
void a_ec_add(a_ec_point *R,const a_ec_point *P,const a_ec_point *Q,const fq_t A,const root_extension *E){jac p,q,r;j_init(&p,E);j_init(&q,E);j_init(&r,E);j_from(&p,P,E);j_from(&q,Q,E);j_add(&r,&p,&q,A,E);j_to(R,&r,E);j_clear(&r,E);j_clear(&q,E);j_clear(&p,E);}
void a_ec_mul(a_ec_point *R,const a_ec_point *P,const fmpz_t n,const fq_t A,const root_extension *E){
    fmpz_t m;fmpz_init(m);fmpz_abs(m,n);jac base,acc;j_init(&base,E);j_init(&acc,E);j_from(&base,P,E);j_inf(&acc,E);
    for(slong i=fmpz_bits(m);i>0;i--){j_double(&acc,&acc,A,E);if(fmpz_tstbit(m,i-1))j_add(&acc,&acc,&base,A,E);}
    j_to(R,&acc,E);if(fmpz_sgn(n)<0)a_ec_neg(R,R,E);j_clear(&acc,E);j_clear(&base,E);fmpz_clear(m);
}
void a_ec_frobenius(a_ec_point *R,const a_ec_point *P,slong power,const root_extension *E){fq_frobenius(R->x,P->x,power,E->field);fq_frobenius(R->y,P->y,power,E->field);R->infinity=P->infinity;}

static int random_below(fmpz_t out,const fmpz_t p,a_ec_rng rng,void *context){
    slong bits=fmpz_bits(p),words=(bits+FLINT_BITS-1)/FLINT_BITS;ulong *buf=malloc((size_t)words*sizeof(ulong));if(!buf)return A_EC_INVALID;
    int status=A_EC_OK;do{if(!rng(context,(unsigned char*)buf,(size_t)words*sizeof(ulong))){status=A_EC_RNG_ERROR;break;}
        unsigned rem=(unsigned)(bits%FLINT_BITS);if(rem)buf[words-1]&=(((ulong)1<<rem)-1);fmpz_set_ui_array(out,buf,words);
    }while(fmpz_cmp(out,p)>=0);memset(buf,0,(size_t)words*sizeof(ulong));free(buf);return status;
}
int a_ec_random_field(fq_t out,a_ec_rng rng,void *context,const root_extension *E){
    if(!rng)return A_EC_INVALID;fmpz_poly_t poly;fmpz_poly_init(poly);fmpz_t c;fmpz_init(c);int status=A_EC_OK;
    for(slong j=0;j<fq_ctx_degree(E->field);j++){status=random_below(c,fq_ctx_prime(E->field),rng,context);if(status!=A_EC_OK)break;fmpz_poly_set_coeff_fmpz(poly,j,c);}
    if(status==A_EC_OK)fq_set_fmpz_poly(out,poly,E->field);fmpz_clear(c);fmpz_poly_clear(poly);return status;
}
int a_ec_random_point(a_ec_point *out,const fq_t A,a_ec_rng rng,void *context,unsigned long budget,unsigned long *trials,const root_extension *E){
    if(!a_ec_curve_valid(A,E)||!rng)return A_EC_INVALID;a_ec_point P;a_ec_init(&P,E);fq_t r;fq_init(r,E->field);int status=A_EC_RETRY;
    for(unsigned long k=0;!budget||k<budget;k++){
        if(trials)++*trials;status=a_ec_random_field(P.x,rng,context,E);if(status!=A_EC_OK)break;rhs(r,P.x,A,E);
        if(!root_extension_sqrt(P.y,r,E)){status=A_EC_RETRY;continue;}
        unsigned char sign;if(!rng(context,&sign,1)){status=A_EC_RNG_ERROR;break;}sign&=1;
        if(fq_is_zero(P.y,E->field)&&sign){status=A_EC_RETRY;continue;}if(sign)fq_neg(P.y,P.y,E->field);
        P.infinity=0;a_ec_set(out,&P,E);status=A_EC_OK;break;
    }fq_clear(r,E->field);a_ec_clear(&P,E);return status;
}
static int is_small_prime(unsigned long n){if(n<2)return 0;for(unsigned long j=2;j<=n/j;j++)if(n%j==0)return 0;return 1;}
int a_ec_torsion_project(a_ec_point *out,const a_ec_point *P,const fmpz_t r,const fq_t A,const root_extension *E){
    if(fmpz_cmp_ui(r,2)<0||!a_ec_curve_valid(A,E)||!a_ec_on_curve(P,A,E))return A_EC_INVALID;
    fmpz_t S,cofactor,a,q,c;fmpz_init(S);fmpz_init(cofactor);fmpz_init(a);fmpz_init(q);fmpz_init(c);int status=A_EC_INVALID;
    fmpz_pow_ui(S,fq_ctx_prime(E->field),(ulong)E->relative_degree);if(E->relative_degree&1)fmpz_add_ui(S,S,1);else fmpz_sub_ui(S,S,1);
    if(!fmpz_divisible(S,r))goto done;fmpz_divexact(cofactor,S,r);
    if(E->relative_degree==1||fmpz_bits(r)*(ulong)E->relative_degree>=fmpz_bits(cofactor)){
        a_ec_mul(out,P,cofactor,A,E);status=A_EC_OK;goto done;
    }
    /* m=-p, a=m mod r, q=(m-a)/r, c=(a^k-1)/r. The identity
     * (m^k-1)/r=q*(m^(k-1)+a*m^(k-2)+...+a^(k-1))+c lets us
     * replace long base-p scalar blocks by the known p^2-Frobenius. */
    fmpz_neg(q,fq_ctx_prime(E->field));fmpz_mod(a,q,r);fmpz_sub(q,q,a);fmpz_divexact(q,q,r);
    fmpz_pow_ui(c,a,(ulong)E->relative_degree);fmpz_sub_ui(c,c,1);if(!fmpz_divisible(c,r))goto done;fmpz_divexact(c,c,r);
    a_ec_point H,R,T;a_ec_init(&H,E);a_ec_init(&R,E);a_ec_init(&T,E);a_ec_set(&H,P,E);a_ec_set(&R,P,E);
    for(slong j=1;j<E->relative_degree;j++){a_ec_frobenius(&H,&H,2,E);a_ec_mul(&R,&R,a,A,E);a_ec_add(&H,&H,&R,A,E);}
    a_ec_mul(&H,&H,q,A,E);a_ec_mul(&T,P,c,A,E);a_ec_add(&H,&H,&T,A,E);if(E->relative_degree&1)a_ec_neg(&H,&H,E);a_ec_set(out,&H,E);status=A_EC_OK;
    a_ec_clear(&T,E);a_ec_clear(&R,E);a_ec_clear(&H,E);
done:fmpz_clear(c);fmpz_clear(q);fmpz_clear(a);fmpz_clear(cofactor);fmpz_clear(S);return status;
}
int a_ec_torsion_basis(a_ec_point *P,a_ec_point *Q,unsigned long ell,unsigned int e,const fq_t A,a_ec_rng rng,void *context,unsigned long budget,unsigned long *trials,const root_extension *E){
    if(!e||!is_small_prime(ell)||!a_ec_curve_valid(A,E)||!rng)return A_EC_INVALID;
    fmpz_t order,ellpower,cofactor,topscalar,n;fmpz_init(order);fmpz_init(ellpower);fmpz_init(cofactor);fmpz_init(topscalar);fmpz_init(n);
    fq_t check;fq_init(check,E->field);fq_frobenius(check,A,2,E->field);int status=A_EC_INVALID;if(!fq_equal(check,A,E->field))goto end_numbers;
    fmpz_pow_ui(order,fq_ctx_prime(E->field),(ulong)E->relative_degree);if(E->relative_degree&1)fmpz_add_ui(order,order,1);else fmpz_sub_ui(order,order,1);
    fmpz_set_ui(ellpower,ell);fmpz_pow_ui(ellpower,ellpower,e);if(!fmpz_divisible(order,ellpower))goto end_numbers;
    fmpz_divexact(cofactor,order,ellpower);fmpz_divexact_ui(topscalar,ellpower,ell);fmpz_set_ui(n,ell);
    a_ec_point p,q,top_p,top_q,sum,kill,random;a_ec_point *points[]={&p,&q,&top_p,&top_q,&sum,&kill,&random};for(int j=0;j<7;j++)a_ec_init(points[j],E);
    int have_p=0;status=A_EC_RETRY;
    for(unsigned long k=0;!budget||k<budget;k++){
        status=a_ec_random_point(&random,A,rng,context,1,trials,E);if(status<0)break;if(status==A_EC_RETRY)continue;
        status=a_ec_torsion_project(&q,&random,ellpower,A,E);if(status!=A_EC_OK)break;a_ec_mul(&top_q,&q,topscalar,A,E);if(top_q.infinity){status=A_EC_RETRY;continue;}
        a_ec_mul(&kill,&top_q,n,A,E);if(!kill.infinity){status=A_EC_INVALID;break;}
        if(!have_p){a_ec_set(&p,&q,E);a_ec_set(&top_p,&top_q,E);have_p=1;status=A_EC_RETRY;continue;}
        a_ec_infinity(&sum,E);int dependent=0;
        for(unsigned long j=0;j<ell;j++){if(a_ec_equal(&sum,&top_q,E)){dependent=1;break;}a_ec_add(&sum,&sum,&top_p,A,E);}
        if(dependent){status=A_EC_RETRY;continue;}a_ec_set(P,&p,E);a_ec_set(Q,&q,E);status=A_EC_OK;break;
    }
    for(int j=6;j>=0;j--)a_ec_clear(points[j],E);
end_numbers:fq_clear(check,E->field);fmpz_clear(n);fmpz_clear(topscalar);fmpz_clear(cofactor);fmpz_clear(ellpower);fmpz_clear(order);return status;
}
