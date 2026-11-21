#include "a-extension-ec.h"
#include <stdlib.h>
#include <string.h>

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
