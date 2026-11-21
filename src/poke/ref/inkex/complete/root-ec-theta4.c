#include "root-ec-theta4.h"
#include <string.h>

int root_ec_theta4_init(root_ec_theta4 *T,const fq_t A,const root_extension *E){
    memset(T,0,sizeof(*T));T->E=E;
    fq_t *vars[]={&T->A,&T->t,&T->a,&T->b,&T->c,&T->k0,&T->k2};
    for(int j=0;j<7;j++)fq_init(*vars[j],E->field);
    fq_set(T->A,A,E->field);T->null=b_theta4_vector_new(4,E->field);
    if(!T->null||!a_ec_curve_valid(A,E))return 0;
    fq_t r,s,u,v,D,one,two;fq_t *tmp[]={&r,&s,&u,&v,&D,&one,&two};
    for(int j=0;j<7;j++)fq_init(*tmp[j],E->field);fq_one(one,E->field);fq_set_ui(two,2,E->field);
    fq_add(r,A,two,E->field);fq_sub(s,A,two,E->field);fq_inv(s,s,E->field);fq_mul(r,r,s,E->field);
    if(!root_extension_sqrt(s,r,E)||!root_extension_sqrt(T->t,s,E)||!root_extension_sqrt(T->b,T->t,E))goto done;
    /* a+c=sqrt(2+2t^2), a-c=sqrt(2-2t^2). */
    fq_sqr(r,T->t,E->field);fq_mul_ui(r,r,2,E->field);fq_add(s,two,r,E->field);fq_sub(u,two,r,E->field);
    if(!root_extension_sqrt(s,s,E)||!root_extension_sqrt(u,u,E))goto done;
    fq_add(T->a,s,u,E->field);fq_sub(T->c,s,u,E->field);fq_mul_fmpz(T->a,T->a,E->half,E->field);fq_mul_fmpz(T->c,T->c,E->half,E->field);
    fq_sqr(r,T->a,E->field);fq_sqr(s,T->c,E->field);fq_sub(D,r,s,E->field);
    if(fq_is_zero(D,E->field)||fq_is_zero(T->b,E->field))goto done;
    fq_inv(D,D,E->field);fq_add(s,s,one,E->field);fq_mul(s,s,T->a,E->field);fq_mul_ui(s,s,4,E->field);fq_mul(T->k0,s,D,E->field);
    fq_add(r,r,one,E->field);fq_mul(r,r,T->c,E->field);fq_mul_ui(r,r,4,E->field);fq_mul(T->k2,r,D,E->field);
    /* Check the defining identities and the recovered Montgomery model. */
    fq_mul(r,T->a,T->c,E->field);fq_sqr(s,T->t,E->field);if(!fq_equal(r,s,E->field))goto done;
    fq_sqr(r,T->a,E->field);fq_sqr(s,T->c,E->field);fq_add(r,r,s,E->field);if(!fq_equal(r,two,E->field))goto done;
    fq_sqr(r,T->b,E->field);if(!fq_equal(r,T->t,E->field))goto done;
    fq_sqr(r,T->t,E->field);fq_sqr(r,r,E->field);fq_add(s,one,r,E->field);fq_mul_ui(s,s,2,E->field);fq_sub(r,one,r,E->field);fq_mul(r,r,A,E->field);fq_add(r,r,s,E->field);if(!fq_is_zero(r,E->field))goto done;
    fq_set(T->null,T->a,E->field);fq_set(T->null+1,T->b,E->field);fq_set(T->null+2,T->c,E->field);fq_set(T->null+3,T->b,E->field);T->ready=1;
done:for(int j=6;j>=0;j--)fq_clear(*tmp[j],E->field);return T->ready;
}
void root_ec_theta4_clear(root_ec_theta4 *T){
    if(!T->E)return;const fq_ctx_struct *f=T->E->field;b_theta4_vector_free(T->null,4,f);
    fq_t *vars[]={&T->A,&T->t,&T->a,&T->b,&T->c,&T->k0,&T->k2};for(int j=6;j>=0;j--)fq_clear(*vars[j],f);memset(T,0,sizeof(*T));
}
int root_ec_theta4_from_null(root_ec_theta4 *T,const fq_struct *null,const root_extension *E){
    memset(T,0,sizeof(*T));T->E=E;
    fq_t *vars[]={&T->A,&T->t,&T->a,&T->b,&T->c,&T->k0,&T->k2};for(int j=0;j<7;j++)fq_init(*vars[j],E->field);
    T->null=b_theta4_vector_new(4,E->field);if(!T->null)return 0;
    const fq_ctx_struct *f=E->field;if(!fq_equal(null+1,null+3,f)||fq_is_zero(null,f)||fq_is_zero(null+1,f)||fq_is_zero(null+2,f))return 0;
    fq_set(T->a,null,f);fq_set(T->b,null+1,f);fq_set(T->c,null+2,f);
    fq_t a2,b2,c2,S,D,w,v,one;fq_t *tmp[]={&a2,&b2,&c2,&S,&D,&w,&v,&one};for(int j=0;j<8;j++)fq_init(*tmp[j],f);
    fq_one(one,f);fq_sqr(a2,T->a,f);fq_sqr(b2,T->b,f);fq_sqr(c2,T->c,f);fq_add(S,a2,c2,f);fq_sub(D,a2,c2,f);
    if(fq_is_zero(S,f)||fq_is_zero(D,f))goto done;
    fq_mul(w,T->a,T->c,f);fq_mul(w,w,S,f);fq_sqr(v,b2,f);fq_mul_ui(v,v,2,f);if(!fq_equal(w,v,f))goto done;
    fq_inv(w,S,f);fq_mul(T->t,b2,w,f);fq_mul_ui(T->t,T->t,2,f);
    fq_sqr(w,T->t,f);fq_sqr(w,w,f);fq_sub(v,one,w,f);if(fq_is_zero(v,f))goto done;
    fq_add(w,one,w,f);fq_mul_ui(w,w,2,f);fq_neg(w,w,f);fq_inv(v,v,f);fq_mul(T->A,w,v,f);if(!a_ec_curve_valid(T->A,E))goto done;
    fq_inv(D,D,f);fq_mul_ui(w,c2,3,f);fq_add(w,w,a2,f);fq_mul(w,w,T->a,f);fq_mul_ui(w,w,2,f);fq_mul(T->k0,w,D,f);
    fq_mul_ui(w,a2,3,f);fq_add(w,w,c2,f);fq_mul(w,w,T->c,f);fq_mul_ui(w,w,2,f);fq_mul(T->k2,w,D,f);
    for(int j=0;j<4;j++)fq_set(T->null+j,null+j,f);T->ready=1;
done:for(int j=7;j>=0;j--)fq_clear(*tmp[j],f);return T->ready;
}
int root_ec_theta4_base_change(root_ec_theta4 *out,const root_ec_theta4 *in,const root_extension *E){
    if(out==in)return 0;memset(out,0,sizeof(*out));out->E=E;
    fq_t *dst[]={&out->A,&out->t,&out->a,&out->b,&out->c,&out->k0,&out->k2};
    const fq_struct *src[]={in->A,in->t,in->a,in->b,in->c,in->k0,in->k2};
    for(int j=0;j<7;j++)fq_init(*dst[j],E->field);out->null=b_theta4_vector_new(4,E->field);
    if(!out->null||!in->ready||fmpz_cmp(fq_ctx_prime(E->field),fq_ctx_prime(in->E->field)))return 0;
    fmpz_t a,b;fmpz_init(a);fmpz_init(b);int ok=1;
    for(int j=0;j<7;j++){if(!root_extension_descend(a,b,src[j],in->E)){ok=0;break;}root_extension_embed(*dst[j],a,b,E);}
    if(ok){fq_set(out->null,out->a,E->field);fq_set(out->null+1,out->b,E->field);fq_set(out->null+2,out->c,E->field);fq_set(out->null+3,out->b,E->field);out->ready=1;}
    fmpz_clear(b);fmpz_clear(a);return ok;
}
int root_ec_theta4_map(fq_struct *out,const a_ec_point *P,const root_ec_theta4 *T){
    if(!T->ready||!a_ec_on_curve(P,T->A,T->E))return 0;const fq_ctx_struct *f=T->E->field;
    if(P->infinity){for(int j=0;j<4;j++)fq_set(out+j,T->null+j,f);return 1;}
    fq_t x2,plus,minus,y2,w;fq_t *v[]={&x2,&plus,&minus,&y2,&w};for(int j=0;j<5;j++)fq_init(*v[j],f);
    fq_one(w,f);fq_sqr(x2,P->x,f);fq_add(plus,x2,w,f);fq_sub(minus,x2,w,f);fq_mul_ui(y2,P->y,2,f);
    fq_mul(out,T->a,plus,f);fq_mul(w,T->k0,P->x,f);fq_sub(out,out,w,f);
    fq_mul(out+2,T->c,plus,f);fq_mul(w,T->k2,P->x,f);fq_add(out+2,out+2,w,f);
    fq_add(w,minus,y2,f);fq_mul(out+1,T->b,w,f);fq_sub(w,minus,y2,f);fq_mul(out+3,T->b,w,f);
    for(int j=4;j>=0;j--)fq_clear(*v[j],f);return 1;
}
int root_ec_theta4_unmap(a_ec_point *out,const fq_struct *P,const root_ec_theta4 *T){
    if(!T->ready)return 0;const fq_ctx_struct *f=T->E->field;int ok=0;
    fq_t q00,q01,q11,U,V,num,den,w,scale;fq_t *v[]={&q00,&q01,&q11,&U,&V,&num,&den,&w,&scale};for(int j=0;j<9;j++)fq_init(*v[j],f);
    fq_mul(q00,P,T->a,f);fq_mul(w,P+2,T->c,f);fq_add(q00,q00,w,f);
    fq_add(q01,P+1,P+3,f);fq_mul(q01,q01,T->b,f);
    fq_mul(q11,P,T->c,f);fq_mul(w,P+2,T->a,f);fq_add(q11,q11,w,f);
    if(!fq_is_zero(q00,f)){fq_set(U,q00,f);fq_set(V,q01,f);}else if(!fq_is_zero(q11,f)){fq_set(U,q01,f);fq_set(V,q11,f);}else goto done;
    fq_mul(w,T->t,U,f);fq_add(num,V,w,f);fq_sub(den,V,w,f);
    a_ec_point R;a_ec_init(&R,T->E);fq_struct *check=b_theta4_vector_new(4,f);if(!check){a_ec_clear(&R,T->E);goto done;}
    if(fq_is_zero(den,f)){a_ec_infinity(&R,T->E);for(int j=0;j<4;j++)fq_set(check+j,T->null+j,f);}
    else{
        fq_inv(den,den,f);fq_mul(R.x,num,den,f);R.infinity=0;
        /* Determine the projective scale from the even coordinate pair. */
        fq_sqr(w,R.x,f);fq_one(num,f);fq_add(w,w,num,f);fq_mul(U,T->a,w,f);fq_mul(num,T->k0,R.x,f);fq_sub(U,U,num,f);
        fq_mul(V,T->c,w,f);fq_mul(num,T->k2,R.x,f);fq_add(V,V,num,f);
        if(!fq_is_zero(U,f)){fq_inv(scale,U,f);fq_mul(scale,scale,P,f);}else if(!fq_is_zero(V,f)){fq_inv(scale,V,f);fq_mul(scale,scale,P+2,f);}else goto cleanup;
        if(fq_is_zero(scale,f))goto cleanup;fq_mul(scale,scale,T->b,f);fq_mul_ui(scale,scale,4,f);fq_inv(scale,scale,f);
        fq_sub(R.y,P+1,P+3,f);fq_mul(R.y,R.y,scale,f);if(!root_ec_theta4_map(check,&R,T))goto cleanup;
    }
    {int pivot=-1;for(int j=0;j<4;j++)if(!fq_is_zero(check+j,f)){pivot=j;break;}if(pivot<0||fq_is_zero(P+pivot,f))goto cleanup;
    for(int j=0;j<4;j++){fq_mul(U,check+j,P+pivot,f);fq_mul(V,P+j,check+pivot,f);if(!fq_equal(U,V,f))goto cleanup;}}
    a_ec_set(out,&R,T->E);ok=1;
cleanup:b_theta4_vector_free(check,4,f);a_ec_clear(&R,T->E);
done:for(int j=8;j>=0;j--)fq_clear(*v[j],f);return ok;
}
int root_ec_theta4_label(unsigned *alpha,unsigned *beta,const a_ec_point *P,const fq_t zeta,const root_ec_theta4 *T){
 if(!T->ready)return 0;const fq_ctx_struct *f=T->E->field;int ok=0;fq_t w;fq_init(w,f);fq_sqr(w,zeta,f);fq_neg(w,w,f);if(!fq_is_one(w,f)){fq_clear(w,f);return 0;}fq_clear(w,f);
 fq_struct *point=b_theta4_vector_new(4,f),*orbit=b_theta4_vector_new(4,f);if(!point||!orbit)goto done;
 if(!root_ec_theta4_map(point,P,T))goto done;b_theta4_ctx B;if(!b_theta4_init(&B,1,T->null,f))goto done;
 for(unsigned a=0;a<4&&!ok;a++)for(unsigned b=0;b<4;b++){b_theta4_action(orbit,T->null,a,b,zeta,&B);if(b_theta4_equal_projective(point,orbit,&B)){*alpha=a;*beta=b;ok=1;break;}}
 b_theta4_clear(&B);
done:b_theta4_vector_free(orbit,4,f);b_theta4_vector_free(point,4,f);return ok;
}
