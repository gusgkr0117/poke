/* Pollack--Trevino, Integers18A(2018), section4, unconditional algorithm.
 * No primality oracle or factorization of the input is used.
 * Internal Hurwitz elements store twice their four coordinates.
 */
#include <stdio.h>
#include <stdlib.h>
#include <intbig.h>
#define FS_CHECK(c) do {if(!(c)){fprintf(stderr,"four-squares invariant failed line%d: %s\n",__LINE__,#c);abort();}}while(0)
typedef struct {ibz_t c[4];} a_hurwitz;
typedef struct {ibz_t c[2];} a_gaussian;
static void h_init(a_hurwitz *h){for(int i=0;i<4;i++)ibz_init(&h->c[i]);}
static void h_clear(a_hurwitz *h){for(int i=0;i<4;i++)ibz_finalize(&h->c[i]);}
static void h_copy(a_hurwitz *h,const a_hurwitz *x){for(int i=0;i<4;i++)mpz_set(h->c[i],x->c[i]);}
static void h_scalar(a_hurwitz *h,const ibz_t *n){mpz_mul_2exp(h->c[0],*n,1);for(int i=1;i<4;i++)mpz_set_ui(h->c[i],0);}
static void h_one(a_hurwitz *h){mpz_set_ui(h->c[0],2);for(int i=1;i<4;i++)mpz_set_ui(h->c[i],0);}
static int h_zero(const a_hurwitz *h){return !mpz_sgn(h->c[0])&&!mpz_sgn(h->c[1])&&!mpz_sgn(h->c[2])&&!mpz_sgn(h->c[3]);}
static int h_integral(const a_hurwitz *h){return mpz_even_p(h->c[0])&&mpz_even_p(h->c[1])&&mpz_even_p(h->c[2])&&mpz_even_p(h->c[3]);}
static void h_conj(a_hurwitz *out,const a_hurwitz *x){mpz_set(out->c[0],x->c[0]);for(int i=1;i<4;i++)mpz_neg(out->c[i],x->c[i]);}
static void h_norm(ibz_t *n,const a_hurwitz *h){mpz_set_ui(*n,0);for(int i=0;i<4;i++)mpz_addmul(*n,h->c[i],h->c[i]);FS_CHECK(mpz_divisible_2exp_p(*n,2));mpz_fdiv_q_2exp(*n,*n,2);}
/* Raw Hamilton product of numerator vectors, without denominator reduction. */
static void h_raw(a_hurwitz *out,const a_hurwitz *a,const a_hurwitz *b){
    a_hurwitz r;h_init(&r);
    mpz_mul(r.c[0],a->c[0],b->c[0]);for(int i=1;i<4;i++)mpz_submul(r.c[0],a->c[i],b->c[i]);
    mpz_mul(r.c[1],a->c[0],b->c[1]);mpz_addmul(r.c[1],a->c[1],b->c[0]);mpz_addmul(r.c[1],a->c[2],b->c[3]);mpz_submul(r.c[1],a->c[3],b->c[2]);
    mpz_mul(r.c[2],a->c[0],b->c[2]);mpz_submul(r.c[2],a->c[1],b->c[3]);mpz_addmul(r.c[2],a->c[2],b->c[0]);mpz_addmul(r.c[2],a->c[3],b->c[1]);
    mpz_mul(r.c[3],a->c[0],b->c[3]);mpz_addmul(r.c[3],a->c[1],b->c[2]);mpz_submul(r.c[3],a->c[2],b->c[1]);mpz_addmul(r.c[3],a->c[3],b->c[0]);
    h_copy(out,&r);h_clear(&r);
}
static void h_mul(a_hurwitz *out,const a_hurwitz *a,const a_hurwitz *b){h_raw(out,a,b);for(int i=0;i<4;i++){FS_CHECK(mpz_even_p(out->c[i]));mpz_fdiv_q_2exp(out->c[i],out->c[i],1);}}
static void fs_nearest_integer(ibz_t *q,const ibz_t *num,const ibz_t *den){
    ibz_t a,b;ibz_init(&a);ibz_init(&b);mpz_mul_2exp(a,*num,1);mpz_add(a,a,*den);mpz_mul_2exp(b,*den,1);mpz_fdiv_q(*q,a,b);ibz_finalize(&b);ibz_finalize(&a);
}
/* alpha=q beta+r if left_quotient is1; alpha=beta q+r otherwise. */
static void h_remainder(a_hurwitz *r,const a_hurwitz *alpha,const a_hurwitz *beta,int left_quotient){
    a_hurwitz cb,raw,qi,qh,q,prod;h_init(&cb);h_init(&raw);h_init(&qi);h_init(&qh);h_init(&q);h_init(&prod);
    ibz_t den,err,di,dh,tmp,nr,nb;ibz_init(&den);ibz_init(&err);ibz_init(&di);ibz_init(&dh);ibz_init(&tmp);ibz_init(&nr);ibz_init(&nb);
    h_conj(&cb,beta);if(left_quotient)h_raw(&raw,alpha,&cb);else h_raw(&raw,&cb,alpha);
    mpz_set_ui(den,0);for(int i=0;i<4;i++)mpz_addmul(den,beta->c[i],beta->c[i]);FS_CHECK(mpz_sgn(den)>0);
    mpz_set_ui(di,0);mpz_set_ui(dh,0);
    for(int i=0;i<4;i++){
        fs_nearest_integer(&qi.c[i],&raw.c[i],&den);mpz_mul_2exp(qi.c[i],qi.c[i],1);
        mpz_fdiv_q(qh.c[i],raw.c[i],den);mpz_mul_2exp(qh.c[i],qh.c[i],1);mpz_add_ui(qh.c[i],qh.c[i],1);
        mpz_mul_2exp(tmp,raw.c[i],1);mpz_mul(err,qi.c[i],den);mpz_sub(err,tmp,err);mpz_addmul(di,err,err);
        mpz_mul(err,qh.c[i],den);mpz_sub(err,tmp,err);mpz_addmul(dh,err,err);
    }
    h_copy(&q,mpz_cmp(di,dh)<=0?&qi:&qh);
    if(left_quotient)h_mul(&prod,&q,beta);else h_mul(&prod,beta,&q);
    for(int i=0;i<4;i++)mpz_sub(r->c[i],alpha->c[i],prod.c[i]);
    h_norm(&nr,r);h_norm(&nb,beta);mpz_mul_2exp(nr,nr,1);FS_CHECK(mpz_cmp(nr,nb)<=0);
    ibz_finalize(&nb);ibz_finalize(&nr);ibz_finalize(&tmp);ibz_finalize(&dh);ibz_finalize(&di);ibz_finalize(&err);ibz_finalize(&den);
    h_clear(&prod);h_clear(&q);h_clear(&qh);h_clear(&qi);h_clear(&raw);h_clear(&cb);
}
static void h_gcd(a_hurwitz *out,const a_hurwitz *x,const a_hurwitz *y,int right_divisor){
    a_hurwitz a,b,r;h_init(&a);h_init(&b);h_init(&r);h_copy(&a,x);h_copy(&b,y);
    while(!h_zero(&b)){h_remainder(&r,&a,&b,right_divisor);h_copy(&a,&b);h_copy(&b,&r);}
    h_copy(out,&a);h_clear(&r);h_clear(&b);h_clear(&a);
}
static void h_make_integral(a_hurwitz *h){
    if(h_integral(h))return;
    a_hurwitz unit,r;h_init(&unit);h_init(&r);int done=0;
    for(unsigned mask=0;mask<16&&!done;mask++){
        for(int i=0;i<4;i++)mpz_set_si(unit.c[i],(mask>>i)&1?1:-1);
        h_mul(&r,&unit,h);if(h_integral(&r)){h_copy(h,&r);done=1;}
    }
    FS_CHECK(done);h_clear(&r);h_clear(&unit);
}
static void g_init(a_gaussian *g){ibz_init(&g->c[0]);ibz_init(&g->c[1]);}
static void g_clear(a_gaussian *g){ibz_finalize(&g->c[1]);ibz_finalize(&g->c[0]);}
static void g_copy(a_gaussian *g,const a_gaussian *h){mpz_set(g->c[0],h->c[0]);mpz_set(g->c[1],h->c[1]);}
static void g_one(a_gaussian *g){mpz_set_ui(g->c[0],1);mpz_set_ui(g->c[1],0);}
static void g_mul(a_gaussian *out,const a_gaussian *a,const a_gaussian *b){
    a_gaussian r;g_init(&r);mpz_mul(r.c[0],a->c[0],b->c[0]);mpz_submul(r.c[0],a->c[1],b->c[1]);mpz_mul(r.c[1],a->c[0],b->c[1]);mpz_addmul(r.c[1],a->c[1],b->c[0]);g_copy(out,&r);g_clear(&r);
}
static void g_norm(ibz_t *n,const a_gaussian *g){mpz_mul(*n,g->c[0],g->c[0]);mpz_addmul(*n,g->c[1],g->c[1]);}
static void g_gcd(a_gaussian *out,const a_gaussian *x,const a_gaussian *y){
    a_gaussian a,b,q,r,prod;g_init(&a);g_init(&b);g_init(&q);g_init(&r);g_init(&prod);g_copy(&a,x);g_copy(&b,y);
    ibz_t den,num;ibz_init(&den);ibz_init(&num);
    while(mpz_sgn(b.c[0])||mpz_sgn(b.c[1])){
        g_norm(&den,&b);mpz_mul(num,a.c[0],b.c[0]);mpz_addmul(num,a.c[1],b.c[1]);fs_nearest_integer(&q.c[0],&num,&den);
        mpz_mul(num,a.c[1],b.c[0]);mpz_submul(num,a.c[0],b.c[1]);fs_nearest_integer(&q.c[1],&num,&den);
        g_mul(&prod,&q,&b);for(int i=0;i<2;i++)mpz_sub(r.c[i],a.c[i],prod.c[i]);g_copy(&a,&b);g_copy(&b,&r);
    }
    g_copy(out,&a);ibz_finalize(&num);ibz_finalize(&den);g_clear(&prod);g_clear(&r);g_clear(&q);g_clear(&b);g_clear(&a);
}
static int small_prime(unsigned long n){if(n<2)return 0;for(unsigned long j=2;j<=n/j;j++)if(n%j==0)return 0;return 1;}
static void small_two_squares(a_gaussian *g,unsigned long n){
    for(unsigned long x=0;x*x<=n;x++)for(unsigned long y=0;x*x+y*y<=n;y++)if(x*x+y*y==n){mpz_set_ui(g->c[0],x);mpz_set_ui(g->c[1],y);return;}
    FS_CHECK(0);
}
static void small_four_squares(a_hurwitz *h,unsigned long n){
    for(unsigned long a=0;a*a<=n;a++)for(unsigned long b=0;a*a+b*b<=n;b++)for(unsigned long c=0;a*a+b*b+c*c<=n;c++)for(unsigned long d=0;a*a+b*b+c*c+d*d<=n;d++)if(a*a+b*b+c*c+d*d==n){mpz_set_ui(h->c[0],2*a);mpz_set_ui(h->c[1],2*b);mpz_set_ui(h->c[2],2*c);mpz_set_ui(h->c[3],2*d);return;}
    FS_CHECK(0);
}

/* Return1 success,0 explicit trial-budget exhaustion,-1 RNG failure.
 * budget=0 means the proved expected-polynomial Las Vegas algorithm.
 */
int a_four_squares(a_hurwitz *out,const ibz_t *input,unsigned long budget,unsigned long *trials){
    if(mpz_sgn(*input)<0)return 0;
    if(mpz_cmp_ui(*input,20)<=0){small_four_squares(out,mpz_get_ui(*input));return 1;}
    ibz_t n,N,P,g,x,y,r,rem,s,w,hi,e,tmp,norm;
    ibz_init(&n);ibz_init(&N);ibz_init(&P);ibz_init(&g);ibz_init(&x);ibz_init(&y);ibz_init(&r);ibz_init(&rem);ibz_init(&s);ibz_init(&w);ibz_init(&hi);ibz_init(&e);ibz_init(&tmp);ibz_init(&norm);
    a_gaussian outer,small,part,ga,gb,large,zw;g_init(&outer);g_init(&small);g_init(&part);g_init(&ga);g_init(&gb);g_init(&large);g_init(&zw);g_one(&outer);
    a_hurwitz h,scalar,result,outerh;h_init(&h);h_init(&scalar);h_init(&result);h_init(&outerh);
    mpz_set(n,*input);unsigned long T=mpz_sizeinbase(n,2);
    for(unsigned long ell=2;ell<=T;ell++)if((ell==2||ell%4==1)&&small_prime(ell)){
        if(!mpz_divisible_ui_p(n,ell))continue;small_two_squares(&small,ell);
        do{mpz_divexact_ui(n,n,ell);g_mul(&outer,&outer,&small);}while(mpz_divisible_ui_p(n,ell));
    }
    int status=0;
    if(mpz_cmp_ui(n,20)<=0){small_four_squares(&result,mpz_get_ui(n));status=1;goto finish;}
    /* Recompute the small-prime threshold after removing known factors.
     * bitlength differs from natural log by fixed constants only; section4's
     * counting proof and polynomial bound are unchanged (documented). */
    T=mpz_sizeinbase(n,2);mpz_set_ui(P,1);
    for(unsigned long ell=3;ell<=T;ell+=4)if(small_prime(ell))mpz_mul_ui(P,P,ell);
    mpz_gcd(g,P,n);mpz_divexact(P,P,g);mpz_mul(N,n,P);
    for(unsigned long k=0;!budget||k<budget;k++){
        if(trials)++*trials;
        if(!ibz_rand_interval(&x,&ibz_const_one,&N)||!ibz_rand_interval(&y,&ibz_const_one,&N)){status=-1;goto finish;}
        mpz_mul(r,x,x);mpz_addmul(r,y,y);mpz_neg(r,r);mpz_mod(r,r,N);
        if(mpz_fdiv_ui(r,4)!=1)continue;mpz_gcd(g,r,N);if(mpz_cmp_ui(g,1)!=0)continue;
        mpz_set(rem,r);g_one(&part);
        for(unsigned long ell=5;ell<=T;ell+=4)if(small_prime(ell)&&mpz_divisible_ui_p(rem,ell)){
            small_two_squares(&small,ell);
            do{mpz_divexact_ui(rem,rem,ell);g_mul(&part,&part,&small);}while(mpz_divisible_ui_p(rem,ell));
        }
        if(mpz_cmp_ui(rem,1)==0)g_one(&large);
        else {
            mpz_sub_ui(hi,rem,1);FS_CHECK(mpz_divisible_2exp_p(hi,2));mpz_fdiv_q_2exp(e,hi,2);
            if(!ibz_rand_interval(&w,&ibz_const_one,&hi)){status=-1;goto finish;}
            mpz_powm(s,w,e,rem);mpz_mul(tmp,s,s);mpz_mod(tmp,tmp,rem);if(mpz_cmp(tmp,hi)!=0)continue;
            mpz_set(ga.c[0],rem);mpz_set_ui(ga.c[1],0);mpz_set(gb.c[0],s);mpz_set_ui(gb.c[1],1);
            g_gcd(&large,&ga,&gb);g_norm(&norm,&large);FS_CHECK(mpz_cmp(norm,rem)==0);
        }
        g_mul(&zw,&part,&large);g_norm(&norm,&zw);FS_CHECK(mpz_cmp(norm,r)==0);
        mpz_gcd(g,zw.c[0],zw.c[1]);FS_CHECK(mpz_cmp_ui(g,1)==0);
        mpz_mul_2exp(h.c[0],x,1);mpz_mul_2exp(h.c[1],y,1);mpz_mul_2exp(h.c[2],zw.c[0],1);mpz_mul_2exp(h.c[3],zw.c[1],1);
        h_norm(&norm,&h);FS_CHECK(mpz_divisible_p(norm,n));h_scalar(&scalar,&n);
        h_gcd(&result,&scalar,&h,1);h_norm(&norm,&result);FS_CHECK(mpz_cmp(norm,n)==0);h_make_integral(&result);status=1;break;
    }
finish:
    if(status==1){
        mpz_mul_2exp(outerh.c[0],outer.c[0],1);mpz_mul_2exp(outerh.c[1],outer.c[1],1);mpz_set_ui(outerh.c[2],0);mpz_set_ui(outerh.c[3],0);
        h_mul(out,&result,&outerh);FS_CHECK(h_integral(out));h_norm(&norm,out);FS_CHECK(mpz_cmp(norm,*input)==0);
    }
    h_clear(&outerh);h_clear(&result);h_clear(&scalar);h_clear(&h);
    g_clear(&zw);g_clear(&large);g_clear(&gb);g_clear(&ga);g_clear(&part);g_clear(&small);g_clear(&outer);
    ibz_finalize(&norm);ibz_finalize(&tmp);ibz_finalize(&e);ibz_finalize(&hi);ibz_finalize(&w);ibz_finalize(&s);ibz_finalize(&rem);ibz_finalize(&r);ibz_finalize(&y);ibz_finalize(&x);ibz_finalize(&g);ibz_finalize(&P);ibz_finalize(&N);ibz_finalize(&n);return status;
}

/* Matrix for left (right=0) or right (right=1) Hamilton multiplication.
 * Input must have integer coordinates; result is an integer4x4 matrix. */
void a_hamilton_matrix(ibz_t mat[4][4],const a_hurwitz *x,int right){
    a_hurwitz e,col;h_init(&e);h_init(&col);FS_CHECK(h_integral(x));
    for(int j=0;j<4;j++){
        for(int i=0;i<4;i++)mpz_set_ui(e.c[i],i==j?2:0);
        if(right)h_mul(&col,&e,x);else h_mul(&col,x,&e);
        for(int i=0;i<4;i++){FS_CHECK(mpz_even_p(col.c[i]));mpz_fdiv_q_2exp(mat[i][j],col.c[i],1);}
    }
    h_clear(&col);h_clear(&e);
}

#ifndef A_FOUR_SQUARES_LIBRARY
int main(void){
    ibz_t n,g,tmp;ibz_init(&n);ibz_init(&g);ibz_init(&tmp);a_hurwitz x,y;h_init(&x);h_init(&y);
    unsigned long total=0;
    for(unsigned long k=0;k<=512;k++){mpz_set_ui(n,k);FS_CHECK(a_four_squares(&x,&n,100000,&total)==1);}
    printf("small_inputs=513 exact_norms=PASS trials=%lu\n",total);fflush(stdout);
    ibz_t X[4][4],Y[4][4];for(int i=0;i<4;i++)for(int j=0;j<4;j++){ibz_init(&X[i][j]);ibz_init(&Y[i][j]);}
    for(unsigned k=0;k<12;k++){
        unsigned bits=64+96*k;mpz_set_ui(n,1);mpz_mul_2exp(n,n,bits);mpz_add_ui(n,n,12345+2*k);
        unsigned long trials=0;FS_CHECK(a_four_squares(&x,&n,100000,&trials)==1);a_hamilton_matrix(X,&x,0);
        mpz_add_ui(n,n,2);FS_CHECK(a_four_squares(&y,&n,100000,&trials)==1);a_hamilton_matrix(Y,&y,1);
        for(int i=0;i<4;i++)for(int j=0;j<4;j++){
            mpz_set_ui(g,0);mpz_set_ui(tmp,0);for(int l=0;l<4;l++){mpz_addmul(g,X[i][l],Y[l][j]);mpz_addmul(tmp,Y[i][l],X[l][j]);}FS_CHECK(mpz_cmp(g,tmp)==0);
        }
        h_norm(&n,&x);for(int i=0;i<4;i++)for(int j=0;j<4;j++){mpz_set_ui(g,0);for(int l=0;l<4;l++)mpz_addmul(g,X[l][i],X[l][j]);FS_CHECK(i==j?mpz_cmp(g,n)==0:mpz_sgn(g)==0);}
        h_norm(&n,&y);for(int i=0;i<4;i++)for(int j=0;j<4;j++){mpz_set_ui(g,0);for(int l=0;l<4;l++)mpz_addmul(g,Y[l][i],Y[l][j]);FS_CHECK(i==j?mpz_cmp(g,n)==0:mpz_sgn(g)==0);}
        printf("bits=%u four_squares_pair=PASS norm_matrices=PASS left_right_commute=PASS trials=%lu\n",bits,trials);fflush(stdout);
    }
    for(int i=0;i<4;i++)for(int j=0;j<4;j++){ibz_finalize(&Y[i][j]);ibz_finalize(&X[i][j]);}
    h_clear(&y);h_clear(&x);ibz_finalize(&tmp);ibz_finalize(&g);ibz_finalize(&n);return 0;
}
#endif
