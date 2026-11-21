#include "root-extension.h"
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_poly.h>
#include "root-field-random.inc"
typedef struct {root_random_bytes rng;void *context;int failed;} root_rng_monitor;
static int monitored_random(void *data,unsigned char *out,size_t count){
    root_rng_monitor *M=data;int ok=M->rng(M->context,out,count);if(!ok)M->failed=1;return ok;
}
int root_extension_init(root_extension *E,const fmpz_t p,slong k){return root_extension_init_rng(E,p,k,root_system_random,NULL);}
int root_extension_init_rng(root_extension *E,const fmpz_t p,slong k,root_random_bytes rng,void *context) {
    memset(E,0,sizeof(*E));E->last_error=-1;
    if(!rng || k<1 || k>WORD_MAX/8 || fmpz_cmp_ui(p,3)<0 || fmpz_fdiv_ui(p,4)!=3)return 0;
    root_rng_monitor monitor={rng,context,0};root_random_bytes original_rng=rng;void *original_context=context;
    rng=monitored_random;context=&monitor;slong n=2*k;
    fmpz_mod_ctx_t modp;fmpz_mod_ctx_init(modp,p);fmpz_mod_poly_t modulus;fmpz_mod_poly_init(modulus,modp);
    fmpz_t c,index,bound;fmpz_init(c);fmpz_init(index);fmpz_init(bound);fmpz_set_si(bound,n-1);int ok=0;
    /* A bounded sparse prefix improves arithmetic. Its success probability is
     * not used: after 4n attempts every coefficient is independent uniform. */
    for(slong sparse_left=4*n;;){
        fmpz_mod_poly_zero(modulus,modp);fmpz_mod_poly_set_coeff_ui(modulus,n,1,modp);
        if(k==1)fmpz_mod_poly_set_coeff_ui(modulus,0,1,modp); /* X^2+1 */
        else if(sparse_left){
            sparse_left--;
            if(!root_random_below(c,p,rng,context))break;fmpz_mod_poly_set_coeff_fmpz(modulus,0,c,modp);
            if(!root_random_below(c,p,rng,context)||!root_random_below(index,bound,rng,context))break;
            fmpz_mod_poly_set_coeff_fmpz(modulus,1+fmpz_get_si(index),c,modp);
        }else{
            int sampled=1;for(slong j=0;j<n;j++){
                if(!root_random_below(c,p,rng,context)){sampled=0;break;}
                fmpz_mod_poly_set_coeff_fmpz(modulus,j,c,modp);
            }if(!sampled)break;
        }
        if(fmpz_mod_poly_is_irreducible(modulus,modp)){ok=1;break;}
    }
    if(ok)fq_ctx_init_modulus(E->field,modulus,modp,"z");
    fmpz_clear(bound);fmpz_clear(index);fmpz_clear(c);fmpz_mod_poly_clear(modulus,modp);fmpz_mod_ctx_clear(modp);
    if(!ok){if(monitor.failed)E->last_error=-2;return 0;}
    fq_init(E->i,E->field);fmpz_init(E->half);E->relative_degree=k;E->rng=original_rng;E->rng_context=original_context;
    if(!root_sqrt_cache_init(&E->square_roots,E->field,rng,context)){root_extension_clear(E);if(monitor.failed)E->last_error=-2;return 0;}
    fmpz_add_ui(E->half,p,1);fmpz_fdiv_q_2exp(E->half,E->half,1);
    fq_t minus_one,check;fq_init(minus_one,E->field);fq_init(check,E->field);
    fq_set_si(minus_one,-1,E->field);
    ok=root_extension_sqrt(E->i,minus_one,E);
    if(ok){fq_sqr(check,E->i,E->field);ok=fq_equal(check,minus_one,E->field);}
    if(ok){fq_frobenius(check,E->i,1,E->field);fq_add(check,check,E->i,E->field);ok=fq_is_zero(check,E->field);}
    fq_clear(check,E->field);fq_clear(minus_one,E->field);
    if(!ok)root_extension_clear(E);else E->last_error=0;
    return ok;
}
void root_extension_clear(root_extension *E) {
    root_sqrt_cache_clear(&E->square_roots,E->field);fq_clear(E->i,E->field);fmpz_clear(E->half);fq_ctx_clear(E->field);E->relative_degree=0;
}
void root_extension_embed(fq_t out,const fmpz_t real,const fmpz_t imag,const root_extension *E) {
    fq_t t;fq_init(t,E->field);fq_mul_fmpz(t,E->i,imag,E->field);
    fq_set_fmpz(out,real,E->field);fq_add(out,out,t,E->field);fq_clear(t,E->field);
}
int root_extension_descend(fmpz_t real,fmpz_t imag,const fq_t x,const root_extension *E) {
    int ok=0;fq_t conjugate,a,b,check;fmpz_t ar,br;
    fq_init(conjugate,E->field);fq_init(a,E->field);fq_init(b,E->field);fq_init(check,E->field);
    fmpz_init(ar);fmpz_init(br);
    fq_frobenius(check,x,2,E->field);
    if(!fq_equal(check,x,E->field))goto done;
    fq_frobenius(conjugate,x,1,E->field);
    fq_add(a,x,conjugate,E->field);fq_mul_fmpz(a,a,E->half,E->field);
    fq_sub(b,x,conjugate,E->field);fq_mul_fmpz(b,b,E->half,E->field);
    fq_mul(b,b,E->i,E->field);fq_neg(b,b,E->field); /* divide by i */
    if(!fq_get_fmpz(ar,a,E->field) || !fq_get_fmpz(br,b,E->field))goto done;
    root_extension_embed(check,ar,br,E);
    if(!fq_equal(check,x,E->field))goto done;
    fmpz_set(real,ar);fmpz_set(imag,br);ok=1;
done:
    fmpz_clear(br);fmpz_clear(ar);fq_clear(check,E->field);fq_clear(b,E->field);
    fq_clear(a,E->field);fq_clear(conjugate,E->field);return ok;
}
