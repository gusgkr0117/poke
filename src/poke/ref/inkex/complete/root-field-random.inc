/* Included once by root-extension.c so every existing build gets the same
 * randomness contract. The arithmetic routines themselves accept callbacks. */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <bcrypt.h>
#else
#include <errno.h>
#include <sys/random.h>
#endif

int root_system_random(void *unused,unsigned char *out,size_t count){
    (void)unused;
#ifdef _WIN32
    /* Resolve only the system copy; no extra linker dependency is imposed on
     * standalone arithmetic tests. No seed or mutable global PRNG is used. */
    HMODULE module=LoadLibraryExW(L"bcrypt.dll",NULL,LOAD_LIBRARY_SEARCH_SYSTEM32);
    if(!module)return 0;
    typedef NTSTATUS (WINAPI *generate_type)(BCRYPT_ALG_HANDLE,PUCHAR,ULONG,ULONG);
    generate_type generate=(generate_type)(void*)GetProcAddress(module,"BCryptGenRandom");
    int ok=generate!=NULL;
    while(ok&&count){ULONG take=count>ULONG_MAX?ULONG_MAX:(ULONG)count;
        ok=generate(NULL,out,take,BCRYPT_USE_SYSTEM_PREFERRED_RNG)>=0;
        out+=take;count-=take;
    }
    FreeLibrary(module);return ok;
#else
    while(count){size_t take=count>1048576?1048576:count;
        ssize_t got=getrandom(out,take,0);
        if(got<0){if(errno==EINTR)continue;return 0;}
        if(!got)return 0;out+=(size_t)got;count-=(size_t)got;
    }return 1;
#endif
}
int root_random_below(fmpz_t out,const fmpz_t bound,root_random_bytes rng,void *context){
    if(!rng||fmpz_sgn(bound)<=0)return 0;
    slong bits=fmpz_bits(bound),words=(bits+FLINT_BITS-1)/FLINT_BITS;
    if((size_t)words>SIZE_MAX/sizeof(ulong))return 0;
    ulong *buf=malloc((size_t)words*sizeof(ulong));if(!buf)return 0;
    fmpz_t value;fmpz_init(value);int ok=0;
    for(;;){if(!rng(context,(unsigned char*)buf,(size_t)words*sizeof(ulong)))break;
        unsigned rem=(unsigned)(bits%FLINT_BITS);
        if(rem)buf[words-1]&=(((ulong)1<<rem)-1);
        fmpz_set_ui_array(value,buf,words);
        if(fmpz_cmp(value,bound)<0){fmpz_set(out,value);ok=1;break;}
    }
    memset(buf,0,(size_t)words*sizeof(ulong));free(buf);fmpz_clear(value);return ok;
}
int root_random_field(fq_t out,const fq_ctx_t field,root_random_bytes rng,void *context){
    fmpz_poly_t value;fmpz_poly_init(value);fmpz_t c;fmpz_init(c);int ok=1;
    for(slong j=0;j<fq_ctx_degree(field);j++){
        if(!root_random_below(c,fq_ctx_prime(field),rng,context)){ok=0;break;}
        fmpz_poly_set_coeff_fmpz(value,j,c);
    }
    if(ok)fq_set_fmpz_poly(out,value,field);
    fmpz_clear(c);fmpz_poly_clear(value);return ok;
}
void root_sqrt_cache_clear(root_sqrt_cache *S,const fq_ctx_t field){
    if(S->ready){fq_clear(S->two_primary_generator,field);fmpz_clear(S->root_exponent);fmpz_clear(S->odd_part);}
    memset(S,0,sizeof(*S));
}
int root_sqrt_cache_init(root_sqrt_cache *S,const fq_ctx_t field,root_random_bytes rng,void *context){
    memset(S,0,sizeof(*S));if(!rng||fmpz_is_even(fq_ctx_prime(field)))return 0;
    fmpz_init(S->odd_part);fmpz_init(S->root_exponent);fq_init(S->two_primary_generator,field);S->ready=1;
    fmpz_t order,euler;fmpz_init(order);fmpz_init(euler);fq_ctx_order(order,field);fmpz_sub_ui(order,order,1);
    S->valuation=fmpz_val2(order);fmpz_tdiv_q_2exp(S->odd_part,order,S->valuation);
    fmpz_add_ui(S->root_exponent,S->odd_part,1);fmpz_tdiv_q_2exp(S->root_exponent,S->root_exponent,1);
    fmpz_tdiv_q_2exp(euler,order,1);
    fq_t z,test,minus;fq_init(z,field);fq_init(test,field);fq_init(minus,field);fq_set_si(minus,-1,field);int ok=0;
    for(;;){if(!root_random_field(z,field,rng,context))break;
        if(fq_is_zero(z,field))continue;
        fq_pow(test,z,euler,field);if(!fq_equal(test,minus,field))continue;
        fq_pow(S->two_primary_generator,z,S->odd_part,field);ok=1;break;
    }
    fq_clear(minus,field);fq_clear(test,field);fq_clear(z,field);fmpz_clear(euler);fmpz_clear(order);
    if(!ok)root_sqrt_cache_clear(S,field);return ok;
}
int root_sqrt_cached(fq_t out,const fq_t in,const root_sqrt_cache *S,const fq_ctx_t field){
    if(!S->ready)return 0;
    if(fq_is_zero(in,field)||fq_is_one(in,field)){fq_set(out,in,field);return 1;}
    fq_t c,t,x,b,test;fq_init(c,field);fq_init(t,field);fq_init(x,field);fq_init(b,field);fq_init(test,field);
    fq_set(c,S->two_primary_generator,field);fq_pow(t,in,S->odd_part,field);fq_pow(x,in,S->root_exponent,field);
    flint_bitcnt_t m=S->valuation;int ok=1;
    while(!fq_is_one(t,field)){
        fq_set(test,t,field);flint_bitcnt_t i=0;
        while(i<m&&!fq_is_one(test,field)){fq_sqr(test,test,field);i++;}
        if(i>=m){ok=0;break;}
        fq_set(b,c,field);for(flint_bitcnt_t j=0;j<m-i-1;j++)fq_sqr(b,b,field);
        fq_mul(x,x,b,field);fq_sqr(c,b,field);fq_mul(t,t,c,field);m=i;
    }
    if(ok){fq_sqr(test,x,field);ok=fq_equal(test,in,field);}
    if(ok)fq_set(out,x,field);
    fq_clear(test,field);fq_clear(b,field);fq_clear(x,field);fq_clear(t,field);fq_clear(c,field);return ok;
}
int root_fq_sqrt(fq_t out,const fq_t in,const fq_ctx_t field){
    if(fq_is_zero(in,field)||fq_is_one(in,field)){fq_set(out,in,field);return 1;}
    root_sqrt_cache S;if(!root_sqrt_cache_init(&S,field,root_system_random,NULL))return 0;
    int ok=root_sqrt_cached(out,in,&S,field);root_sqrt_cache_clear(&S,field);return ok;
}
int root_extension_sqrt(fq_t out,const fq_t in,const root_extension *E){return root_sqrt_cached(out,in,&E->square_roots,E->field);}
