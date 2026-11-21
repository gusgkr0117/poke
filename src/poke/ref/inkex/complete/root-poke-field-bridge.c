#include "root-poke-field-bridge.h"
#include <encoded_sizes.h>
static void read_le(fmpz_t x,const unsigned char *p,size_t n) {
    fmpz_zero(x);while(n){fmpz_mul_2exp(x,x,8);fmpz_add_ui(x,x,p[--n]);}
}
static void write_le(unsigned char *p,size_t n,const fmpz_t x) {
    fmpz_t t;fmpz_init(t);fmpz_set(t,x);
    for(size_t j=0;j<n;j++){p[j]=(unsigned char)fmpz_fdiv_ui(t,256);fmpz_fdiv_q_2exp(t,t,8);}fmpz_clear(t);
}
void root_from_poke_fp2(fq_t out,const fp2_t *x,const root_extension *E) {
    unsigned char bytes[FP2_ENCODED_BYTES];fmpz_t a,b;fmpz_init(a);fmpz_init(b);
    fp2_encode(bytes,x);read_le(a,bytes,FP_ENCODED_BYTES);read_le(b,bytes+FP_ENCODED_BYTES,FP_ENCODED_BYTES);
    root_extension_embed(out,a,b,E);fmpz_clear(b);fmpz_clear(a);
}
int root_to_poke_fp2(fp2_t *out,const fq_t x,const root_extension *E) {
    unsigned char bytes[FP2_ENCODED_BYTES];fmpz_t a,b;fmpz_init(a);fmpz_init(b);
    int ok=root_extension_descend(a,b,x,E);
    if(ok){write_le(bytes,FP_ENCODED_BYTES,a);write_le(bytes+FP_ENCODED_BYTES,FP_ENCODED_BYTES,b);fp2_t value;fp2_decode(&value,bytes);*out=value;}
    fmpz_clear(b);fmpz_clear(a);return ok;
}
