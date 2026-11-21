// SPDX-License-Identifier: Apache-2.0
/* Portable 64-bit Montgomery arithmetic for the manuscript INKE primes.
 * Fixed loop bounds and one conditional subtraction; no external arithmetic
 * runtime is used by field operations. Requires 2*p < R and p = -1 mod 2^64.
 * This backend is checked against GMP independently of the elliptic tests. */
#include <fp.h>
#include <inke_field_constants.h>

#if RADIX != 64
#error This Montgomery backend requires 64-bit limbs.
#endif
typedef unsigned __int128 wide_t;
const uint64_t p[NWORDS_FIELD]=INKE_FIELD_MODULUS;
const uint64_t ZERO[NWORDS_FIELD]={0};
const uint64_t ONE[NWORDS_FIELD]=INKE_FIELD_R;
static const fp_t R2=INKE_FIELD_R2;

static void reduce_once(fp_t *out,const fp_t *a)
{
    fp_t difference; uint64_t borrow=0;
    for (int k=0;k<NWORDS_FIELD;++k) {
        wide_t sub=(wide_t)p[k]+borrow;
        difference[k]=(uint64_t)((wide_t)(*a)[k]-sub);
        borrow=((wide_t)(*a)[k]<sub);
    }
    uint64_t mask=(uint64_t)0-borrow;
    for (int k=0;k<NWORDS_FIELD;++k)
        (*out)[k]=(difference[k]&~mask)|((*a)[k]&mask);
}
void fp_add(fp_t *out,const fp_t *a,const fp_t *b)
{
    fp_t sum; uint64_t carry=0;
    for (int k=0;k<NWORDS_FIELD;++k) {
        wide_t t=(wide_t)(*a)[k]+(*b)[k]+carry;
        sum[k]=(uint64_t)t; carry=(uint64_t)(t>>64);
    }
    reduce_once(out,&sum);
}
void fp_sub(fp_t *out,const fp_t *a,const fp_t *b)
{
    fp_t difference; uint64_t borrow=0;
    for (int k=0;k<NWORDS_FIELD;++k) {
        wide_t sub=(wide_t)(*b)[k]+borrow;
        difference[k]=(uint64_t)((wide_t)(*a)[k]-sub);
        borrow=((wide_t)(*a)[k]<sub);
    }
    uint64_t mask=(uint64_t)0-borrow,carry=0;
    for (int k=0;k<NWORDS_FIELD;++k) {
        wide_t t=(wide_t)difference[k]+(p[k]&mask)+carry;
        (*out)[k]=(uint64_t)t; carry=(uint64_t)(t>>64);
    }
}
void fp_mul(fp_t *out,const fp_t *a,const fp_t *b)
{
    uint64_t t[2*NWORDS_FIELD+1]={0};
    for (int i=0;i<NWORDS_FIELD;++i) {
        uint64_t carry=0;
        for (int j=0;j<NWORDS_FIELD;++j) {
            wide_t v=(wide_t)(*a)[i]*(*b)[j]+t[i+j]+carry;
            t[i+j]=(uint64_t)v; carry=(uint64_t)(v>>64);
        }
        t[i+NWORDS_FIELD]=carry;
    }
    for (int i=0;i<NWORDS_FIELD;++i) {
        uint64_t m=t[i],carry=0; /* -p[0]^{-1} = 1 mod 2^64. */
        for (int j=0;j<NWORDS_FIELD;++j) {
            wide_t v=(wide_t)m*p[j]+t[i+j]+carry;
            t[i+j]=(uint64_t)v; carry=(uint64_t)(v>>64);
        }
        for (int j=i+NWORDS_FIELD;j<=2*NWORDS_FIELD;++j) {
            wide_t v=(wide_t)t[j]+carry;
            t[j]=(uint64_t)v; carry=(uint64_t)(v>>64);
        }
    }
    fp_t result;
    for (int k=0;k<NWORDS_FIELD;++k) result[k]=t[k+NWORDS_FIELD];
    reduce_once(out,&result);
}
void fp_sqr(fp_t *out,const fp_t *a) { fp_mul(out,a,a); }
void fp_tomont(fp_t *out,const fp_t *a) { fp_mul(out,a,&R2); }
void fp_frommont(fp_t *out,const fp_t *a)
{ fp_t one={1}; fp_mul(out,a,&one); }
void fp_mont_setone(fp_t *out) { memcpy(*out,ONE,sizeof(fp_t)); }
