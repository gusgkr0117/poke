// SPDX-License-Identifier: Apache-2.0
/* Independent GMP oracle for the new portable Montgomery backend. */
#include <fp.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <inke_field_constants.h>
extern const uint64_t p[NWORDS_FIELD];
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"field oracle failed at line %d\n",__LINE__); exit(1); } } while(0)
static void as_integer(mpz_t z,const fp_t *a)
{ mpz_import(z,NWORDS_FIELD,-1,sizeof(uint64_t),0,0,*a); }
static void as_field(fp_t *a,const mpz_t z)
{ memset(*a,0,sizeof(fp_t)); size_t written; mpz_export(*a,&written,-1,sizeof(uint64_t),0,0,z); }
static void equal_integer(const fp_t *a,const mpz_t z,const mpz_t prime)
{
    fp_t ordinary; mpz_t got,reduced; mpz_inits(got,reduced,NULL);
    fp_frommont(&ordinary,a); as_integer(got,&ordinary); mpz_mod(reduced,z,prime);
    CHECK(mpz_cmp(got,reduced)==0); as_integer(got,a); CHECK(mpz_cmp(got,prime)<0);
    mpz_clears(got,reduced,NULL);
}
int main(void)
{
    mpz_t prime,x,y,expected; mpz_inits(prime,x,y,expected,NULL);
    mpz_import(prime,NWORDS_FIELD,-1,sizeof(uint64_t),0,0,p);
    CHECK(mpz_probab_prime_p(prime,40)>0);
    CHECK(p[0]==UINT64_MAX && mpz_sizeinbase(prime,2)<64*NWORDS_FIELD);
    gmp_randstate_t random; gmp_randinit_mt(random); gmp_randseed_ui(random,0x494e4b45);
    for (int k=0;k<1000;++k) {
        mpz_urandomb(x,random,64*NWORDS_FIELD); mpz_urandomm(y,random,prime);
        if (k<4) { mpz_set_ui(x,k); mpz_sub_ui(y,prime,k+1); }
        fp_t a,b,ordinary,result,alias;
        as_field(&ordinary,x); fp_tomont(&a,&ordinary); equal_integer(&a,x,prime);
        as_field(&ordinary,y); fp_tomont(&b,&ordinary); equal_integer(&b,y,prime);
        fp_add(&result,&a,&b); mpz_add(expected,x,y); equal_integer(&result,expected,prime);
        fp_sub(&result,&a,&b); mpz_sub(expected,x,y); equal_integer(&result,expected,prime);
        fp_mul(&result,&a,&b); mpz_mul(expected,x,y); equal_integer(&result,expected,prime);
        fp_copy(&alias,&a); fp_mul(&alias,&alias,&b); CHECK(fp_is_equal(&alias,&result));
        fp_copy(&alias,&b); fp_mul(&alias,&a,&alias); CHECK(fp_is_equal(&alias,&result));
        fp_sqr(&result,&a); mpz_mul(expected,x,x); equal_integer(&result,expected,prime);
        fp_neg(&result,&a); mpz_neg(expected,x); equal_integer(&result,expected,prime);
        if (k<16 && mpz_invert(expected,x,prime)) {
            fp_copy(&result,&a); fp_inv(&result); equal_integer(&result,expected,prime);
            fp_sqr(&result,&a); fp_sqrt(&result); fp_sqr(&result,&result);
            mpz_mul(expected,x,x); equal_integer(&result,expected,prime);
        }
    }
    gmp_randclear(random); mpz_clears(prime,x,y,expected,NULL);
    printf("GMP oracle: 1000 cases passed for %d-bit INKE field.\n",BITS); return 0;
}
