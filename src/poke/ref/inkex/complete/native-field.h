#ifndef INKE_NATIVE_FIELD_H
#define INKE_NATIVE_FIELD_H

#include <fp2.h>
#include <gmp.h>

/* All level-four coordinates use the native Montgomery representation of Fp².
 * The array-of-one notation lets the quotient formulas pass them by address. */
typedef fp2_t fq_struct;
typedef fp2_t fq_t[1];
typedef struct {
    mpz_t prime;
    fp_t inverse_four;
} fq_ctx_struct;
typedef fq_ctx_struct fq_ctx_t[1];
typedef unsigned long ulong;
typedef long slong;

typedef mpz_t fmpz_t;
#define fmpz_abs mpz_abs
#define fmpz_add mpz_add
#define fmpz_clear mpz_clear
#define fmpz_cmp mpz_cmp
#define fmpz_equal(x,y) (mpz_cmp((x),(y)) == 0)
#define fmpz_init mpz_init
#define fmpz_mul_2exp mpz_mul_2exp
#define fmpz_mul_ui mpz_mul_ui
#define fmpz_neg mpz_neg
#define fmpz_one(x) mpz_set_ui((x), 1)
#define fmpz_set_mpz mpz_set
#define fmpz_sgn mpz_sgn
#define fmpz_tstbit mpz_tstbit
static inline size_t fmpz_bits(const mpz_t x)
{
    return mpz_sgn(x) ? mpz_sizeinbase(x, 2) : 0;
}

static inline mpz_srcptr fq_ctx_prime(const fq_ctx_t f) { return f->prime; }
static inline slong fq_ctx_degree(const fq_ctx_t f) { (void)f; return 2; }

/* Native coordinates are stored inline and have no per-element allocation. */
static inline void fq_init(fq_t x, const fq_ctx_t f) { (void)f; fp2_set_zero(x); }
static inline void fq_clear(fq_t x, const fq_ctx_t f) { (void)x; (void)f; }
static inline void fq_zero(fq_t x, const fq_ctx_t f) { (void)f; fp2_set_zero(x); }
static inline void fq_one(fq_t x, const fq_ctx_t f) { (void)f; fp2_set_one(x); }
static inline void fq_set(fq_t r, const fq_t x, const fq_ctx_t f) { (void)f; *r = *x; }
static inline void fq_set_ui(fq_t r, ulong n, const fq_ctx_t f) { (void)f; fp2_set_small(r, n); }
static inline int fq_equal(const fq_t x, const fq_t y, const fq_ctx_t f) { (void)f; return fp2_is_equal(x, y) != 0; }
static inline int fq_is_zero(const fq_t x, const fq_ctx_t f) { (void)f; return fp2_is_zero(x) != 0; }
static inline int fq_is_one(const fq_t x, const fq_ctx_t f) { (void)f; return fp2_is_one(x) != 0; }
static inline void fq_swap(fq_t x, fq_t y, const fq_ctx_t f) { (void)f; fp2_t t = *x; *x = *y; *y = t; }
static inline void fq_add(fq_t r, const fq_t x, const fq_t y, const fq_ctx_t f) { (void)f; fp2_add(r, x, y); }
static inline void fq_sub(fq_t r, const fq_t x, const fq_t y, const fq_ctx_t f) { (void)f; fp2_sub(r, x, y); }
static inline void fq_mul(fq_t r, const fq_t x, const fq_t y, const fq_ctx_t f) { (void)f; fp2_mul(r, x, y); }
static inline void fq_neg(fq_t r, const fq_t x, const fq_ctx_t f) { (void)f; fp2_neg(r, x); }
static inline void fq_sqr(fq_t r, const fq_t x, const fq_ctx_t f) { (void)f; fp2_sqr(r, x); }
static inline void fq_inv(fq_t r, const fq_t x, const fq_ctx_t f) { (void)f; *r = *x; fp2_inv(r); }
static inline void fq_div(fq_t r, const fq_t x, const fq_t y, const fq_ctx_t f)
{
    fp2_t inverse;
    fq_inv(&inverse, y, f);
    fp2_mul(r, x, &inverse);
}

void fq_mul_ui(fq_t out, const fq_t x, ulong n, const fq_ctx_t f);
void fq_mul_fmpz(fq_t out, const fq_t x, const mpz_t n, const fq_ctx_t f);
void fq_pow_ui(fq_t out, const fq_t x, ulong n, const fq_ctx_t f);
int root_fq_sqrt(fq_t out, const fq_t x, const fq_ctx_t f);

#endif
