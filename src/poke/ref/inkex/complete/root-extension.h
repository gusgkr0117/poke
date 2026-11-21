#ifndef INKE_RESEARCH_EXTENSION_H
#define INKE_RESEARCH_EXTENSION_H
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <stddef.h>
/* A successful callback supplies independent uniform random bytes in the
 * mathematical algorithm. A failure is an error, never a new ideal draw. */
typedef int (*root_random_bytes)(void *context,unsigned char *out,size_t count);
typedef struct {
    fmpz_t odd_part,root_exponent;
    fq_t two_primary_generator;
    flint_bitcnt_t valuation;
    int ready;
} root_sqrt_cache;
/* C research backend for F_(p^(2k)), with an explicit F_p[i] embedding.
 * The caller supplies a verified public prime p == 3 mod 4. */
typedef struct {
    fq_ctx_t field;
    fq_t i;
    fmpz_t half;
    slong relative_degree;
    root_random_bytes rng;
    void *rng_context;
    root_sqrt_cache square_roots;
    /* Failure of the initializer: -2 means its byte callback failed; -1
     * means invalid input/allocation. Successful initialization sets 0. */
    int last_error;
} root_extension;
int root_system_random(void *unused,unsigned char *out,size_t count);
int root_random_below(fmpz_t out,const fmpz_t bound,root_random_bytes rng,void *context);
int root_random_field(fq_t out,const fq_ctx_t field,root_random_bytes rng,void *context);
int root_sqrt_cache_init(root_sqrt_cache *S,const fq_ctx_t field,root_random_bytes rng,void *context);
void root_sqrt_cache_clear(root_sqrt_cache *S,const fq_ctx_t field);
/* Deterministic after cache initialization. Nonsquare/invalid inputs preserve
 * out; successful answers are checked by squaring, including aliased calls. */
int root_sqrt_cached(fq_t out,const fq_t in,const root_sqrt_cache *S,const fq_ctx_t field);
/* Generic theta contexts lack an owning root_extension. This uncached wrapper
 * samples its nonsquare explicitly; it never calls FLINT's fixed-seed sqrt. */
int root_fq_sqrt(fq_t out,const fq_t in,const fq_ctx_t field);
int root_extension_sqrt(fq_t out,const fq_t in,const root_extension *E);
int root_extension_init_rng(root_extension *E,const fmpz_t p,slong k,root_random_bytes rng,void *context);
int root_extension_init(root_extension *E,const fmpz_t p,slong k);
void root_extension_clear(root_extension *E);
void root_extension_embed(fq_t out,const fmpz_t real,const fmpz_t imag,const root_extension *E);
/* Rejects elements outside the embedded F_(p^2), without altering outputs. */
int root_extension_descend(fmpz_t real,fmpz_t imag,const fq_t x,const root_extension *E);
#endif
