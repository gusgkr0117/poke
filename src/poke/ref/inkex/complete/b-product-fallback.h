#ifndef INKE_RESEARCH_PRODUCT_FALLBACK_H
#define INKE_RESEARCH_PRODUCT_FALLBACK_H
#include "b-poke-basis.h"
#include "b-theta-balanced.h"
#include <inke_keygen.h>
typedef void (*a_kani_progress)(void *context,size_t prime_index,unsigned long ell,
                               unsigned completed,unsigned total);
/* Public field setup is reused between finish/decrypt calls. The borrowed
 * base and the RNG context must outlive this zero-initialized context. */
typedef struct {
 const root_extension *base;
 a_ec_rng rng;void *rng_context;
 unsigned long point_budget;
 a_kani_progress progress;void *progress_context;
 unsigned long calls;
 int ready;
} b_product_context;
int b_product_context_init(b_product_context *context,const root_extension *base,
 a_ec_rng rng,void *rng_context,unsigned long point_budget);
void b_product_context_clear(b_product_context *context);
/* Implements the inke_product_fallback contract on maximal Fp2 curves.
 * No source isogeny representation beyond the supplied graph is required.
 * The two supplied bases have order 2^(exponent+2). The point on the first
 * lifted graph generator labels the target, including for isomorphic factors.
 * All theta coordinates and evaluations stay over the base field. */
int b_product_fallback(ec_curve_t *target,ec_basis_t *left_images,
 ec_basis_t *right_images,const ec_curve_t *left_curve,
 const ec_curve_t *right_curve,const ec_basis_t *left_two,
 const ec_basis_t *right_two,unsigned int exponent,const ibz_t *q,
 const ec_basis_t *left_points,const ec_basis_t *right_points,
 void *context,int *stage);
#endif
