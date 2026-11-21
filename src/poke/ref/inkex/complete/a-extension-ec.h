#ifndef INKE_RESEARCH_EXTENSION_EC_H
#define INKE_RESEARCH_EXTENSION_EC_H
#include <stddef.h>
#include "root-extension.h"
/* The caller supplies cryptographically uniform bytes. Return 1 on success,
 * 0 on failure; the core does not use FLINT's deterministic test PRNG. */
typedef int (*a_ec_rng)(void *context,unsigned char *out,size_t count);
enum { A_EC_OK=1,A_EC_RETRY=0,A_EC_INVALID=-1,A_EC_RNG_ERROR=-2 };
typedef struct { fq_t x,y; int infinity; } a_ec_point;
void a_ec_init(a_ec_point *P,const root_extension *E);
void a_ec_clear(a_ec_point *P,const root_extension *E);
void a_ec_infinity(a_ec_point *P,const root_extension *E);
void a_ec_set(a_ec_point *P,const a_ec_point *Q,const root_extension *E);
int a_ec_equal(const a_ec_point *P,const a_ec_point *Q,const root_extension *E);
int a_ec_curve_valid(const fq_t A,const root_extension *E);
int a_ec_on_curve(const a_ec_point *P,const fq_t A,const root_extension *E);
void a_ec_neg(a_ec_point *R,const a_ec_point *P,const root_extension *E);
/* All point operations support output/input aliasing. */
void a_ec_add(a_ec_point *R,const a_ec_point *P,const a_ec_point *Q,
              const fq_t A,const root_extension *E);
void a_ec_mul(a_ec_point *R,const a_ec_point *P,const fmpz_t n,
              const fq_t A,const root_extension *E);
void a_ec_frobenius(a_ec_point *R,const a_ec_point *P,slong power,
                    const root_extension *E);
/* Exact uniform field sampling; RNG error is distinct. */
int a_ec_random_field(fq_t out,a_ec_rng rng,void *context,const root_extension *E);
/* Uniform over the affine rational points (infinity is excluded).
 * budget=0 means no candidate cap. Failure leaves out unchanged. */
int a_ec_random_point(a_ec_point *out,const fq_t A,a_ec_rng rng,void *context,
                       unsigned long budget,unsigned long *trials,const root_extension *E);
/* Exact [|(-p)^k-1|/r] projection under the same maximal-model promise.
 * Uses the Frobenius polynomial when r is small, ordinary scalar
 * multiplication otherwise. Failure preserves out. */
int a_ec_torsion_project(a_ec_point *out,const a_ec_point *P,const fmpz_t r,
                          const fq_t A,const root_extension *E);
/* Caller promises A is a maximal F_(p^2) model with Frobenius [-p].
 * ell must be a small prime, e>0, and ell^e divides |(-p)^k-1|.
 * Returns an exact basis of E[ell^e], verified by order-ell independence.
 * budget limits random affine-point candidates; zero means unbounded.
 * RETRY/RNG/invalid returns preserve both outputs. */
int a_ec_torsion_basis(a_ec_point *P,a_ec_point *Q,unsigned long ell,unsigned int e,
                        const fq_t A,a_ec_rng rng,void *context,unsigned long budget,
                        unsigned long *trials,const root_extension *E);
#endif
