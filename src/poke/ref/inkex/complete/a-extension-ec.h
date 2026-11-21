#ifndef INKE_RESEARCH_EXTENSION_EC_H
#define INKE_RESEARCH_EXTENSION_EC_H
#include <stddef.h>
#include "root-extension.h"
/* Return 1 if the callback supplied the requested random bytes. */
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
#endif
