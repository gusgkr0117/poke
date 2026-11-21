#ifndef INKE_RESEARCH_EXTENSION_PAIRING_H
#define INKE_RESEARCH_EXTENSION_PAIRING_H
#include "a-extension-ec.h"
/* Independent rational-function Weil pairing. Random auxiliary divisors
 * avoid zeros and poles; no final-exponentiation/Tate-pairing shortcut.
 * Inputs must be killed by n, prime to the characteristic. A budget of
 * zero allows unbounded auxiliary-point retries. Failure preserves out. */
int a_ec_weil(fq_t out,const a_ec_point *P,const a_ec_point *Q,const fmpz_t n,
               const fq_t A,a_ec_rng rng,void *context,unsigned long budget,
               const root_extension *E);
#endif
