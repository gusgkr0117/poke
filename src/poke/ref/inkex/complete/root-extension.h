#ifndef INKE_RESEARCH_EXTENSION_H
#define INKE_RESEARCH_EXTENSION_H

#include "native-field.h"

typedef int (*root_random_bytes)(void *context, unsigned char *out, size_t count);

/* The product quotient takes values in Fp[i], where i² = -1, with the
 * same i and Montgomery representation as the public curve coordinates. */
typedef struct {
    fq_ctx_t field;
    fq_t i;
    mpz_t half;
    slong relative_degree;
    root_random_bytes rng;
    void *rng_context;
} root_extension;

int root_system_random(void *unused, unsigned char *out, size_t count);
int root_extension_init_rng(root_extension *E, const mpz_t prime, slong k,
                            root_random_bytes rng, void *context);
int root_extension_init(root_extension *E, const mpz_t prime, slong k);
void root_extension_clear(root_extension *E);
int root_extension_sqrt(fq_t out, const fq_t x, const root_extension *E);

#endif
