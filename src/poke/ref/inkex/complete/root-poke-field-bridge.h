#ifndef INKE_RESEARCH_POKE_FIELD_BRIDGE_H
#define INKE_RESEARCH_POKE_FIELD_BRIDGE_H
#include "root-extension.h"

/* The complete quotient and the public curves use identical coordinates. */
void root_from_poke_fp2(fq_t out, const fp2_t *x, const root_extension *E);
int root_to_poke_fp2(fp2_t *out, const fq_t x, const root_extension *E);
#endif
