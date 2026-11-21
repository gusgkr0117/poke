#ifndef INKE_RESEARCH_POKE_FIELD_BRIDGE_H
#define INKE_RESEARCH_POKE_FIELD_BRIDGE_H
#include "root-extension.h"
#include <fp2.h>
/* E must use the same prime as this compiled POKE field backend. */
void root_from_poke_fp2(fq_t out,const fp2_t *x,const root_extension *E);
/* Rejects extension elements outside F_(p^2); output stays unchanged. */
int root_to_poke_fp2(fp2_t *out,const fq_t x,const root_extension *E);
#endif
