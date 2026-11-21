#include "root-poke-field-bridge.h"

void root_from_poke_fp2(fq_t out, const fp2_t *x, const root_extension *E)
{
    (void)E;
    *out = *x;
}

int root_to_poke_fp2(fp2_t *out, const fq_t x, const root_extension *E)
{
    if (!E || E->relative_degree != 1) return 0;
    *out = *x;
    return 1;
}
