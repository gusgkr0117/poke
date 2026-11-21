// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_RESEARCH_POKE_BASIS_H
#define INKE_RESEARCH_POKE_BASIS_H
#include "a-extension-ec.h"
#include <inke_conversion.h>
/* The field prime must equal the linked POKE parameter prime. All fq and
 * a_ec_point outputs are initialized by the caller. Failures preserve outputs.
 * Lift P,Q with the relative sign fixed by the supplied P-Q. Only a common
 * sign remains; no basis change, restriction, or random sampling occurs. */
int b_poke_basis_lift(fq_t A,a_ec_point out[3],const ec_curve_t *curve,
                     const ec_basis_t *basis,const root_extension *E);
/* Export ONE model and the signed P,Q,P-Q triple, requiring F_(p^2) descent.
 * The curve is normalized to C=1 and its entire cache/storage is reset. */
int b_poke_basis_export(ec_curve_t *curve,ec_basis_t *basis,const fq_t A,
                       const a_ec_point points[3],const root_extension *E);
#endif
