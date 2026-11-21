#ifndef INKE_RESEARCH_EC_HALVING_H
#define INKE_RESEARCH_EC_HALVING_H
#include "root-ec-tate.h"
/* One half of a point on E0:y²=x³+x, in the SAME supplied field.
 * Returns0 if the point is invalid or no rational half exists; no field
 * extension or random search is silently performed. Exact membership and
 * doubling are checked. Output unchanged on failure, aliasing allowed. */
int root_ec_half_e0(a_ec_point *out,const a_ec_point *P,const root_extension *E);
/* Refine previously known maximal-order action matrices from low_bits to
 * bits in a compatible lifted elliptic basis. The old basis must equal
 * [2^(bits-low_bits)]basis, so its action is the low_bits residue.
 * Only bits-low_bits discrete-log bits are solved. Wrong supplied residues
 * fail the exact residual reconstruction. Outputs unchanged on failure. */
int root_ec_refine_actions(fmpz_mat_struct out[4],const fmpz_mat_struct known[4],
 const a_ec_point basis[2],const a_ec_point half[2],unsigned low_bits,
 unsigned bits,const root_extension *E);
#endif
