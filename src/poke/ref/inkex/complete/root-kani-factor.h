#ifndef INKE_RESEARCH_KANI_FACTOR_H
#define INKE_RESEARCH_KANI_FACTOR_H
#include "root-tate-projector.h"
#include "root-theta-factor.h"
#include "root-field-embedding.h"
/* Isolate the designated elliptic factor of an already evaluated quotient.
 * lattice_points contains 2g consecutive level4 vectors, in graph-lattice
 * B/4 order. If lattice_embedding is non-NULL, they are over its sup field
 * and are checked/projectively descended to E. The payload vectors and
 * target theta are over E. Pnum/eta use one coherently aligned Tate action.
 * All signed points share one target Montgomery model. No target-factor
 * search by j-invariant or independent change of basis is performed.
 * Output curve A and initialized point array are unchanged on failure.
 * stage is optional: 1=input,2=labels,3=projector/frame,4=factor,5=points.
 */
int root_kani_extract(fq_t A,a_ec_point *points,size_t count,
 unsigned g,b_theta4_ctx *target,const fq_struct *lattice_points,
 const root_field_embedding *lattice_embedding,const fq_struct *payload,
 const fmpz_mat_t Pnum,const fmpz_mat_t eta,const fmpz_t n1,const fmpz_t N,
 const fq_t zeta4,const fq_t zeta8,const root_extension *E,unsigned *stage);
#endif
