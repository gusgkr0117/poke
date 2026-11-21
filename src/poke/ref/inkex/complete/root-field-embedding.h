#ifndef INKE_RESEARCH_FIELD_EMBEDDING_H
#define INKE_RESEARCH_FIELD_EMBEDDING_H
#include "root-extension.h"
#include <flint/fmpz_mod_mat.h>
typedef struct {
 const root_extension *sub,*sup;
 fmpz_mod_mat_t embed,project;
 int ready,last_error;
} root_field_embedding;
/* Construct a compatible finite-field embedding when subdegree|supdegree.
 * It MUST carry sub.i to sup.i, fixing the common named Fp2 constants.
 * Bool result; last_error is 0 on success, -2 on byte-RNG failure,
 * -1 on invalid input/allocation/invariant failure, retained after cleanup. */
int root_field_embedding_init(root_field_embedding *L,const root_extension *sub,
                                const root_extension *sup);
void root_field_embedding_clear(root_field_embedding *L);
void root_field_up(fq_t out,const fq_t in,const root_field_embedding *L);
/* Exact subfield membership is checked by re-embedding the proposed result.
 * Output remains unchanged on failure. */
int root_field_down(fq_t out,const fq_t in,const root_field_embedding *L);
/* Normalize the entire nonzero projective vector by one common pivot,
 * then descend all coordinates. This never changes coordinates separately.
 * Outputs remain unchanged on failure. Caller proves the theta marking
 * descends; this function checks the proposed projective descent exactly. */
int root_field_projective_down(fq_struct *out,const fq_struct *in,size_t count,
                                 const root_field_embedding *L);
#endif
