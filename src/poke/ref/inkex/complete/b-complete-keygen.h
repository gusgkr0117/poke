#ifndef INKE_RESEARCH_COMPLETE_KEYGEN_H
#define INKE_RESEARCH_COMPLETE_KEYGEN_H
#include "b-product-fallback.h"
#include "b-fast-keygen.h"
/* Public, per-parameter setup; zero initialize before init. It can be reused
 * serially across keys. The RNG context must remain alive until clear. */
typedef struct {
 root_extension field;
 b_product_context product;
 int ready;
} b_complete_public;
int b_complete_public_init(b_complete_public *public_context,a_ec_rng rng,
 void *rng_context,unsigned long point_budget);
void b_complete_public_clear(b_complete_public *public_context);
/* A single selected receiver instance; zero initialize before init.
 * Clearing means abandoning this key,
 * so do not clear/reinitialize on RETRY. Public context outlives the state.
 * The verified map is retained on each retry of the final construction. */
typedef struct {
 struct {
     ibz_t q;
     unsigned int a;
     int has_q,has_conversion,complete;
     int phase,last_stage;
     unsigned long finish_failures;
     inke_conversion_result conversion;
 } key;
 b_complete_public *public_context;
 b_fast_counts fast_counts;
 unsigned long fast_trials;
 int ready;
} b_complete_keygen;
int b_complete_keygen_init(b_complete_keygen *state,b_complete_public *public_context,
 unsigned long mixed_trials,unsigned long prime_trials,
 unsigned long norm_trials,unsigned long unit_trials);
/* fast_trials is a positive per-step norm-candidate work batch. Exhaustion
 * returns RETRY; subsequent steps continue at the same q without an ideal
 * fallback. Legacy mixed/prime/norm/unit arguments are ignored. */
int b_complete_keygen_init_with_fast_trials(b_complete_keygen *state,
 b_complete_public *public_context,unsigned long mixed_trials,
 unsigned long prime_trials,unsigned long norm_trials,unsigned long unit_trials,
 unsigned long fast_trials);
void b_complete_keygen_clear(b_complete_keygen *state);
int b_complete_keygen_step(inke_sk_t *sk,inke_pk_t *pk,b_complete_keygen *state,
 a_kani_progress progress,void *progress_context);
/* Return 1 on success and 0 on failure; message is unchanged on failure.
 * Internally repeat only product RETRY on the same ciphertext and secret.
 * RNG_ERROR and INVALID return immediately (stage -2 and -3 respectively).
 * Other pre-product validation failures retain the PKE API's stage code. */
int b_complete_decrypt(unsigned char message[32],const inke_ct_t *ciphertext,
 const inke_sk_t *secret,b_complete_public *public_context,int *stage);
#endif
