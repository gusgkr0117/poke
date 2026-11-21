#ifndef INKE_RESEARCH_KANI_FULLCHAIN_H
#define INKE_RESEARCH_KANI_FULLCHAIN_H
#include "a-theta-oddbalanced.h"
#include "b-theta-balanced.h"
typedef struct {
 ulong ell;unsigned exponent;const root_extension *field;
 const fq_struct *const *generators,*const *complements;
 const fq_struct *pairing4;
} a_kani_prime_input;
typedef struct {const root_extension *field;fq_struct *points;size_t count;} a_kani_group;
typedef struct {ulong ell;unsigned exponent;const root_extension *field;fq_struct *pairing4;} a_kani_prime;
typedef struct {
 const root_extension *base;b_theta4_ctx theta;fq_t zeta4,zeta8;
 a_kani_group *groups;size_t group_count,payload_count;
 a_kani_prime *primes;size_t prime_count,completed;
 fmpz_t polarization_multiplier;int ready,last_error;
} a_kani_fullchain;
typedef void (*a_kani_progress)(void *context,size_t prime_index,ulong ell,unsigned completed,unsigned total);
/* Distinct prime components; 2, if present, must be first. Binary input K has
 * order 2^(e+2), and L/pairing4 are unused. Odd K,L have order ell^e and their
 * initial e_{ell^e}(K_i,L_j)^(4ell^(e-1)) is pairing4. Every component and
 * payload begins in the same base-changed source marking. Original payload
 * group indices stay unchanged; internal full frames follow them. */
int a_kani_fullchain_init(a_kani_fullchain *S,unsigned g,const root_extension *base,
 const fq_struct *null,const fq_t zeta4,const fq_t zeta8,
 const a_kani_prime_input *primes,size_t prime_count,
 const a_theta_payload_input *payloads,size_t payload_count);
/* Executes one full prime-power component with stage progress notifications.
 * Atomic across that component: failure preserves the previous committed
 * target and every payload/frame. It never samples a new ideal or kernel.
 * last_error=-2 means RNG failure in an auxiliary field/embedding;
 * -1 means another failure or no next component; 0 means success. */
int a_kani_fullchain_next(a_kani_fullchain *S,a_kani_progress progress,void *context);
void a_kani_fullchain_clear(a_kani_fullchain *S);
#endif
