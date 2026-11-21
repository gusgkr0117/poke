#ifndef INKE_RESEARCH_CLAPOTI_PREPARATION_H
#define INKE_RESEARCH_CLAPOTI_PREPARATION_H
#include <quaternion.h>
#include <stddef.h>
typedef struct {
 quat_alg_elem_t beta1,beta2,theta;
 ibz_t d1,d2,u,v,N,n1,n2;ibz_mat_4x4_t X,Y;
 unsigned long *primes;unsigned *exponents,*extensions;size_t prime_count;
 int verified;
} a_clapoti_certificate;
typedef struct {unsigned long second_candidates,four_square_trials;int used_complete_residue_box;} a_clapoti_counts;
enum {A_CLAPOTI_OK=1,A_CLAPOTI_INVALID=-1,A_CLAPOTI_RNG_ERROR=-2,A_CLAPOTI_NO_MEMORY=-3};
void a_clapoti_certificate_init(a_clapoti_certificate *C);
void a_clapoti_certificate_clear(a_clapoti_certificate *C);
/* Caller guarantees I is a valid integral left O0 ideal. The boundary checks
 * integrality/index/norm but does not independently test left-order closure.
 * Complete expected-polynomial algebraic preparation on precisely I.
 * The first element is found by finite local avoidance; the second uses
 * the complete residue box after a 32-trial optional prefix. Four-square
 * calls use budget0 (unconditional Las Vegas mode). No q/I is resampled.
 * The linked POKE system RNG supplies all random bytes. RNG failure is
 * explicit, not rejection of the selected ideal. out is unchanged on failure. */
int a_clapoti_prepare(a_clapoti_certificate *out,const quat_left_ideal_t *I,a_clapoti_counts *counts);
/* Checks exact ideal transport, every Kani identity, exact N factorization,
 * and sufficient torsion-field containment. No particular planner profile
 * is required: a smaller valid certificate proves the same selected map. */
int a_clapoti_certificate_verify(const a_clapoti_certificate *C,const quat_left_ideal_t *I);
/* Additional constraints used by the complete polynomial-extension planner:
 * retained full planned 2/3 part and distinct exponent-one auxiliary primes.
 * Production prepare checks this stronger predicate before committing. */
int a_clapoti_complete_profile_verify(const a_clapoti_certificate *C,const quat_left_ideal_t *I);
#endif
