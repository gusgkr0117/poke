#ifndef INKE_RESEARCH_GAUSSIAN_PREPARATION_H
#define INKE_RESEARCH_GAUSSIAN_PREPARATION_H
#include <quaternion.h>
/* Fixed-parameter plans can be prepared once and shared read-only. Their
 * construction is public precomputation, whose cost must be accounted for
 * consistently in comparisons. Plans MUST be zero-initialized before their
 * first init call. certificate_init initializes certificate storage. */
typedef struct {void *plans;unsigned count;int ready;} b_gaussian_plans;
typedef struct {
    quat_alg_elem_t beta1,beta2,x,y,eta,transport;
    quat_left_ideal_t Ia;
    ibz_t d1,d2,u,v,N,n1,n2;
    unsigned primes[11],exponents[11],extensions[11],common_extension;
    int verified;
} b_gaussian_certificate;
typedef struct {
    unsigned candidate_norms,selected_plan;
    unsigned long two_square_tests,first_square_successes;
} b_gaussian_counts;
enum {B_GAUSSIAN_INVALID=-1,B_GAUSSIAN_NO_MEMORY=-2,
      B_GAUSSIAN_MISS=0,B_GAUSSIAN_OK=1};
int b_gaussian_plans_init(b_gaussian_plans *plans);
void b_gaussian_plans_clear(b_gaussian_plans *plans);
void b_gaussian_certificate_init(b_gaussian_certificate *certificate);
void b_gaussian_certificate_clear(b_gaussian_certificate *certificate);
/* Read-only input I must be an integral ideal of the linked E0 order.
 * Bounded exact search: no random resampling, no input substitution.
 * On MISS use the complete converter on precisely this same I.
 * All failures preserve out. counts may be NULL; otherwise overwritten.
 * The initialized plan table may be shared by concurrent calls. */
int b_gaussian_prepare(b_gaussian_certificate *out,const quat_left_ideal_t *I,
    const b_gaussian_plans *plans,b_gaussian_counts *counts);
/* Exact certificate check against the originally selected I, including
 * Ia*transport=n1*I, not merely quaternion norms or isomorphism class. */
int b_gaussian_certificate_verify(const b_gaussian_certificate *certificate,
    const quat_left_ideal_t *I);
#endif
