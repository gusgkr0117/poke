#ifndef INKE_RESEARCH_KANI_GEOMETRY_H
#define INKE_RESEARCH_KANI_GEOMETRY_H
#include "a-kani-fullchain.h"
#include "a-kani-two-kernel.h"
#include "root-ec-theta4.h"
#include "root-tate-projector.h"
typedef struct {
 unsigned dimension;const fmpz *n1,*n2,*N,*endomorphism;unsigned denominator;
 /* For dimension8: theta=endomorphism, n1=u*d1,n2=v*d2, eta=YX theta.
  * For dimension2: eta=endomorphism and these four pointers are unused. */
 const fmpz_mat_struct *X,*Y;const fmpz *d1,*d2;
 size_t prime_count;const ulong *primes;const unsigned *exponents,*extensions;
} a_kani_geometry_spec;
typedef struct {
 root_extension field;root_field_embedding embedding;root_ec_theta4 model;
 fq_struct *generators,*complements,*pairing4;
 const fq_struct *generator_ptrs[8],*complement_ptrs[8];
 int field_ready;
} a_kani_geometry_component;
typedef struct {
 const root_ec_theta4 *base_model;unsigned dimension;size_t count,prime_count;
 a_kani_geometry_component *components;a_kani_prime_input *prime_inputs;
 fq_struct *null,*lattice;const fq_struct *lattice_ptrs[16];
 fmpz_mat_t Pnum,eta;fmpz_t n1,N;fq_t zeta4,zeta8;
 a_kani_fullchain chain;int ready,started,last_error;
} a_kani_geometry;
/* Caller supplies a verified same-I certificate. The geometry additionally
 * checks its norm/similitude and exact degree-factorization identities.
 * Fields are owned here; the caller's base_model on E0/Fp2 must outlive G.
 * The first prime is2. Its extension must contain E0[2^(e2+3)] even when the
 * quaternion denominator is1, for coherent maximal-order action matrices.
 * budget0 gives unbounded uniform torsion sampling; finite exhaustion is
 * RETRY, not a reason to change I. RNG errors remain distinct.
 * For bool start/next calls, last_error is -2 for RNG failure, -1 for
 * other failure, and 0 on success. Init returns A_EC_* directly. */
int a_kani_geometry_init(a_kani_geometry *G,const root_ec_theta4 *base_model,
 const a_kani_geometry_spec *spec,a_ec_rng rng,void *context,unsigned long budget);
/* Tensor product of g copies of the SAME signed elliptic theta marking. */
int a_kani_geometry_product(fq_struct *out,const a_ec_point *points,unsigned g,const root_ec_theta4 *model);
/* Registers the internal 4/16 target-lattice points as chain group0.
 * Caller product payload groups follow at1, preserving their order. */
int a_kani_geometry_start(a_kani_geometry *G,const a_theta_payload_input *payloads,size_t groups);
int a_kani_geometry_next(a_kani_geometry *G,a_kani_progress progress,void *context);
void a_kani_geometry_clear(a_kani_geometry *G);
#endif
