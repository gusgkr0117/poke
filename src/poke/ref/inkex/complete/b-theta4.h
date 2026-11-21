#ifndef INKE_RESEARCH_THETA4_H
#define INKE_RESEARCH_THETA4_H
#include <flint/fq.h>
#include <stddef.h>
#include <stdint.h>

/* Research implementation of the level-4 Riemann equations.
 * Coordinates are indexed by little-endian base-4 digits. Affine lifts are
 * retained: routines never normalize the coordinate vector implicitly. */
typedef struct {
  uint32_t k, l;
  int ready;
  fq_t inverse_denominator;
} b_theta4_chart;
typedef struct {
  const fq_ctx_struct *field;
  unsigned g;
  size_t count, characters;
  fq_struct *null;
  b_theta4_chart *charts;
} b_theta4_ctx;

fq_struct *b_theta4_vector_new(size_t n, const fq_ctx_t field);
void b_theta4_vector_free(fq_struct *v, size_t n, const fq_ctx_t field);
int b_theta4_init(b_theta4_ctx *T, unsigned g, const fq_struct *null,
                  const fq_ctx_t field);
void b_theta4_clear(b_theta4_ctx *T);
void b_theta4_copy(fq_struct *out, const fq_struct *in, const b_theta4_ctx *T);
int b_theta4_equal_projective(const fq_struct *a, const fq_struct *b,
                              const b_theta4_ctx *T);
int b_theta4_nonzero(const fq_struct *v, const b_theta4_ctx *T);
void b_theta4_inverse(fq_struct *out, const fq_struct *in,
                      const b_theta4_ctx *T);
/* Action convention: out_i = zeta^(beta dot i) in_(i+alpha).
 * zeta is a fixed primitive fourth root in the field. */
void b_theta4_action(fq_struct *out, const fq_struct *in, uint32_t alpha,
                     uint32_t beta, const fq_t zeta, const b_theta4_ctx *T);
/* Complete finite chart search. Returns 0 if the input theta null does not
 * provide a nonzero Riemann chart; no division by zero is attempted. */
/* Complete projective addition from full level-four coordinates. The output
 * scale is unspecified; no difference point, roots or earlier isogeny
 * preimages are needed. Aliasing is allowed; failure preserves out. */
int b_theta4_add(fq_struct *out, const fq_struct *P, const fq_struct *Q,
                 b_theta4_ctx *T);
int b_theta4_diff_add(fq_struct *out, const fq_struct *P, const fq_struct *Q,
                      const fq_struct *PmQ, b_theta4_ctx *T);
int b_theta4_three_way_add(fq_struct *out, const fq_struct *P,
                           const fq_struct *Q, const fq_struct *R,
                           const fq_struct *PQ, const fq_struct *QR,
                           const fq_struct *PR, b_theta4_ctx *T);
int b_theta4_scalar(fq_struct *out, const fq_struct *P, ulong scalar,
                    b_theta4_ctx *T);
/* Returns nP+Q and nP from consistent projective inputs P+Q,P,Q.
 * Retains the affine scaling law of DL Lemma7 Eq65. */
int b_theta4_scalar_add(fq_struct *nPQ, fq_struct *nP, ulong n,
                        const fq_struct *PQ, const fq_struct *P,
                        const fq_struct *Q, b_theta4_ctx *T);
/* Computes c with target=c*source, checking every coordinate. */
int b_theta4_scale_factor(fq_t c, const fq_struct *source,
                          const fq_struct *target, const b_theta4_ctx *T);
/* DL d=2 good lifts. alpha is the standard first-Lagrangian coordinate of
 * 2P. The basis routine returns -1 if symmetric compatibility fails, -2
 * if the requested fourth root is absent, 0 for invalid input/arithmetic. */
int b_theta4_good_basis2(fq_struct *out, const fq_struct *P, uint32_t alpha,
                         const fq_t zeta, b_theta4_ctx *T);
/* Given good P, arbitrary Q and P+Q, rescale P+Q so that ScalarMult(2)
 * is ThetaAct(Q,alpha). For Q another good basis point this yields its
 * pair lift; for Q a general evaluation point it yields the coset lift. */
int b_theta4_good_coset2(fq_struct *out, const fq_struct *PQ,
                         const fq_struct *P, const fq_struct *Q, uint32_t alpha,
                         const fq_t zeta, b_theta4_ctx *T);
/* Same normalization without extracting a square root. If P represents a
 * good lift multiplied by r_P^-1, pass sP=r_P^2. Q is kept as supplied.
 * Returns s with sqrt(s)*PQ a good coset lift relative to the good P.
 * The returned squared scale is sufficient for degree-2 evaluation. */
int b_theta4_good_coset2_squared_scale(fq_t s, const fq_struct *PQ,
    const fq_struct *P,const fq_t sP,const fq_struct *Q,uint32_t alpha,
    const fq_t zeta,b_theta4_ctx *T);
/* DL Algorithm8 m=4,d=2, standard first-Lagrangian kernel. shifted is a
 * good lift of x+sum(parity_i*T_i), with 2*T_i standard4torsion. Writes
 * only indices of the requested parity. Call for every parity [0,2^g).
 * The caller can stream shifts instead of storing all2^g vectors. */
void b_theta4_two_standard_block(fq_struct *out,unsigned parity,
    const fq_struct *shifted,const b_theta4_ctx *T);
void b_theta4_two_standard_scaled_block(fq_struct *out,unsigned parity,
    const fq_struct *shifted,const fq_t squared_scale,const b_theta4_ctx *T);
/* Complete trace-section charts: seed is a binary vector in [0,2^g).
 * Seed i uses theta_i*theta_-i instead of theta_0^2. */
void b_theta4_two_trace_block(fq_struct *out,unsigned parity,unsigned seed,
    const fq_struct *shifted,const fq_t squared_scale,const b_theta4_ctx *T);
/* Natural genus-1 level-4 -> level-2 Kummer projection. */
int b_theta4_to_kummer(fq_t U, fq_t V, const fq_struct *P,
                       const b_theta4_ctx *T);
/* Linear Fourier change on one tensor axis. Same change must be applied to
 * null and all points. This helper does not update torsion labels. */
void b_theta4_fourier_axis(fq_struct *out, const fq_struct *in, unsigned axis,
                           const fq_t zeta, const b_theta4_ctx *T);
#endif
