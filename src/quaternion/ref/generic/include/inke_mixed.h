// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_MIXED_H
#define INKE_MIXED_H
#include <inke_uniform.h>
#include <stddef.h>

/* Whole 3-byte draws: retained only by this initialized, non-copyable state. */
#define INKE_MIXED_RANDOM_BUFFER_BYTES 12288u

/* Arithmetic stage of the mixed 2/3 conversion. No curve/isogeny evaluation
 * is performed here. A successful witness is checked by exact arithmetic;
 * a bounded search does not assert a success probability or runtime theorem. */
typedef struct {
    quat_left_ideal_t input, reduced;
    quat_alg_elem_t beta;
    ibz_mat_4x4_t basis;
    ibz_t target;
    /* Nrd(B*c/D) = sum_{i<=j} norm_form[i][j]*c_i*c_j / norm_denom.
     * Off-diagonal entries include the factor two; no integrality of the
     * individual polynomial coefficients is assumed. */
    ibz_mat_4x4_t norm_form;
    ibz_t norm_denom;
    unsigned char random_buffer[INKE_MIXED_RANDOM_BUFFER_BYTES];
    size_t random_position;
    const quat_alg_t *algebra;
    unsigned long attempts, generators, positive_targets, square_solutions;
    int ready;
} inke_mixed_state;

typedef struct {
    quat_alg_elem_t mu1, mu2, theta, transport;
    ibz_t d1, d2;
} inke_mixed_witness;

void inke_mixed_state_init(inke_mixed_state *state);
void inke_mixed_state_clear(inke_mixed_state *state);
void inke_mixed_witness_init(inke_mixed_witness *witness);
void inke_mixed_witness_clear(inke_mixed_witness *witness);

/* Snapshot I, compute J=I*conj(beta)/Nrd(I), and reduce its basis.
 * target must be positive and divisible by 6. order/algebra must outlive
 * state. prepare() is not called again when a search budget is exhausted. */
int inke_mixed_prepare(inke_mixed_state *state, const quat_left_ideal_t *ideal,
                       const ibz_t *target, const quat_alg_t *algebra);

/* Resume the same state for at most budget generator candidates. Neither
 * input nor reduced ideal is replaced on retry. out is initialized and is
 * meaningful only on INKE_SAMPLE_OK. Buffered OS randomness survives a
 * bounded retry, but is erased on clear, reprepare, or an error. */
int inke_mixed_search(inke_mixed_witness *out, inke_mixed_state *state,
                      unsigned long budget);

/* Independently check norm, membership, and
 * (J*conj(mu1)/Nrd(J))*transport = d1*I.
 * Auxiliary J*conj(mu1)/Nrd(J) need not be primitive. */
int inke_mixed_verify(const inke_mixed_state *state,
                      const inke_mixed_witness *witness);
#endif
