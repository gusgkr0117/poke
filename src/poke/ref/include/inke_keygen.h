// SPDX-License-Identifier: Apache-2.0
#ifndef INKE_UNIFORM_KEYGEN_H
#define INKE_UNIFORM_KEYGEN_H
#include <inke.h>
#include <inke_mixed.h>
#include <inke_conversion.h>

/* Stateful C integration of uniform ideal sampling and the mixed conversion.
 * The public exponent comes from the selected parameter build. These APIs do not
 * implement the general dimension-eight fallback or its runtime guarantee.
 * RETRY means resume this context: never discard q or a sampled ideal.
 * sk/pk are written only on success. Use the matching encryption and
 * decryption APIs for these parameters, bases and key encodings. */
typedef struct {
    ibz_t q;
    inke_mixed_state mixed;
    inke_sample_counts sampling;
    unsigned long witnesses,conversion_failures,kani_failures;
    int product_steps;
    unsigned int a;
    int has_q,has_ideal,complete,last_stage;
} inke_keygen_context;

void inke_keygen_init(inke_keygen_context *ctx);
void inke_keygen_clear(inke_keygen_context *ctx);
int inke_keygen_uniform(inke_sk_t *sk,inke_pk_t *pk,
    inke_keygen_context *ctx,unsigned long norm_budget);

/* Complete Keygen from a previously verified conversion of the ORIGINAL
 * selected ideal of norm q*(2^a-q). conversion->two must be the image of
 * exact BASIS_EVEN at full TORSION_PLUS_EVEN_POWER, under one target map.
 * Checks q and full-order degree pairing, then uses the existing factor,
 * masking and canonical-model code. Outputs are atomic. On RETRY keep both
 * q and this conversion; never resample a receiver input. stage is optional.
 * This boundary does not assert completeness of the caller's converter. */
int inke_keygen_from_conversion(inke_sk_t *sk,inke_pk_t *pk,const ibz_t *q,
    const inke_conversion_result *conversion,int *stage);

/* Optional complete product-isogeny implementation. left_two and right_two
 * are corresponding bases of order 2^(exponent+2). Multiplying each of the
 * two lifts (-q*left_two, right_two) by 4 gives the kernel generators. Return the
 * elliptic factor on which the first lift has order 4; its other projection
 * is zero. Optional left/right_points are
 * mapped from (P,0)/(0,P), respectively, using ONE target identification for
 * every image. A NULL points input has a NULL corresponding output.
 * Return INKE_SAMPLE_OK/RETRY/INVALID/RNG_ERROR and preserve outputs unless
 * successful. No receiver ideal or secret is sampled by this callback.
 * This interface alone does not assert that a complete implementation is
 * linked. The ordinary APIs below retain the existing fast implementation. */
typedef int (*inke_product_fallback)(ec_curve_t *target,
    ec_basis_t *left_images,ec_basis_t *right_images,
    const ec_curve_t *left_curve,const ec_curve_t *right_curve,
    const ec_basis_t *left_two,const ec_basis_t *right_two,
    unsigned int exponent,const ibz_t *q,
    const ec_basis_t *left_points,const ec_basis_t *right_points,
    void *context,int *stage);
int inke_keygen_from_conversion_with_fallback(inke_sk_t *sk,inke_pk_t *pk,
    const ibz_t *q,const inke_conversion_result *conversion,
    inke_product_fallback fallback,void *fallback_context,int *stage);

/* Matching PKE integration for 32-byte messages, including all cyclic
 * 3^b-subgroups and selection of the degree-q factor during decryption. */
int inke_encrypt_uniform(inke_ct_t *ct,const inke_pk_t *pk,
    const unsigned char message[32],const unsigned char *seed,size_t seed_len);
int inke_decrypt_uniform(unsigned char message[32],const inke_ct_t *ct,
    const inke_sk_t *sk);
int inke_decrypt_uniform_with_fallback(unsigned char message[32],
    const inke_ct_t *ct,const inke_sk_t *sk,
    inke_product_fallback fallback,void *fallback_context,int *stage);
#endif
