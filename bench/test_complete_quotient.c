// SPDX-License-Identifier: Apache-2.0
#include "b-complete-keygen.h"
#include "inke_curve_model.h"
#include <curve_extras.h>
#include <endomorphism_action.h>
#include <id2iso.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <stdio.h>
#include <string.h>

static int equal_basis(const ec_basis_t *a, const ec_basis_t *b)
{
    return ec_is_equal(&a->P, &b->P) && ec_is_equal(&a->Q, &b->Q) &&
           ec_is_equal(&a->PmQ, &b->PmQ);
}

/* The elliptic model can differ. One isomorphism must match both triples. */
static int aligned_payloads(const ec_curve_t *target, const ec_basis_t actual[2],
                            const ec_curve_t *expected_curve, const ec_basis_t expected[2])
{
    inke_model_map a[6], b[6];
    unsigned an, bn;
    if (!inke_model_candidates(a, &an, target) ||
        !inke_model_candidates(b, &bn, expected_curve)) return 0;
    for (unsigned i = 0; i < an; ++i) {
        for (unsigned j = 0; j < bn; ++j) {
            if (!fp2_is_equal(&a[i].curve.A, &b[j].curve.A)) continue;
            ec_basis_t p[2], q[2];
            for (unsigned k = 0; k < 2; ++k) {
                if (!inke_model_apply_basis(&p[k], &actual[k], &a[i]) ||
                    !inke_model_apply_basis(&q[k], &expected[k], &b[j])) return 0;
            }
            if (equal_basis(&p[0], &q[0]) && equal_basis(&p[1], &q[1])) return 1;
        }
    }
    return 0;
}

/* For a known 3-isogeny E -> F and exponent two, the degree labels
 * q=1 and q=3 select E and F respectively. Compare both torsion payloads. */
static int check_factors(b_product_context *context)
{
    ec_curve_t source = CURVE_E0, right, target;
    ec_basis_t lift = BASIS_EVEN;
    unsigned shift = TORSION_PLUS_EVEN_POWER - 4;
    ec_dbl_iter(&lift.P, shift, &source, &lift.P);
    ec_dbl_iter(&lift.Q, shift, &source, &lift.Q);
    ec_dbl_iter(&lift.PmQ, shift, &source, &lift.PmQ);
    ec_basis_t images[2] = {lift, BASIS_THREE};
    ec_isog_odd_t phi = {0};
    phi.curve = source;
    phi.degree[0] = 1;
    ec_set_zero(&phi.ker_minus);
    ibz_t q;
    ibz_init(&q);
    mpz_divexact_ui(q, TORSION_PLUS_3POWER, 3);
    ec_mul_ibz(&phi.ker_plus, &source, &q, &BASIS_THREE.P);
    ec_eval_three(&right, &phi, (ec_point_t *)images, 6);
    ec_basis_t triple = BASIS_THREE;
    mpz_set_ui(q, 3);
    ec_mul_ibz(&triple.P, &source, &q, &triple.P);
    ec_mul_ibz(&triple.Q, &source, &q, &triple.Q);
    ec_mul_ibz(&triple.PmQ, &source, &q, &triple.PmQ);
    int ok = 0, stage;
    for (unsigned degree = 1; degree <= 3; degree += 2) {
        ec_basis_t actual[2];
        mpz_set_ui(q, degree);
        if (b_product_fallback(&target, &actual[0], &actual[1], &source, &right,
                &lift, &images[0], 2, &q, &BASIS_THREE, &images[1], context,
                &stage) != INKE_SAMPLE_OK) goto done;
        const ec_basis_t expected[2] = {degree == 1 ? BASIS_THREE : images[1],
                                       degree == 1 ? triple : images[1]};
        if (!aligned_payloads(&target, actual, degree == 1 ? &source : &right,
                              expected)) goto done;
        ec_curve_t curve_only;
        if (b_product_fallback(&curve_only, NULL, NULL, &source, &right, &lift,
                &images[0], 2, &q, NULL, NULL, context, &stage) != INKE_SAMPLE_OK)
            goto done;
        fp2_t lhs, rhs;
        ec_j_inv(&lhs, &target);
        ec_j_inv(&rhs, &curve_only);
        if (!fp2_is_equal(&lhs, &rhs)) goto done;
    }
    ec_curve_t sentinel;
    memset(&sentinel, 0xa5, sizeof sentinel);
    target = sentinel;
    mpz_set_ui(q, 2);
    ok = b_product_fallback(&target, NULL, NULL, &source, &right, &lift, &images[0],
                            2, &q, NULL, NULL, context, &stage) == INKE_SAMPLE_INVALID
         && !memcmp(&target, &sentinel, sizeof target);
done:
    ibz_finalize(&q);
    return ok;
}

static void report_quotient(void *operation, size_t index, unsigned long prime,
                            unsigned completed, unsigned total)
{
    (void)index;
    (void)prime;
    if (completed % 128 == 0 || completed == total) {
        printf("%s: %u/%u quotient steps\n", (const char *)operation, completed, total);
        fflush(stdout);
    }
}

static int forced_round_trip(b_complete_public *params)
{
    b_complete_keygen state = {0};
    inke_pk_t pk;
    inke_sk_t sk;
    inke_ct_t ciphertext;
    const unsigned char message[32] = "complete quotient round trip";
    const unsigned char seed[32] = "complete quotient test coins";
    unsigned char recovered[32];
    unsigned long calls = params->product.calls;
    int stage = 0, status;
    if (b_complete_keygen_init(&state, params, 0, 0, 0, 0) != INKE_SAMPLE_OK) return 0;
    do {
        status = b_complete_keygen_step(&sk, &pk, &state, report_quotient, "Key generation");
    } while (status == INKE_SAMPLE_RETRY);
    b_complete_keygen_clear(&state);
    if (status != INKE_SAMPLE_OK || params->product.calls <= calls) return 0;
    if (!inke_encrypt_uniform(&ciphertext, &pk, message, seed, sizeof seed)) return 0;
    calls = params->product.calls;
    params->product.progress = report_quotient;
    params->product.progress_context = "Decryption";
    return b_complete_decrypt(recovered, &ciphertext, &sk, params, &stage) &&
           params->product.calls > calls && !memcmp(message, recovered, sizeof message);
}

int main(void)
{
    b_complete_public params = {0};
    if (b_complete_public_init(&params, root_system_random, NULL, 10000) != INKE_SAMPLE_OK)
        return 1;
    int ok = check_factors(&params.product);
    if (ok) {
        puts("Verified both degree labels, torsion images and invalid-input handling.");
        fflush(stdout);
        ok = forced_round_trip(&params);
    }
    b_complete_public_clear(&params);
    if (!ok) {
        fputs("Complete quotient test failed.\n", stderr);
        return 1;
    }
    puts("Verified a forced complete quotient in key generation and decryption.");
    return 0;
}
