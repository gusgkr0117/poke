// SPDX-License-Identifier: Apache-2.0
#include <inke_uniform.h>
#include "../internal.h"
#include <stdio.h>
#include <stdlib.h>

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); exit(1); \
} } while (0)

static void init_order(quat_order_t *order)
{
    quat_lattice_init(order);
    mpz_set_ui(order->denom, 2);
    mpz_set_ui(order->basis[0][0], 2);
    mpz_set_ui(order->basis[1][1], 2);
    mpz_set_ui(order->basis[1][2], 1);
    mpz_set_ui(order->basis[2][2], 1);
    mpz_set_ui(order->basis[0][3], 1);
    mpz_set_ui(order->basis[3][3], 1);
}

static void check_ideal(const quat_left_ideal_t *ideal, const ibz_t *d,
                        const quat_order_t *order, const quat_alg_t *algebra)
{
    CHECK(mpz_cmp(ideal->norm, *d) == 0);
    ibz_t index, square, content;
    ibz_init(&index); ibz_init(&square); ibz_init(&content);
    quat_alg_coord_t coords;
    quat_alg_coord_init(&coords);
    quat_alg_elem_t x, o, product;
    quat_alg_elem_init(&x); quat_alg_elem_init(&o); quat_alg_elem_init(&product);
    mpz_set(x.denom, ideal->lattice.denom);
    mpz_set(o.denom, order->denom);
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) mpz_set(x.coord[i], ideal->lattice.basis[i][j]);
        CHECK(quat_lattice_contains(&coords, order, &x, algebra));
        for (int i = 0; i < 4; ++i) mpz_gcd(content, content, coords[i]);
        for (int k = 0; k < 4; ++k) {
            for (int i = 0; i < 4; ++i) mpz_set(o.coord[i], order->basis[i][k]);
            quat_alg_mul(&product, &o, &x, algebra);
            CHECK(quat_lattice_contains(NULL, &ideal->lattice, &product, algebra));
        }
    }
    CHECK(mpz_cmp_ui(content, 1) == 0);
    quat_lattice_index(&index, &ideal->lattice, order);
    mpz_mul(square, *d, *d);
    CHECK(mpz_cmp(index, square) == 0);
    quat_alg_elem_finalize(&product); quat_alg_elem_finalize(&o); quat_alg_elem_finalize(&x);
    quat_alg_coord_finalize(&coords);
    ibz_finalize(&content); ibz_finalize(&square); ibz_finalize(&index);
}

/* Exhaust all residues, including nonunits, so this checks equal fiber
 * cardinalities rather than relying on a statistical uniformity test. */
static void test_fibers(unsigned int modulus, unsigned int expected_ideals,
                         unsigned int expected_fiber)
{
    ibz_t p, d;
    ibz_init(&p); ibz_init(&d);
    mpz_set_ui(p, 7); mpz_set_ui(d, modulus);
    quat_alg_t algebra;
    quat_alg_init_set(&algebra, &p);
    quat_order_t order;
    init_order(&order);
    quat_alg_elem_t gamma, unit;
    quat_alg_elem_init(&gamma); quat_alg_elem_init(&unit);
    if (modulus == 5) {
        mpz_set_ui(gamma.coord[0], 2); mpz_set_ui(gamma.coord[1], 1);
    } else {
        mpz_set_ui(gamma.coord[0], 1); mpz_set_ui(gamma.coord[1], 1);
        mpz_set_ui(gamma.coord[2], 1);
    }
    quat_left_ideal_t ideal, representatives[12];
    quat_left_ideal_init(&ideal);
    for (int i = 0; i < 12; ++i) quat_left_ideal_init(&representatives[i]);
    unsigned int counts[12] = {0}, distinct = 0, units = 0;
    quat_alg_coord_t coords;
    quat_alg_coord_init(&coords);
    ibq_t norm;
    ibq_init(&norm);
    ibz_t n, gcd;
    ibz_init(&n); ibz_init(&gcd);
    unsigned int total = modulus * modulus * modulus * modulus;
    for (unsigned int value = 0; value < total; ++value) {
        unsigned int digits = value;
        for (int i = 0; i < 4; ++i) {
            mpz_set_ui(coords[i], digits % modulus);
            digits /= modulus;
        }
        ibz_mat_4x4_eval(&unit.coord, &order.basis, &coords);
        mpz_set(unit.denom, order.denom);
        quat_alg_norm(&norm, &unit, &algebra);
        CHECK(ibq_to_ibz(&n, &norm));
        mpz_gcd(gcd, n, d);
        int status = inke_ideal_from_unit(&ideal, &gamma, &unit, &d, &order, &algebra);
        if (mpz_cmp_ui(gcd, 1) != 0) {
            CHECK(status == INKE_SAMPLE_INVALID);
            continue;
        }
        CHECK(status == INKE_SAMPLE_OK);
        check_ideal(&ideal, &d, &order, &algebra);
        ++units;
        unsigned int j = 0;
        while (j < distinct && !quat_lideal_equals(&ideal, &representatives[j], &algebra)) ++j;
        if (j == distinct) {
            CHECK(distinct < 12);
            ibz_mat_4x4_copy(&representatives[j].lattice.basis, &ideal.lattice.basis);
            mpz_set(representatives[j].lattice.denom, ideal.lattice.denom);
            mpz_set(representatives[j].norm, ideal.norm);
            representatives[j].parent_order = &order;
            ++distinct;
        }
        ++counts[j];
    }
    CHECK(distinct == expected_ideals);
    CHECK(units == expected_ideals * expected_fiber);
    for (unsigned int i = 0; i < distinct; ++i) CHECK(counts[i] == expected_fiber);
    printf("d=%u: %u units, %u ideals, exactly %u preimages each\n",
           modulus, units, distinct, expected_fiber);
    ibz_finalize(&gcd); ibz_finalize(&n); ibq_finalize(&norm);
    quat_alg_coord_finalize(&coords);
    for (int i = 0; i < 12; ++i) quat_left_ideal_finalize(&representatives[i]);
    quat_left_ideal_finalize(&ideal);
    quat_alg_elem_finalize(&unit); quat_alg_elem_finalize(&gamma);
    quat_lattice_finalize(&order); quat_alg_finalize(&algebra);
    ibz_finalize(&d); ibz_finalize(&p);
}

static void test_rng_boundaries(void)
{
    const unsigned int sizes[] = {1, 7, 63, 64, 65, 127, 128, 129, 255, 256, 257};
    const unsigned char seed[] = "INKE interval boundary regression";
    shake256ctx state;
    shake256_absorb(&state, seed, sizeof(seed));
    ibz_t lo, hi, x;
    ibz_init(&lo); ibz_init(&hi); ibz_init(&x);
    const int64_t signed_values[] = {0, 1, -1, INT64_C(140227657289781369),
        -INT64_C(140227657289781369), INT64_MAX, INT64_MIN};
    for (unsigned int i = 0; i < sizeof(signed_values) / sizeof(signed_values[0]); ++i) {
        ibz_set(&x, signed_values[i]);
        CHECK(ibz_get(&x) == signed_values[i]);
    }
    CHECK(ibz_rand_interval_i(&x, INT64_MAX, INT64_MAX));
    CHECK(ibz_get(&x) == INT64_MAX);
    CHECK(ibz_rand_interval_i(&x, INT64_MIN, INT64_MIN));
    CHECK(ibz_get(&x) == INT64_MIN);
    for (int i = 0; i < 64; ++i) {
        CHECK(ibz_rand_interval_minm_m(&x, INT64_MAX));
        ibz_set(&hi, INT64_MAX);
        CHECK(mpz_cmpabs(x, hi) <= 0);
    }
    mpz_set_si(lo, -19);
    for (unsigned int i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i) {
        mpz_set_ui(hi, 1); mpz_mul_2exp(hi, hi, sizes[i]);
        mpz_sub_ui(hi, hi, 1);
        int high_seen[2] = {0};
        for (int j = 0; j < 64; ++j) {
            CHECK(ibz_rand_interval(&x, &ibz_const_zero, &hi));
            CHECK(mpz_sgn(x) >= 0 && mpz_cmp(x, hi) <= 0);
            high_seen[0] |= mpz_tstbit(x, sizes[i] - 1);
            CHECK(ibz_rand_interval_with_state(&x, &ibz_const_zero, &hi, &state));
            CHECK(mpz_sgn(x) >= 0 && mpz_cmp(x, hi) <= 0);
            high_seen[1] |= mpz_tstbit(x, sizes[i] - 1);
        }
        CHECK(high_seen[0] && high_seen[1]);
    }
    CHECK(ibz_rand_interval(&x, &lo, &lo) && mpz_cmp(x, lo) == 0);
    CHECK(ibz_rand_interval_with_state(&x, &lo, &lo, &state) && mpz_cmp(x, lo) == 0);
    CHECK(!ibz_rand_interval(&x, &hi, &lo));
    CHECK(!ibz_rand_interval_with_state(&x, &hi, &lo, &state));
    ibz_finalize(&x); ibz_finalize(&hi); ibz_finalize(&lo);
    puts("RNG interval boundaries: passed");
}

static void test_manuscript_parameters(void)
{
    const unsigned int a[] = {256, 384, 512};
    const unsigned int f2[] = {257, 386, 512};
    const unsigned int f3[] = {163, 244, 324};
    const unsigned int cofactor[] = {11, 23, 517};
    ibz_t p, power, q, saved_q, d, saved_norm;
    ibz_init(&p); ibz_init(&power); ibz_init(&q); ibz_init(&saved_q);
    ibz_init(&d); ibz_init(&saved_norm);
    for (int level = 0; level < 3; ++level) {
        mpz_ui_pow_ui(p, 3, f3[level]);
        mpz_mul_2exp(p, p, f2[level]);
        mpz_mul_ui(p, p, cofactor[level]); mpz_sub_ui(p, p, 1);
        CHECK(mpz_probab_prime_p(p, 32));
        quat_alg_t algebra;
        quat_alg_init_set(&algebra, &p);
        quat_order_t order;
        init_order(&order);
        quat_left_ideal_t ideal;
        quat_left_ideal_init(&ideal);
        mpz_set_ui(power, 1); mpz_mul_2exp(power, power, a[level]);
        for (int sample = 0; sample < 3; ++sample) {
            if (sample == 0) mpz_set_ui(q, 5);
            else if (sample == 1) mpz_sub_ui(q, power, 5);
            else CHECK(inke_sample_q(&q, a[level]) == INKE_SAMPLE_OK);
            CHECK(mpz_fdiv_ui(q, 6) == 5 && mpz_cmp(q, power) < 0);
            mpz_set(saved_q, q);
            mpz_sub(d, power, q); mpz_mul(d, d, q);
            inke_sample_counts counts = {0};
            CHECK(inke_sample_ideal(&ideal, &q, a[level], &order, &algebra,
                                    200000, 20000, 1024, &counts) == INKE_SAMPLE_OK);
            CHECK(mpz_cmp(q, saved_q) == 0);
            check_ideal(&ideal, &d, &order, &algebra);
            mpz_set(saved_norm, ideal.norm);
            CHECK(inke_sample_ideal(&ideal, &q, a[level], &order, &algebra,
                                    0, 0, 0, NULL) == INKE_SAMPLE_RETRY);
            CHECK(mpz_cmp(q, saved_q) == 0 && mpz_cmp(ideal.norm, saved_norm) == 0);
            printf("a=%u sample=%d: ideal verified; prime/norm/unit candidates=%lu/%lu/%lu\n",
                   a[level], sample, counts.prime_candidates,
                   counts.norm_candidates, counts.unit_candidates);
        }
        quat_left_ideal_finalize(&ideal);
        quat_lattice_finalize(&order); quat_alg_finalize(&algebra);
    }
    ibz_finalize(&saved_norm); ibz_finalize(&d); ibz_finalize(&saved_q);
    ibz_finalize(&q); ibz_finalize(&power); ibz_finalize(&p);
}

int main(void)
{
    test_rng_boundaries();
    test_fibers(5, 6, 80);
    test_fibers(9, 12, 324);
    test_manuscript_parameters();
    puts("All INKE uniform ideal sampling checks passed (not a full Keygen test).");
    return 0;
}
