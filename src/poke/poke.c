#include <hd.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <klpt.h>
#include <quaternion.h>
#include <gmp.h>
#include <intbig.h>
#include <ec.h>
#include <rng.h>
#include <id2iso.h>

int ibz_random_unit(ibz_t *q, const ibz_t *modulus) {
    ibz_t gcd;
    ibz_init(&gcd);
    while (1) {
        ibz_rand_interval(q, &ibz_const_zero, modulus);
        ibz_gcd(&gcd, q, modulus);
        if (ibz_is_one(&gcd)) {
            break;
        }
    }
    return 1;
}

int test() {
    ibz_t q, alpha, beta, gamma, delta, rhs, deg;
    ibz_t A;
    
    ibz_init(&q); ibz_init(&alpha); ibz_init(&beta);
    ibz_init(&gamma); ibz_init(&delta); ibz_init(&rhs); ibz_init(&deg);
    ibz_init(&A);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    // Set q to a random value in the range [0, TORSION_PLUS_2POWER)
    for(int i = 0; i < 1000; i++) {
        ibz_rand_interval(&q, &ibz_const_zero, &A);
        ibz_sub(&deg, &A, &q);
        if (ibz_divides(&q, &ibz_const_two) == 0 && ibz_divides(&q, &ibz_const_three) == 0 && ibz_divides(&q, &ibz_const_five) == 0
            && ibz_divides(&deg, &ibz_const_three) == 0 && ibz_divides(&deg, &ibz_const_five) == 0) {
            break;
        }
    }
    ibz_mul(&rhs, &deg, &q);
    ibz_mul(&rhs, &rhs, &TORSION_PLUS_3POWER);
    ibz_random_unit(&alpha, &A);
    ibz_random_unit(&beta, &A);
    ibz_random_unit(&gamma, &TORSION_PLUS_3POWER);
    ibz_random_unit(&delta, &TORSION_PLUS_5POWER);

    gmp_printf("q = %Zd\n", q);

    gmp_printf("alpha = %Zd\n", alpha);
    gmp_printf("beta = %Zd\n", beta);
    gmp_printf("gamma = %Zd\n", gamma);
    gmp_printf("delta = %Zd\n", delta);
    ec_curve_t curve = CURVE_E0;
    fp2_t j_inv;
    ec_j_inv(&j_inv, &curve);
    fp2_print("j_invariant",&j_inv);

    quat_alg_elem_t tau;
    quat_alg_elem_init(&tau);
    represent_integer(&tau, &rhs, &QUATALG_PINFTY);

    quat_alg_elem_print(&tau);

    ec_basis_t E0_two, E0_three;
    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);

    endomorphism_application_even_basis(&E0_two, &curve, &tau, TORSION_PLUS_EVEN_POWER);
    endomorphism_application_even_basis(&E0_three, &curve, &tau, TORSION_PLUS_ODD_POWERS[0]);

    point_print("tau(P2)", E0_two.P);
    point_print("tau(Q2)", E0_two.Q);
    point_print("tau(P3)", E0_three.P);
    point_print("tau(Q3)", E0_three.Q);

    ec_point_t kernel_point;
    ec_isog_odd_t isog;
    isog.curve = curve;

    ibz_t three_m1_order, remainder;
    ibz_init(&three_m1_order);
    ibz_init(&remainder);
    ibz_div(&three_m1_order, &remainder, &TORSION_PLUS_3POWER, &ibz_const_three);
    ec_point_t test_point;
    ec_mul_ibz(&test_point, &curve, &three_m1_order, &E0_three.P);
    if (ec_is_zero(&test_point)) {
        copy_point(&isog.ker_plus, &E0_three.Q);
    } else {
        copy_point(&isog.ker_plus, &E0_three.P);
    }

    ec_mul_ibz(&test_point, &curve, &three_m1_order, &isog.ker_plus);
    if (ec_is_zero(&test_point)){
        printf("the order of the kernel point is not 3^m\n");
        return 1;
    }
    
    isog.degree[0] = TORSION_PLUS_ODD_POWERS[0];
    ec_curve_t E1;

    ec_eval_odd_basis(&E1, &isog, &E0_two, 2);

    curve_print("Image curve", E1);

    // Evaluating the theta-based 2-dim isogeny
    theta_couple_curve_t E01;
    theta_couple_point_t T1, T2, T1m2;
    theta_chain_t hd_isog;

    E01.E1 = curve; E01.E2 = E1;

    copy_point(&T1.P1, &BASIS_EVEN.P);
    copy_point(&T2.P1, &BASIS_EVEN.Q);
    copy_point(&T1m2.P1, &BASIS_EVEN.PmQ);
    ec_mul_ibz(&T1.P1, &curve, &q, &T1.P1);
    ec_neg(&T1.P1, &T1.P1);
    ec_mul_ibz(&T2.P1, &curve, &q, &T2.P1);
    ec_neg(&T2.P1, &T2.P1);
    ec_mul_ibz(&T1m2.P1, &curve, &q, &T1m2.P1);
    ec_neg(&T1m2.P1, &T1m2.P1);
    ibz_t inverse;
    ibz_init(&inverse);
    ibz_invmod(&inverse, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);

    ec_mul_ibz(&E0_two.P, &E1, &inverse, &E0_two.P);
    ec_mul_ibz(&E0_two.Q, &E1, &inverse, &E0_two.Q);
    ec_mul_ibz(&E0_two.PmQ, &E1, &inverse, &E0_two.PmQ);

    T1.P2 = E0_two.P;
    T2.P2 = E0_two.Q;
    T1m2.P2 = E0_two.PmQ;
    // Check the order of the T1.P2
    ibz_t two_m1_order;
    ibz_init(&two_m1_order);
    ibz_div(&two_m1_order, &remainder, &TORSION_PLUS_2POWER, &ibz_const_two);
    ec_mul_ibz(&test_point, &curve, &three_m1_order, &T1.P1);
    if (ec_is_zero(&test_point)){
        printf("the order of the kernel point is not 3^m\n");
        return 1;
    }

    ec_mul_ibz(&test_point, &curve, &three_m1_order, &T1.P2);
    if (ec_is_zero(&test_point)){
        printf("the order of the kernel point is not 3^m\n");
        return 1;
    }

    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER-2, &E01, &T1, &T2, &T1m2, strategies[2], 1);

    return 0; 
}

int main() {
    int res = 1;

    test();

    return res;
}