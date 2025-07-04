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
#include <biextension.h>

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

int ec_dlog_6(digit_t *scalarP, digit_t *scalarQ, ec_basis_t *base, ec_point_t *R, ec_curve_t *E) {
    ec_basis_t three_base, two_base;
    ec_point_t R3, R2;
    digit_t scalarP3[NWORDS_ORDER], scalarQ3[NWORDS_ORDER];
    digit_t scalarP2[NWORDS_ORDER], scalarQ2[NWORDS_ORDER];
    ibz_t iP3, iQ3, iP2, iQ2, t1, t2, t3;

    ibz_init(&iP3);
    ibz_init(&iQ3);
    ibz_init(&iP2);
    ibz_init(&iQ2);
    ibz_init(&t1);
    ibz_init(&t2);
    ibz_init(&t3);

    copy_point(&R3, R);
    copy_point(&three_base.P, &base->P);
    copy_point(&three_base.Q, &base->Q);
    copy_point(&three_base.PmQ, &base->PmQ);
    for(int i = 0; i < POWER_OF_2; i++) {
        ec_dbl(&R3, E, &R3);
        ec_dbl(&three_base.P, E, &three_base.P);
        ec_dbl(&three_base.Q, E, &three_base.Q);
        ec_dbl(&three_base.PmQ, E, &three_base.PmQ);
    }

    ec_dlog_3(scalarP3, scalarQ3, &three_base, &R3, E);
    ec_point_t test_point;
    xDBLMUL(&test_point, &three_base.P, scalarP3, &three_base.Q, scalarQ3, &three_base.PmQ, E);
    assert(ec_is_equal(&test_point, &R3));

    copy_point(&R2, R);
    copy_point(&two_base.P, &base->P);
    copy_point(&two_base.Q, &base->Q);
    copy_point(&two_base.PmQ, &base->PmQ);
    for(int i = 0; i < POWER_OF_3; i++) {
        ec_mul_ibz(&R2, E, &ibz_const_three, &R2);
        ec_mul_ibz(&two_base.P, E, &ibz_const_three, &two_base.P);
        ec_mul_ibz(&two_base.Q, E, &ibz_const_three, &two_base.Q);
        ec_mul_ibz(&two_base.PmQ, E, &ibz_const_three, &two_base.PmQ);
    }

    ec_dlog_2(scalarP2, scalarQ2, &two_base, &R2, E);
    xDBLMUL(&test_point, &two_base.P, scalarP2, &two_base.Q, scalarQ2, &two_base.PmQ, E);
    assert(ec_is_equal(&test_point, &R2));

    // Chinese Remainder Theorem
    ibz_copy_digits(&iP3, scalarP3, NWORDS_ORDER);
    ibz_copy_digits(&iQ3, scalarQ3, NWORDS_ORDER);
    ibz_copy_digits(&iP2, scalarP2, NWORDS_ORDER);
    ibz_copy_digits(&iQ2, scalarQ2, NWORDS_ORDER);

    ibz_crt(&iP3, &iP3, &iP2, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);
    ibz_crt(&iQ3, &iQ3, &iQ2, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);
    ibz_to_digits(scalarP, &iP3);
    ibz_to_digits(scalarQ, &iQ3);

    xDBLMUL(&test_point, &base->P, scalarP, &base->Q, scalarQ, &base->PmQ, E);
    if(!ec_is_equal(&test_point, R)){
        ibz_copy_digits(&iP3, scalarP3, NWORDS_ORDER);
        ibz_copy_digits(&iQ3, scalarQ3, NWORDS_ORDER);
        ibz_mod(&iP3, &iP3, &TORSION_PLUS_3POWER);
        ibz_sub(&iP3, &TORSION_PLUS_3POWER, &iP3);
        ibz_mod(&iQ3, &iQ3, &TORSION_PLUS_3POWER);
        ibz_sub(&iQ3, &TORSION_PLUS_3POWER, &iQ3);
        ibz_crt(&iP3, &iP3, &iP2, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);
        ibz_crt(&iQ3, &iQ3, &iQ2, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);
        ibz_to_digits(scalarP, &iP3);
        ibz_to_digits(scalarQ, &iQ3);
    }
    // Test if the computed scalars are correct
    xDBLMUL(&test_point, &base->P, scalarP, &base->Q, scalarQ, &base->PmQ, E);
    assert(ec_is_equal(&test_point, R));

    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&t3);
    ibz_finalize(&iP3);
    ibz_finalize(&iQ3);
    ibz_finalize(&iP2);
    ibz_finalize(&iQ2);

    return 1;
}

int eval_dimtwo_isog(theta_chain_t *phi, ec_basis_t *evalPQ, ec_basis_t *PQ, theta_couple_curve_t *E01) {
    theta_couple_point_t output_points, input_points;
    ec_point_t imP, imQ, imPQ, imR, imS, imRS;

    input_points.P1 = PQ->P;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&imP, &output_points.P1);

    input_points.P1 = PQ->Q;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&imQ, &output_points.P1);

    input_points.P1 = PQ->PmQ;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&imPQ, &output_points.P1);

    ec_basis_t RS;
    ec_curve_to_basis_6(&RS, &E01->E2);

    input_points.P2 = RS.P;
    ec_set_zero(&input_points.P1);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&imR, &output_points.P1);

    input_points.P2 = RS.Q;
    ec_set_zero(&input_points.P1);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&imS, &output_points.P1);

    input_points.P2 = RS.PmQ;
    ec_set_zero(&input_points.P1);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&imRS, &output_points.P1);

    ec_basis_t imRS_basis;
    imRS_basis.P = imR;
    imRS_basis.Q = imS;
    imRS_basis.PmQ = imRS;
    digit_t x[NWORDS_FIELD], y[NWORDS_FIELD];
    jac_point_t jacR, jacS;
    jac_point_t evalP, evalQ, evalPmQ;
    lift_basis(&jacR, &jacS, &RS, &E01->E2);

    ec_dlog_6(x, y, &imRS_basis, &imP, &phi->codomain.E1);
    ec_point_t test_point;
    xDBLMUL(&test_point, &imRS_basis.P, x, &imRS_basis.Q, y, &imRS_basis.PmQ, &phi->codomain.E1);
    if (!ec_is_equal(&test_point, &imP)) {
        printf("x*R + y*S != imP\n");
        return 0;
    }

    DBLMUL_generic(&evalP, &jacR, x, &jacS, y, &E01->E2, NWORDS_FIELD);
    ec_dlog_6(x, y, &imRS_basis, &imQ, &phi->codomain.E1);
    xDBLMUL(&test_point, &imRS_basis.P, x, &imRS_basis.Q, y, &imRS_basis.PmQ, &phi->codomain.E1);
    if (!ec_is_equal(&test_point, &imQ)) {
        printf("x*R + y*S != imQ\n");
        return 0;
    }

    DBLMUL_generic(&evalQ, &jacR, x, &jacS, y, &E01->E2, NWORDS_FIELD);
    ADD(&evalPmQ, &evalP, &evalQ, &E01->E2);
    
    jac_to_xz(&evalPQ->P, &evalP);
    jac_to_xz(&evalPQ->Q, &evalQ);
    jac_to_xz(&evalPQ->PmQ, &evalPmQ);

    return 1;
}

int test() {
    ibz_t q, alpha, beta, gamma, delta, rhs, deg;
    ibz_t A, q_bound;
    
    ibz_init(&q); ibz_init(&alpha); ibz_init(&beta);
    ibz_init(&gamma); ibz_init(&delta); ibz_init(&rhs); ibz_init(&deg);
    ibz_init(&A); ibz_init(&q_bound);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_div_2exp(&q_bound, &TORSION_PLUS_2POWER, 4);
    // Set q to a random value in the range [0, TORSION_PLUS_2POWER)
    for(int i = 0; i < 1000; i++) {
        ibz_rand_interval(&q, &ibz_const_zero, &q_bound);
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
    // fp2_print("j_invariant",&j_inv);

    quat_alg_elem_t tau;
    quat_alg_elem_init(&tau);
    if(represent_integer(&tau, &rhs, &QUATALG_PINFTY) == 0) {
        printf("Failed to represent integer in non-diagonal form\n");
        return 1;
    }

    quat_alg_elem_print(&tau);
    // check if the norm of tau equals to rhs
    ibq_t tau_norm;
    ibq_init(&tau_norm);
    quat_alg_norm(&tau_norm, &tau, &QUATALG_PINFTY);
    // gmp_printf("nrd(tau) : %Qd\n", tau_norm);
    // gmp_printf("rhs : %Zd\n", rhs);

    ec_basis_t E0_two, E0_three;
    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);

    endomorphism_application_three_basis(&E0_three, &curve, &tau, TORSION_PLUS_ODD_POWERS[0]);
    endomorphism_application_even_basis(&E0_two, &curve, &tau, TORSION_PLUS_EVEN_POWER);

    // point_print("tau(P2): ", E0_two.P);
    // point_print("tau(Q2): ", E0_two.Q);
    // point_print("tau(P3): ", E0_three.P);
    // point_print("tau(Q3): ", E0_three.Q);
    // ec_point_t tmpP, tmpQ;
    // ec_mul_ibz(&tmpP, &curve, &ibz_const_five, &BASIS_THREE.P);
    // ec_mul_ibz(&tmpQ, &curve, &ibz_const_five, &BASIS_THREE.Q);
    // point_print("5*P3: ", tmpP);
    // point_print("5*Q3: ", tmpQ);

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
        ec_mul_ibz(&test_point, &curve, &three_m1_order, &E0_three.Q);
        if (ec_is_zero(&test_point)) {
            printf("invalid 3^b-torsion points\n");
            return 1;
        }
        copy_point(&isog.ker_plus, &E0_three.Q);
    } else {
        ec_mul_ibz(&test_point, &curve, &three_m1_order, &E0_three.P);
        if (ec_is_zero(&test_point)) {
            printf("invalid 3^b-torsion points\n");
            return 1;
        }
        copy_point(&isog.ker_plus, &E0_three.P);
    }

    ec_set_zero(&isog.ker_minus);

    ec_mul_ibz(&test_point, &curve, &three_m1_order, &isog.ker_plus);
    if (ec_is_zero(&test_point)){
        printf("the order of the kernel point is not 3^m\n");
        return 1;
    }
    
    isog.degree[0] = TORSION_PLUS_ODD_POWERS[0];
    isog.degree[1] = 0;
    ec_curve_t E1;

    ec_eval_odd_basis(&E1, &isog, &E0_two, 1);

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
    ec_mul_ibz(&T2.P1, &curve, &q, &T2.P1);
    ec_mul_ibz(&T1m2.P1, &curve, &q, &T1m2.P1);

    fp2_t w0, w1, w0tw1;
    ec_point_t AC, A24;
    copy_point(&AC, &CURVE_E0_A24);
    A24_from_AC(&A24, &AC);
    weil(&w0, TORSION_PLUS_EVEN_POWER, &T1.P1, &T2.P1, &T1m2.P1, &A24);

    ibz_t inverse;
    ibz_init(&inverse);
    ibz_invmod(&inverse, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);

    ec_mul_ibz(&E0_two.P, &E1, &inverse, &E0_two.P);
    ec_mul_ibz(&E0_two.Q, &E1, &inverse, &E0_two.Q);
    ec_mul_ibz(&E0_two.PmQ, &E1, &inverse, &E0_two.PmQ);

    T1.P2 = E0_two.P;
    T2.P2 = E0_two.Q;
    T1m2.P2 = E0_two.PmQ;

    ec_curve_normalize_A24(&E1);
    copy_point(&A24, &E1.A24);
    weil(&w1, TORSION_PLUS_EVEN_POWER, &T1.P2, &T2.P2, &T1m2.P2, &A24);

    fp2_mul(&w0tw1, &w0, &w1);
    fp2_mul(&w0tw1, &w0tw1, &w0tw1);
    fp2_mul(&w0tw1, &w0tw1, &w0tw1);

    assert(fp2_is_one(&w0tw1));

    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &E01, &T1, &T2, &T1m2, strategies[2], 1);

    // Point evaluation via 2-dim isogeny
    jac_point_t P, Q, R, S, X, Y;
    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);

    lift_basis(&P, &Q, &E0_two, &curve);
    lift_basis(&R, &S, &E0_three, &curve);

    // lift_basis(&X, &Y, &BASIS_FIVE, &curve);
    // P23x = P + R + X and Q23x = Q + S + Y
    jac_point_t P23x, Q23x, PmQ23x;
    ADD(&P23x, &P, &R, &curve);
    // ADD(&P23x, &P23x, &X, &curve);
    ADD(&Q23x, &Q, &S, &curve);
    // ADD(&Q23x, &Q23x, &Y, &curve);
    ADD(&PmQ23x, &P23x, &Q23x, &curve);

    ec_basis_t eval_points, imPQ23x;
    jac_to_xz(&eval_points.P, &P23x);
    jac_to_xz(&eval_points.Q, &Q23x);
    jac_to_xz(&eval_points.PmQ, &PmQ23x);
    
    eval_dimtwo_isog(&hd_isog, &imPQ23x, &eval_points, &E01);

    printf("Succeed\n");

    ibz_finalize(&inverse);
    ibz_finalize(&three_m1_order);
    ibz_finalize(&remainder);
    ibz_finalize(&q);
    return 0; 
}

int main() {
    int res = 1;

    test();

    return res;
}