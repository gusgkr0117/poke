#include "poke.h"
#include <hd.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <klpt.h>
#include <quaternion.h>
#include <intbig.h>
#include <ec.h>
#include <rng.h>
#include <id2iso.h>
#include <biextension.h>
#include <fips202.h>

int ibz_random_matrix(ibz_mat_2x2_t mat22, const ibz_t *modulus) {
    ibz_t gcd, det;
    ibz_init(&gcd);
    ibz_init(&det);
    while (1) {
        ibz_rand_interval(&mat22[0][0], &ibz_const_zero, modulus);
        ibz_rand_interval(&mat22[0][1], &ibz_const_zero, modulus);
        ibz_rand_interval(&mat22[1][0], &ibz_const_zero, modulus);
        ibz_rand_interval(&mat22[1][1], &ibz_const_zero, modulus);
        ibz_mul(&det, &mat22[0][0], &mat22[1][1]);
        ibz_mul(&gcd, &mat22[0][1], &mat22[1][0]);
        ibz_sub(&det, &det, &gcd);
        ibz_gcd(&gcd, &det, modulus);
        if (ibz_is_one(&gcd)) {
            break;
        }
    }
    ibz_finalize(&gcd);
    ibz_finalize(&det);
    return 1;
}

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
    ibz_finalize(&gcd);
    return 1;
}

int eval_dimtwo_isog(theta_chain_t *phi, ibz_t *q, ec_basis_t *evalPQ, ec_basis_t *PQ, theta_couple_curve_t *E01) {
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
    ec_curve_to_basis_35(&RS, &E01->E2);
    // point_print("RS.P : ", RS.P);
    // point_print("RS.Q : ", RS.Q);
    // point_print("RS.PmQ : ", RS.PmQ);
    // curve_print("E01->E2 : ", E01->E2);

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
    digit_t x[NWORDS_FIELD] = {0}, y[NWORDS_FIELD] = {0};
    jac_point_t jacR, jacS;
    jac_point_t evalP, evalQ, evalPmQ, pointT;
    ibz_t t1,t2,t3,t4;
    ibz_init(&t1);
    ibz_init(&t2);
    ibz_init(&t3);
    ibz_init(&t4);

    lift_basis(&jacR, &jacS, &RS, &E01->E2);

    ec_dlog_35(x, y, &imRS_basis, &imP, &phi->codomain.E1);
    ec_point_t test_point;
    xDBLMUL(&test_point, &imRS_basis.P, x, &imRS_basis.Q, y, &imRS_basis.PmQ, &phi->codomain.E1);
    if (!ec_is_equal(&test_point, &imP)) {
        printf("x*R + y*S != imP\n");
        return 0;
    }
    ibz_copy_digits(&t1, x, NWORDS_ORDER);
    ibz_copy_digits(&t2, y, NWORDS_ORDER);

    DBLMUL_generic(&evalP, &jacR, x, &jacS, y, &E01->E2, NWORDS_ORDER);

    ec_dlog_35(x, y, &imRS_basis, &imQ, &phi->codomain.E1);
    xDBLMUL(&test_point, &imRS_basis.P, x, &imRS_basis.Q, y, &imRS_basis.PmQ, &phi->codomain.E1);
    if (!ec_is_equal(&test_point, &imQ)) {
        printf("x*R + y*S != imQ\n");
        return 0;
    }
    ibz_copy_digits(&t3, x, NWORDS_ORDER);
    ibz_copy_digits(&t4, y, NWORDS_ORDER);

    DBLMUL_generic(&evalQ, &jacR, x, &jacS, y, &E01->E2, NWORDS_ORDER);

    // Test if (t1 - t3) * imR + (t2 - t4) * imS = imP - imQ
    // If not, evalQ = -evalQ
    ibz_sub(&t1, &t1, &t3);
    ibz_sub(&t2, &t2, &t4);
    ibz_mod(&t1, &t1, &TORSION_PLUS_23CPOWER);
    ibz_mod(&t2, &t2, &TORSION_PLUS_23CPOWER);
    memset(x, 0, NWORDS_ORDER * RADIX / 8);
    memset(y, 0, NWORDS_ORDER * RADIX / 8);
    ibz_to_digits(x, &t1);
    ibz_to_digits(y, &t2);
    xDBLMUL(&test_point, &imRS_basis.P, x, &imRS_basis.Q, y, &imRS_basis.PmQ, &phi->codomain.E1);
    if (!ec_is_equal(&test_point, &imPQ)) {
        // test otherwise
        ibz_add(&t1, &t1, &t3);
        ibz_add(&t1, &t1, &t3);
        ibz_add(&t2, &t2, &t4);
        ibz_add(&t2, &t2, &t4);
        ibz_mod(&t1, &t1, &TORSION_PLUS_23CPOWER);
        ibz_mod(&t2, &t2, &TORSION_PLUS_23CPOWER);
        memset(x, 0, NWORDS_ORDER * RADIX / 8);
        memset(y, 0, NWORDS_ORDER * RADIX / 8);
        ibz_to_digits(x, &t1);
        ibz_to_digits(y, &t2);
        xDBLMUL(&test_point, &imRS_basis.P, x, &imRS_basis.Q, y, &imRS_basis.PmQ, &phi->codomain.E1);
        assert(ec_is_equal(&test_point, &imPQ));
        jac_neg(&evalQ, &evalQ);
    }


    jac_neg(&pointT, &evalQ);
    ADD(&evalPmQ, &evalP, &pointT, &E01->E2);
    
    jac_to_xz(&evalPQ->P, &evalP);
    jac_to_xz(&evalPQ->Q, &evalQ);
    jac_to_xz(&evalPQ->PmQ, &evalPmQ);
    // ec_mul_ibz(&evalPQ->P, &E01->E2, q, &evalPQ->P);
    // ec_mul_ibz(&evalPQ->Q, &E01->E2, q, &evalPQ->Q);
    // ec_mul_ibz(&evalPQ->PmQ, &E01->E2, q, &evalPQ->PmQ);

    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&t3);
    ibz_finalize(&t4);
    return 1;
}

int keygen(poke_sk_t *sk, poke_pk_t *pk) {
    ibz_t q, alpha, beta, gamma, delta, rhs, deg;
    ibz_t A, q_bound;
    ec_point_t pointT;
    
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
    ibz_random_unit(&delta, &TORSION_PLUS_CPOWER);

    memset(&sk->deg, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->alpha, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->beta, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->delta, 0, NWORDS_ORDER * RADIX / 8);

    ibz_to_digits(sk->deg, &q);
    ibz_to_digits(sk->alpha, &alpha);
    ibz_to_digits(sk->beta, &beta);
    ibz_to_digits(sk->delta, &delta);

    // gmp_printf("q = %Zd\n", q);

    // gmp_printf("alpha = %Zd\n", alpha);
    // gmp_printf("beta = %Zd\n", beta);
    // gmp_printf("gamma = %Zd\n", gamma);
    // gmp_printf("delta = %Zd\n", delta);
    ec_curve_t curve = CURVE_E0;
    // fp2_t j_inv;
    // ec_j_inv(&j_inv, &curve);

    quat_alg_elem_t tau;
    quat_alg_elem_init(&tau);
    if(represent_integer(&tau, &rhs, &QUATALG_PINFTY) == 0) {
        printf("Failed to represent integer in non-diagonal form\n");
        return 1;
    }

    // quat_alg_elem_print(&tau);
    // check if the norm of tau equals to rhs
    ibq_t tau_norm;
    ibq_init(&tau_norm);
    quat_alg_norm(&tau_norm, &tau, &QUATALG_PINFTY);

    ec_basis_t E0_two, E0_three, E0_five;
    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);

    endomorphism_application_three_basis(&E0_three, &curve, &tau, TORSION_PLUS_ODD_POWERS[0]);
    endomorphism_application_even_basis(&E0_two, &curve, &tau, TORSION_PLUS_EVEN_POWER);

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
    isog.degree[2] = 0;
    ec_curve_t E1;

    ec_eval_odd_basis(&E1, &isog, &E0_two, 1);

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

    // fp2_t w0, w1, w0tw1;
    // ec_point_t AC, A24;
    // copy_point(&AC, &CURVE_E0_A24);
    // A24_from_AC(&A24, &AC);
    // weil(&w0, TORSION_PLUS_EVEN_POWER, &T1.P1, &T2.P1, &T1m2.P1, &A24);

    ibz_t inverse;
    ibz_init(&inverse);
    ibz_invmod(&inverse, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);

    ec_mul_ibz(&E0_two.P, &E1, &inverse, &E0_two.P);
    ec_mul_ibz(&E0_two.Q, &E1, &inverse, &E0_two.Q);
    ec_mul_ibz(&E0_two.PmQ, &E1, &inverse, &E0_two.PmQ);

    T1.P2 = E0_two.P;
    T2.P2 = E0_two.Q;
    T1m2.P2 = E0_two.PmQ;

    // ec_curve_normalize_A24(&E1);
    // copy_point(&A24, &E1.A24);
    // weil(&w1, TORSION_PLUS_EVEN_POWER, &T1.P2, &T2.P2, &T1m2.P2, &A24);

    // fp2_mul(&w0tw1, &w0, &w1);
    // fp2_mul(&w0tw1, &w0tw1, &w0tw1);
    // fp2_mul(&w0tw1, &w0tw1, &w0tw1);

    // assert(fp2_is_one(&w0tw1));

    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &E01, &T1, &T2, &T1m2, strategies[2], 1);

    // Point evaluation via 2-dim isogeny
    jac_point_t P, Q, R, S, X, Y, T;
    ec_basis_t tmp_basis;
    copy_point(&tmp_basis.P, &E0_two.P);
    copy_point(&tmp_basis.Q, &E0_two.Q);
    copy_point(&tmp_basis.PmQ, &E0_two.PmQ);

    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);
    copy_point(&E0_five.P, &BASIS_C.P);
    copy_point(&E0_five.Q, &BASIS_C.Q);
    copy_point(&E0_five.PmQ, &BASIS_C.PmQ);

    // lift_basis(&P, &Q, &E0_two, &curve);
    lift_basis(&R, &S, &E0_three, &curve);
    lift_basis(&X, &Y, &E0_five, &curve);

    // P23x = R + X and Q23x = S + Y
    jac_point_t P23x, Q23x, PmQ23x;
    ADD(&P23x, &R, &X, &curve);
    ADD(&Q23x, &S, &Y, &curve);
    jac_neg(&T, &Q23x);
    ADD(&PmQ23x, &P23x, &T, &curve);

    ec_basis_t eval_points, imPQ23x;
    jac_to_xz(&eval_points.P, &P23x);
    jac_to_xz(&eval_points.Q, &Q23x);
    jac_to_xz(&eval_points.PmQ, &PmQ23x);
    
    eval_dimtwo_isog(&hd_isog, &deg, &imPQ23x, &eval_points, &E01);

    ibz_t cofactor;
    // Compute X_A, Y_B
    ibz_init(&cofactor);
    ibz_invmod(&cofactor, &TORSION_PLUS_23POWER, &TORSION_PLUS_CPOWER);
    ibz_mul(&cofactor, &cofactor, &TORSION_PLUS_23POWER);
    ec_mul_ibz(&pk->PQxy.P, &E01.E2, &cofactor, &imPQ23x.P);
    ec_mul_ibz(&pk->PQxy.P, &E01.E2, &delta, &pk->PQxy.P);
    ec_mul_ibz(&pk->PQxy.Q, &E01.E2, &cofactor, &imPQ23x.Q);
    ec_mul_ibz(&pk->PQxy.Q, &E01.E2, &delta, &pk->PQxy.Q);
    ec_mul_ibz(&pk->PQxy.PmQ, &E01.E2, &cofactor, &imPQ23x.PmQ);
    ec_mul_ibz(&pk->PQxy.PmQ, &E01.E2, &delta, &pk->PQxy.PmQ);

    // Compute P2, Q2
    ec_mul(&pk->PQ2.P, &E01.E2, sk->alpha, &tmp_basis.P);
    ec_mul(&pk->PQ2.Q, &E01.E2, sk->beta, &tmp_basis.Q);
    xADD(&pointT, &tmp_basis.P, &tmp_basis.Q, &tmp_basis.PmQ);
    xDBLMUL(&pk->PQ2.PmQ, &tmp_basis.P, sk->alpha, &tmp_basis.Q, sk->beta, &pointT, &E01.E2);  

    // Compute P3, Q3
    ibz_invmod(&cofactor, &TORSION_PLUS_2CPOWER, &TORSION_PLUS_3POWER);
    ibz_mul(&cofactor, &cofactor, &TORSION_PLUS_2CPOWER);
    ec_mul_ibz(&pk->PQ3.P, &E01.E2, &cofactor, &imPQ23x.P);
    ec_mul_ibz(&pk->PQ3.P, &E01.E2, &gamma, &pk->PQ3.P);
    ec_mul_ibz(&pk->PQ3.Q, &E01.E2, &cofactor, &imPQ23x.Q);
    ec_mul_ibz(&pk->PQ3.Q, &E01.E2, &gamma, &pk->PQ3.Q);
    ec_mul_ibz(&pk->PQ3.PmQ, &E01.E2, &cofactor, &imPQ23x.PmQ);
    ec_mul_ibz(&pk->PQ3.PmQ, &E01.E2, &gamma, &pk->PQ3.PmQ);

    pk->EA = E01.E2;

    ibz_finalize(&cofactor);
    ibz_finalize(&remainder);
    ibz_finalize(&inverse);
    ibz_finalize(&q);
    ibz_finalize(&alpha);
    ibz_finalize(&beta);
    ibz_finalize(&delta);
    ibz_finalize(&gamma);
    ibz_finalize(&rhs);
    ibz_finalize(&deg);
    ibz_finalize(&A);
    ibz_finalize(&q_bound);
    return 1; 
}

int encrypt(poke_ct_t *ct, const poke_pk_t *pk, const unsigned char *m, const size_t m_len) {
    ibz_mat_2x2_t mask_xy;
    ec_isog_odd_t isogB, isogB_prime;
    ibz_t beta, omega, omega_inv, TT, A;
    ec_curve_t EB, EAB;
    ec_basis_t E0_two, E0_xy, EA_two, EA_xy, EAB_xy, eval_basis[2];
    ec_point_t pointT;

    digit_t beta_scalar[NWORDS_ORDER] = {0}, omega_scalar[NWORDS_ORDER] = {0}, omega_inv_scalar[NWORDS_ORDER] = {0}, one_scalar[NWORDS_ORDER] = {1};
    digit_t mask_xy_scalar[6][NWORDS_ORDER] = {0};

    ibz_init(&beta);
    ibz_init(&omega);
    ibz_init(&omega_inv);
    ibz_init(&TT);
    ibz_init(&A);
    ibz_mat_2x2_init(&mask_xy);

    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_random_unit(&beta, &TORSION_PLUS_3POWER);
    ibz_random_unit(&omega, &A);
    ibz_invmod(&omega_inv, &omega, &A);
    ibz_to_digits(omega_scalar, &omega);
    ibz_to_digits(omega_inv_scalar, &omega_inv);
    ibz_to_digits(beta_scalar, &beta);
    ibz_random_matrix(mask_xy, &TORSION_PLUS_CPOWER);
    for(int i = 0; i < 4; i++) {
        ibz_to_digits(mask_xy_scalar[i], &mask_xy[i / 2][i % 2]);
    }
    ibz_sub(&TT, &mask_xy[0][0], &mask_xy[1][0]);
    ibz_mod(&TT, &TT, &TORSION_PLUS_CPOWER);
    ibz_to_digits(mask_xy_scalar[4], &TT);
    ibz_sub(&TT, &mask_xy[0][1], &mask_xy[1][1]);
    ibz_mod(&TT, &TT, &TORSION_PLUS_CPOWER);
    ibz_to_digits(mask_xy_scalar[5], &TT);

    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_xy.P, &BASIS_C.P);
    copy_point(&E0_xy.Q, &BASIS_C.Q);
    copy_point(&E0_xy.PmQ, &BASIS_C.PmQ);
    copy_point(&EA_two.P, &pk->PQ2.P);
    copy_point(&EA_two.Q, &pk->PQ2.Q);
    copy_point(&EA_two.PmQ, &pk->PQ2.PmQ);
    copy_point(&EA_xy.P, &pk->PQxy.P);
    copy_point(&EA_xy.Q, &pk->PQxy.Q);
    copy_point(&EA_xy.PmQ, &pk->PQxy.PmQ);


    // Compute the isogeny E0 -> EB
    isogB.curve = CURVE_E0;
    isogB.degree[0] = POWER_OF_3;
    isogB.degree[1] = 0;
    isogB.degree[2] = 0;
    ec_set_zero(&isogB.ker_minus);
    // kernel = P + beta * Q
    xDBLMUL(&isogB.ker_plus, &BASIS_THREE.P, one_scalar, &BASIS_THREE.Q, beta_scalar, &BASIS_THREE.PmQ, &isogB.curve);
    
    eval_basis[0] = E0_two;
    eval_basis[1] = E0_xy;
    ec_eval_odd_basis(&EB, &isogB, eval_basis, 2);
    E0_two = eval_basis[0];
    E0_xy = eval_basis[1];

    ct->EB = EB;

    // Masking evaluated basis points
    xMUL(&ct->PQ2_B.P, &E0_two.P, omega_scalar, &EB);
    xMUL(&ct->PQ2_B.Q, &E0_two.Q, omega_inv_scalar, &EB);
    xADD(&pointT, &E0_two.P, &E0_two.Q, &E0_two.PmQ);
    xDBLMUL(&ct->PQ2_B.PmQ, &E0_two.P, omega_scalar, &E0_two.Q, omega_inv_scalar, &pointT, &EB);

    xDBLMUL(&ct->PQxy_B.P, &E0_xy.P, mask_xy_scalar[0], &E0_xy.Q, mask_xy_scalar[1], &E0_xy.PmQ, &EB);
    xDBLMUL(&ct->PQxy_B.Q, &E0_xy.P, mask_xy_scalar[2], &E0_xy.Q, mask_xy_scalar[3], &E0_xy.PmQ, &EB);
    xDBLMUL(&ct->PQxy_B.PmQ, &E0_xy.P, mask_xy_scalar[4], &E0_xy.Q, mask_xy_scalar[5], &E0_xy.PmQ, &EB);

    // Compute the isogeny EA -> EAB
    isogB_prime.curve = pk->EA;
    isogB_prime.degree[0] = TORSION_PLUS_ODD_POWERS[0];
    isogB_prime.degree[1] = 0;
    isogB_prime.degree[2] = 0;
    ec_set_zero(&isogB_prime.ker_minus);
    // kernel = P + beta * Q
    xDBLMUL(&isogB_prime.ker_plus, &pk->PQ3.P, one_scalar, &pk->PQ3.Q, beta_scalar, &pk->PQ3.PmQ, &isogB_prime.curve);

    eval_basis[0] = EA_two;
    eval_basis[1] = EA_xy;
    ec_eval_odd_basis(&EAB, &isogB_prime, eval_basis, 2);
    EA_two = eval_basis[0];
    EA_xy = eval_basis[1];

    ct->EAB = EAB;

    // Masking evaluated basis points
    xMUL(&ct->PQ2_AB.P, &EA_two.P, omega_scalar, &EAB);
    xMUL(&ct->PQ2_AB.Q, &EA_two.Q, omega_inv_scalar, &EAB);
    xADD(&pointT, &EA_two.P, &EA_two.Q, &EA_two.PmQ);
    xDBLMUL(&ct->PQ2_AB.PmQ, &EA_two.P, omega_scalar, &EA_two.Q, omega_inv_scalar, &pointT, &EAB);

    xDBLMUL(&EAB_xy.P, &EA_xy.P, mask_xy_scalar[0], &EA_xy.Q, mask_xy_scalar[1], &EA_xy.PmQ, &EAB);
    xDBLMUL(&EAB_xy.Q, &EA_xy.P, mask_xy_scalar[2], &EA_xy.Q, mask_xy_scalar[3], &EA_xy.PmQ, &EAB);
    xDBLMUL(&EAB_xy.PmQ, &EA_xy.P, mask_xy_scalar[4], &EA_xy.Q, mask_xy_scalar[5], &EA_xy.PmQ, &EAB);

    // TODO : ct <- SHA256(EAB_xy.P || EAB_xy.Q) xor m
    unsigned char hash_input[4 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};

    ec_normalize_point(&EAB_xy.P);
    ec_normalize_point(&EAB_xy.Q);

    memcpy(hash_input, &EAB_xy.P.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &EAB_xy.P.x.im[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + 2 * NWORDS_FIELD * RADIX / 8, &EAB_xy.Q.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + 3 * NWORDS_FIELD * RADIX / 8, &EAB_xy.Q.x.im[0], NWORDS_FIELD * RADIX / 8);
    

    SHAKE256(hash_output, sizeof(hash_output), hash_input, sizeof(hash_input));

    // ct->ct = m xor hash_output
    memset(ct->ct, 0, sizeof(ct->ct));
    for (size_t i = 0; i < 32; i++) {
        if (i >= m_len) {
            ct->ct[i] = hash_output[i];
        } else {
            ct->ct[i] = m[i] ^ hash_output[i];
        }
    }

    ibz_mat_2x2_finalize(&mask_xy);
    ibz_finalize(&TT);
    ibz_finalize(&A);
    ibz_finalize(&beta);
    ibz_finalize(&omega);
    ibz_finalize(&omega_inv);
    return 1;
}

int decrypt(unsigned char *m, size_t *m_len, const poke_ct_t *ct, const poke_sk_t *sk) {
    unsigned char hash_input[4 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};
    digit_t T1_scalar[NWORDS_ORDER] = {0}, T2_scalar[NWORDS_ORDER] = {0};
    theta_chain_t hd_isog;
    theta_couple_curve_t EBAB;
    theta_couple_point_t T1, T2, T1m2;
    ec_basis_t eval_points;
    ibz_t alpha_inv, beta_inv, deg, A;
    ec_point_t pointT;

    ibz_init(&alpha_inv);
    ibz_init(&beta_inv);
    ibz_init(&deg);
    ibz_init(&A);

    ibz_copy_digits(&deg, sk->deg, NWORDS_ORDER);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_copy_digits(&alpha_inv, sk->alpha, NWORDS_ORDER);
    ibz_copy_digits(&beta_inv, sk->beta, NWORDS_ORDER);
    ibz_invmod(&alpha_inv, &alpha_inv, &TORSION_PLUS_2POWER);
    ibz_invmod(&beta_inv, &beta_inv, &TORSION_PLUS_2POWER);

    EBAB.E1 = ct->EB;
    EBAB.E2 = ct->EAB;

    ec_mul_ibz(&T1.P1, &EBAB.E1, &deg, &ct->PQ2_B.P);
    ec_mul_ibz(&T2.P1, &EBAB.E1, &deg, &ct->PQ2_B.Q);
    ec_mul_ibz(&T1m2.P1, &EBAB.E1, &deg, &ct->PQ2_B.PmQ);
    ibz_to_digits(T1_scalar, &alpha_inv);
    ibz_to_digits(T2_scalar, &beta_inv);
    xMUL(&T1.P2, &ct->PQ2_AB.P, T1_scalar, &EBAB.E2);
    xMUL(&T2.P2, &ct->PQ2_AB.Q, T2_scalar, &EBAB.E2);
    xADD(&pointT, &ct->PQ2_AB.P, &ct->PQ2_AB.Q, &ct->PQ2_AB.PmQ);
    xDBLMUL(&T1m2.P2, &ct->PQ2_AB.P, T1_scalar, &ct->PQ2_AB.Q, T2_scalar, &pointT, &EBAB.E2);

    fp2_t w0, w1, w0tw1;
    ec_point_t A24;
    ec_curve_normalize_A24(&EBAB.E1);
    copy_point(&A24, &EBAB.E1.A24);
    weil(&w0, TORSION_PLUS_EVEN_POWER, &T1.P1, &T2.P1, &T1m2.P1, &A24);
    
    ec_curve_normalize_A24(&EBAB.E2);
    copy_point(&A24, &EBAB.E2.A24);
    weil(&w1, TORSION_PLUS_EVEN_POWER, &T1.P2, &T2.P2, &T1m2.P2, &A24);

    fp2_mul(&w0tw1, &w0, &w1);
    fp2_mul(&w0tw1, &w0tw1, &w0tw1);
    fp2_mul(&w0tw1, &w0tw1, &w0tw1);

    if(!fp2_is_one(&w0tw1)) {
        printf("The Weil pairing is not equal to 1\n");
        fp2_print("w0 : ", &w0);
        fp2_print("w1 : ", &w1);
        fp2_print("w0tw1 : ", &w0tw1);
        ibz_finalize(&alpha_inv);
        ibz_finalize(&beta_inv);
        ibz_finalize(&deg);
        return 0;
    }

    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &EBAB, &T1, &T2, &T1m2, strategies[2], 1);

    eval_points.P = ct->PQxy_B.P;
    eval_points.Q = ct->PQxy_B.Q;
    eval_points.PmQ = ct->PQxy_B.PmQ;
    ibz_sub(&deg, &A, &deg);
    if (!eval_dimtwo_isog(&hd_isog, &deg, &eval_points, &eval_points, &EBAB)) {
        printf("Failed to evaluate the 2-dim isogeny\n");
        return 0;
    }

    xMUL(&eval_points.P, &eval_points.P, sk->delta, &EBAB.E2);
    xMUL(&eval_points.Q, &eval_points.Q, sk->delta, &EBAB.E2);
    xMUL(&eval_points.PmQ, &eval_points.PmQ, sk->delta, &EBAB.E2);

    ec_normalize_point(&eval_points.P);
    ec_normalize_point(&eval_points.Q);

    memcpy(hash_input, &eval_points.P.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &eval_points.P.x.im[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + 2 * NWORDS_FIELD * RADIX / 8, &eval_points.Q.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + 3 * NWORDS_FIELD * RADIX / 8, &eval_points.Q.x.im[0], NWORDS_FIELD * RADIX / 8);

    SHAKE256(hash_output, sizeof(hash_output), hash_input, sizeof(hash_input));

    // m = ct->ct xor hash_output
    memset(m, 0, sizeof(hash_output));
    for (size_t i = 0; i < sizeof(hash_output); i++) {
        m[i] = ct->ct[i] ^ hash_output[i];
    }
    *m_len = sizeof(hash_output);

    ibz_finalize(&alpha_inv);
    ibz_finalize(&beta_inv);
    ibz_finalize(&deg);
    ibz_finalize(&A);
    return 1;
}