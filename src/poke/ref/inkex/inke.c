#include <inke.h>
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

int ibz_random_unit(ibz_t *q, const ibz_t *modulus, shake256ctx *state) {
    ibz_t gcd;
    ibz_init(&gcd);
    while (1) {
        if (state != NULL) {
            ibz_rand_interval_with_state(q, &ibz_const_zero, modulus, state);
        } else {
            ibz_rand_interval(q, &ibz_const_zero, modulus);
        }
        ibz_gcd(&gcd, q, modulus);
        if (ibz_is_one(&gcd)) {
            break;
        }
    }
    ibz_finalize(&gcd);
    return 1;
}

int eval_dimtwo_isog_with_middle(theta_chain_t *phi, ibz_t *q, ec_basis_t *evalPQmid, ec_basis_t *evalPQ, ec_basis_t *PQ, theta_couple_curve_t *E01) {
    theta_couple_point_t output_points, input_points;
    ec_point_t imP, imQ, imPQ, imR, imS, imRS;

    input_points.P1 = PQ->P;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&evalPQmid->P, &output_points.P1);

    input_points.P1 = PQ->Q;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&evalPQmid->Q, &output_points.P1);

    input_points.P1 = PQ->PmQ;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&evalPQmid->PmQ, &output_points.P1);

    ec_basis_t RS;
    ec_curve_to_basis_3(&RS, &E01->E2);

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

    ec_dlog_3(x, y, &imRS_basis, &evalPQmid->P, &phi->codomain.E1);
    ec_point_t test_point;
    
    ibz_copy_digits(&t1, x, NWORDS_ORDER);
    ibz_copy_digits(&t2, y, NWORDS_ORDER);

    DBLMUL_generic(&evalP, &jacR, x, &jacS, y, &E01->E2, NWORDS_ORDER);

    ec_dlog_3(x, y, &imRS_basis, &evalPQmid->Q, &phi->codomain.E1);
    
    ibz_copy_digits(&t3, x, NWORDS_ORDER);
    ibz_copy_digits(&t4, y, NWORDS_ORDER);

    DBLMUL_generic(&evalQ, &jacR, x, &jacS, y, &E01->E2, NWORDS_ORDER);

    // Test if (t1 - t3) * imR + (t2 - t4) * imS = imP - imQ
    // If not, evalQ = -evalQ
    ibz_sub(&t1, &t1, &t3);
    ibz_sub(&t2, &t2, &t4);
    ibz_mod(&t1, &t1, &TORSION_PLUS_3POWER);
    ibz_mod(&t2, &t2, &TORSION_PLUS_3POWER);
    memset(x, 0, NWORDS_ORDER * RADIX / 8);
    memset(y, 0, NWORDS_ORDER * RADIX / 8);
    ibz_to_digits(x, &t1);
    ibz_to_digits(y, &t2);
    xDBLMUL_bounded(&test_point, &imRS_basis.P, x, &imRS_basis.Q, y, &imRS_basis.PmQ, &phi->codomain.E1, TORSION_3POWER_BYTES * 8);
    if (!ec_is_equal(&test_point, &evalPQmid->PmQ)) {
        jac_neg(&evalQ, &evalQ);
    }


    jac_neg(&pointT, &evalQ);
    ADD(&evalPmQ, &evalP, &pointT, &E01->E2);
    
    jac_to_xz(&evalPQ->P, &evalP);
    jac_to_xz(&evalPQ->Q, &evalQ);
    jac_to_xz(&evalPQ->PmQ, &evalPmQ);

    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&t3);
    ibz_finalize(&t4);
    return 1;
}

int keygen(inke_sk_t *sk, inke_pk_t *pk) {
    ibz_t q, alpha, beta, gamma, gamma1, rhs, deg;
    ibz_t A, q_bound;
    ec_point_t pointT;
    
    ibz_init(&q); ibz_init(&alpha); ibz_init(&beta);
    ibz_init(&gamma); ibz_init(&gamma1); ibz_init(&rhs); ibz_init(&deg);
    ibz_init(&A); ibz_init(&q_bound);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_div_2exp(&q_bound, &TORSION_PLUS_2POWER, 4);
    // Set q to a random value in the range [0, TORSION_PLUS_2POWER)
    for(int i = 0; i < 1000; i++) {
        ibz_rand_interval(&q, &ibz_const_zero, &q_bound);
        ibz_sub(&deg, &A, &q);
        if (ibz_divides(&q, &ibz_const_two) == 0 && ibz_divides(&q, &ibz_const_three) == 0
            && ibz_divides(&deg, &ibz_const_three) == 0) {
            break;
        }
    }
    ibz_mul(&rhs, &deg, &q);
    ibz_mul(&rhs, &rhs, &TORSION_PLUS_3POWER);
    ibz_random_unit(&alpha, &A, NULL);
    ibz_random_unit(&beta, &A, NULL);
    ibz_random_unit(&gamma, &TORSION_PLUS_3POWER, NULL);
    ibz_random_unit(&gamma1, &TORSION_PLUS_3POWER, NULL);

    memset(&sk->deg, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->alpha, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->beta, 0, NWORDS_ORDER * RADIX / 8);

    ibz_to_digits(sk->deg, &q);
    ibz_to_digits(sk->alpha, &alpha);
    ibz_to_digits(sk->beta, &beta);

    ec_curve_t curve = CURVE_E0;
    quat_alg_elem_t tau;
    quat_alg_elem_init(&tau);
    if(represent_integer(&tau, &rhs, &QUATALG_PINFTY) == 0) {
        printf("Failed to represent integer in non-diagonal form\n");
        return 1;
    }

    ec_basis_t E0_two, E0_three, E0_five;
    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);
    
    endomorphism_application_three_basis(&E0_three, &curve, &tau, TORSION_PLUS_ODD_POWERS[0]);
    endomorphism_application_even_basis(&E0_two, &curve, &tau, TORSION_PLUS_EVEN_POWER);
    quat_alg_elem_finalize(&tau);

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
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isog.degree[i] = 0;
    }
    ec_curve_t E1;

    ec_eval_three(&E1, &isog, (ec_point_t*)&E0_two, 3);

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

    ibz_t inverse;
    ibz_init(&inverse);
    ibz_invmod(&inverse, &TORSION_PLUS_3POWER, &TORSION_PLUS_2POWER);

    ec_mul_ibz(&E0_two.P, &E1, &inverse, &E0_two.P);
    ec_mul_ibz(&E0_two.Q, &E1, &inverse, &E0_two.Q);
    ec_mul_ibz(&E0_two.PmQ, &E1, &inverse, &E0_two.PmQ);

    T1.P2 = E0_two.P;
    T2.P2 = E0_two.Q;
    T1m2.P2 = E0_two.PmQ;

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

    ec_basis_t eval_points, imPQ3, imPQxy, imPQ31;

    copy_point(&eval_points.P, &E0_three.P);
    copy_point(&eval_points.Q, &E0_three.Q);
    copy_point(&eval_points.PmQ, &E0_three.PmQ);
    
    eval_dimtwo_isog_with_middle(&hd_isog, &deg, &imPQ31, &imPQ3, &eval_points, &E01);

    // Compute P2, Q2
    ec_mul(&pk->PQ2.P, &E01.E2, sk->alpha, &tmp_basis.P);
    ec_mul(&pk->PQ2.Q, &E01.E2, sk->beta, &tmp_basis.Q);
    xADD(&pointT, &tmp_basis.P, &tmp_basis.Q, &tmp_basis.PmQ);
    xDBLMUL_bounded(&pk->PQ2.PmQ, &tmp_basis.P, sk->alpha, &tmp_basis.Q, sk->beta, &pointT, &E01.E2, TORSION_2POWER_BYTES * 8);  

    // Compute P3, Q3
    ec_mul_ibz(&pk->PQ3.P, &E01.E2, &gamma, &imPQ3.P);
    ec_mul_ibz(&pk->PQ3.Q, &E01.E2, &gamma, &imPQ3.Q);
    ec_mul_ibz(&pk->PQ3.PmQ, &E01.E2, &gamma, &imPQ3.PmQ);

    ec_mul_ibz(&pk->PQA13.P, &hd_isog.codomain.E1, &gamma1, &imPQ31.P);
    ec_mul_ibz(&pk->PQA13.Q, &hd_isog.codomain.E1, &gamma1, &imPQ31.Q);
    ec_mul_ibz(&pk->PQA13.PmQ, &hd_isog.codomain.E1, &gamma1, &imPQ31.PmQ);

    pk->EA1 = hd_isog.codomain.E1;
    pk->EA = E01.E2;

    ibz_finalize(&three_m1_order);
    ibz_finalize(&remainder);
    ibz_finalize(&inverse);
    ibz_finalize(&q);
    ibz_finalize(&alpha);
    ibz_finalize(&beta);
    ibz_finalize(&gamma);
    ibz_finalize(&gamma1);
    ibz_finalize(&rhs);
    ibz_finalize(&deg);
    ibz_finalize(&A);
    ibz_finalize(&q_bound);
    return 1; 
}

int encrypt(inke_ct_t *ct, const inke_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len) {
    ibz_mat_2x2_t mask_xy;
    ec_isog_odd_t isogB, isogB_prime1, isogB_prime;
    ibz_t beta, omega, omega_inv, TT, A;
    ec_curve_t EB, EA1B, EAB;
    ec_basis_t E0_two, E0_xy, EA_two, EA1_xy, EA1B_xy, eval_basis[2];
    ec_point_t pointT;
    shake256ctx state;

    digit_t beta_scalar[NWORDS_ORDER] = {0}, omega_scalar[NWORDS_ORDER] = {0}, omega_inv_scalar[NWORDS_ORDER] = {0}, one_scalar[NWORDS_ORDER] = {1};
    digit_t mask_xy_scalar[6][NWORDS_ORDER] = {0};

    ibz_init(&beta);
    ibz_init(&omega);
    ibz_init(&omega_inv);
    ibz_init(&TT);
    ibz_init(&A);
    ibz_mat_2x2_init(&mask_xy);

    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    if (seed != NULL) {
        shake256_absorb(&state, seed, seed_len);
        ibz_random_unit(&beta, &TORSION_PLUS_3POWER, &state);
        ibz_random_unit(&omega, &A, &state);
        shake256_ctx_release(&state);
    }
    else {
        ibz_random_unit(&beta, &TORSION_PLUS_3POWER, NULL);
        ibz_random_unit(&omega, &A, NULL);
    }
    ibz_invmod(&omega_inv, &omega, &A);
    ibz_to_digits(omega_scalar, &omega);
    ibz_to_digits(omega_inv_scalar, &omega_inv);
    ibz_to_digits(beta_scalar, &beta);

    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&EA_two.P, &pk->PQ2.P);
    copy_point(&EA_two.Q, &pk->PQ2.Q);
    copy_point(&EA_two.PmQ, &pk->PQ2.PmQ);

    // Compute the isogeny E0 -> EB
    isogB.curve = CURVE_E0;
    isogB.degree[0] = POWER_OF_3;
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isogB.degree[i] = 0;
    }
    ec_set_zero(&isogB.ker_minus);
    // kernel = P + beta * Q
    xDBLMUL_bounded(&isogB.ker_plus, &BASIS_THREE.P, one_scalar, &BASIS_THREE.Q, beta_scalar, &BASIS_THREE.PmQ, &isogB.curve, TORSION_3POWER_BYTES * 8);
    
    eval_basis[0] = E0_two;
    ec_eval_three(&EB, &isogB, (ec_point_t*)eval_basis, 3);
    E0_two = eval_basis[0];

    ct->EB = EB;

    // Masking evaluated basis points
    xMUL(&ct->PQ2_B.P, &E0_two.P, omega_scalar, &EB);
    xMUL(&ct->PQ2_B.Q, &E0_two.Q, omega_inv_scalar, &EB);
    xADD(&pointT, &E0_two.P, &E0_two.Q, &E0_two.PmQ);
    xDBLMUL_bounded(&ct->PQ2_B.PmQ, &E0_two.P, omega_scalar, &E0_two.Q, omega_inv_scalar, &pointT, &EB, TORSION_2POWER_BYTES * 8);

    // Compute the isogeny EA1 -> EA1B
    isogB_prime1.curve = pk->EA1;
    isogB_prime1.degree[0] = POWER_OF_3;
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isogB_prime1.degree[i] = 0;
    }
    ec_set_zero(&isogB_prime1.ker_minus);
    // kernel = P + beta * Q
    xDBLMUL_bounded(&isogB_prime1.ker_plus, &pk->PQA13.P, one_scalar, &pk->PQA13.Q, beta_scalar, &pk->PQA13.PmQ, &isogB_prime1.curve, TORSION_3POWER_BYTES * 8);
    
    ec_eval_three(&EA1B, &isogB_prime1, (ec_point_t*)&eval_basis[0], 0);

    // Compute the isogeny EA -> EAB
    isogB_prime.curve = pk->EA;
    isogB_prime.degree[0] = TORSION_PLUS_ODD_POWERS[0];
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isogB_prime.degree[i] = 0;
    }
    ec_set_zero(&isogB_prime.ker_minus);

    // kernel = P + beta * Q
    xDBLMUL_bounded(&isogB_prime.ker_plus, &pk->PQ3.P, one_scalar, &pk->PQ3.Q, beta_scalar, &pk->PQ3.PmQ, &isogB_prime.curve, TORSION_3POWER_BYTES * 8);

    eval_basis[0] = EA_two;
    ec_eval_three(&EAB, &isogB_prime, (ec_point_t*)eval_basis, 3);
    EA_two = eval_basis[0];

    ct->EAB = EAB;

    // Masking evaluated basis points
    xMUL(&ct->PQ2_AB.P, &EA_two.P, omega_scalar, &EAB);
    xMUL(&ct->PQ2_AB.Q, &EA_two.Q, omega_inv_scalar, &EAB);
    xADD(&pointT, &EA_two.P, &EA_two.Q, &EA_two.PmQ);
    xDBLMUL_bounded(&ct->PQ2_AB.PmQ, &EA_two.P, omega_scalar, &EA_two.Q, omega_inv_scalar, &pointT, &EAB, TORSION_2POWER_BYTES * 8);

    unsigned char hash_input[2 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};

    fp2_t EA1B_j_inv;
    ec_j_inv(&EA1B_j_inv, &EA1B);

    memcpy(hash_input, &EA1B_j_inv.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &EA1B_j_inv.im[0], NWORDS_FIELD * RADIX / 8);

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

int decrypt(unsigned char *m, size_t *m_len, const inke_ct_t *ct, const inke_sk_t *sk) {
    unsigned char hash_input[2 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};
    digit_t T1_scalar[NWORDS_ORDER] = {0}, T2_scalar[NWORDS_ORDER] = {0};
    theta_chain_t hd_isog;
    theta_couple_curve_t EBAB;
    theta_couple_point_t T1, T2, T1m2;
    ec_basis_t eval_points;
    ibz_t alpha_inv, beta_inv, deg;
    ec_point_t pointT;

    ibz_init(&alpha_inv);
    ibz_init(&beta_inv);
    ibz_init(&deg);

    ibz_copy_digits(&deg, sk->deg, NWORDS_ORDER);
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
    xDBLMUL_bounded(&T1m2.P2, &ct->PQ2_AB.P, T1_scalar, &ct->PQ2_AB.Q, T2_scalar, &pointT, &EBAB.E2, TORSION_2POWER_BYTES * 8);

    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &EBAB, &T1, &T2, &T1m2, strategies[2], 1);

    fp2_t EA1B_j_inv;
    ec_j_inv(&EA1B_j_inv, &hd_isog.codomain.E1);

    memcpy(hash_input, &EA1B_j_inv.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &EA1B_j_inv.im[0], NWORDS_FIELD * RADIX / 8);

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
    return 1;
}

////
//// Key Encapsulation Mechanism using Fujisaki-Okamoto transform
////

const unsigned char G_hash_str[9] = "encrypt_";
const size_t G_hash_str_len = 8;

int ct_encode(unsigned char *encoded_ct, inke_ct_t *ct) {
    jac_point_t P2, Q2, Px, Qx;
    jac_point_t R, S, RmS;
    // total_len += NWORDS_FIELD * 2; // EB
    // total_len += NWORDS_ORDER * 6; // PQ2_B + PQxy_B -> 4/3 * lambda
    // total_len += NWORDS_FIELD * 2; // EAB
    // total_len += NWORDS_ORDER * 6; // PQ2_AB -> lambda

    fp2_encode(encoded_ct, &ct->EB.A);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 2 / 8, &ct->PQ2_B.P.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 4 / 8, &ct->PQ2_B.Q.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 6 / 8, &ct->PQ2_B.PmQ.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 8 / 8, &ct->EAB.A);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 10 / 8, &ct->PQ2_AB.P.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 12 / 8, &ct->PQ2_AB.Q.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 14 / 8, &ct->PQ2_AB.PmQ.x);

    return 1;
}

int ct_decode(inke_ct_t *ct, const unsigned char *encoded_ct) {
    ec_basis_t added_basis;

    fp2_decode(&ct->EB.A, encoded_ct);
    fp2_decode(&ct->PQ2_B.P.x, encoded_ct + NWORDS_FIELD * RADIX * 2 / 8);
    fp2_decode(&ct->PQ2_B.Q.x, encoded_ct + NWORDS_FIELD * RADIX * 4 / 8);
    fp2_decode(&ct->PQ2_B.PmQ.x, encoded_ct + NWORDS_FIELD * RADIX * 6 / 8);
    fp2_set_one(&ct->PQ2_B.P.z);
    fp2_set_one(&ct->PQ2_B.Q.z);
    fp2_set_one(&ct->PQ2_B.PmQ.z);
    fp2_decode(&ct->EAB.A, encoded_ct + NWORDS_FIELD * RADIX * 8 / 8);
    fp2_decode(&ct->PQ2_AB.P.x, encoded_ct + NWORDS_FIELD * RADIX * 10 / 8);
    fp2_decode(&ct->PQ2_AB.Q.x, encoded_ct + NWORDS_FIELD * RADIX * 12 / 8);
    fp2_decode(&ct->PQ2_AB.PmQ.x, encoded_ct + NWORDS_FIELD * RADIX * 14 / 8);
    fp2_set_one(&ct->PQ2_AB.P.z);
    fp2_set_one(&ct->PQ2_AB.Q.z);
    fp2_set_one(&ct->PQ2_AB.PmQ.z);

    return 1;
}

int encaps(unsigned char *key, inke_ct_t *ct, const inke_pk_t *pk) {
    unsigned char m[32];
    unsigned char tt[32 + G_hash_str_len];
    unsigned char gm[32];
    unsigned char encoded_ct[32 + NWORDS_FIELD * 16 * RADIX / 8];

    randombytes(m, 32);
    memcpy(tt, G_hash_str, G_hash_str_len);
    memcpy(tt, m, 32);
    SHAKE256(gm, 32, tt, 32);
    encrypt(ct, pk, m, 32, gm, 32);        // ct <- Enc(pk, m; G(m))
    memcpy(encoded_ct, m, 32);
    ct_encode(encoded_ct + 32, ct);
    SHAKE256(key, 32, encoded_ct, 32 + NWORDS_FIELD * 16 * RADIX / 8);    // K <- H(m, ct)
    return 1;
}

int decaps(unsigned char *key, inke_ct_t *ct, const inke_pk_t *pk, const inke_sk_t *sk, unsigned char *dummy_m) {
    unsigned char m[32];
    unsigned char tt[32 + G_hash_str_len];
    unsigned char gm[32];
    unsigned char test_ct_bytes[NWORDS_FIELD * 16 * RADIX / 8];
    unsigned char ct_bytes[32 + NWORDS_FIELD * 16 * RADIX / 8];
    size_t m_len;
    inke_ct_t test_ct;

    decrypt(m, &m_len, ct, sk);
    memcpy(tt, G_hash_str, G_hash_str_len);
    memcpy(tt, m, 32);
    SHAKE256(gm, 32, tt, 32);
    encrypt(&test_ct, pk, m, m_len, gm, 32);    // ct <- Enc(pk, m; G(m))
    ct_encode(test_ct_bytes, &test_ct);
    ct_encode(ct_bytes + 32, ct);
    if (memcmp(ct_bytes + 32, test_ct_bytes, NWORDS_FIELD * 16 * RADIX / 8) != 0) {
        memcpy(ct_bytes, dummy_m, 32);  // K <- H(s, ct)
    } else {
        memcpy(ct_bytes, m, 32);        // K <- H(m, ct)
    }
    SHAKE256(key, 32, ct_bytes, 32 + NWORDS_FIELD * 16 * RADIX / 8);

    return 1;
}