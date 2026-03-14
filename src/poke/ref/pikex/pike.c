#include <pike.h>
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

int ibz_random_matrix(ibz_mat_2x2_t mat22, const ibz_t *modulus, shake256ctx *state) {
    ibz_t gcd, det;
    ibz_init(&gcd);
    ibz_init(&det);
    while (1) {
        if (state != NULL) {
            ibz_rand_interval_with_state(&mat22[0][0], &ibz_const_zero, modulus, state);
            ibz_rand_interval_with_state(&mat22[0][1], &ibz_const_zero, modulus, state);
            ibz_rand_interval_with_state(&mat22[1][0], &ibz_const_zero, modulus, state);
            ibz_rand_interval_with_state(&mat22[1][1], &ibz_const_zero, modulus, state);
        } else {
            ibz_rand_interval(&mat22[0][0], &ibz_const_zero, modulus);
            ibz_rand_interval(&mat22[0][1], &ibz_const_zero, modulus);
            ibz_rand_interval(&mat22[1][0], &ibz_const_zero, modulus);
            ibz_rand_interval(&mat22[1][1], &ibz_const_zero, modulus);
        }
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

int eval_dimtwo_isog(theta_chain_t *phi, ibz_t *q, ec_basis_t *evalPQ, ec_basis_t *PQ, theta_couple_curve_t *E01, bool is_five) {
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
    if (is_five) ec_curve_to_basis_5(&RS, &E01->E2);
    else ec_curve_to_basis_3(&RS, &E01->E2);

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
    ec_basis_t imPQ_basis;
    imRS_basis.P = imR;
    imRS_basis.Q = imS;
    imRS_basis.PmQ = imRS;
    
    imPQ_basis.P = imP;
    imPQ_basis.Q = imQ;
    imPQ_basis.PmQ = imPQ;

    digit_t x1[NWORDS_ORDER] = {0}, x2[NWORDS_ORDER] = {0};
    digit_t x3[NWORDS_ORDER] = {0}, x4[NWORDS_ORDER] = {0};
    digit_t x5[NWORDS_ORDER] = {0}, x6[NWORDS_ORDER] = {0};
    digit_t t[NWORDS_ORDER] = {0};

    if (is_five) {
        const int nwords = (THREE_FIVE_bitlen + RADIX) / RADIX;
        ec_dlog_weil_5(x1, x2, x3, x4, &imRS_basis, &imPQ_basis, &phi->codomain.E1);
        ec_biscalar_mul_bounded(&evalPQ->P, &E01->E2, x1, x2, &RS, THREE_FIVE_bitlen);
        ec_biscalar_mul_bounded(&evalPQ->Q, &E01->E2, x3, x4, &RS, THREE_FIVE_bitlen);
        mp_add(x5, x1, FIVEpF, nwords);
        mp_sub(t, x5, x3, nwords);
        if (mp_compare(t, FIVEpF, nwords) != -1) {
            mp_sub(x5, t, FIVEpF, nwords);
        } else {
            memcpy(x5, t, NWORDS_ORDER * RADIX / 8);
        }

        mp_add(x6, x2, FIVEpF, nwords);
        mp_sub(t, x6, x4, nwords);
        if (mp_compare(t, FIVEpF, nwords) != -1) {
            mp_sub(x6, t, FIVEpF, nwords);
        } else {
            memcpy(x6, t, NWORDS_ORDER * RADIX / 8);
        }
        ec_biscalar_mul_bounded(&evalPQ->PmQ, &E01->E2, x5, x6, &RS, THREE_FIVE_bitlen);

    } else {
        const int nwords = (THREEpF_bitlen + RADIX) / RADIX;
        ec_dlog_weil_3(x1, x2, x3, x4, &imRS_basis, &imPQ_basis, &phi->codomain.E1);
        ec_biscalar_mul_bounded(&evalPQ->P, &E01->E2, x1, x2, &RS, THREEpF_bitlen);
        ec_biscalar_mul_bounded(&evalPQ->Q, &E01->E2, x3, x4, &RS, THREEpF_bitlen);
        mp_add(x5, x1, THREEpF, nwords);
        mp_sub(t, x5, x3, nwords);
        
        if (mp_compare(t, THREEpF, nwords) != -1) {
            mp_sub(x5, t, THREEpF, nwords);
        } else {
            memcpy(x5, t, NWORDS_ORDER * RADIX / 8);
        }

        mp_add(x6, x2, THREEpF, nwords);
        mp_sub(t, x6, x4, nwords);
        if (mp_compare(t, THREEpF, nwords) != -1) {
            mp_sub(x6, t, THREEpF, nwords);
        } else {
            memcpy(x6, t, NWORDS_ORDER * RADIX / 8);
        }

        ec_biscalar_mul_bounded(&evalPQ->PmQ, &E01->E2, x5, x6, &RS, THREEpF_bitlen);
    }

    return 1;
}

int eval_dimtwo_isog_single(theta_chain_t *phi, ibz_t *q, ec_point_t *evalP, ec_point_t *P, theta_couple_curve_t *E01, bool is_five) {
    theta_couple_point_t output_points, input_points;
    ec_point_t imP, imR, imS, imRS;

    copy_point(&input_points.P1, P);
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case(&output_points, phi, &input_points, E01);
    copy_point(&imP, &output_points.P1);

    ec_basis_t RS;
    if (is_five) ec_curve_to_basis_5(&RS, &E01->E2);
    else ec_curve_to_basis_3(&RS, &E01->E2);

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

    if (is_five) {
        ;; // Not used
        //ec_dlog_5(x, y, &imRS_basis, &imP, &phi->codomain.E1);
    }
    else
    {
        ec_dlog_weil_3_single(x, y, &imRS_basis, &imP, &phi->codomain.E1);
    }

    ec_curve_normalize_A24(&E01->E2);
    xDBLMUL_bounded(evalP, &RS.P, x, &RS.Q, y, &RS.PmQ, &E01->E2, THREEpF_bitlen);
    return 1;
}

int keygen(pike_sk_t *sk, pike_pk_t *pk) {
    ibz_t q, alpha, beta, iota, gamma1, gamma2, rhs, deg, tmp1, tmp2;
    ibz_t A, q_bound;
    ec_point_t pointT;
    jac_point_t tmpjS;
    
    ibz_init(&q); ibz_init(&alpha); ibz_init(&beta); ibz_init(&iota);
    ibz_init(&gamma1); ibz_init(&gamma2); ibz_init(&rhs); ibz_init(&deg);
    ibz_init(&tmp1); ibz_init(&tmp2);
    ibz_init(&A); ibz_init(&q_bound);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_div_2exp(&q_bound, &TORSION_PLUS_2POWER, 4);

    // Set q to a random value in the range [0, TORSION_PLUS_2POWER)
    for(int i = 0; i < 1000; i++) {
        ibz_rand_interval(&q, &ibz_const_zero, &q_bound);
        ibz_sub(&deg, &A, &q);
        ibz_mul(&tmp1, &TORSION_D, &ibz_const_two);
        ibz_mul(&tmp1, &tmp1, &ibz_const_three);
        ibz_mul(&tmp1, &tmp1, &ibz_const_five);
        ibz_mul(&tmp2, &q, &deg);
        ibz_gcd(&tmp1, &tmp1, &tmp2);
        if (ibz_is_one(&tmp1)) {
            break;
        }
    }
    ibz_mul(&rhs, &deg, &q);
    ibz_mul(&rhs, &rhs, &TORSION_PLUS_3POWER);
    ibz_mul(&rhs, &rhs, &TORSION_PLUS_3POWER);
    ibz_random_unit(&alpha, &A, NULL);
    ibz_random_unit(&beta, &A, NULL);
    ibz_random_unit(&gamma1, &TORSION_PLUS_3POWER, NULL);
    ibz_random_unit(&gamma2, &TORSION_ODD_MINUS, NULL);
    ibz_random_unit(&iota, &TORSION_D, NULL);

    memset(&sk->deg, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->alpha, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->beta, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->iota, 0, NWORDS_ORDER * RADIX / 8);

    ibz_to_digits(sk->deg, &q);
    ibz_to_digits(sk->alpha, &alpha);
    ibz_to_digits(sk->beta, &beta);
    ibz_to_digits(sk->iota, &iota);

    ec_curve_t curve = CURVE_E0;

    quat_alg_elem_t tau, tau_conj;
    quat_alg_elem_init(&tau);
    quat_alg_elem_init(&tau_conj);
    if(represent_integer(&tau, &rhs, &QUATALG_PINFTY) == 0) {
        printf("Failed to represent integer in non-diagonal form\n");
        return 1;
    }
    quat_alg_conj(&tau_conj, &tau);

    ec_basis_t E0_two, E0_three, E0_five, E0_three_conj;
    ec_point_t E0_point_S;
    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);
    copy_point(&E0_five.P, &BASIS_FIVE.P);
    copy_point(&E0_five.Q, &BASIS_FIVE.Q);
    copy_point(&E0_five.PmQ, &BASIS_FIVE.PmQ);
    copy_point(&E0_three_conj.P, &BASIS_THREE.P);
    copy_point(&E0_three_conj.Q, &BASIS_THREE.Q);
    copy_point(&E0_three_conj.PmQ, &BASIS_THREE.PmQ);
    copy_point(&E0_point_S, &BASIS_S.P);

    endomorphism_application_three_basis(&E0_three, &curve, &tau, TORSION_PLUS_ODD_POWERS[0]);
    endomorphism_application_three_basis(&E0_three_conj, &curve, &tau_conj, TORSION_PLUS_ODD_POWERS[0]);
    endomorphism_application_five_basis(&E0_five, &curve, &tau, TORSION_MINUS_ODD_POWERS[0]);
    endomorphism_application_even_basis(&E0_two, &curve, &tau, TORSION_PLUS_EVEN_POWER);
    endomorphism_application_single_point(&E0_point_S, &curve, &tau, &TORSION_D);

    ec_point_t K1, K2;
    ec_point_init(&K1);
    ec_point_init(&K2);
    ibz_t scalar, remainder;
    ibz_init(&scalar); ibz_init(&remainder);
    ibz_div(&scalar, &remainder, &TORSION_PLUS_3POWER, &ibz_const_three);
    digit_t uscl[NWORDS_ORDER] = {0};
    ibz_to_digits(uscl, &scalar);
    xMUL(&K1, &E0_three_conj.P, uscl, &curve);
    if(ec_is_zero(&K1)){
        copy_point(&K1, &E0_three_conj.Q);
    }      
    else{
        copy_point(&K1, &E0_three_conj.P);        
    }
    xMUL(&K2, &E0_three.P, uscl, &curve);
    if(ec_is_zero(&K2)){
        copy_point(&K2, &E0_three.Q);
    }      
    else{
        copy_point(&K2, &E0_three.P);        
    }

    ec_isog_odd_t isog1, isog2, isog3;
    ec_curve_t E1, E2, E3;
    isog1.curve = curve;
    copy_point(&isog1.ker_plus, &K1);
    ec_set_zero(&isog1.ker_minus);
    isog1.degree[0] = TORSION_PLUS_ODD_POWERS[0];
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isog1.degree[i] = 0;
    }
    ec_point_t eval_points[7];
    for(int i = 0; i < 7; i++){
        ec_point_init(&eval_points[i]);
    }
    copy_point(&eval_points[0], &E0_two.P);
    copy_point(&eval_points[1], &E0_two.Q);
    copy_point(&eval_points[2], &E0_two.PmQ);
    copy_point(&eval_points[3], &E0_five.P);
    copy_point(&eval_points[4], &E0_five.Q);
    copy_point(&eval_points[5], &E0_five.PmQ);

    copy_point(&eval_points[6], &E0_point_S);


    ec_point_t eval_points1[5];
    for(int i = 0; i < 5; i++){
        ec_point_init(&eval_points1[i]);
    }
    copy_point(&eval_points1[0], &BASIS_EVEN.P);
    copy_point(&eval_points1[1], &BASIS_EVEN.Q);
    copy_point(&eval_points1[2], &BASIS_EVEN.PmQ);
    copy_point(&eval_points1[3], &BASIS_THREE.P);
    copy_point(&eval_points1[4], &BASIS_THREE.Q);
    ec_eval_three(&E1, &isog1, (ec_point_t *)&eval_points1, 5);
    ec_point_t Kdual;
    ec_point_init(&Kdual);
    xMUL(&Kdual, &eval_points1[3], uscl, &E1);
    if(ec_is_zero(&Kdual)){
        copy_point(&Kdual, &eval_points1[4]);
    }
    else{
        copy_point(&Kdual, &eval_points1[3]);
    }
    isog2.curve = curve;
    copy_point(&isog2.ker_plus, &K2);
    ec_set_zero(&isog2.ker_minus);
    isog2.degree[0] = TORSION_PLUS_ODD_POWERS[0];
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isog2.degree[i] = 0;
    }
    ec_eval_three(&E2, &isog2, (ec_point_t *)&eval_points, 7);

    ibz_mul(&scalar, &TORSION_PLUS_3POWER, &q);
    ibz_mod(&scalar, &scalar, &TORSION_PLUS_2POWER);

    // Evaluating the theta-based 2-dim isogeny
    theta_couple_curve_t E01;
    theta_couple_point_t T1, T2, T1m2;
    theta_chain_t hd_isog;

    E01.E1 = E1; E01.E2 = E2;
    T1.P1 = eval_points1[0];
    T2.P1 = eval_points1[1];
    T1m2.P1 = eval_points1[2];
    ec_mul_ibz(&T1.P1, &E01.E1, &scalar, &T1.P1);
    ec_mul_ibz(&T2.P1, &E01.E1, &scalar, &T2.P1);
    ec_mul_ibz(&T1m2.P1, &E01.E1, &scalar, &T1m2.P1);
    copy_point(&T1.P2, &eval_points[0]);
    copy_point(&T2.P2, &eval_points[1]);
    copy_point(&T1m2.P2, &eval_points[2]);
    
    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &E01, &T1, &T2, &T1m2, strategies[2], 1);    
    eval_dimtwo_isog_single(&hd_isog, &q, &Kdual, &Kdual, &E01, 0);
    isog3.curve = E2;
    copy_point(&isog3.ker_plus, &Kdual);
    ec_set_zero(&isog3.ker_minus);
    isog3.degree[0] = TORSION_PLUS_ODD_POWERS[0];
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isog3.degree[i] = 0;
    }
    ec_eval_three(&E3, &isog3, (ec_point_t *)&eval_points, 7);
    
    pk->EA = E3;
    ec_mul(&pk->PQ2.P, &E3, sk->alpha, &eval_points[0]);
    ec_mul(&pk->PQ2.Q, &E3, sk->beta, &eval_points[1]);
    xADD(&pointT, &eval_points[0], &eval_points[1], &eval_points[2]);
    ec_curve_normalize_A24(&E3);
    xDBLMUL_bounded(&pk->PQ2.PmQ, &eval_points[0], sk->alpha, &eval_points[1], sk->beta, &pointT, &E3, POWER_OF_2);  
    ec_mul_ibz(&pk->PQ5.P, &E3, &gamma2, &eval_points[3]);
    ec_mul_ibz(&pk->PQ5.Q, &E3, &gamma2, &eval_points[4]);
    ec_mul_ibz(&pk->PQ5.PmQ, &E3, &gamma2, &eval_points[5]);
    E01.E1 = curve; E01.E2 = E3;
    ibz_mul(&scalar, &TORSION_PLUS_3POWER, &scalar);
    ibz_mod(&scalar, &scalar, &TORSION_PLUS_2POWER);
    copy_point(&T1.P2, &eval_points[0]);
    copy_point(&T2.P2, &eval_points[1]);
    copy_point(&T1m2.P2, &eval_points[2]);
    T1.P1 = BASIS_EVEN.P;
    T2.P1 = BASIS_EVEN.Q;
    T1m2.P1 = BASIS_EVEN.PmQ;
    ec_mul_ibz(&T1.P1, &curve, &scalar, &T1.P1);
    ec_mul_ibz(&T2.P1, &curve, &scalar, &T2.P1);
    ec_mul_ibz(&T1m2.P1, &curve, &scalar, &T1m2.P1);
    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &E01, &T1, &T2, &T1m2, strategies[2], 1); 
    copy_point(&pk->PQ3.P, &BASIS_THREE.P);
    copy_point(&pk->PQ3.Q, &BASIS_THREE.Q);
    copy_point(&pk->PQ3.PmQ, &BASIS_THREE.PmQ);
    eval_dimtwo_isog(&hd_isog, &deg, &pk->PQ3, &pk->PQ3, &E01, false);
    ec_mul_ibz(&pk->PQ3.P, &E3, &gamma1, &pk->PQ3.P);
    ec_mul_ibz(&pk->PQ3.Q, &E3, &gamma1, &pk->PQ3.Q);
    ec_mul_ibz(&pk->PQ3.PmQ, &E3, &gamma1, &pk->PQ3.PmQ);

    // Evaluating a D-torsion single point X
    ec_mul(&pk->imPs, &pk->EA, sk->iota, &eval_points[6]);
    
    ibz_finalize(&scalar);
    ibz_finalize(&remainder);
    ibz_finalize(&q);
    ibz_finalize(&alpha);
    ibz_finalize(&beta);
    ibz_finalize(&iota);
    ibz_finalize(&gamma1);
    ibz_finalize(&gamma2);
    ibz_finalize(&rhs);
    ibz_finalize(&deg);
    ibz_finalize(&tmp1);
    ibz_finalize(&tmp2);
    ibz_finalize(&A);
    ibz_finalize(&q_bound);
    quat_alg_elem_finalize(&tau);
    quat_alg_elem_finalize(&tau_conj);
    return 1; 
}

int encrypt(pike_ct_t *ct, const pike_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len) {
    ec_isog_odd_t isogB1, isogB2;
    ibz_t beta1, beta2, omega, omega_inv, A, t1, t2, t1t2;
    ec_curve_t EB, EAB;
    ec_basis_t E0_two, EA_two;
    ec_point_t pointT, eval_points[5];
    shake256ctx state;

    digit_t beta1_scalar[NWORDS_ORDER] = {0}, beta2_scalar[NWORDS_ORDER] = {0}, omega_scalar[NWORDS_ORDER] = {0}, omega_inv_scalar[NWORDS_ORDER] = {0}, t1_scalar[NWORDS_ORDER] = {0}, t2_scalar[NWORDS_ORDER] = {0}, one_scalar[NWORDS_ORDER] = {1};

    ibz_init(&beta1);
    ibz_init(&beta2);
    ibz_init(&omega);
    ibz_init(&omega_inv);
    ibz_init(&A);
    ibz_init(&t1);
    ibz_init(&t2);
    ibz_init(&t1t2);


    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    if (seed != NULL) {
        shake256_absorb(&state, seed, seed_len);
        ibz_random_unit(&beta1, &TORSION_PLUS_3POWER, &state);
        ibz_random_unit(&beta2, &TORSION_ODD_MINUS, &state);
        ibz_random_unit(&omega, &A, &state);
        ibz_random_unit(&t1, &TORSION_D, &state);
        ibz_random_unit(&t2, &TORSION_D, &state);
        ibz_mul(&t1t2, &t1, &t2);
        ibz_mod(&t1t2, &t1t2, &TORSION_D);
        shake256_ctx_release(&state);
    }
    else {
        ibz_random_unit(&beta1, &TORSION_PLUS_3POWER, NULL);
        ibz_random_unit(&beta2, &TORSION_ODD_MINUS, NULL);
        ibz_random_unit(&omega, &A, NULL);
        ibz_random_unit(&t1, &TORSION_D, NULL);
        ibz_random_unit(&t2, &TORSION_D, NULL);
        ibz_mul(&t1t2, &t1, &t2);
        ibz_mod(&t1t2, &t1t2, &TORSION_D);
    }
    ibz_invmod(&omega_inv, &omega, &A);
    ibz_to_digits(omega_scalar, &omega);
    ibz_to_digits(omega_inv_scalar, &omega_inv);
    ibz_to_digits(beta1_scalar, &beta1);
    ibz_to_digits(beta2_scalar, &beta2);
    ibz_to_digits(t1_scalar, &t1);
    ibz_to_digits(t2_scalar, &t2);

    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&EA_two.P, &pk->PQ2.P);
    copy_point(&EA_two.Q, &pk->PQ2.Q);
    copy_point(&EA_two.PmQ, &pk->PQ2.PmQ);

    // Compute the isogeny E0 -> EB
    isogB1.curve = CURVE_E0;
    isogB1.degree[0] = POWER_OF_3;
    isogB1.degree[1] = 0;
    ec_biscalar_mul_bounded(&isogB1.ker_plus, &CURVE_E0, one_scalar, beta1_scalar, &BASIS_THREE, TORSION_3POWER_BYTES * 8);
    ec_set_zero(&isogB1.ker_minus);
    ec_biscalar_mul_bounded(&pointT, &CURVE_E0, one_scalar, beta2_scalar, &BASIS_FIVE, FIVEpF_bitlen);
    copy_point(&eval_points[0], &BASIS_EVEN.P);
    copy_point(&eval_points[1], &BASIS_EVEN.Q);
    copy_point(&eval_points[2], &BASIS_EVEN.PmQ);
    copy_point(&eval_points[3], &BASIS_S.Q);
    copy_point(&eval_points[4], &pointT);
    ec_eval_three(&EB, &isogB1, (ec_point_t*)eval_points, 5);
    
    isogB2.curve = EB;
    isogB2.degree[0] = 0;
    isogB2.degree[1] = POWER_OF_5;
    ec_set_zero(&isogB2.ker_plus);
    copy_point(&isogB2.ker_minus, &eval_points[4]);
    ec_eval_five(&EB, &isogB2, (ec_point_t*)eval_points, 4);

    ct->EB = EB;
    // Masking evaluated points
    xMUL(&ct->PQ2_B.P, &eval_points[0], omega_scalar, &EB);
    xMUL(&ct->PQ2_B.Q, &eval_points[1], omega_inv_scalar, &EB);
    // necessary?
    xADD(&pointT, &eval_points[0], &eval_points[1], &eval_points[2]);

    ec_curve_normalize_A24(&EB);
    xDBLMUL_bounded(&ct->PQ2_B.PmQ, &eval_points[0], omega_scalar, &eval_points[1], omega_inv_scalar, &pointT, &EB, POWER_OF_2);
    xMUL(&ct->QsB, &eval_points[3], t1_scalar, &EB);

    isogB1.curve = pk->EA;
    isogB1.degree[0] = POWER_OF_3;
    isogB1.degree[1] = 0;
    ec_biscalar_mul_bounded(&isogB1.ker_plus, &pk->EA, one_scalar, beta1_scalar, &pk->PQ3, THREEpF_bitlen);
    ec_set_zero(&isogB1.ker_minus);
    ec_biscalar_mul_bounded(&pointT, &pk->EA, one_scalar, beta2_scalar, &pk->PQ5, FIVEpF_bitlen);
    copy_point(&eval_points[0], &pk->PQ2.P);
    copy_point(&eval_points[1], &pk->PQ2.Q);
    copy_point(&eval_points[2], &pk->PQ2.PmQ);
    copy_point(&eval_points[3], &pk->imPs);
    copy_point(&eval_points[4], &pointT);
    ec_eval_three(&EAB, &isogB1, (ec_point_t*)eval_points, 5);

    isogB2.curve = EAB;
    isogB2.degree[0] = 0;
    isogB2.degree[1] = POWER_OF_5;
    ec_set_zero(&isogB2.ker_plus);
    copy_point(&isogB2.ker_minus, &eval_points[4]);
    ec_eval_five(&EAB, &isogB2, (ec_point_t*)eval_points, 4);

    ct->EAB = EAB;
    // Masking evaluated points
    xMUL(&ct->PQ2_AB.P, &eval_points[0], omega_scalar, &EAB);
    xMUL(&ct->PQ2_AB.Q, &eval_points[1], omega_inv_scalar, &EAB);
    // necessary?
    xADD(&pointT, &eval_points[0], &eval_points[1], &eval_points[2]);

    ec_curve_normalize_A24(&EAB);
    xDBLMUL_bounded(&ct->PQ2_AB.PmQ, &eval_points[0], omega_scalar, &eval_points[1], omega_inv_scalar, &pointT, &EAB, POWER_OF_2);
    
    xMUL(&ct->PsAB, &eval_points[3], t2_scalar, &EAB);

    unsigned char hash_input[NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};

    fp2_t shared_sec;
    digit_t exp[NWORDS_ORDER] = {0};
    ibz_to_digits(exp, &t1t2);
    fp2_pow_vartime(&shared_sec, &Pairing_value, exp, TORSION_D->_mp_size);
    // fp_add(&shared_sec.re, &shared_sec.re, &shared_sec.re);

    memcpy(hash_input, &shared_sec.re, NWORDS_FIELD * RADIX / 8);

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


    ibz_finalize(&A);
    ibz_finalize(&beta1);
    ibz_finalize(&beta2);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&omega);
    ibz_finalize(&omega_inv);
    return 1;
}

int decrypt(unsigned char *m, size_t *m_len, const pike_ct_t *ct, const pike_sk_t *sk) {

    unsigned char hash_input[NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};
    digit_t T1_scalar[NWORDS_ORDER] = {0}, T2_scalar[NWORDS_ORDER] = {0};
    theta_chain_t hd_isog;
    theta_couple_curve_t EBAB;
    theta_couple_point_t T1, T2, T1m2, tmp;
    ibz_t alpha_inv, beta_inv, deg, iota, scalar, A, d;
    ec_point_t pointT, Psmid, Qsmid, PmQsmid;
    jac_point_t jPsmid, jQsmid, jPmQsmid;

    ibz_init(&alpha_inv);
    ibz_init(&beta_inv);
    ibz_init(&deg);
    ibz_init(&iota);
    ibz_init(&scalar);
    ibz_init(&A);
    ibz_init(&d);

    ibz_copy_digits(&deg, sk->deg, NWORDS_ORDER);
    ibz_mul(&scalar, &deg, &TORSION_ODD_PLUS);
    ibz_mod(&scalar, &scalar, &TORSION_PLUS_2POWER);
    ibz_mul(&scalar, &TORSION_ODD_PLUS, &scalar);
    ibz_mod(&scalar, &scalar, &TORSION_PLUS_2POWER);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_copy_digits(&alpha_inv, sk->alpha, NWORDS_ORDER);
    ibz_copy_digits(&beta_inv, sk->beta, NWORDS_ORDER);
    ibz_invmod(&alpha_inv, &alpha_inv, &TORSION_PLUS_2POWER);
    ibz_invmod(&beta_inv, &beta_inv, &TORSION_PLUS_2POWER);

    ibz_sub(&d, &A, &deg);
    ibz_mul(&d, &deg, &d); // d = q(2^a-q)
    ibz_mod(&d, &d, &TORSION_D);
    ibz_copy_digits(&iota, sk->iota, NWORDS_ORDER);
    ibz_mul(&d, &iota, &d); // d = q(2^a-q) * iota
    ibz_mod(&d, &d, &TORSION_D);
    ibz_mul(&d, &d, &TORSION_ODD_PLUS); // d = q(2^a-q) * iota * C1
    ibz_mod(&d, &d, &TORSION_D);
    ibz_mul(&d, &d, &TORSION_ODD_PLUS); // d = q(2^a-q) * iota * C1^2
    ibz_mod(&d, &d, &TORSION_D);
    ibz_mul(&d, &d, &TORSION_ODD_PLUS); // d = q(2^a-q) * iota * C1^3
    ibz_mod(&d, &d, &TORSION_D);
    ibz_mul(&d, &d, &TORSION_ODD_MINUS);  // d = q(2^a-q) * iota * C1^3 * C2
    ibz_mod(&d, &d, &TORSION_D);
    ibz_invmod(&d, &d, &TORSION_D);
    // ibz_mod(&d, &d, &TORSION_D);

    EBAB.E1 = ct->EB;
    EBAB.E2 = ct->EAB;

    ec_mul_ibz(&T1.P1, &EBAB.E1, &scalar, &ct->PQ2_B.P);
    ec_mul_ibz(&T2.P1, &EBAB.E1, &scalar, &ct->PQ2_B.Q);
    ec_mul_ibz(&T1m2.P1, &EBAB.E1, &scalar, &ct->PQ2_B.PmQ);
    ibz_to_digits(T1_scalar, &alpha_inv);
    ibz_to_digits(T2_scalar, &beta_inv);
    xMUL(&T1.P2, &ct->PQ2_AB.P, T1_scalar, &EBAB.E2);
    xMUL(&T2.P2, &ct->PQ2_AB.Q, T2_scalar, &EBAB.E2);
    xADD(&pointT, &ct->PQ2_AB.P, &ct->PQ2_AB.Q, &ct->PQ2_AB.PmQ);
    ec_curve_normalize_A24(&EBAB.E2);
    xDBLMUL_bounded(&T1m2.P2, &ct->PQ2_AB.P, T1_scalar, &ct->PQ2_AB.Q, T2_scalar, &pointT, &EBAB.E2, POWER_OF_2);
    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &EBAB, &T1, &T2, &T1m2, strategies[2], 1);

    tmp.P1 = ct->QsB;
    ec_set_zero(&tmp.P2);
    theta_chain_eval_special_case(&tmp, &hd_isog, &tmp, &EBAB);
    copy_point(&Psmid, &tmp.P1);
    lift_point(&jPsmid, &Psmid, &hd_isog.codomain.E1);

    tmp.P2 = ct->PsAB;
    ec_set_zero(&tmp.P1);
    theta_chain_eval_special_case(&tmp, &hd_isog, &tmp, &EBAB);
    copy_point(&Qsmid, &tmp.P1);
    lift_point(&jQsmid, &Qsmid, &hd_isog.codomain.E1);

    ADD(&jPmQsmid, &jPsmid, &jQsmid, &hd_isog.codomain.E1);
    jac_to_xz(&PmQsmid, &jPmQsmid);

    fp2_t shared_sec;
    digit_t n[NWORDS_ORDER] = {0};
    digit_t exp[NWORDS_ORDER] = {0};
    ibz_to_digits(n, &TORSION_D);
    ibz_to_digits(exp, &d);
    ec_curve_normalize_A24(&hd_isog.codomain.E1);
    tate_odd(&shared_sec, n, P_COFACTOR_FOR_6FG_BITLENGTH, &Qsmid, &Psmid, &PmQsmid, &hd_isog.codomain.E1.A24);
    // TODO: constant-time implementation (use Lucas sequence)
    fp2_pow_vartime(&shared_sec, &shared_sec, exp, TORSION_D->_mp_size);
    // fp_add(&shared_sec.re, &shared_sec.re, &shared_sec.re);

    memcpy(hash_input, &shared_sec.re, NWORDS_FIELD * RADIX / 8);

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
    ibz_finalize(&iota);
    ibz_finalize(&scalar);
    ibz_finalize(&A);
    ibz_finalize(&d);
    return 1;
}

////
//// Key Encapsulation Mechanism using Fujisaki-Okamoto transform
////

const unsigned char G_hash_str[9] = "encrypt_";
const size_t G_hash_str_len = 8;

int ct_encode(unsigned char *encoded_ct, pike_ct_t *ct) {
    ec_basis_t added_basis;
    jac_point_t P2, Q2, Px, Qx;
    jac_point_t R, S, RmS;
    // total_len += NWORDS_FIELD * 2; // EB
    // total_len += NWORDS_ORDER * 6; // PQ2_B + PQxy_B -> 4/3 * lambda
    // total_len += NWORDS_FIELD * 2; // EAB
    // total_len += NWORDS_ORDER * 6; // PQ2_AB -> lambda
    lift_basis(&P2, &Q2, &ct->PQ2_B, &ct->EB);
    // lift_basis(&Px, &Qx, &ct->PQxy_B, &ct->EB);
    ADD(&R, &P2, &Px, &ct->EB);
    ADD(&S, &Q2, &Qx, &ct->EB);
    jac_neg(&RmS, &S);
    ADD(&RmS, &R, &RmS, &ct->EB);
    jac_to_xz(&added_basis.P, &R);
    jac_to_xz(&added_basis.Q, &S);
    jac_to_xz(&added_basis.PmQ, &RmS);

    fp2_encode(encoded_ct, &ct->EB.A);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 2 / 8, &added_basis.P.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 4 / 8, &added_basis.Q.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 6 / 8, &added_basis.PmQ.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 8 / 8, &ct->EAB.A);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 10 / 8, &ct->PQ2_AB.P.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 12 / 8, &ct->PQ2_AB.Q.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 14 / 8, &ct->PQ2_AB.PmQ.x);

    return 1;
}

int ct_decode(pike_ct_t *ct, const unsigned char *encoded_ct) {
    ec_basis_t added_basis;
    ibz_t five_inv, two_inv;
    ibz_init(&five_inv); ibz_init(&two_inv);

    fp2_decode(&ct->EB.A, encoded_ct);
    fp2_decode(&added_basis.P.x, encoded_ct + NWORDS_FIELD * RADIX * 2 / 8);
    fp2_decode(&added_basis.Q.x, encoded_ct + NWORDS_FIELD * RADIX * 4 / 8);
    fp2_decode(&added_basis.PmQ.x, encoded_ct + NWORDS_FIELD * RADIX * 6 / 8);
    fp2_set_one(&added_basis.P.z);
    fp2_set_one(&added_basis.Q.z);
    fp2_set_one(&added_basis.PmQ.z);
    fp2_decode(&ct->EAB.A, encoded_ct + NWORDS_FIELD * RADIX * 8 / 8);
    fp2_decode(&ct->PQ2_AB.P.x, encoded_ct + NWORDS_FIELD * RADIX * 10 / 8);
    fp2_decode(&ct->PQ2_AB.Q.x, encoded_ct + NWORDS_FIELD * RADIX * 12 / 8);
    fp2_decode(&ct->PQ2_AB.PmQ.x, encoded_ct + NWORDS_FIELD * RADIX * 14 / 8);
    fp2_set_one(&ct->PQ2_AB.P.z);
    fp2_set_one(&ct->PQ2_AB.Q.z);
    fp2_set_one(&ct->PQ2_AB.PmQ.z);

    ibz_invmod(&five_inv, &TORSION_PLUS_CPOWER, &TORSION_PLUS_2POWER);
    ibz_invmod(&two_inv, &TORSION_PLUS_2POWER, &TORSION_PLUS_CPOWER);

    ec_mul_ibz(&ct->PQ2_B.P, &ct->EB, &TORSION_PLUS_CPOWER, &added_basis.P);
    ec_mul_ibz(&ct->PQ2_B.Q, &ct->EB, &TORSION_PLUS_CPOWER, &added_basis.Q);
    ec_mul_ibz(&ct->PQ2_B.PmQ, &ct->EB, &TORSION_PLUS_CPOWER, &added_basis.PmQ);
    ec_mul_ibz(&ct->PQ2_B.P, &ct->EB, &five_inv, &added_basis.P);
    ec_mul_ibz(&ct->PQ2_B.Q, &ct->EB, &five_inv, &added_basis.Q);
    ec_mul_ibz(&ct->PQ2_B.PmQ, &ct->EB, &five_inv, &added_basis.PmQ);

    ibz_finalize(&five_inv);
    ibz_finalize(&two_inv);

    return 1;
}

int encaps(unsigned char *key, pike_ct_t *ct, const pike_pk_t *pk) {
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

int decaps(unsigned char *key, pike_ct_t *ct, const pike_pk_t *pk, const pike_sk_t *sk, unsigned char *dummy_m) {
    unsigned char m[32];
    unsigned char tt[32 + G_hash_str_len];
    unsigned char gm[32];
    unsigned char test_ct_bytes[NWORDS_FIELD * 16 * RADIX / 8];
    unsigned char ct_bytes[32 + NWORDS_FIELD * 16 * RADIX / 8];
    size_t m_len;
    pike_ct_t test_ct;

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