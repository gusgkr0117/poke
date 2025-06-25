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
// helper function to apply a matrix to a basis of E[3^f]
// works in place
void
matrix_application_three_basis(ec_basis_t *bas, ec_curve_t *E, ibz_mat_2x2_t *mat, int f)
{
    digit_t scalars[2][NWORDS_FIELD] = { 0 };

    ibz_t tmp, pow_three;
    ibz_init(&tmp);
    ibz_init(&pow_three);
    ibz_pow(&pow_three, &ibz_const_three, f);

    ec_basis_t tmp_bas;
    copy_point(&tmp_bas.P, &bas->P);
    copy_point(&tmp_bas.Q, &bas->Q);
    copy_point(&tmp_bas.PmQ, &bas->PmQ);

    // reduction mod 3f
    ibz_mod(&(*mat)[0][0], &(*mat)[0][0], &pow_three);
    ibz_mod(&(*mat)[0][1], &(*mat)[0][1], &pow_three);
    ibz_mod(&(*mat)[1][0], &(*mat)[1][0], &pow_three);
    ibz_mod(&(*mat)[1][1], &(*mat)[1][1], &pow_three);

    
    ibz_mat_2x2_print(mat);

    jac_point_t P, Q, R;
    

    // first basis element
    ibz_to_digit_array(scalars[0], &(*mat)[0][0]);
    // ibz_set(&mat[0][1],0);
    ibz_to_digit_array(scalars[1], &(*mat)[1][0]);
    lift_basis(&P, &Q, &tmp_bas, E);
    DBLMUL_generic(&R, &P, scalars[0], &Q, scalars[1], E, NWORDS_FIELD);
    jac_to_xz(&bas->P, &R);
    // ec_biscalar_mul(&bas->P, E, scalars[0], scalars[1], &tmp_bas);
    point_print("bas->P: ", bas->P);
    ibz_to_digit_array(scalars[0], &(*mat)[0][1]);
    ibz_to_digit_array(scalars[1], &(*mat)[1][1]);
    lift_basis(&P, &Q, &tmp_bas, E);
    DBLMUL_generic(&R, &P, scalars[0], &Q, scalars[1], E, NWORDS_FIELD);
    jac_to_xz(&bas->Q, &R);
    // ec_biscalar_mul(&bas->Q, E, scalars[0], scalars[1], &tmp_bas);

    ibz_sub(&tmp, &(*mat)[0][0], &(*mat)[0][1]);
    ibz_mod(&tmp, &tmp, &pow_three);
    ibz_to_digit_array(scalars[0], &tmp);
    ibz_sub(&tmp, &(*mat)[1][0], &(*mat)[1][1]);
    ibz_mod(&tmp, &tmp, &pow_three);
    ibz_to_digit_array(scalars[1], &tmp);
    lift_basis(&P, &Q, &tmp_bas, E);
    DBLMUL_generic(&R, &P, scalars[0], &Q, scalars[1], E, NWORDS_FIELD);
    jac_to_xz(&bas->PmQ, &R);
    // ec_biscalar_mul(&bas->PmQ, E, scalars[0], scalars[1], &tmp_bas);

    ibz_finalize(&tmp);
    ibz_finalize(&pow_three);
}

void endomorphism_application_three_basis(ec_basis_t *bas, ec_curve_t *E, quat_alg_elem_t *theta, int f) {
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_vec_4_t coeffs;
    ibz_vec_4_init(&coeffs);
    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);

    ibz_t content;
    ibz_init(&content);

    // // decomposing theta on the basis
    quat_alg_make_primitive(&coeffs, &content, theta, &MAXORD_O0, &QUATALG_PINFTY);
    assert(ibz_get(&content) % 2 == 1);

    ibz_vec_4_print(&coeffs);

    ibz_set(&mat[0][0], 0);
    ibz_set(&mat[0][1], 0);
    ibz_set(&mat[1][0], 0);
    ibz_set(&mat[1][1], 0);

    // computing the matrix
    for (unsigned i = 0; i < 2; ++i) {
        ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
        for (unsigned j = 0; j < 2; ++j) {
            ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&mat[i][j], &mat[i][j], &content);
            // ibz_mod(&mat[i][j],&mat[i][j],&twopow);
        }
    }

    // and now we apply it
    matrix_application_three_basis(bas, E, &mat, f);

    ibz_vec_4_finalize(&coeffs);
    ibz_mat_2x2_finalize(&mat);
    ibz_finalize(&content);
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
    fp2_print("j_invariant",&j_inv);

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
    gmp_printf("nrd(tau) : %Qd\n", tau_norm);
    gmp_printf("rhs : %Zd\n", rhs);

    ec_basis_t E0_two, E0_three;
    copy_point(&E0_two.P, &BASIS_EVEN.P);
    copy_point(&E0_two.Q, &BASIS_EVEN.Q);
    copy_point(&E0_two.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&E0_three.P, &BASIS_THREE.P);
    copy_point(&E0_three.Q, &BASIS_THREE.Q);
    copy_point(&E0_three.PmQ, &BASIS_THREE.PmQ);

    point_print("P2", E0_two.P);
    point_print("Q2", E0_two.Q);
    point_print("P3", E0_three.P);
    point_print("Q3", E0_three.Q);

    // ibz_set(&(tau.coord[0]), 0);
    // ibz_set(&(tau.coord[1]), 0);
    // ibz_set(&(tau.coord[2]), 1);
    // ibz_set(&(tau.coord[3]), 0);

    // ibz_set(&(tau.denom), 1);

    endomorphism_application_three_basis(&E0_three, &curve, &tau, TORSION_PLUS_ODD_POWERS[0]);
    endomorphism_application_even_basis(&E0_two, &curve, &tau, TORSION_PLUS_EVEN_POWER);

    point_print("tau(P2): ", E0_two.P);
    point_print("tau(Q2): ", E0_two.Q);
    point_print("tau(P3): ", E0_three.P);
    point_print("tau(Q3): ", E0_three.Q);
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
    fp2_print("e(P,Q) : ", &w0);

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
    fp2_print("e(phi(P),phi(Q)) : ", &w1);
    fp2_mul(&w0tw1, &w0, &w1);
    fp2_mul(&w0tw1, &w0tw1, &w0tw1);
    fp2_mul(&w0tw1, &w0tw1, &w0tw1);

    fp2_print("(e(P, Q) x e(phi(P),phi(Q)))^4 = ", &w0tw1);
    assert(fp2_is_one(&w0tw1));

    theta_chain_comput_strategy(&hd_isog, TORSION_PLUS_EVEN_POWER - 2, &E01, &T1, &T2, &T1m2, strategies[2], 1);

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