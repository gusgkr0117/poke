#include <hd.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <klpt.h>
#include <quaternion.h>
#include <gmp.h>
#include <intbig.h>
#include <ec.h>
#include <rng.h>

// Solve x^2+y^2=p using cornacchia algorithm
// This function assumes that p is a non-negative integer and returns 1 if a solution is found,
// and 0 if no solution exists. The solution is stored in x and y.
int ibz_qf_solve(ibz_t *x, ibz_t *y, ibz_t *p){
    ibz_t d, p_minus_one;
    ibz_init(&d);
    ibz_init(&p_minus_one);
    ibz_sub(&p_minus_one, p, &ibz_const_one); // p - 1

    ibz_sqrt_mod_p(&d, &p_minus_one, p);

    // Euclidean algorithm for d and p
    ibz_t a, b, r;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&r);
    ibz_copy(&a, &d);
    ibz_copy(&b, p);
    while (!ibz_is_zero(&b)) {
        ibz_mod(&r, &a, &b);
        ibz_copy(&a, &b);
        ibz_copy(&b, &r);
    }

    return 1; // Solution found
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
    
    ibz_init(&q); ibz_init(&alpha); ibz_init(&beta);
    ibz_init(&gamma); ibz_init(&delta); ibz_init(&rhs); ibz_init(&deg);
    // Set q to a random value in the range [0, TORSION_PLUS_2POWER)
    ibz_rand_interval(&q, &ibz_const_zero, &TORSION_PLUS_2POWER);
    ibz_sub(&deg, &TORSION_PLUS_2POWER, &q);
    ibz_mul(&rhs, &deg, &q);
    ibz_mul(&rhs, &rhs, &TORSION_PLUS_3POWER);
    ibz_random_unit(&alpha, &TORSION_PLUS_2POWER);
    ibz_random_unit(&beta, &TORSION_PLUS_2POWER);
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
    return 0; 
}

int main() {
    int res = 1;

    test();

    return res;
}