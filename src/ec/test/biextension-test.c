
#include <tutil.h>
#include <ec_params.h>
#include "ec.h"
#include "fp.h"
#include "intbig.h"
#include <stdio.h>
#include <curve_extras.h>
#include <biextension.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>

int
compare_words(digit_t *a, digit_t *b, unsigned int nwords)
{ // Comparing "nword" elements, a=b? : (1) a>b, (0) a=b, (-1) a<b
  // SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING
  // ONLY.
    int i;

    for (i = nwords - 1; i >= 0; i--) {
        if (a[i] > b[i])
            return 1;
        else if (a[i] < b[i])
            return -1;
    }

    return 0;
}

int dlog_2_test() {
    ec_curve_t E0 = CURVE_E0;
    ec_basis_t even_basis = BASIS_EVEN;
    ibz_t a, b, c, d;
    digit_t x[NWORDS_ORDER] = {0}, y[NWORDS_ORDER] = {0};
    ibz_init(&a); ibz_init(&b);
    ibz_init(&c); ibz_init(&d);
    ibz_rand_interval(&a, &ibz_const_zero, &TORSION_PLUS_2POWER);
    ibz_rand_interval(&b, &ibz_const_zero, &TORSION_PLUS_2POWER);

    ibz_to_digits(x, &a);
    ibz_to_digits(y, &b);

    ibz_copy_digits(&a, x, NWORDS_ORDER);
    ibz_copy_digits(&b, y, NWORDS_ORDER);

    ec_point_t R = {0}, S = {0};
    digit_t z[NWORDS_ORDER] = {0}, w[NWORDS_ORDER] = {0};
    ec_biscalar_mul(&R, &E0, x, y, &even_basis);

    ec_dlog_2(z ,w, &even_basis, &R, &E0);
    ec_biscalar_mul(&S, &E0, z, w, &even_basis);
    
    if (!is_point_equal(&S, &R)) {
        point_print("R : ", R);
        point_print("S : ", S);
        printf("dlog_2_test failed\n");
        return 1;
    }

    ibz_finalize(&a); ibz_finalize(&b);
    return 0;
}

int dlog_3_test() {
    ec_curve_t E0 = CURVE_E0;
    ec_basis_t three_basis = BASIS_THREE;
    ibz_t a, b, c, d;
    digit_t x[NWORDS_ORDER] = {0}, y[NWORDS_ORDER] = {0};
    ibz_init(&a); ibz_init(&b);
    ibz_init(&c); ibz_init(&d);
    ibz_rand_interval(&a, &ibz_const_zero, &TORSION_PLUS_3POWER);
    ibz_rand_interval(&b, &ibz_const_zero, &TORSION_PLUS_3POWER);
    // ibz_set(&a, 1);
    // ibz_set(&b, 1);

    ibz_to_digits(x, &a);
    ibz_to_digits(y, &b);

    ibz_copy_digits(&a, x, NWORDS_ORDER);
    ibz_copy_digits(&b, y, NWORDS_ORDER);

    ec_point_t R = {0}, S = {0};
    digit_t z[NWORDS_ORDER] = {0}, w[NWORDS_ORDER] = {0};
    ec_biscalar_mul(&R, &E0, x, y, &three_basis);

    ec_dlog_3(z ,w, &three_basis, &R, &E0);
    ec_biscalar_mul(&S, &E0, z, w, &three_basis);
    
    if (!is_point_equal(&S, &R)) {
        // point_print("R : ", R);
        // point_print("S : ", S);
        ibz_copy_digits(&c, z, NWORDS_ORDER);
        ibz_copy_digits(&d, w, NWORDS_ORDER);
        gmp_printf("x, y: %Zx, %Zx\n", a, b);
        gmp_printf("z, w: %Zx, %Zx\n", c, d);

        printf("dlog_3_test failed\n");
        return 1;
    }

    ibz_finalize(&a); ibz_finalize(&b);
    return 0;
}

int dlog_5_test() {
    ec_curve_t E0 = {{{0x3940dda05af233c7, 0xe0c71dc8daee1591, 0x35e393cf98c157e9, 0xd2252bee0da49865, 0xcc4a01728e08dada, 0xbffac66fba663d66, 0x000003c75e8d09a6}, {0x9a8b8a94847e2c26, 0x58fdcc22f9ad5a13, 0x195169d7cd6ca624, 0xbd48d79e4afd5257, 0xa8c29a3da26b5655, 0x982604249a9ddd98, 0x00004b889895687f}}, {{0x1}}};

    // ec_basis_t five_basis = BASIS_FIVE;
    ec_basis_t five_basis;
    five_basis.P = (ec_point_t){{{0x0b3c5b941eddf0fd, 0x063705eb9a1d3819, 0xbe4ef1db81778f65, 0xffb19c771359a1b2, 0x9c0e8404bd0430e5, 0x3839b5b65c3cd848, 0x00000a6b32d2b48e}, {0x2c1bc37438e53431, 0xdef9b132930390c6, 0x4c95e6a74a997e75, 0xfa9c74d0e44c4476, 0xe05ec3d08b6ffce9, 0x630a5c4c63be60c0, 0x000065d841e1c3a7}}, {1}};
    five_basis.Q = (ec_point_t){{{0xd2221aded7a4cea8, 0x563a24bf3290766c, 0x99f459ff7ecb9995, 0x2d9476c453e3533d, 0x55323e929d876000, 0x229ea3b7155ee93d, 0x0000690822912c3c}, {0xd4b60aa2bdb7c1d4, 0xbd69a2702760b7d7, 0xdf617f7ddeb27a9e, 0x4a9562f07010f402, 0x0d6acce79c00c8a6, 0x258e7e125d7983ab, 0x000036fe85655e85}}, {1}};
    five_basis.PmQ = (ec_point_t){{{0x8181e636d919c900, 0xab296d10f2ce4855, 0x7309ebfd1a12e3f7, 0xbf41705725367344, 0x8a510665947547ce, 0x035116a08999bf43, 0x0000453f10e6f8f5}, {0xee587fe43c4082dc, 0x0935443949ece58f, 0xeb3082e610c83b16, 0xa2194820d49bf4cf, 0x113941e341617040, 0x4682ff797f946ee9, 0x00004f082479b59a}}, {1}};

    ibz_t a, b, c, d;
    digit_t x[NWORDS_ORDER] = {0}, y[NWORDS_ORDER] = {0};
    ibz_init(&a); ibz_init(&b);
    ibz_init(&c); ibz_init(&d);
    ibz_rand_interval(&a, &ibz_const_zero, &TORSION_PLUS_5POWER);
    ibz_rand_interval(&b, &ibz_const_zero, &TORSION_PLUS_5POWER);

    ibz_to_digits(x, &a);
    ibz_to_digits(y, &b);

    ibz_copy_digits(&a, x, NWORDS_ORDER);
    ibz_copy_digits(&b, y, NWORDS_ORDER);

    ec_point_t R = {0}, S = {0};
    digit_t z[NWORDS_ORDER] = {0}, w[NWORDS_ORDER] = {0};
    ec_biscalar_mul(&R, &E0, x, y, &five_basis);

    ec_dlog_5(z ,w, &five_basis, &R, &E0);
    ec_biscalar_mul(&S, &E0, z, w, &five_basis);
    
    if (!is_point_equal(&S, &R)) {
        point_print("R : ", R);
        point_print("five_basis.P : ", five_basis.P);
        point_print("five_basis.Q : ", five_basis.Q);
        ibz_copy_digits(&c, z, NWORDS_ORDER);
        ibz_copy_digits(&d, w, NWORDS_ORDER);
        curve_print("E0 : ", E0);
        gmp_printf("x, y: %Zx, %Zx\n", a, b);
        gmp_printf("z, w: %Zx, %Zx\n", c, d);

        printf("dlog_5_test failed\n");
        return 1;
    }

    ibz_finalize(&a); ibz_finalize(&b);
    return 0;
}

int main (int argc, char* argv[]) {
    int test_num = 1;
    if (argc > 1) {
        test_num = atoi(argv[1]);
    }
    // for (int i = 0; i < test_num; i++){
    //     if (dlog_2_test())
    //         return 1;
    // }

    // for (int i = 0; i < test_num; i++){
    //     if (dlog_3_test())
    //         return 1;
    // }

    for (int i = 0; i < test_num; i++){
        if (dlog_5_test())
            return 1;
    }
    return 0;
}