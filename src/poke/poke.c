#include <hd.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <gmp.h>
#include <intbig.h>

int test() {
    // ibz_t q, low, high;
    // ibz_init(&q);
    // ibz_init(&low);
    // ibz_init(&high);
    // ibz_pow(&high, &ibz_const_two, TORSION_PLUS_EVEN_POWER);
    // ibz_rand_interval(&q, &low, &high);
    // ibz_printf("q = %Zd\n", &q);
    // ibz_finalize(&q);
    // ibz_finalize(&low);
    // ibz_finalize(&high);
    ec_curve_t curve = CURVE_E0;
    fp2_t j_inv;
    ec_j_inv(&j_inv, &curve);
    fp2_print("j_invariant",&j_inv);
    
    return 0; 
}

int main() {
    int res = 1;

    test();

    return res;
}