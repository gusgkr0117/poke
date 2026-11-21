// SPDX-License-Identifier: Apache-2.0
/* Independent checks of the new parameter constants at their full orders.
 * In particular, comparisons modulo 2^f do not discard the two high bits. */
#include <inke_conversion.h>
#include <endomorphism_action.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <biextension.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c); exit(1); } } while (0)
extern const uint64_t p[NWORDS_FIELD];
void weil_odd(fp2_t *r,const digit_t *n,int n_bitlen,ec_point_t *P,
    ec_point_t *Q,ec_point_t *PQ,const ec_point_t *A24);

static void element(quat_alg_elem_t *q,long a,long b,long c,long d,unsigned long denominator)
{
    mpz_set_si(q->coord[0],a); mpz_set_si(q->coord[1],b);
    mpz_set_si(q->coord[2],c); mpz_set_si(q->coord[3],d);
    mpz_set_ui(q->denom,denominator);
}
static int matrix_equal(const ibz_mat_2x2_t *a,const ibz_mat_2x2_t *b,const ibz_t *m)
{
    ibz_t difference; ibz_init(&difference); int equal=1;
    for (int i=0;i<2;++i) for (int j=0;j<2;++j) {
        mpz_sub(difference,(*a)[i][j],(*b)[i][j]);
        equal=equal && mpz_divisible_p(difference,*m);
    }
    ibz_finalize(&difference); return equal;
}
static void matrix_for_element(ibz_mat_2x2_t *out,const quat_alg_elem_t *q,const ibz_t *m)
{
    ibz_vec_4_t coefficients; ibz_vec_4_init(&coefficients);
    CHECK(quat_lattice_contains(&coefficients,&MAXORD_O0,q,&QUATALG_PINFTY));
    const ibz_mat_2x2_t *generators[3]={&ACTION_GEN2,&ACTION_GEN3,&ACTION_GEN4};
    for (int i=0;i<2;++i) for (int j=0;j<2;++j) {
        mpz_set_ui((*out)[i][j],0);
        if (i==j) mpz_set((*out)[i][j],coefficients[0]);
        for (int k=0;k<3;++k) mpz_addmul((*out)[i][j],coefficients[k+1],(*generators[k])[i][j]);
        mpz_mod((*out)[i][j],(*out)[i][j],*m);
    }
    ibz_vec_4_finalize(&coefficients);
}
static void check_norm_trace(const ibz_mat_2x2_t *matrix,const quat_alg_elem_t *q,const ibz_t *modulus)
{
    ibq_t rational; ibq_init(&rational);
    ibz_t reference,result; ibz_init(&reference); ibz_init(&result);
    quat_alg_norm(&rational,q,&QUATALG_PINFTY); CHECK(ibq_to_ibz(&reference,&rational));
    mpz_mul(result,(*matrix)[0][0],(*matrix)[1][1]);
    mpz_submul(result,(*matrix)[0][1],(*matrix)[1][0]);
    mpz_sub(result,result,reference); CHECK(mpz_divisible_p(result,*modulus));
    quat_alg_trace(&rational,q); CHECK(ibq_to_ibz(&reference,&rational));
    mpz_add(result,(*matrix)[0][0],(*matrix)[1][1]);
    mpz_sub(result,result,reference); CHECK(mpz_divisible_p(result,*modulus));
    ibz_finalize(&result); ibz_finalize(&reference); ibq_finalize(&rational);
}
static void check_ring(const ibz_t *modulus)
{
    quat_alg_elem_t basis[4],product,q; quat_alg_elem_init(&product); quat_alg_elem_init(&q);
    for (int i=0;i<4;++i) quat_alg_elem_init(&basis[i]);
    element(&basis[0],1,0,0,0,1); element(&basis[1],0,1,0,0,1);
    element(&basis[2],0,1,1,0,2); element(&basis[3],1,0,0,1,2);
    ibz_mat_2x2_t identity,left,right,actual,expected;
    ibz_mat_2x2_init(&identity); ibz_mat_2x2_init(&left); ibz_mat_2x2_init(&right);
    ibz_mat_2x2_init(&actual); ibz_mat_2x2_init(&expected);
    ibz_mat_2x2_set(&identity,1,0,0,1);
    const ibz_mat_2x2_t *matrices[4]={&identity,&ACTION_GEN2,&ACTION_GEN3,&ACTION_GEN4};
    CHECK(matrix_equal(&ACTION_I,&ACTION_GEN2,modulus));
    for (int i=0;i<4;++i) {
        matrix_for_element(&left,&basis[i],modulus);
        CHECK(matrix_equal(&left,matrices[i],modulus));
        check_norm_trace(matrices[i],&basis[i],modulus);
        for (int j=0;j<4;++j) {
            quat_alg_mul(&product,&basis[i],&basis[j],&QUATALG_PINFTY);
            matrix_for_element(&expected,&product,modulus);
            ibz_2x2_mul_mod(&actual,matrices[i],matrices[j],modulus);
            CHECK(matrix_equal(&actual,&expected,modulus));
        }
    }
    element(&q,0,0,1,0,1); matrix_for_element(&actual,&q,modulus);
    CHECK(matrix_equal(&actual,&ACTION_J,modulus)); check_norm_trace(&ACTION_J,&q,modulus);
    element(&q,0,0,0,1,1); matrix_for_element(&actual,&q,modulus);
    CHECK(matrix_equal(&actual,&ACTION_K,modulus)); check_norm_trace(&ACTION_K,&q,modulus);
    ibz_2x2_mul_mod(&actual,&ACTION_I,&ACTION_J,modulus);
    CHECK(matrix_equal(&actual,&ACTION_K,modulus));
    /* 2 G3 = I+J and G4 = 1 + I G3 hold at the full modulus. */
    ibz_mat_2x2_add(&actual,&ACTION_GEN3,&ACTION_GEN3);
    ibz_mat_2x2_add(&expected,&ACTION_I,&ACTION_J);
    CHECK(matrix_equal(&actual,&expected,modulus));
    ibz_2x2_mul_mod(&actual,&ACTION_I,&ACTION_GEN3,modulus);
    ibz_mat_2x2_add(&expected,&actual,&identity);
    CHECK(matrix_equal(&expected,&ACTION_GEN4,modulus));
    ibz_mat_2x2_add(&actual,&ACTION_GEN4,&ACTION_GEN4);
    ibz_mat_2x2_add(&expected,&identity,&ACTION_K);
    CHECK(matrix_equal(&actual,&expected,modulus));
    /* Integral combinations exercise the complete maximal-order basis. */
    for (long t=-7;t<=7;++t) {
        element(&q,2*t+3,t*t+2*t-1,t*t-1,3,2);
        matrix_for_element(&actual,&q,modulus);
        check_norm_trace(&actual,&q,modulus);
    }
    ibz_mat_2x2_finalize(&expected); ibz_mat_2x2_finalize(&actual);
    ibz_mat_2x2_finalize(&right); ibz_mat_2x2_finalize(&left); ibz_mat_2x2_finalize(&identity);
    for (int i=0;i<4;++i) quat_alg_elem_finalize(&basis[i]);
    quat_alg_elem_finalize(&q); quat_alg_elem_finalize(&product);
}

static void matrix_image(jac_point_t *out,const jac_point_t points[2],
    const ibz_mat_2x2_t *matrix,int column,const ibz_t *modulus,const ec_curve_t *curve)
{
    digit_t coefficients[2][NWORDS_ORDER]={{0}};
    ibz_t reduced; ibz_init(&reduced);
    for (int i=0;i<2;++i) {
        mpz_mod(reduced,(*matrix)[i][column],*modulus);
        ibz_to_digits(coefficients[i],&reduced);
    }
    DBLMUL_generic(out,&points[0],coefficients[0],&points[1],coefficients[1],curve,NWORDS_ORDER);
    ibz_finalize(&reduced);
}
static void direct_i(jac_point_t *out,const jac_point_t *P)
{
    fp2_t imaginary; fp2_set_zero(&imaginary); fp_set_one(&imaginary.im);
    *out=*P; fp2_neg(&out->x,&P->x); fp2_mul(&out->y,&P->y,&imaginary);
}
static void direct_j(jac_point_t *out,const jac_point_t *P)
{
    fp2_conj(&out->x,&P->x); fp2_conj(&out->y,&P->y); fp2_conj(&out->z,&P->z);
}
static void to_basis(ec_basis_t *out,const jac_point_t points[2],const ec_curve_t *curve)
{
    jac_point_t negative=points[1],difference;
    fp2_neg(&negative.y,&negative.y); ADD(&difference,&points[0],&negative,curve);
    jac_to_xz(&out->P,&points[0]); jac_to_xz(&out->Q,&points[1]); jac_to_xz(&out->PmQ,&difference);
}
static int same_basis(const ec_basis_t *a,const ec_basis_t *b)
{ return ec_is_equal(&a->P,&b->P) && ec_is_equal(&a->Q,&b->Q) && ec_is_equal(&a->PmQ,&b->PmQ); }

static void check_points(const ec_basis_t *constant_basis,const ibz_t *modulus,int ell,int exponent)
{
    ec_curve_t curve=CURVE_E0; ec_basis_t basis=*constant_basis,reconstructed,converted;
    const ec_point_t *points[3]={&basis.P,&basis.Q,&basis.PmQ};
    for (int k=0;k<3;++k) {
        CHECK(ec_is_on_curve(&curve,points[k]));
        CHECK(ell==2?test_point_order_twof(points[k],&curve,exponent):
                      test_point_order_threef(points[k],&curve,exponent));
    }
    jac_point_t lifted[2],direct[2],mapped[2];
    lift_basis(&lifted[0],&lifted[1],&basis,&curve);
    to_basis(&reconstructed,lifted,&curve); CHECK(same_basis(&basis,&reconstructed));
    ec_point_t ac={curve.A,curve.C},a24; A24_from_AC(&a24,&ac);
    fp2_t pairing,power,tmp;
    if (ell==2) weil(&pairing,exponent,&basis.P,&basis.Q,&basis.PmQ,&a24);
    else weil_odd(&pairing,THREEpF,THREEpF_bitlen,&basis.P,&basis.Q,&basis.PmQ,&a24);
    CHECK(!fp2_is_zero(&pairing)); power=pairing;
    for (int k=1;k<exponent;++k) {
        fp2_sqr(&tmp,&power);
        if (ell==2) power=tmp; else fp2_mul(&power,&tmp,&power);
    }
    CHECK(!fp2_is_one(&power)); fp2_sqr(&tmp,&power);
    if (ell==2) power=tmp; else fp2_mul(&power,&tmp,&power);
    CHECK(fp2_is_one(&power));
    quat_alg_elem_t q; quat_alg_elem_init(&q);
    const ibz_mat_2x2_t *matrices[4]={&ACTION_I,&ACTION_J,&ACTION_GEN3,&ACTION_GEN4};
    for (int action=0;action<4;++action) {
        if (action==0) element(&q,0,1,0,0,1);
        if (action==1) element(&q,0,0,1,0,1);
        if (action==2) element(&q,0,1,1,0,2);
        if (action==3) element(&q,1,0,0,1,2);
        for (int k=0;k<2;++k) {
            matrix_image(&mapped[k],lifted,matrices[action],k,modulus,&curve);
            if (action<2) {
                if (action==0) direct_i(&direct[k],&lifted[k]); else direct_j(&direct[k],&lifted[k]);
                CHECK(is_jac_equal(&mapped[k],&direct[k]));
            } else if (ell==3) {
                jac_point_t other;
                if (action==2) { direct_i(&direct[k],&lifted[k]); direct_j(&other,&lifted[k]); }
                else { direct_j(&other,&lifted[k]); direct_i(&direct[k],&other); other=lifted[k]; }
                ADD(&direct[k],&direct[k],&other,&curve);
                ibz_t half; ibz_init(&half); mpz_add_ui(half,*modulus,1); mpz_divexact_ui(half,half,2);
                digit_t scalar[NWORDS_ORDER]={0}; ibz_to_digits(scalar,&half);
                MUL_generic(&direct[k],&direct[k],scalar,NWORDS_ORDER,&curve); ibz_finalize(&half);
                CHECK(is_jac_equal(&mapped[k],&direct[k]));
            }
        }
        to_basis(&reconstructed,mapped,&curve);
        CHECK(inke_endomorphism_basis(&converted,&basis,&q,modulus));
        CHECK(same_basis(&reconstructed,&converted));
    }
    quat_alg_elem_finalize(&q);
    printf("Full %d^%d basis, primitive pairing, signed point actions, and P-Q checked.\n",ell,exponent);
}

int main(void)
{
    CHECK(fp2_is_zero(&CURVE_E0.A)); CHECK(fp2_is_one(&CURVE_E0.C));
    ibz_t modulus,prime,three; ibz_init(&modulus); ibz_init(&prime); ibz_init(&three);
    mpz_import(prime,NWORDS_FIELD,-1,sizeof(uint64_t),0,0,p);
    CHECK(mpz_cmp(prime,QUATALG_PINFTY.p)==0);
    CHECK(TORSION_PLUS_EVEN_POWER >= INKE_PUBLIC_TWO_EXPONENT + 2);
    mpz_set_ui(modulus,1); mpz_mul_2exp(modulus,modulus,TORSION_PLUS_EVEN_POWER);
    CHECK(mpz_cmp(modulus,TORSION_PLUS_2POWER)==0);
    mpz_ui_pow_ui(three,3,POWER_OF_3); CHECK(mpz_cmp(three,TORSION_PLUS_3POWER)==0);
    mpz_mul(modulus,modulus,three);
    check_ring(&modulus);
    puts("Full-modulus order multiplication, trace, and determinant checks passed."); fflush(stdout);
    check_points(&BASIS_EVEN,&TORSION_PLUS_2POWER,2,TORSION_PLUS_EVEN_POWER);
    check_points(&BASIS_THREE,&TORSION_PLUS_3POWER,3,POWER_OF_3);
    ibz_finalize(&three); ibz_finalize(&prime); ibz_finalize(&modulus);
    puts("Manuscript precomputation arithmetic and full-order action checks passed.");
    return 0;
}
