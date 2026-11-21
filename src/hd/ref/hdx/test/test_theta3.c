// SPDX-License-Identifier: Apache-2.0
#include "theta_isogeny3.h"
#include "theta_isogenies.h"
#include <curve_extras.h>
#include <endomorphism_action.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c); exit(1); \
} } while (0)

static void glue_one(theta_point_t *out, const ec_point_t *P, int side,
    const theta_gluing_t *glue)
{
    theta_couple_point_t pair;
    ec_set_zero(&pair.P1); ec_set_zero(&pair.P2);
    if (side==0) pair.P1=*P; else pair.P2=*P;
    gluing_eval_point_special_case(out,&pair,glue);
}

static int full_order(theta_point_t P, theta_structure_t A, int f)
{
    double_iter(&P,&A,&P,f-1);
    if (inke_theta_equal(&P,&A.null_point)) return 0;
    double_point(&P,&A,&P);
    return inke_theta_equal(&P,&A.null_point);
}

static void check_chain(theta_structure_t A, const theta_gluing_t *glue,
    theta_couple_curve_t *E, int f)
{
    /* Three successive (3,3) steps, from two independent 27-torsion
     * generators on different elliptic factors. All four nonzero kernel
     * classes are retained so a different generator pair can use the same
     * subgroup if the first coefficient chart vanishes. */
    ec_point_t P=BASIS_THREE.P,Q=BASIS_THREE.Q;
    digit_t three_scalar[NWORDS_ORDER]={0}; three_scalar[0]=3;
    for (int k=3;k<(int)TORSION_PLUS_ODD_POWERS[0];++k) {
        ec_mul(&P,&E->E1,three_scalar,&P);
        ec_mul(&Q,&E->E2,three_scalar,&Q);
    }
    theta_point_t gen[4],points[6];
    glue_one(&gen[0],&P,0,glue); glue_one(&gen[1],&Q,1,glue);
    theta_couple_jac_point_t both;
    lift_point(&both.P1,&P,&E->E1); lift_point(&both.P2,&Q,&E->E2);
    gluing_eval_point_no_help(&gen[2],&both,E,glue);
    fp2_neg(&both.P2.y,&both.P2.y);
    gluing_eval_point_no_help(&gen[3],&both,E,glue);
    const ec_point_t *two[3]={&BASIS_EVEN.P,&BASIS_EVEN.Q,&BASIS_EVEN.PmQ};
    for (int side=0;side<2;++side) for (int k=0;k<3;++k)
        glue_one(&points[3*side+k],two[k],side,glue);
    for (int remaining=3;remaining>0;--remaining) {
        theta_point_t kernel[4];
        for (int k=0;k<4;++k) {
            kernel[k]=gen[k];
            for (int t=1;t<remaining;++t)
                CHECK(inke_theta_triple(&kernel[k],&A,&kernel[k]));
        }
        inke_theta_isogeny3 phi;
        int found=0;
        for (int k=0;k<4 && !found;++k) for (int l=k+1;l<4 && !found;++l)
            found=inke_theta_isogeny3_compute(&phi,&A,&kernel[k],&kernel[l]);
        CHECK(found);
        for (int k=0;k<4;++k) {
            theta_point_t image;
            CHECK(inke_theta_isogeny3_eval(&image,&phi,&kernel[k]));
            CHECK(inke_theta_equal(&image,&phi.codomain.null_point));
            CHECK(inke_theta_isogeny3_eval(&gen[k],&phi,&gen[k]));
            CHECK(inke_theta_equal(&gen[k],&phi.codomain.null_point)==(remaining==1));
        }
        for (int k=0;k<6;++k) {
            theta_point_t twice,left,right;
            double_point(&twice,&A,&points[k]);
            CHECK(inke_theta_isogeny3_eval(&left,&phi,&twice));
            CHECK(inke_theta_isogeny3_eval(&points[k],&phi,&points[k]));
            double_point(&right,&phi.codomain,&points[k]);
            CHECK(inke_theta_equal(&left,&right));
            CHECK(full_order(points[k],phi.codomain,f));
        }
        A=phi.codomain;
    }
}

int main(void)
{
    theta_couple_curve_t E={CURVE_E0,CURVE_E0};
    ec_basis_t eight=BASIS_EVEN;
    int f=(int)TORSION_PLUS_EVEN_POWER;
    ec_dbl_iter(&eight.P,f-3,&E.E1,&eight.P);
    ec_dbl_iter(&eight.Q,f-3,&E.E1,&eight.Q);
    ec_dbl_iter(&eight.PmQ,f-3,&E.E1,&eight.PmQ);
    jac_point_t basis[3];
    lift_basis(&basis[0],&basis[1],&eight,&E.E1);
    ADD(&basis[2],&basis[0],&basis[1],&E.E1);
    ec_point_t three=BASIS_THREE.P, three2=BASIS_THREE.Q;
    digit_t scalar[NWORDS_ORDER]={0}; scalar[0]=3;
    for (int i=1;i<(int)TORSION_PLUS_ODD_POWERS[0];++i) {
        ec_mul(&three,&E.E1,scalar,&three);
        ec_mul(&three2,&E.E2,scalar,&three2);
    }
    CHECK(test_point_order_threef(&three,&E.E1,1));
    jac_point_t three_jac,three_jac2;
    lift_point(&three_jac,&three,&E.E1);
    lift_point(&three_jac2,&three2,&E.E2);

    unsigned int passed=0,unsupported=0;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) if (i!=j) {
        theta_couple_jac_point_t T1={basis[0],basis[i]},T2={basis[1],basis[j]};
        theta_gluing_t glue;
        memset(&glue,0,sizeof(glue));
        gluing_comput(&glue,&E,&T1,&T2);
        theta_structure_t A={0}; A.null_point=glue.codomain;
        theta_point_t gen[4];
        glue_one(&gen[0],&three,0,&glue);
        glue_one(&gen[1],&three2,1,&glue);
        theta_couple_jac_point_t both={three_jac,three_jac2};
        gluing_eval_point_no_help(&gen[2],&both,&E,&glue);
        fp2_neg(&both.P2.y,&both.P2.y);
        gluing_eval_point_no_help(&gen[3],&both,&E,&glue);
        inke_theta_isogeny3 phi;
        int found=0;
        for (int k=0;k<4 && !found;++k) for (int l=k+1;l<4 && !found;++l) {
            found=inke_theta_isogeny3_compute(&phi,&A,&gen[k],&gen[l]);
        }
        if (!found) {
            ++unsupported; continue;
        }
        theta_point_t image,doubled,left,right;
        for (int k=0;k<4;++k) {
            CHECK(inke_theta_isogeny3_eval(&image,&phi,&gen[k]));
            CHECK(inke_theta_equal(&image,&phi.codomain.null_point));
        }
        /* Independent elliptic arithmetic supplies P and 2P before gluing.
         * The odd-degree map must preserve the exact 2^f order. */
        const ec_point_t *two[3]={&BASIS_EVEN.P,&BASIS_EVEN.Q,&BASIS_EVEN.PmQ};
        for (int side=0;side<2;++side) for (int k=0;k<3;++k) {
            theta_point_t P;
            glue_one(&P,two[k],side,&glue);
            CHECK(full_order(P,A,f));
            ec_point_t twice;
            ec_dbl(&twice,&E.E1,two[k]);
            glue_one(&doubled,&twice,side,&glue);
            double_point(&left,&A,&P);
            CHECK(inke_theta_equal(&left,&doubled));
            ec_mul(&twice,&E.E1,scalar,two[k]);
            glue_one(&doubled,&twice,side,&glue);
            CHECK(inke_theta_triple(&left,&A,&P));
            CHECK(inke_theta_equal(&left,&doubled));
            CHECK(inke_theta_isogeny3_eval(&image,&phi,&P));
            CHECK(full_order(image,phi.codomain,f));
            CHECK(inke_theta_isogeny3_eval(&left,&phi,&doubled));
            CHECK(inke_theta_triple(&right,&phi.codomain,&image));
            CHECK(inke_theta_equal(&left,&right));
            double_point(&doubled,&A,&P);
            CHECK(inke_theta_isogeny3_eval(&left,&phi,&doubled));
            double_point(&right,&phi.codomain,&image);
            CHECK(inke_theta_equal(&left,&right));
            /* Aliasing is required by the eventual chain evaluator. */
            right=P;
            CHECK(inke_theta_isogeny3_eval(&right,&phi,&right));
            CHECK(inke_theta_equal(&right,&image));
        }
        inke_theta_isogeny3 rejected;
        CHECK(!inke_theta_isogeny3_compute(&rejected,&A,&gen[0],&gen[0]));
        CHECK(!inke_theta_isogeny3_eval(&image,&rejected,&gen[0]));
        CHECK(!inke_theta_isogeny3_compute(&rejected,&A,&A.null_point,&gen[1]));
        check_chain(A,&glue,&E,f);
        ++passed;
    }
    CHECK(passed==4 && unsupported==2);
    printf("C (3,3) step and length-3 chain: %u fixtures passed, %u unsupported charts; full 2^%d torsion retained.\n",
        passed,unsupported,f);
    puts("Uses upstream field parameters; not a revised Keygen benchmark.");
    return 0;
}
