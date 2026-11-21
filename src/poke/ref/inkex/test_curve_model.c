// SPDX-License-Identifier: Apache-2.0
/* Finite-model invariance and signed point-map checks. These tests exercise
 * the exceptional automorphism groups as well as a generic j-invariant. */
#include <inke_curve_model.h>
#include <curve_extras.h>
#include <encoded_sizes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c); exit(1); } } while (0)

static int same_basis(const ec_basis_t *a,const ec_basis_t *b)
{
    return ec_is_equal(&a->P,&b->P) && ec_is_equal(&a->Q,&b->Q) &&
        ec_is_equal(&a->PmQ,&b->PmQ);
}

static void check_normalized(const ec_basis_t *basis)
{
    const ec_point_t *points[3]={&basis->P,&basis->Q,&basis->PmQ};
    for (int i=0;i<3;++i) {
        if (fp2_is_zero(&points[i]->z)) CHECK(fp2_is_one(&points[i]->x));
        else CHECK(fp2_is_one(&points[i]->z));
    }
}

static void check_on_curve(const jac_point_t *P,const ec_curve_t *curve)
{
    CHECK(fp2_is_one(&curve->C));
    if (fp2_is_zero(&P->z)) return;
    fp2_t x2,z2,z4,right,term,left;
    fp2_sqr(&x2,&P->x); fp2_sqr(&z2,&P->z); fp2_sqr(&z4,&z2);
    fp2_mul(&right,&x2,&P->x);
    fp2_mul(&term,&curve->A,&x2); fp2_mul(&term,&term,&z2);
    fp2_add(&right,&right,&term); fp2_mul(&term,&P->x,&z4);
    fp2_add(&right,&right,&term); fp2_sqr(&left,&P->y);
    CHECK(fp2_is_equal(&left,&right));
}

static void find_point(jac_point_t *P,const ec_curve_t *curve,unsigned int start)
{
    for (unsigned int x=start;x<start+1000;++x) {
        fp2_t rhs,x2;
        fp2_set_small(&P->x,x); fp2_sqr(&x2,&P->x);
        fp2_mul(&rhs,&curve->A,&x2); fp2_mul(&x2,&x2,&P->x);
        fp2_add(&rhs,&rhs,&x2); fp2_add(&rhs,&rhs,&P->x);
        if (fp2_is_zero(&rhs) || !fp2_is_square(&rhs)) continue;
        P->y=rhs; fp2_sqrt(&P->y); fp2_set_one(&P->z);
        check_on_curve(P,curve); return;
    }
    CHECK(0);
}

static void basis_from_points(ec_basis_t *basis,const jac_point_t *P,
    const jac_point_t *Q,const ec_curve_t *curve)
{
    jac_point_t negative=*Q,difference;
    fp2_neg(&negative.y,&negative.y); ADD(&difference,P,&negative,curve);
    CHECK(!fp2_is_zero(&difference.z));
    jac_to_xz(&basis->P,P); jac_to_xz(&basis->Q,Q);
    jac_to_xz(&basis->PmQ,&difference);
}

/* An independent signed lift of x'=(x-r)/t: the Jacobian coordinates are
 * (X-r Z^2:Y:u Z), where u^2=t. The negative u gives the other y sign. */
static void map_jacobian(jac_point_t *out,const jac_point_t *P,
    const inke_model_map *map)
{
    fp2_t root=map->xmap.D,z2,offset,check;
    CHECK(fp2_is_one(&map->xmap.Nx)); CHECK(fp2_is_square(&root));
    fp2_sqrt(&root); fp2_sqr(&check,&root);
    CHECK(fp2_is_equal(&check,&map->xmap.D));
    *out=*P; fp2_sqr(&z2,&P->z); fp2_mul(&offset,&map->xmap.Nz,&z2);
    fp2_sub(&out->x,&P->x,&offset); fp2_mul(&out->z,&root,&P->z);
}

static void check_polynomial(const ec_curve_t *source,const inke_model_map *map)
{
    const fp2_t *r=&map->xmap.Nz,*t=&map->xmap.D;
    fp2_t r2,r3,Ar,root_value,derivative,t2,one,coefficient;
    CHECK(fp2_is_one(&map->xmap.Nx)); CHECK(fp2_is_one(&map->curve.C));
    CHECK(!map->curve.is_A24_computed_and_normalized);
    CHECK(fp2_is_square(t)); CHECK(!fp2_is_zero(t));
    fp2_sqr(&r2,r); fp2_mul(&r3,&r2,r);
    fp2_mul(&root_value,&source->A,&r2); fp2_add(&root_value,&root_value,&r3);
    fp2_add(&root_value,&root_value,r); CHECK(fp2_is_zero(&root_value));
    fp2_mul(&Ar,&source->A,r); fp2_add(&derivative,&r2,&r2);
    fp2_add(&derivative,&derivative,&r2); fp2_add(&derivative,&derivative,&Ar);
    fp2_add(&derivative,&derivative,&Ar); fp2_set_one(&one);
    fp2_add(&derivative,&derivative,&one); fp2_sqr(&t2,t);
    CHECK(fp2_is_equal(&derivative,&t2));
    fp2_add(&coefficient,r,r); fp2_add(&coefficient,&coefficient,r);
    fp2_add(&coefficient,&coefficient,&source->A);
    fp2_mul(&t2,&map->curve.A,t); CHECK(fp2_is_equal(&coefficient,&t2));
    fp2_t js,jt; ec_j_inv(&js,source); ec_j_inv(&jt,&map->curve);
    CHECK(fp2_is_equal(&js,&jt));
}

static void check_fixture(const ec_curve_t *source,unsigned int automorphism_count,
    const char *name)
{
    inke_model_map all[INKE_MODEL_MAX_MAPS],canonical[INKE_MODEL_MAX_MAPS];
    unsigned int total,count;
    CHECK(inke_model_enumerate(all,&total,source));
    CHECK(inke_model_candidates(canonical,&count,source));
    CHECK(count==automorphism_count);
    unsigned char least[FP2_ENCODED_BYTES],encoded[FP2_ENCODED_BYTES];
    fp2_encode(least,&canonical[0].curve.A);
    jac_point_t P,Q,twice;
    find_point(&P,source,1); find_point(&Q,source,1001); DBL(&twice,&P,source);
    ec_basis_t basis,canonical_images[INKE_MODEL_MAX_MAPS];
    basis_from_points(&basis,&P,&Q,source);
    for (unsigned int k=0;k<count;++k) {
        CHECK(fp2_is_equal(&canonical[k].curve.A,&canonical[0].curve.A));
        CHECK(inke_model_apply_basis(&canonical_images[k],&basis,&canonical[k]));
        check_normalized(&canonical_images[k]);
    }
    for (unsigned int k=0;k<total;++k) {
        check_polynomial(source,&all[k]);
        fp2_encode(encoded,&all[k].curve.A); CHECK(memcmp(least,encoded,sizeof(least))<=0);
        for (unsigned int l=0;l<k;++l)
            CHECK(!fp2_is_equal(&all[l].xmap.Nz,&all[k].xmap.Nz) ||
                  !fp2_is_equal(&all[l].xmap.D,&all[k].xmap.D));
        jac_point_t imageP,imageQ,imageTwice,doubleImage;
        map_jacobian(&imageP,&P,&all[k]); map_jacobian(&imageQ,&Q,&all[k]);
        map_jacobian(&imageTwice,&twice,&all[k]);
        check_on_curve(&imageP,&all[k].curve); check_on_curve(&imageQ,&all[k].curve);
        DBL(&doubleImage,&imageP,&all[k].curve);
        CHECK(is_jac_equal(&imageTwice,&doubleImage));
        ec_basis_t image,from_jacobian,inplace=basis;
        CHECK(inke_model_apply_basis(&image,&basis,&all[k]));
        CHECK(inke_model_apply_basis(&inplace,&inplace,&all[k]));
        basis_from_points(&from_jacobian,&imageP,&imageQ,&all[k].curve);
        CHECK(same_basis(&image,&from_jacobian)); CHECK(same_basis(&image,&inplace));
        check_normalized(&image);
        /* The inverse gives exactly the original normalized affine x values. */
        inke_model_map inverse=all[k]; inverse.curve=*source;
        inverse.xmap.Nx=all[k].xmap.D;
        fp2_neg(&inverse.xmap.Nz,&all[k].xmap.Nz); fp2_set_one(&inverse.xmap.D);
        ec_basis_t recovered,normalized;
        CHECK(inke_model_apply_basis(&recovered,&image,&inverse));
        inke_model_map identity={0}; identity.curve=*source;
        fp2_set_one(&identity.xmap.Nx); fp2_set_one(&identity.xmap.D);
        CHECK(inke_model_apply_basis(&normalized,&basis,&identity));
        CHECK(same_basis(&recovered,&normalized)); check_normalized(&recovered);
        const ec_point_t *a[3]={&recovered.P,&recovered.Q,&recovered.PmQ};
        const ec_point_t *b[3]={&normalized.P,&normalized.Q,&normalized.PmQ};
        for (int l=0;l<3;++l) {
            CHECK(fp2_is_equal(&a[l]->x,&b[l]->x)); CHECK(fp2_is_equal(&a[l]->z,&b[l]->z));
        }
        /* Changing the source model permutes the complete set of canonical
         * basis images, with multiplicities. This tests the uniform map law. */
        inke_model_map again[INKE_MODEL_MAX_MAPS]; unsigned int again_count;
        CHECK(inke_model_candidates(again,&again_count,&all[k].curve));
        CHECK(again_count==count); unsigned int seen[INKE_MODEL_MAX_MAPS]={0};
        for (unsigned int l=0;l<again_count;++l) {
            CHECK(fp2_is_equal(&again[l].curve.A,&canonical[0].curve.A));
            ec_basis_t combined; CHECK(inke_model_apply_basis(&combined,&image,&again[l]));
            int found=0;
            for (unsigned int m=0;m<count;++m) if (!seen[m] && same_basis(&combined,&canonical_images[m])) {
                seen[m]=1; found=1; break;
            }
            CHECK(found);
        }
    }
    /* Arbitrary projective representatives must not affect the model or
     * published affine point coordinates. Include a populated stale cache. */
    ec_curve_t rescaled_curve=*source; fp2_t scale; fp2_set_small(&scale,7);
    fp2_mul(&rescaled_curve.A,&rescaled_curve.A,&scale); fp2_mul(&rescaled_curve.C,&rescaled_curve.C,&scale);
    rescaled_curve.is_A24_computed_and_normalized=true;
    fp2_set_small(&rescaled_curve.A24.x,83); fp2_set_small(&rescaled_curve.A24.z,41);
    inke_model_map from_scaled[INKE_MODEL_MAX_MAPS]; unsigned int scaled_count;
    CHECK(inke_model_candidates(from_scaled,&scaled_count,&rescaled_curve)); CHECK(scaled_count==count);
    for (unsigned int k=0;k<count;++k) {
        CHECK(fp2_is_equal(&from_scaled[k].curve.A,&canonical[0].curve.A));
        CHECK(!from_scaled[k].curve.is_A24_computed_and_normalized);
    }
    ec_basis_t scaled_basis=basis,scaled_image;
    ec_point_t *points[3]={&scaled_basis.P,&scaled_basis.Q,&scaled_basis.PmQ};
    for (int k=0;k<3;++k) {
        fp2_set_small(&scale,11+k); fp2_mul(&points[k]->x,&points[k]->x,&scale);
        fp2_mul(&points[k]->z,&points[k]->z,&scale);
    }
    CHECK(inke_model_apply_basis(&scaled_image,&scaled_basis,&canonical[0]));
    CHECK(same_basis(&scaled_image,&canonical_images[0])); check_normalized(&scaled_image);
    printf("%s: %u complete maps, %u canonical maps; model and point invariance passed.\n",name,total,count);
}

static void check_invalid_and_infinity(void)
{
    ec_curve_t invalid; memset(&invalid,0,sizeof(invalid)); ec_curve_init(&invalid);
    inke_model_map output[INKE_MODEL_MAX_MAPS],before[INKE_MODEL_MAX_MAPS];
    memset(output,0xa5,sizeof(output)); memcpy(before,output,sizeof(output)); unsigned int count=91;
    fp2_set_small(&invalid.A,2);
    CHECK(!inke_model_candidates(output,&count,&invalid));
    CHECK(count==91 && memcmp(output,before,sizeof(output))==0);
    fp2_set_zero(&invalid.C);
    CHECK(!inke_model_enumerate(output,&count,&invalid));
    CHECK(count==91 && memcmp(output,before,sizeof(output))==0);
    ec_curve_init(&invalid); CHECK(inke_model_candidates(output,&count,&invalid));
    ec_basis_t basis,result,saved; memset(&basis,0,sizeof(basis));
    ec_set_zero(&basis.P); ec_set_zero(&basis.Q); ec_set_zero(&basis.PmQ);
    CHECK(inke_model_apply_basis(&result,&basis,&output[0])); check_normalized(&result);
    CHECK(ec_is_zero(&result.P) && ec_is_zero(&result.Q) && ec_is_zero(&result.PmQ));
    saved=result; fp2_set_zero(&basis.Q.x); /* Invalid (0:0), after a valid P. */
    CHECK(!inke_model_apply_basis(&result,&basis,&output[0]));
    CHECK(memcmp(&result,&saved,sizeof(result))==0);
}

int main(void)
{
    ec_curve_t fixture; memset(&fixture,0,sizeof(fixture)); ec_curve_init(&fixture);
    check_fixture(&fixture,2,"j=1728");
    fp2_set_one(&fixture.A); check_fixture(&fixture,1,"generic j");
    fp2_set_small(&fixture.A,3); CHECK(fp2_is_square(&fixture.A)); fp2_sqrt(&fixture.A);
    fp2_t j; ec_j_inv(&j,&fixture); CHECK(fp2_is_zero(&j));
    check_fixture(&fixture,3,"j=0");
    check_invalid_and_infinity();
    puts("Canonical Montgomery models, automorphism multiplicities, signed maps and atomic failures passed.");
    return 0;
}
