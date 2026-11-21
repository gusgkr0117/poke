// SPDX-License-Identifier: Apache-2.0
#include <inke_curve_model.h>
#include <curve_extras.h>
#include <encoded_sizes.h>
#include <string.h>

static int checked_square_root(fp2_t *root,const fp2_t *value)
{
    if (!fp2_is_square(value)) return 0;
    *root=*value; fp2_sqrt(root);
    fp2_t check; fp2_sqr(&check,root);
    return fp2_is_equal(&check,value)!=0;
}

int inke_model_enumerate(inke_model_map maps[INKE_MODEL_MAX_MAPS],
    unsigned int *count,const ec_curve_t *source)
{
    if (!maps || !count || !source || fp2_is_zero(&source->C)) return 0;
    inke_model_map result[INKE_MODEL_MAX_MAPS]; unsigned int used=0;
    memset(result,0,sizeof(result));
    fp2_t A=source->C,discriminant,four,root,roots[3];
    fp2_inv(&A); fp2_mul(&A,&source->A,&A);
    fp2_set_small(&four,4); fp2_sqr(&discriminant,&A);
    fp2_sub(&discriminant,&discriminant,&four);
    if (fp2_is_zero(&discriminant)) return 0;
    fp2_set_zero(&roots[0]); unsigned int root_count=1;
    if (checked_square_root(&root,&discriminant)) {
        fp2_sub(&roots[1],&root,&A); fp2_half(&roots[1],&roots[1]);
        fp2_neg(&roots[2],&root); fp2_sub(&roots[2],&roots[2],&A);
        fp2_half(&roots[2],&roots[2]); root_count=3;
    }
    for (unsigned int k=0;k<root_count;++k) {
        fp2_t r2,Ar,derivative,one,t;
        fp2_sqr(&r2,&roots[k]); fp2_mul(&Ar,&A,&roots[k]);
        fp2_add(&derivative,&r2,&r2); fp2_add(&derivative,&derivative,&r2);
        fp2_add(&derivative,&derivative,&Ar); fp2_add(&derivative,&derivative,&Ar);
        fp2_set_one(&one); fp2_add(&derivative,&derivative,&one);
        if (fp2_is_zero(&derivative)) return 0;
        if (!checked_square_root(&t,&derivative)) continue;
        for (unsigned int sign=0;sign<2;++sign) {
            if (sign) fp2_neg(&t,&t);
            /* x=(r+t*x'), y=u^3*y', u^2=t. Requiring t to be a
             * square excludes maps defined only over an extension. */
            if (!fp2_is_square(&t)) continue;
            if (used>=INKE_MODEL_MAX_MAPS) return 0;
            inke_model_map *candidate=&result[used++];
            ec_curve_init(&candidate->curve);
            fp2_add(&candidate->curve.A,&roots[k],&roots[k]);
            fp2_add(&candidate->curve.A,&candidate->curve.A,&roots[k]);
            fp2_add(&candidate->curve.A,&candidate->curve.A,&A);
            fp2_t inverse=t; fp2_inv(&inverse);
            fp2_mul(&candidate->curve.A,&candidate->curve.A,&inverse);
            fp2_set_one(&candidate->xmap.Nx);
            candidate->xmap.Nz=roots[k]; candidate->xmap.D=t;
        }
    }
    /* r=0,t=1 always supplies the identity map for a valid input. */
    if (!used) return 0;
    memcpy(maps,result,used*sizeof(*maps)); *count=used; return 1;
}

int inke_model_candidates(inke_model_map maps[INKE_MODEL_MAX_MAPS],
    unsigned int *count,const ec_curve_t *source)
{
    if (!maps || !count) return 0;
    inke_model_map all[INKE_MODEL_MAX_MAPS],result[INKE_MODEL_MAX_MAPS];
    unsigned int total,used=0;
    if (!inke_model_enumerate(all,&total,source)) return 0;
    unsigned char smallest[FP2_ENCODED_BYTES],encoding[FP2_ENCODED_BYTES];
    for (unsigned int i=0;i<total;++i) {
        fp2_encode(encoding,&all[i].curve.A);
        int comparison=used?memcmp(encoding,smallest,sizeof(smallest)):-1;
        if (comparison<0) { used=0; memcpy(smallest,encoding,sizeof(smallest)); }
        if (comparison<=0) result[used++]=all[i];
    }
    if (!used) return 0;
    memcpy(maps,result,used*sizeof(*maps)); *count=used; return 1;
}

int inke_model_apply_basis(ec_basis_t *out,const ec_basis_t *basis,
    const inke_model_map *map)
{
    if (!out || !basis || !map || fp2_is_zero(&map->xmap.D)) return 0;
    ec_basis_t result=*basis;
    ec_point_t *points[3]={&result.P,&result.Q,&result.PmQ};
    ec_isom_t xmap=map->xmap;
    for (int i=0;i<3;++i) {
        if (fp2_is_zero(&points[i]->z)) {
            if (fp2_is_zero(&points[i]->x)) return 0;
            ec_set_zero(points[i]); continue;
        }
        ec_iso_eval(points[i],&xmap);
        if (fp2_is_zero(&points[i]->z)) return 0;
        fp2_t inverse=points[i]->z; fp2_inv(&inverse);
        fp2_mul(&points[i]->x,&points[i]->x,&inverse); fp2_set_one(&points[i]->z);
    }
    *out=result; return 1;
}
