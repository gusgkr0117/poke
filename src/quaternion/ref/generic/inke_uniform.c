// SPDX-License-Identifier: Apache-2.0
#include <inke_uniform.h>

static int integral_norm(ibz_t *out, const quat_alg_elem_t *x,
                         const quat_alg_t *algebra)
{
    ibq_t norm;
    ibq_init(&norm);
    quat_alg_norm(&norm, x, algebra);
    int ok = ibq_to_ibz(out, &norm);
    ibq_finalize(&norm);
    return ok;
}

int inke_sample_q(ibz_t *q, unsigned int a)
{
    if (a < 4 || (a & 1)) return INKE_SAMPLE_INVALID;
    ibz_t minimum, maximum, k;
    ibz_init(&minimum);
    ibz_init(&maximum);
    ibz_init(&k);
    /* q=6k+5, with 2^(a-2) <= q <= 3*2^(a-2). */
    mpz_set_ui(minimum, 1);
    mpz_mul_2exp(minimum, minimum, a-2);
    mpz_mul_ui(maximum, minimum, 3);
    mpz_sub_ui(minimum, minimum, 5);
    mpz_cdiv_q_ui(minimum, minimum, 6);
    mpz_sub_ui(maximum, maximum, 5);
    mpz_fdiv_q_ui(maximum, maximum, 6);
    int ok = ibz_rand_interval(&k, &minimum, &maximum);
    if (ok) {
        mpz_mul_ui(*q, k, 6);
        mpz_add_ui(*q, *q, 5);
    }
    ibz_finalize(&k);
    ibz_finalize(&maximum);
    ibz_finalize(&minimum);
    return ok ? INKE_SAMPLE_OK : INKE_SAMPLE_RNG_ERROR;
}

int inke_q_is_admissible(const ibz_t *q, unsigned int a)
{
    if (a < 4 || (a & 1) || mpz_fdiv_ui(*q,6) != 5) return 0;
    ibz_t bound; ibz_init(&bound);
    mpz_set_ui(bound,1); mpz_mul_2exp(bound,bound,a-2);
    int ok = mpz_cmp(*q,bound) >= 0;
    mpz_mul_ui(bound,bound,3); ok = ok && mpz_cmp(*q,bound) <= 0;
    ibz_finalize(&bound); return ok;
}

int inke_ideal_from_unit(quat_left_ideal_t *out,
                         const quat_alg_elem_t *gamma,
                         const quat_alg_elem_t *unit,
                         const ibz_t *d,
                         const quat_order_t *order,
                         const quat_alg_t *algebra)
{
    if (mpz_cmp_ui(*d, 1) <= 0 || mpz_sgn(gamma->denom) == 0 ||
        mpz_sgn(unit->denom) == 0 || mpz_sgn(order->denom) <= 0)
        return INKE_SAMPLE_INVALID;
    ibz_t n, gcd;
    ibz_init(&n);
    ibz_init(&gcd);
    quat_alg_elem_t product;
    quat_alg_elem_init(&product);
    quat_left_ideal_t candidate;
    quat_left_ideal_init(&candidate);
    int status = INKE_SAMPLE_INVALID;

    if (!quat_lattice_contains(NULL, order, gamma, algebra) ||
        !quat_lattice_contains(NULL, order, unit, algebra) ||
        !quat_alg_is_primitive(gamma, order, algebra) ||
        !integral_norm(&n, gamma, algebra) || mpz_sgn(n) <= 0 ||
        !mpz_divisible_p(n, *d)) goto done;
    mpz_gcd(gcd, *d, algebra->p);
    if (mpz_cmp_ui(gcd, 1) != 0) goto done;
    if (!integral_norm(&n, unit, algebra)) goto done;
    mpz_gcd(gcd, n, *d);
    if (mpz_cmp_ui(gcd, 1) != 0) goto done;

    quat_alg_mul(&product, gamma, unit, algebra);
    /* The product can have scalar content prime to d. Removing that
     * content leaves the ideal modulo d unchanged. This constructor also
     * checks the content against d; verify that the prescribed norm stays d. */
    quat_lideal_make_primitive_then_create(&candidate, &product, d, order, algebra);
    if (mpz_cmp(candidate.norm, *d) != 0) goto done;
    ibz_mat_4x4_copy(&out->lattice.basis, &candidate.lattice.basis);
    ibz_copy(&out->lattice.denom, &candidate.lattice.denom);
    ibz_copy(&out->norm, &candidate.norm);
    out->parent_order = order;
    status = INKE_SAMPLE_OK;
done:
    quat_left_ideal_finalize(&candidate);
    quat_alg_elem_finalize(&product);
    ibz_finalize(&gcd);
    ibz_finalize(&n);
    return status;
}

/* SQIsign v2 Algorithm 3.12 at q_form=1, isogenyCond=false, O0.
 * Its bounded failure changes only the auxiliary m; q stays fixed.
 * Uniformity comes from the subsequent unit action, not the seed law.
 * Candidate budgets count polynomial-cost trials, not machine steps. */
/* Cornacchia for x^2+y^2=r. The quadratic nonresidue is sampled,
 * not found by an unproved bounded scan through small integers.
 * A finite inner budget is part of the advertised solver budget.
 * RNG failure is propagated separately. */
static int seed_two_squares(ibz_t *x,ibz_t *y,const ibz_t *r) {
    if(mpz_cmp_ui(*r,2)==0) {mpz_set_ui(*x,1);mpz_set_ui(*y,1);return 1;}
    if(mpz_cmp_ui(*r,2)<0 || mpz_fdiv_ui(*r,4)!=1) return 0;
    ibz_t hi,e,w,t,u,v,s;
    ibz_init(&hi);ibz_init(&e);ibz_init(&w);ibz_init(&t);
    ibz_init(&u);ibz_init(&v);ibz_init(&s);
    mpz_sub_ui(hi,*r,1);mpz_fdiv_q_2exp(e,hi,2);
    int ok=0;
    unsigned long budget=mpz_sizeinbase(*r,2)+128;
    for(unsigned long j=0;j<budget;j++) {
        if(!ibz_rand_interval(&w,&ibz_const_one,&hi)) {ok=-1;goto done;}
        mpz_powm(t,w,e,*r);mpz_mul(s,t,t);mpz_mod(s,s,*r);
        if(mpz_cmp(s,hi)==0) {ok=1;break;}
    }
    if(!ok)goto done;
    mpz_set(u,*r);mpz_set(v,t);
    for(;;) {
        mpz_mul(s,v,v);
        if(mpz_cmp(s,*r)<=0)break;
        mpz_mod(t,u,v);mpz_set(u,v);mpz_set(v,t);
    }
    mpz_sub(s,*r,s);
    ok=mpz_perfect_square_p(s);
    if(ok) {mpz_set(*x,v);mpz_sqrt(*y,s);}
done:
    ibz_finalize(&s);ibz_finalize(&v);ibz_finalize(&u);ibz_finalize(&t);
    ibz_finalize(&w);ibz_finalize(&e);ibz_finalize(&hi);return ok;
}

static int seed_norm(quat_alg_elem_t *gamma,const ibz_t *target,
                       const quat_order_t *order,const quat_alg_t *algebra,
                       unsigned long trial_budget,inke_sample_counts *counts) {
    if(mpz_even_p(*target) || mpz_cmp(*target,algebra->p)<=0 ||
       mpz_fdiv_ui(algebra->p,4)!=3) return INKE_SAMPLE_INVALID;
    ibz_t four_m,sq_bound,z_bound,t_bound,lo,rem,tmp,n;
    ibz_init(&four_m);ibz_init(&sq_bound);ibz_init(&z_bound);
    ibz_init(&t_bound);ibz_init(&lo);ibz_init(&rem);ibz_init(&tmp);ibz_init(&n);
    mpz_mul_2exp(four_m,*target,2);mpz_fdiv_q(sq_bound,four_m,algebra->p);
    /* Source line4 subtracts q_form=1 inside the square root. */
    mpz_sub_ui(z_bound,sq_bound,1);mpz_sqrt(z_bound,z_bound);
    /* Algorithm3.12's ceiling(4M/p) cap, truncated by caller's T. */
    mpz_cdiv_q(tmp,four_m,algebra->p);
    if(mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp)<trial_budget)trial_budget=mpz_get_ui(tmp);
    int status=INKE_SAMPLE_RETRY;
    for(unsigned long k=0;k<trial_budget;k++) {
        if(counts)++counts->norm_candidates;
        if(!ibz_rand_interval(&gamma->coord[2],&ibz_const_one,&z_bound)) {
            status=INKE_SAMPLE_RNG_ERROR;break;
        }
        mpz_mul(tmp,gamma->coord[2],gamma->coord[2]);
        mpz_sub(t_bound,sq_bound,tmp);mpz_sqrt(t_bound,t_bound);mpz_neg(lo,t_bound);
        if(!ibz_rand_interval(&gamma->coord[3],&lo,&t_bound)) {
            status=INKE_SAMPLE_RNG_ERROR;break;
        }
        mpz_mul(rem,gamma->coord[3],gamma->coord[3]);mpz_add(rem,rem,tmp);
        mpz_mul(rem,rem,algebra->p);mpz_sub(rem,four_m,rem);
        if(mpz_cmp_ui(rem,2)<0 || !mpz_probab_prime_p(rem,64))continue;
        int sq=seed_two_squares(&gamma->coord[0],&gamma->coord[1],&rem);
        if(sq<0) {status=INKE_SAMPLE_RNG_ERROR;break;}
        if(!sq)continue;
        mpz_set_ui(gamma->denom,2);
        /* The spec requires content of the numerator to be exactly2.
         * Equivalently the half-integral candidate lies in O and is primitive.
         * No swapping of x,y is done: isogenyCond=false follows the spec. */
        if(!quat_lattice_contains(NULL,order,gamma,algebra) ||
           !quat_alg_is_primitive(gamma,order,algebra))continue;
        ibq_t norm;ibq_init(&norm);quat_alg_norm(&norm,gamma,algebra);
        int exact=ibq_to_ibz(&n,&norm) && mpz_cmp(n,*target)==0;ibq_finalize(&norm);
        if(!exact)continue;
        quat_alg_normalize(gamma);status=INKE_SAMPLE_OK;break;
    }
    ibz_finalize(&n);ibz_finalize(&tmp);ibz_finalize(&rem);ibz_finalize(&lo);
    ibz_finalize(&t_bound);ibz_finalize(&z_bound);ibz_finalize(&sq_bound);ibz_finalize(&four_m);
    return status;
}

int inke_sample_ideal(quat_left_ideal_t *out,
                      const ibz_t *q, unsigned int a,
                      const quat_order_t *order, const quat_alg_t *algebra,
                      unsigned long prime_trials, unsigned long norm_trials,
                      unsigned long unit_trials, inke_sample_counts *counts)
{
    if (a < 4 || (a & 1) || mpz_sgn(*q) <= 0 || mpz_fdiv_ui(*q, 6) != 5 ||
        mpz_cmp_ui(algebra->p, 3) <= 0 || mpz_fdiv_ui(algebra->p, 4) != 3)
        return INKE_SAMPLE_INVALID;
    ibz_t d, power, lo, hi, m, target, maximum, norm, gcd;
    ibz_init(&d); ibz_init(&power); ibz_init(&lo); ibz_init(&hi);
    ibz_init(&m); ibz_init(&target); ibz_init(&maximum);
    ibz_init(&norm); ibz_init(&gcd);
    quat_alg_elem_t gamma, unit;
    quat_alg_elem_init(&gamma); quat_alg_elem_init(&unit);
    quat_alg_coord_t coords;
    quat_alg_coord_init(&coords);
    int status = INKE_SAMPLE_INVALID;
    mpz_set_ui(power, 1);
    mpz_mul_2exp(power, power, a);
    if (mpz_cmp(*q, power) >= 0) goto done;
    mpz_sub(d, power, *q);
    mpz_mul(d, d, *q);
    if (mpz_cmp(d, algebra->p) >= 0) goto done;
    mpz_add_ui(lo, algebra->p, 1);
    mpz_mul_ui(hi, algebra->p, 2);
    mpz_sub_ui(hi, hi, 1);
    status = INKE_SAMPLE_RETRY;
    for (unsigned long k = 0; k < prime_trials; ++k) {
        if (counts) ++counts->prime_candidates;
        if (!ibz_rand_interval(&m, &lo, &hi)) {
            status = INKE_SAMPLE_RNG_ERROR;
            goto done;
        }
        if (!mpz_probab_prime_p(m, 32)) continue;
        mpz_gcd(gcd, m, d);
        if (mpz_cmp_ui(gcd, 1) != 0) continue;
        mpz_mul(target, m, d);
        status = seed_norm(&gamma, &target, order, algebra, norm_trials, counts);
        if (status == INKE_SAMPLE_RNG_ERROR) goto done;
        if (status == INKE_SAMPLE_OK) break;
    }
    if (status != INKE_SAMPLE_OK) goto done;
    status = INKE_SAMPLE_RETRY;
    mpz_sub_ui(maximum, d, 1);
    for (unsigned long k = 0; k < unit_trials; ++k) {
        if (counts) ++counts->unit_candidates;
        for (int i = 0; i < 4; ++i) {
            if (!ibz_rand_interval(&coords[i], &ibz_const_zero, &maximum)) {
                status = INKE_SAMPLE_RNG_ERROR;
                goto done;
            }
        }
        ibz_mat_4x4_eval(&unit.coord, &order->basis, &coords);
        ibz_copy(&unit.denom, &order->denom);
        if (!integral_norm(&norm, &unit, algebra)) {
            status = INKE_SAMPLE_INVALID;
            goto done;
        }
        mpz_gcd(gcd, norm, d);
        if (mpz_cmp_ui(gcd, 1) != 0) continue;
        status = inke_ideal_from_unit(out, &gamma, &unit, &d, order, algebra);
        break;
    }
done:
    quat_alg_coord_finalize(&coords);
    quat_alg_elem_finalize(&unit); quat_alg_elem_finalize(&gamma);
    ibz_finalize(&gcd); ibz_finalize(&norm); ibz_finalize(&maximum);
    ibz_finalize(&target); ibz_finalize(&m); ibz_finalize(&hi); ibz_finalize(&lo);
    ibz_finalize(&power); ibz_finalize(&d);
    return status;
}
