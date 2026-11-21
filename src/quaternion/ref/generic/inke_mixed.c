// SPDX-License-Identifier: Apache-2.0
/* Mixed norm equation following the construction discussed in SPRITZ
 * (ePrint 2026/945, Appendix C.3). Written against POKE's C arithmetic.
 * No Rust/Sage runtime or code is used. See the implementation notes for
 * the equations, exact checks, and the remaining isogeny stage. */
#include <inke_mixed.h>
#include "internal.h"
#include <rng.h>

/* Volatile stores retain the erasure even when the state is about to die. */
static void mixed_erase(void *memory, size_t length)
{
    volatile unsigned char *p = (volatile unsigned char *)memory;
    while (length--) *p++ = 0;
}

static void mixed_random_reset(inke_mixed_state *s)
{
    mixed_erase(s->random_buffer, sizeof(s->random_buffer));
    s->random_position = sizeof(s->random_buffer);
}

/* Exactly the old [1,100000] rejection distribution. Explicit byte decoding
 * agrees with ibz_rand_interval's little-endian byte interpretation. Filling
 * a buffer does not introduce a PRNG or alter the system entropy source. */
static int mixed_random_coefficient(unsigned long *out, inke_mixed_state *s)
{
    for (;;) {
        if (s->random_position >= sizeof(s->random_buffer)) {
            mixed_random_reset(s);
            if (randombytes(s->random_buffer, sizeof(s->random_buffer)) != 0) {
                mixed_random_reset(s);
                return 0;
            }
            s->random_position = 0;
        }
        const unsigned char *p = s->random_buffer + s->random_position;
        unsigned long value = (unsigned long)p[0] |
                              ((unsigned long)p[1] << 8) |
                              (((unsigned long)p[2] & 1u) << 16);
        mixed_erase(s->random_buffer + s->random_position, 3);
        s->random_position += 3;
        if (value < 100000ul) {
            *out = value + 1;
            return 1;
        }
    }
}

/* Cache the numerator form, not a rounded or heuristically integral form. */
static void mixed_prepare_norm(inke_mixed_state *s)
{
    ibz_t term;
    ibz_init(&term);
    mpz_mul(s->norm_denom, s->reduced.lattice.denom, s->reduced.lattice.denom);
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            mpz_set_ui(s->norm_form[i][j], 0);
            if (j<i) continue;
            for (int k=0; k<4; ++k) {
                mpz_mul(term, s->basis[k][i], s->basis[k][j]);
                if (k>=2) mpz_mul(term, term, s->algebra->p);
                mpz_add(s->norm_form[i][j], s->norm_form[i][j], term);
            }
            if (i!=j) mpz_mul_2exp(s->norm_form[i][j], s->norm_form[i][j], 1);
        }
    }
    ibz_finalize(&term);
}

static int mixed_candidate_norm(ibz_t *out, const inke_mixed_state *s,
                                const unsigned long coeff[4], ibz_t *scratch)
{
    mpz_set_ui(*out, 0);
    for (int i=0; i<4; ++i) {
        for (int j=i; j<4; ++j) {
            /* Two small multiplications also work with 32-bit unsigned long:
             * the product of two coefficients need not fit in that type. */
            mpz_mul_ui(*scratch, s->norm_form[i][j], coeff[i]);
            mpz_addmul_ui(*out, *scratch, coeff[j]);
        }
    }
    if (!mpz_divisible_p(*out, s->norm_denom)) return 0;
    mpz_divexact(*out, *out, s->norm_denom);
    return 1;
}

static void elem_copy(quat_alg_elem_t *out, const quat_alg_elem_t *in)
{
    ibz_copy(&out->denom, &in->denom);
    for (int i=0; i<4; ++i) ibz_copy(&out->coord[i], &in->coord[i]);
}

static int elem_equal(const quat_alg_elem_t *a, const quat_alg_elem_t *b)
{
    ibz_t left, right;
    ibz_init(&left); ibz_init(&right);
    int equal = 1;
    for (int i=0; i<4 && equal; ++i) {
        mpz_mul(left,a->coord[i],b->denom);
        mpz_mul(right,b->coord[i],a->denom);
        equal = mpz_cmp(left,right)==0;
    }
    ibz_finalize(&right); ibz_finalize(&left);
    return equal;
}

static void ideal_copy(quat_left_ideal_t *out, const quat_left_ideal_t *in)
{
    ibz_mat_4x4_copy(&out->lattice.basis, &in->lattice.basis);
    ibz_copy(&out->lattice.denom, &in->lattice.denom);
    ibz_copy(&out->norm, &in->norm);
    out->parent_order = in->parent_order;
}

static int norm_integer(ibz_t *norm, const quat_alg_elem_t *x, const quat_alg_t *alg)
{
    ibq_t rational;
    ibq_init(&rational);
    quat_alg_norm(&rational, x, alg);
    int ok = ibq_to_ibz(norm, &rational);
    ibq_finalize(&rational);
    return ok;
}

/* Direct multiplication of a lattice by a quaternion. Avoids the random
 * generator search in the older POKE quat_lideal_mul() helper. */
static void multiply_lattice(quat_lattice_t *out, const quat_lattice_t *in,
                              const quat_alg_elem_t *right, const quat_alg_t *alg)
{
    ibz_mat_4x4_t matrix;
    ibz_mat_4x4_init(&matrix);
    quat_alg_rightmul_mat(&matrix, right, alg);
    ibz_mat_4x4_mul(&out->basis, &matrix, &in->basis);
    ibz_mul(&out->denom, &in->denom, &right->denom);
    quat_lattice_reduce_denom(out, out);
    quat_lattice_hnf(out);
    ibz_mat_4x4_finalize(&matrix);
}

static int integral_lattice(const quat_lattice_t *lattice,
                            const quat_order_t *order, const quat_alg_t *alg)
{
    quat_alg_elem_t element;
    quat_alg_elem_init(&element);
    ibz_copy(&element.denom, &lattice->denom);
    int ok = 1;
    for (int j=0; j<4 && ok; ++j) {
        for (int i=0; i<4; ++i) ibz_copy(&element.coord[i], &lattice->basis[i][j]);
        ok = quat_lattice_contains(NULL, order, &element, alg);
    }
    quat_alg_elem_finalize(&element);
    return ok;
}

void inke_mixed_state_init(inke_mixed_state *s)
{
    quat_left_ideal_init(&s->input); quat_left_ideal_init(&s->reduced);
    quat_alg_elem_init(&s->beta); ibz_mat_4x4_init(&s->basis); ibz_init(&s->target);
    ibz_mat_4x4_init(&s->norm_form); ibz_init(&s->norm_denom);
    mixed_random_reset(s);
    s->algebra = NULL; s->ready = 0;
    s->attempts = s->generators = s->positive_targets = s->square_solutions = 0;
}
void inke_mixed_state_clear(inke_mixed_state *s)
{
    mixed_random_reset(s);
    ibz_finalize(&s->norm_denom); ibz_mat_4x4_finalize(&s->norm_form);
    ibz_finalize(&s->target); ibz_mat_4x4_finalize(&s->basis);
    quat_alg_elem_finalize(&s->beta);
    quat_left_ideal_finalize(&s->reduced); quat_left_ideal_finalize(&s->input);
    s->ready = 0;
}
void inke_mixed_witness_init(inke_mixed_witness *w)
{
    quat_alg_elem_init(&w->mu1); quat_alg_elem_init(&w->mu2);
    quat_alg_elem_init(&w->theta); quat_alg_elem_init(&w->transport);
    ibz_init(&w->d1); ibz_init(&w->d2);
}
void inke_mixed_witness_clear(inke_mixed_witness *w)
{
    ibz_finalize(&w->d2); ibz_finalize(&w->d1);
    quat_alg_elem_finalize(&w->transport); quat_alg_elem_finalize(&w->theta);
    quat_alg_elem_finalize(&w->mu2); quat_alg_elem_finalize(&w->mu1);
}

int inke_mixed_prepare(inke_mixed_state *s, const quat_left_ideal_t *ideal,
                       const ibz_t *target, const quat_alg_t *alg)
{
    s->ready = 0;
    mixed_random_reset(s);
    if (!ideal->parent_order || mpz_sgn(ideal->norm) <= 0 ||
        mpz_sgn(ideal->lattice.denom) <= 0 || mpz_sgn(*target) <= 0 ||
        mpz_fdiv_ui(*target, 6) != 0 ||
        !integral_lattice(&ideal->lattice, ideal->parent_order, alg))
        return INKE_SAMPLE_INVALID;
    ideal_copy(&s->input, ideal);
    ibz_copy(&s->target, target); s->algebra = alg;
    s->attempts = s->generators = s->positive_targets = s->square_solutions = 0;
    ibz_mat_4x4_t reduced;
    ibz_mat_4x4_init(&reduced);
    quat_alg_elem_t conjugate;
    quat_alg_elem_init(&conjugate);
    ibz_t norm, index, square;
    ibz_init(&norm); ibz_init(&index); ibz_init(&square);
    int status = INKE_SAMPLE_INVALID;
    quat_lattice_index(&index, &ideal->lattice, ideal->parent_order);
    mpz_mul(square, ideal->norm, ideal->norm);
    if (mpz_cmp(index, square) != 0) goto done;
    if (quat_lattice_lll(&reduced, &ideal->lattice, &alg->p) != 0) goto done;
    for (int i=0; i<4; ++i) ibz_copy(&s->beta.coord[i], &reduced[i][0]);
    ibz_copy(&s->beta.denom, &ideal->lattice.denom);
    quat_alg_normalize(&s->beta);
    if (!quat_lattice_contains(NULL, &ideal->lattice, &s->beta, alg) ||
        !norm_integer(&norm, &s->beta, alg) || mpz_sgn(norm) <= 0 ||
        !mpz_divisible_p(norm, ideal->norm)) goto done;
    mpz_divexact(s->reduced.norm, norm, ideal->norm);
    s->reduced.parent_order = ideal->parent_order;
    quat_alg_conj(&conjugate, &s->beta);
    mpz_mul(conjugate.denom, conjugate.denom, ideal->norm);
    multiply_lattice(&s->reduced.lattice, &ideal->lattice, &conjugate, alg);
    if (!integral_lattice(&s->reduced.lattice, ideal->parent_order, alg)) goto done;
    quat_lattice_index(&index, &s->reduced.lattice, ideal->parent_order);
    mpz_mul(square, s->reduced.norm, s->reduced.norm);
    if (mpz_cmp(index, square) != 0) goto done;
    if (quat_lattice_lll(&s->basis, &s->reduced.lattice, &alg->p) != 0) goto done;
    mixed_prepare_norm(s);
    s->ready = 1;
    status = INKE_SAMPLE_OK;
done:
    ibz_finalize(&square); ibz_finalize(&index); ibz_finalize(&norm);
    quat_alg_elem_finalize(&conjugate); ibz_mat_4x4_finalize(&reduced);
    return status;
}

/* Euclidean Cornacchia with at most 198 square-root candidates. Even if a
 * probable-prime test accepts a composite, acceptance requires x^2+y^2=n. */
static int prime_squares(ibz_t *x, ibz_t *y, const ibz_t *n)
{
    if (mpz_cmp_ui(*n, 2) == 0) { mpz_set_ui(*x,1); mpz_set_ui(*y,1); return 1; }
    if (mpz_cmp_ui(*n, 2) < 0 || mpz_fdiv_ui(*n, 4) != 1 ||
        !mpz_probab_prime_p(*n, 32)) return 0;
    ibz_t exponent, base, u, v, t;
    ibz_init(&exponent); ibz_init(&base); ibz_init(&u); ibz_init(&v); ibz_init(&t);
    mpz_sub_ui(exponent, *n, 1); mpz_fdiv_q_2exp(exponent, exponent, 2);
    int ok = 0;
    for (unsigned long z=2; z<200; ++z) {
        mpz_set_ui(base,z); mpz_powm(v,base,exponent,*n);
        mpz_mul(t,v,v); mpz_add_ui(t,t,1); mpz_mod(t,t,*n);
        if (mpz_sgn(t) == 0) { ok = 1; break; }
    }
    if (!ok) goto done;
    mpz_set(u,*n); mpz_mul(t,v,v);
    while (mpz_cmp(t,*n) > 0) {
        mpz_mod(t,u,v); mpz_set(u,v); mpz_set(v,t); mpz_mul(t,v,v);
    }
    mpz_sub(t,*n,t); mpz_sqrt(*y,t); mpz_mul(u,*y,*y);
    ok = mpz_cmp(u,t) == 0;
    if (ok) mpz_set(*x,v);
done:
    ibz_finalize(&t); ibz_finalize(&v); ibz_finalize(&u);
    ibz_finalize(&base); ibz_finalize(&exponent);
    return ok;
}

static void complex_multiply(ibz_t *a, ibz_t *b, const ibz_t *x, const ibz_t *y)
{
    ibz_t re, im, tmp;
    ibz_init(&re); ibz_init(&im); ibz_init(&tmp);
    mpz_mul(re,*a,*x); mpz_mul(tmp,*b,*y); mpz_sub(re,re,tmp);
    mpz_mul(im,*a,*y); mpz_mul(tmp,*b,*x); mpz_add(im,im,tmp);
    mpz_set(*a,re); mpz_set(*b,im);
    ibz_finalize(&tmp); ibz_finalize(&im); ibz_finalize(&re);
}

static int two_squares(ibz_t *x, ibz_t *y, const ibz_t *value)
{
    static const unsigned short primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,
        47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,
        151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,
        251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,
        359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,
        461,463,467,479,487,491,499};
    if (mpz_sgn(*value) <= 0) return 0;
    ibz_t n, prime, scale, u, v, check;
    ibz_init(&n); ibz_init(&prime); ibz_init(&scale);
    ibz_init(&u); ibz_init(&v); ibz_init(&check);
    mpz_set(n,*value); mpz_set_ui(*x,1); mpz_set_ui(*y,0);
    int ok = 0;
    for (unsigned int i=0; i<sizeof(primes)/sizeof(primes[0]); ++i) {
        unsigned long ell=primes[i];
        if (!mpz_divisible_ui_p(n,ell)) continue;
        mpz_set_ui(prime,ell);
        mp_bitcnt_t exponent = mpz_remove(n,n,prime);
        if (ell%4 == 3 && (exponent&1)) goto done;
        mpz_pow_ui(scale,prime,(unsigned long)(exponent/2));
        mpz_mul(*x,*x,scale); mpz_mul(*y,*y,scale);
        if (exponent&1) {
            if (!prime_squares(&u,&v,&prime)) goto done;
            complex_multiply(x,y,&u,&v);
        }
    }
    if (mpz_cmp_ui(n,1) != 0) {
        if (!prime_squares(&u,&v,&n)) goto done;
        complex_multiply(x,y,&u,&v);
    }
    mpz_mul(check,*x,*x); mpz_mul(u,*y,*y); mpz_add(check,check,u);
    ok = mpz_cmp(check,*value) == 0;
    mpz_abs(*x,*x); mpz_abs(*y,*y);
done:
    ibz_finalize(&check); ibz_finalize(&v); ibz_finalize(&u);
    ibz_finalize(&scale); ibz_finalize(&prime); ibz_finalize(&n);
    return ok;
}

static void nearest_integer(ibz_t *out, const ibz_t *num, const ibz_t *den)
{
    ibz_t n,d;
    ibz_init(&n); ibz_init(&d);
    mpz_set(n,*num); mpz_set(d,*den);
    if (mpz_sgn(d)<0) { mpz_neg(n,n); mpz_neg(d,d); }
    mpz_mul_2exp(n,n,1); mpz_add(n,n,d); mpz_mul_2exp(d,d,1);
    mpz_fdiv_q(*out,n,d);
    ibz_finalize(&d); ibz_finalize(&n);
}

static void make_mu(quat_alg_elem_t *out, const quat_alg_elem_t *alpha,
                    const ibz_t *n, const ibz_t *A, const ibz_t *B,
                    const ibz_t *X, const ibz_t *Y, int sign)
{
    ibz_t delta;
    ibz_init(&delta); elem_copy(out,alpha);
    for (int i=0; i<2; ++i) {
        const ibz_t *center = i ? B : A, *root = i ? Y : X;
        if (sign>0) mpz_add(delta,*center,*root); else mpz_sub(delta,*center,*root);
        mpz_divexact_ui(delta,delta,2); mpz_mul(delta,delta,*n);
        mpz_mul(delta,delta,out->denom); mpz_add(out->coord[i],out->coord[i],delta);
    }
    quat_alg_normalize(out); ibz_finalize(&delta);
}

int inke_mixed_search(inke_mixed_witness *out, inke_mixed_state *s,
                      unsigned long budget)
{
    if (!s->ready) { mixed_random_reset(s); return INKE_SAMPLE_INVALID; }
    const ibz_t *n = &s->reduced.norm;
    ibz_t na,n2,a0,a1,r,inv,x,t,det,c1,c2,A,B,k,value,X,Y,tmp,gcd;
    ibz_t *ints[]={&na,&n2,&a0,&a1,&r,&inv,&x,&t,&det,&c1,&c2,&A,&B,&k,&value,&X,&Y,&tmp,&gcd};
    for (unsigned int i=0;i<sizeof(ints)/sizeof(ints[0]);++i) ibz_init(ints[i]);
    ibz_mat_2x2_t lattice, reduced;
    ibz_mat_2x2_init(&lattice); ibz_mat_2x2_init(&reduced);
    unsigned long coeff[4];
    quat_alg_elem_t alpha, conjugate;
    quat_alg_elem_init(&alpha); quat_alg_elem_init(&conjugate);
    inke_mixed_witness candidate;
    inke_mixed_witness_init(&candidate);
    int status = INKE_SAMPLE_RETRY;
    mpz_mul(n2,*n,*n);
    for (unsigned long trial=0;trial<budget;++trial) {
        ++s->attempts;
        for (int i=0;i<4;++i) {
            if (!mixed_random_coefficient(&coeff[i],s)) {
                status=INKE_SAMPLE_RNG_ERROR; goto done;
            }
        }
        for (int i=0; i<4; ++i) {
            mpz_set_ui(alpha.coord[i],0);
            for (int j=0; j<4; ++j)
                mpz_addmul_ui(alpha.coord[i],s->basis[i][j],coeff[j]);
        }
        ibz_copy(&alpha.denom,&s->reduced.lattice.denom);
        quat_alg_normalize(&alpha);
        mpz_mul_2exp(a0,alpha.coord[0],1); mpz_mul_2exp(a1,alpha.coord[1],1);
        if (!mpz_divisible_p(a0,alpha.denom) || !mpz_divisible_p(a1,alpha.denom) ||
            !mixed_candidate_norm(&na,s,coeff,&tmp) || !mpz_divisible_p(na,*n)) {
            status=INKE_SAMPLE_INVALID; goto done;
        }
        mpz_divexact(a0,a0,alpha.denom); mpz_divexact(a1,a1,alpha.denom);
        mpz_gcd(gcd,na,n2); if (mpz_cmp(gcd,*n)!=0) continue;
        if (mpz_cmp_ui(*n,1)==0) mpz_set_ui(inv,0);
        else if (!mpz_invert(inv,a0,*n)) continue;
        ++s->generators;
        mpz_divexact(r,na,*n); mpz_mul_2exp(r,r,1); mpz_sub(r,s->target,r);
        if (mpz_sgn(r)<=0) continue;
        mpz_mul(x,a1,inv); mpz_mod(x,x,*n);
        mpz_mul(t,r,inv); mpz_mod(t,t,*n);
        mpz_sub(lattice[0][0],*n,x); mpz_set_ui(lattice[1][0],1);
        mpz_set(lattice[0][1],*n); mpz_set_ui(lattice[1][1],0);
        quat_dim2_lattice_short_basis(&reduced,&lattice,&ibz_const_one);
        mpz_mul(det,reduced[0][0],reduced[1][1]);
        mpz_mul(tmp,reduced[0][1],reduced[1][0]); mpz_sub(det,det,tmp);
        if (mpz_sgn(det)==0) { status=INKE_SAMPLE_INVALID; goto done; }
        mpz_mul(tmp,t,reduced[1][1]); mpz_neg(tmp,tmp); nearest_integer(&c1,&tmp,&det);
        mpz_mul(tmp,t,reduced[1][0]); nearest_integer(&c2,&tmp,&det);
        mpz_mul(A,c1,reduced[0][0]); mpz_addmul(A,c2,reduced[0][1]); mpz_add(A,A,t);
        mpz_mul(B,c1,reduced[1][0]); mpz_addmul(B,c2,reduced[1][1]);
        mpz_set(k,r); mpz_submul(k,a0,A); mpz_submul(k,a1,B);
        if (!mpz_divisible_p(k,*n)) { status=INKE_SAMPLE_INVALID; goto done; }
        mpz_divexact(k,k,*n); mpz_mul_2exp(value,k,1);
        mpz_submul(value,A,A); mpz_submul(value,B,B);
        if (mpz_sgn(value)<=0) continue;
        ++s->positive_targets;
        if (!two_squares(&X,&Y,&value)) continue;
        ++s->square_solutions;
        for (int variant=0;variant<8;++variant) {
            if (variant==4) mpz_swap(X,Y);
            mpz_abs(X,X); mpz_abs(Y,Y);
            if (variant&1) mpz_neg(X,X);
            if (variant&2) mpz_neg(Y,Y);
            if (mpz_fdiv_ui(X,2)!=mpz_fdiv_ui(A,2) ||
                mpz_fdiv_ui(Y,2)!=mpz_fdiv_ui(B,2)) continue;
            make_mu(&candidate.mu1,&alpha,n,&A,&B,&X,&Y,1);
            make_mu(&candidate.mu2,&alpha,n,&A,&B,&X,&Y,-1);
            if (!norm_integer(&candidate.d1,&candidate.mu1,s->algebra) ||
                !norm_integer(&candidate.d2,&candidate.mu2,s->algebra) ||
                !mpz_divisible_p(candidate.d1,*n) || !mpz_divisible_p(candidate.d2,*n)) {
                status=INKE_SAMPLE_INVALID; goto done;
            }
            mpz_divexact(candidate.d1,candidate.d1,*n); mpz_divexact(candidate.d2,candidate.d2,*n);
            if (mpz_sgn(candidate.d1)<=0 || mpz_sgn(candidate.d2)<=0 ||
                mpz_even_p(candidate.d1) || mpz_even_p(candidate.d2) ||
                mpz_divisible_ui_p(candidate.d1,3) || mpz_divisible_ui_p(candidate.d2,3)) continue;
            mpz_gcd(gcd,candidate.d1,candidate.d2); if (mpz_cmp_ui(gcd,1)!=0) continue;
            quat_alg_conj(&conjugate,&candidate.mu1);
            quat_alg_mul(&candidate.theta,&candidate.mu2,&conjugate,s->algebra);
            mpz_mul(candidate.theta.denom,candidate.theta.denom,*n);
            quat_alg_normalize(&candidate.theta);
            quat_alg_mul(&candidate.transport,&candidate.mu1,&s->beta,s->algebra);
            mpz_mul(candidate.transport.denom,candidate.transport.denom,*n);
            quat_alg_normalize(&candidate.transport);
            if (!inke_mixed_verify(s,&candidate)) { status=INKE_SAMPLE_INVALID; goto done; }
            elem_copy(&out->mu1,&candidate.mu1); elem_copy(&out->mu2,&candidate.mu2);
            elem_copy(&out->theta,&candidate.theta); elem_copy(&out->transport,&candidate.transport);
            ibz_copy(&out->d1,&candidate.d1); ibz_copy(&out->d2,&candidate.d2);
            status=INKE_SAMPLE_OK; goto done;
        }
    }
done:
    if (status < 0) mixed_random_reset(s);
    inke_mixed_witness_clear(&candidate);
    quat_alg_elem_finalize(&conjugate); quat_alg_elem_finalize(&alpha);
    mixed_erase(coeff,sizeof(coeff));
    ibz_mat_2x2_finalize(&reduced); ibz_mat_2x2_finalize(&lattice);
    for (unsigned int i=0;i<sizeof(ints)/sizeof(ints[0]);++i) ibz_finalize(ints[i]);
    return status;
}

int inke_mixed_verify(const inke_mixed_state *s, const inke_mixed_witness *w)
{
    if (!s->ready || mpz_sgn(w->d1)<=0 || mpz_sgn(w->d2)<=0 ||
        mpz_sgn(w->mu1.denom)<=0 || mpz_sgn(w->mu2.denom)<=0 ||
        mpz_sgn(w->theta.denom)<=0 || mpz_sgn(w->transport.denom)<=0) return 0;
    const quat_alg_t *alg=s->algebra;
    const quat_order_t *order=s->input.parent_order;
    const ibz_t *n=&s->reduced.norm;
    ibz_t norm,expected,gcd,index;
    ibz_init(&norm); ibz_init(&expected); ibz_init(&gcd); ibz_init(&index);
    quat_alg_elem_t conjugate, computed;
    quat_alg_elem_init(&conjugate); quat_alg_elem_init(&computed);
    quat_lattice_t J1, lhs, rhs;
    quat_lattice_init(&J1); quat_lattice_init(&lhs); quat_lattice_init(&rhs);
    int ok=0;
    mpz_add(expected,w->d1,w->d2); if (mpz_cmp(expected,s->target)!=0) goto done;
    mpz_gcd(gcd,w->d1,w->d2); if (mpz_cmp_ui(gcd,1)!=0) goto done;
    mpz_mul(expected,w->d1,w->d2); mpz_gcd_ui(gcd,expected,6);
    if (mpz_cmp_ui(gcd,1)!=0) goto done;
    if (!quat_lattice_contains(NULL,&s->input.lattice,&s->beta,alg) ||
        !norm_integer(&norm,&s->beta,alg)) goto done;
    mpz_mul(expected,s->input.norm,*n); if (mpz_cmp(norm,expected)!=0) goto done;
    quat_alg_conj(&conjugate,&s->beta); mpz_mul(conjugate.denom,conjugate.denom,s->input.norm);
    multiply_lattice(&lhs,&s->input.lattice,&conjugate,alg);
    if (!quat_lattice_equal(&lhs,&s->reduced.lattice)) goto done;
    if (!quat_lattice_contains(NULL,&s->reduced.lattice,&w->mu1,alg) ||
        !quat_lattice_contains(NULL,&s->reduced.lattice,&w->mu2,alg)) goto done;
    mpz_mul(expected,*n,w->d1);
    if (!norm_integer(&norm,&w->mu1,alg) || mpz_cmp(norm,expected)!=0) goto done;
    mpz_mul(expected,*n,w->d2);
    if (!norm_integer(&norm,&w->mu2,alg) || mpz_cmp(norm,expected)!=0) goto done;
    quat_alg_conj(&conjugate,&w->mu1);
    quat_alg_mul(&computed,&w->mu2,&conjugate,alg); mpz_mul(computed.denom,computed.denom,*n);
    if (!elem_equal(&computed,&w->theta)) goto done;
    quat_alg_mul(&computed,&w->mu1,&s->beta,alg); mpz_mul(computed.denom,computed.denom,*n);
    if (!elem_equal(&computed,&w->transport)) goto done;
    if (!quat_lattice_contains(NULL,order,&w->theta,alg) ||
        !quat_lattice_contains(NULL,order,&w->transport,alg)) goto done;
    mpz_mul(expected,w->d1,w->d2);
    if (!norm_integer(&norm,&w->theta,alg) || mpz_cmp(norm,expected)!=0) goto done;
    mpz_mul(expected,w->d1,s->input.norm);
    if (!norm_integer(&norm,&w->transport,alg) || mpz_cmp(norm,expected)!=0) goto done;
    mpz_mul(conjugate.denom,conjugate.denom,*n);
    multiply_lattice(&J1,&s->reduced.lattice,&conjugate,alg);
    if (!integral_lattice(&J1,order,alg)) goto done;
    quat_lattice_index(&index,&J1,order); mpz_mul(expected,w->d1,w->d1);
    if (mpz_cmp(index,expected)!=0) goto done;
    multiply_lattice(&lhs,&J1,&w->transport,alg);
    ibz_mat_4x4_scalar_mul(&rhs.basis,&w->d1,&s->input.lattice.basis);
    ibz_copy(&rhs.denom,&s->input.lattice.denom);
    ok=quat_lattice_equal(&lhs,&rhs);
done:
    quat_lattice_finalize(&rhs); quat_lattice_finalize(&lhs); quat_lattice_finalize(&J1);
    quat_alg_elem_finalize(&computed); quat_alg_elem_finalize(&conjugate);
    ibz_finalize(&index); ibz_finalize(&gcd); ibz_finalize(&expected); ibz_finalize(&norm);
    return ok;
}
