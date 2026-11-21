#include "root-extension.h"
#include <rng.h>
#include <string.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <bcrypt.h>
#include <limits.h>
#endif

/* Modulus of the linked GF implementation, least-significant word first. */
extern const uint64_t p[NWORDS_FIELD];

int root_system_random(void *unused, unsigned char *out, size_t count)
{
    (void)unused;
#ifdef _WIN32
    HMODULE module = LoadLibraryExW(L"bcrypt.dll", NULL, LOAD_LIBRARY_SEARCH_SYSTEM32);
    if (!module) return 0;
    typedef NTSTATUS (WINAPI *generate_type)(BCRYPT_ALG_HANDLE, PUCHAR, ULONG, ULONG);
    generate_type generate = (generate_type)(void *)GetProcAddress(module, "BCryptGenRandom");
    int ok = generate != NULL;
    while (ok && count) {
        ULONG take = count > ULONG_MAX ? ULONG_MAX : (ULONG)count;
        ok = generate(NULL, out, take, BCRYPT_USE_SYSTEM_PREFERRED_RNG) >= 0;
        out += take;
        count -= take;
    }
    FreeLibrary(module);
    return ok;
#else
    return randombytes(out, count) == 0;
#endif
}

int root_extension_init_rng(root_extension *E, const mpz_t prime, slong k,
                            root_random_bytes rng, void *context)
{
    if (!E || !rng || k != 1) return 0;
    memset(E, 0, sizeof(*E));
    mpz_init(E->field->prime);
    mpz_import(E->field->prime, NWORDS_FIELD, -1, sizeof(*p), 0, 0, p);
    if (mpz_cmp(E->field->prime, prime) != 0 || mpz_fdiv_ui(prime, 4) != 3) {
        mpz_clear(E->field->prime);
        return 0;
    }
    mpz_init(E->half);
    mpz_add_ui(E->half, prime, 1);
    mpz_divexact_ui(E->half, E->half, 2);
    fp_set_small(&E->field->inverse_four, 4);
    fp_inv(&E->field->inverse_four);
    fp_set_zero(&E->i->re);
    fp_set_one(&E->i->im);
    E->relative_degree = 1;
    E->rng = rng;
    E->rng_context = context;
    return 1;
}

int root_extension_init(root_extension *E, const mpz_t prime, slong k)
{
    return root_extension_init_rng(E, prime, k, root_system_random, NULL);
}

void root_extension_clear(root_extension *E)
{
    if (!E || !E->relative_degree) return;
    mpz_clear(E->half);
    mpz_clear(E->field->prime);
    memset(E, 0, sizeof(*E));
}

void fq_mul_ui(fq_t out, const fq_t x, ulong n, const fq_ctx_t f)
{
    (void)f;
    if (n == 2) {
        fp2_add(out, x, x);
        return;
    }
    fp_t scalar;
    fp_set_small(&scalar, n);
    fp_mul(&out->re, &x->re, &scalar);
    fp_mul(&out->im, &x->im, &scalar);
}

void fq_mul_fmpz(fq_t out, const fq_t x, const mpz_t n, const fq_ctx_t f)
{
    unsigned char bytes[sizeof(fp_t)] = {0};
    mpz_t residue;
    mpz_init(residue);
    mpz_mod(residue, n, f->prime);
    mpz_export(bytes, NULL, -1, 1, 0, 0, residue);
    mpz_clear(residue);

    fp_t scalar;
    fp_decode(&scalar, bytes);
    fp_mul(&out->re, &x->re, &scalar);
    fp_mul(&out->im, &x->im, &scalar);
}

void fq_pow_ui(fq_t out, const fq_t x, ulong n, const fq_ctx_t f)
{
    (void)f;
    fp2_t result, power = *x;
    fp2_set_one(&result);
    while (n) {
        if (n & 1) fp2_mul(&result, &result, &power);
        n >>= 1;
        if (n) fp2_sqr(&power, &power);
    }
    *out = result;
}

int root_fq_sqrt(fq_t out, const fq_t x, const fq_ctx_t f)
{
    (void)f;
    fp2_t candidate = *x, check;
    fp2_sqrt(&candidate);
    fp2_sqr(&check, &candidate);
    if (!fp2_is_equal(&check, x)) return 0;
    *out = candidate;
    return 1;
}

int root_extension_sqrt(fq_t out, const fq_t x, const root_extension *E)
{
    return root_fq_sqrt(out, x, E->field);
}
