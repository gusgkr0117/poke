#include <ec.h>

/** @brief POKE secret key
 *
 * @typedef poke_sk_t
 *
 * @struct poke_sk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct poke_sk_t {
    digit_t deg[NWORDS_ORDER];
    digit_t alpha[NWORDS_ORDER];
    digit_t beta[NWORDS_ORDER];
    digit_t delta[NWORDS_ORDER];
} poke_sk_t;

/** @brief POKE public key
 *
 * @typedef poke_pk_t
 *
 * @struct poke_pk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct poke_pk_t {
    ec_curve_t EA;
    ec_basis_t PQ2;
    ec_basis_t PQ3;
    ec_basis_t PQxy;
} poke_pk_t;

/** @brief POKE ciphertext
 *
 * @typedef poke_ct_t
 *
 * @struct poke_ct_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct poke_ct_t {
    ec_curve_t EB;
    ec_curve_t EAB;
    ec_basis_t PQ2_B;
    ec_basis_t PQxy_B;
    ec_basis_t PQ2_AB;
    uint8_t ct[32];
} poke_ct_t;

int decaps(unsigned char *key, poke_ct_t *ct, const poke_pk_t *pk, const poke_sk_t *sk, unsigned char *dummy_m);
int encaps(unsigned char *key, poke_ct_t *ct, const poke_pk_t *pk);
int ct_encode(unsigned char *encoded_ct, poke_ct_t *ct);
int ct_decode(poke_ct_t *ct, const unsigned char *encoded_ct);
int keygen(poke_sk_t *sk, poke_pk_t *pk);
int encrypt(poke_ct_t *ct, const poke_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len);
int decrypt(unsigned char *m, size_t *m_len, const poke_ct_t *ct, const poke_sk_t *sk);