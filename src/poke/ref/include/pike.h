#include <ec.h>

/** @brief PIKE secret key
 *
 * @typedef pike_sk_t
 *
 * @struct pike_sk_t
 *
 */
typedef struct pike_sk_t {
    digit_t deg[NWORDS_ORDER];
    digit_t alpha[NWORDS_ORDER];
    digit_t beta[NWORDS_ORDER];
    digit_t iota[NWORDS_ORDER];
} pike_sk_t;

/** @brief PIKE public key
 *
 * @typedef pike_pk_t
 *
 * @struct pike_pk_t
 *
 */
typedef struct pike_pk_t {
    ec_curve_t EA;
    ec_basis_t PQ2;
    ec_basis_t PQ3;
    ec_basis_t PQ5;
    ec_point_t imPs;
} pike_pk_t;

/** @brief PIKE ciphertext
 *
 * @typedef pike_ct_t
 *
 * @struct pike_ct_t
 *
 */
typedef struct pike_ct_t {
    ec_curve_t EB;
    ec_curve_t EAB;
    ec_basis_t PQ2_B;
    ec_basis_t PQ2_AB;
    ec_point_t QsB;
    ec_point_t PsAB;
    uint8_t ct[32];
} pike_ct_t;

int decaps(unsigned char *key, pike_ct_t *ct, const pike_pk_t *pk, const pike_sk_t *sk, unsigned char *dummy_m);
int encaps(unsigned char *key, pike_ct_t *ct, const pike_pk_t *pk);
int ct_encode(unsigned char *encoded_ct, pike_ct_t *ct);
int ct_decode(pike_ct_t *ct, const unsigned char *encoded_ct);
int keygen(pike_sk_t *sk, pike_pk_t *pk);
int encrypt(pike_ct_t *ct, const pike_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len);
int decrypt(unsigned char *m, size_t *m_len, const pike_ct_t *ct, const pike_sk_t *sk);