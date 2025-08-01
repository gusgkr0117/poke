#include <ec.h>

/** @brief INKE secret key
 *
 * @typedef inke_sk_t
 *
 * @struct inke_sk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct inke_sk_t {
    digit_t deg[NWORDS_ORDER];
    digit_t alpha[NWORDS_ORDER];
    digit_t beta[NWORDS_ORDER];
    digit_t delta[NWORDS_ORDER];
} inke_sk_t;

/** @brief INKE public key
 *
 * @typedef inke_pk_t
 *
 * @struct inke_pk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct inke_pk_t {
    ec_curve_t EA;
    ec_curve_t EA1;
    ec_basis_t PQ2;
    ec_basis_t PQ3;
    ec_basis_t PQA13;
} inke_pk_t;

/** @brief INKE ciphertext
 *
 * @typedef inke_ct_t
 *
 * @struct inke_ct_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct inke_ct_t {
    ec_curve_t EB;
    ec_curve_t EAB;
    ec_basis_t PQ2_B;
    ec_basis_t PQ2_AB;
    uint8_t ct[32];
} inke_ct_t;

int keygen(inke_sk_t *sk, inke_pk_t *pk);
int encrypt(inke_ct_t *ct, const inke_pk_t *pk, const unsigned char *m, const size_t m_len);
int decrypt(unsigned char *m, size_t *m_len, const inke_ct_t *ct, const inke_sk_t *sk);