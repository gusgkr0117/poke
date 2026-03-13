#include <stddef.h>
#include <stdint.h>
#include <torsion_constants.h>
#if 0
#elif 8*DIGIT_LEN == 16
const uint64_t TORSION_PLUS_EVEN_POWER = 0x80;
const uint64_t TORSION_ODD_PRIMES[2] = {0x3, 0x5};
const uint64_t TORSION_ODD_POWERS[2] = {0x51, 0x35};
const uint64_t TORSION_PLUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_PLUS_ODD_POWERS[1] = {0x51};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0x5};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x35};
const size_t DEGREE_COMMITMENT_POWERS[3] = {0x51, 0x35, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 25, ._mp_d = (mp_limb_t[]) {0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xcd4c,0xcb8e,0x4ef0,0x3541,0xa5d2,0x704a,0xceb2,0x805e,0x7f56,0x3a38,0x8780,0x262d,0xaca6,0xafd4,0x83a3,0x7bd,0x4}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x311f,0x91ae,0x463,0x3298,0x113f,0xf4b1,0xfcb7,0xb8b0,0xfd30,0xc7d4,0x35fe,0x57ca,0xff63,0xcb4b,0x550e,0xae2}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x2caf,0xba57,0x85da,0xdee5,0xecd9,0x297c,0xfa4,0x17c0,0x3}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[3] = {{{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x7cc3,0xd56d,0xc69,0xb6bf,0xe834,0xa149,0xd5ce,0x4d98,0x1}}}, {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x6475,0xb7f8,0x6da2,0x147a,0x1f04,0x6eb7,0x3636,0x85a}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x7cc3,0xd56d,0xc69,0xb6bf,0xe834,0xa149,0xd5ce,0x4d98,0x1}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x6475,0xb7f8,0x6da2,0x147a,0x1f04,0x6eb7,0x3636,0x85a}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x1}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x7cc3,0xd56d,0xc69,0xb6bf,0xe834,0xa149,0xd5ce,0x4d98,0x1}}};
const ibz_t TORSION_PLUS_CPOWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 17, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x7cc3,0xd56d,0xc69,0xb6bf,0xe834,0xa149,0xd5ce,0x4d98,0x1}}};
const ibz_t TORSION_PLUS_3CPOWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x6fcf,0x2b23,0x3e11,0x91bb,0x8907,0x2671,0x2d09,0x83fc,0x6}}};
const ibz_t TORSION_PLUS_2CPOWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x5}}};
const ibz_t TORSION_PLUS_23CPOWER = {{._mp_alloc = 0, ._mp_size = 17, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x6fcf,0x2b23,0x3e11,0x91bb,0x8907,0x2671,0x2d09,0x83fc,0x6}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x311f,0x91ae,0x463,0x3298,0x113f,0xf4b1,0xfcb7,0xb8b0,0xfd30,0xc7d4,0x35fe,0x57ca,0xff63,0xcb4b,0x550e,0xae2}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x7cc3,0xd56d,0xc69,0xb6bf,0xe834,0xa149,0xd5ce,0x4d98,0x1}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x6475,0xb7f8,0x6da2,0x147a,0x1f04,0x6eb7,0x3636,0x85a}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x1}}};
#elif 8*DIGIT_LEN == 32
const uint64_t TORSION_PLUS_EVEN_POWER = 0x80;
const uint64_t TORSION_ODD_PRIMES[2] = {0x3, 0x5};
const uint64_t TORSION_ODD_POWERS[2] = {0x51, 0x35};
const uint64_t TORSION_PLUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_PLUS_ODD_POWERS[1] = {0x51};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0x5};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x35};
const size_t DEGREE_COMMITMENT_POWERS[3] = {0x51, 0x35, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 13, ._mp_d = (mp_limb_t[]) {0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xcb8ecd4c,0x35414ef0,0x704aa5d2,0x805eceb2,0x3a387f56,0x262d8780,0xafd4aca6,0x7bd83a3,0x4}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x91ae311f,0x32980463,0xf4b1113f,0xb8b0fcb7,0xc7d4fd30,0x57ca35fe,0xcb4bff63,0xae2550e}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xba572caf,0xdee585da,0x297cecd9,0x17c00fa4,0x3}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[3] = {{{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xd56d7cc3,0xb6bf0c69,0xa149e834,0x4d98d5ce,0x1}}}, {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xb7f86475,0x147a6da2,0x6eb71f04,0x85a3636}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xd56d7cc3,0xb6bf0c69,0xa149e834,0x4d98d5ce,0x1}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xb7f86475,0x147a6da2,0x6eb71f04,0x85a3636}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x1}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xd56d7cc3,0xb6bf0c69,0xa149e834,0x4d98d5ce,0x1}}};
const ibz_t TORSION_PLUS_CPOWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0xd56d7cc3,0xb6bf0c69,0xa149e834,0x4d98d5ce,0x1}}};
const ibz_t TORSION_PLUS_3CPOWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x2b236fcf,0x91bb3e11,0x26718907,0x83fc2d09,0x6}}};
const ibz_t TORSION_PLUS_2CPOWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x5}}};
const ibz_t TORSION_PLUS_23CPOWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x2b236fcf,0x91bb3e11,0x26718907,0x83fc2d09,0x6}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x91ae311f,0x32980463,0xf4b1113f,0xb8b0fcb7,0xc7d4fd30,0x57ca35fe,0xcb4bff63,0xae2550e}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xd56d7cc3,0xb6bf0c69,0xa149e834,0x4d98d5ce,0x1}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xb7f86475,0x147a6da2,0x6eb71f04,0x85a3636}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x1}}};
#elif 8*DIGIT_LEN == 64
const uint64_t TORSION_PLUS_EVEN_POWER = 0x80;
const uint64_t TORSION_ODD_PRIMES[2] = {0x3, 0x5};
const uint64_t TORSION_ODD_POWERS[2] = {0x51, 0x35};
const uint64_t TORSION_PLUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_PLUS_ODD_POWERS[1] = {0x51};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0x5};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x35};
const size_t DEGREE_COMMITMENT_POWERS[3] = {0x51, 0x35, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 7, ._mp_d = (mp_limb_t[]) {0xffffffffffffffff,0xffffffffffffffff,0x35414ef0cb8ecd4c,0x805eceb2704aa5d2,0x262d87803a387f56,0x7bd83a3afd4aca6,0x4}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x3298046391ae311f,0xb8b0fcb7f4b1113f,0x57ca35fec7d4fd30,0xae2550ecb4bff63}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xdee585daba572caf,0x17c00fa4297cecd9,0x3}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[3] = {{{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xb6bf0c69d56d7cc3,0x4d98d5cea149e834,0x1}}}, {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x147a6da2b7f86475,0x85a36366eb71f04}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xb6bf0c69d56d7cc3,0x4d98d5cea149e834,0x1}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x147a6da2b7f86475,0x85a36366eb71f04}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x1}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xb6bf0c69d56d7cc3,0x4d98d5cea149e834,0x1}}};
const ibz_t TORSION_PLUS_CPOWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x0,0x0,0xb6bf0c69d56d7cc3,0x4d98d5cea149e834,0x1}}};
const ibz_t TORSION_PLUS_3CPOWER = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x91bb3e112b236fcf,0x83fc2d0926718907,0x6}}};
const ibz_t TORSION_PLUS_2CPOWER = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x5}}};
const ibz_t TORSION_PLUS_23CPOWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x91bb3e112b236fcf,0x83fc2d0926718907,0x6}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x3298046391ae311f,0xb8b0fcb7f4b1113f,0x57ca35fec7d4fd30,0xae2550ecb4bff63}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xb6bf0c69d56d7cc3,0x4d98d5cea149e834,0x1}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x147a6da2b7f86475,0x85a36366eb71f04}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x1}}};
#endif
