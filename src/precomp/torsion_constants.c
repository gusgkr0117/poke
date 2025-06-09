#include <stddef.h>
#include <stdint.h>
#include <torsion_constants.h>
#if 0
#elif 8*DIGIT_LEN == 16
const uint64_t TORSION_PLUS_EVEN_POWER = 0x81;
const uint64_t TORSION_ODD_PRIMES[4] = {0x3, 0x5, 0x7, 0x1d};
const uint64_t TORSION_ODD_POWERS[4] = {0xa4, 0x12, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[2] = {0x3, 0x5};
const size_t TORSION_PLUS_ODD_POWERS[2] = {0xa4, 0x12};
const uint64_t TORSION_MINUS_ODD_PRIMES[2] = {0x7, 0x1d};
const size_t TORSION_MINUS_ODD_POWERS[2] = {0x1, 0x1};
const size_t DEGREE_COMMITMENT_POWERS[4] = {0xa4, 0x12, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 27, ._mp_d = (mp_limb_t[]) {0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0x9451,0x65f6,0x6b3e,0x3d34,0x99d0,0x1930,0xecff,0x7ef3,0x1477,0xead9,0x93fa,0x56ff,0x17a5,0xc50d,0x4673,0xc612,0x690,0xf418,0x6a0b}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 20, ._mp_d = (mp_limb_t[]) {0xce83,0x6d43,0x853d,0x4448,0xfc10,0x7cc4,0x7724,0xd5b6,0x9d60,0x1d11,0x2bfc,0xfe55,0xe00d,0xa0b0,0xeee6,0x883e,0x9a66,0x4786,0xbbd,0x2a}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[4] = {{{._mp_alloc = 0, ._mp_size = 17, ._mp_d = (mp_limb_t[]) {0x60d1,0xbfb7,0x395,0x3bcc,0x4bc4,0x4453,0x134e,0x42d0,0x628d,0xe45c,0x8d97,0x8ee9,0xb558,0x63e5,0xb34b,0x486e,0xf}}}, {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xe9d9,0x2dac,0x378}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1d}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 19, ._mp_d = (mp_limb_t[]) {0x4a29,0x32fb,0x359f,0x1e9a,0x4ce8,0x8c98,0xf67f,0xbf79,0x8a3b,0x756c,0xc9fd,0xab7f,0x8bd2,0xe286,0x2339,0x6309,0x348,0xfa0c,0x3505}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xcb}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x2}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 17, ._mp_d = (mp_limb_t[]) {0x60d1,0xbfb7,0x395,0x3bcc,0x4bc4,0x4453,0x134e,0x42d0,0x628d,0xe45c,0x8d97,0x8ee9,0xb558,0x63e5,0xb34b,0x486e,0xf}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 25, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0xc1a2,0x7f6e,0x72b,0x7798,0x9788,0x88a6,0x269c,0x85a0,0xc51a,0xc8b8,0x1b2f,0x1dd3,0x6ab1,0xc7cb,0x6696,0x90dd,0x1e}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 20, ._mp_d = (mp_limb_t[]) {0xce83,0x6d43,0x853d,0x4448,0xfc10,0x7cc4,0x7724,0xd5b6,0x9d60,0x1d11,0x2bfc,0xfe55,0xe00d,0xa0b0,0xeee6,0x883e,0x9a66,0x4786,0xbbd,0x2a}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 19, ._mp_d = (mp_limb_t[]) {0x4a29,0x32fb,0x359f,0x1e9a,0x4ce8,0x8c98,0xf67f,0xbf79,0x8a3b,0x756c,0xc9fd,0xab7f,0x8bd2,0xe286,0x2339,0x6309,0x348,0xfa0c,0x3505}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xcb}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x2}}};
#elif 8*DIGIT_LEN == 32
const uint64_t TORSION_PLUS_EVEN_POWER = 0x81;
const uint64_t TORSION_ODD_PRIMES[4] = {0x3, 0x5, 0x7, 0x1d};
const uint64_t TORSION_ODD_POWERS[4] = {0xa4, 0x12, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[2] = {0x3, 0x5};
const size_t TORSION_PLUS_ODD_POWERS[2] = {0xa4, 0x12};
const uint64_t TORSION_MINUS_ODD_PRIMES[2] = {0x7, 0x1d};
const size_t TORSION_MINUS_ODD_POWERS[2] = {0x1, 0x1};
const size_t DEGREE_COMMITMENT_POWERS[4] = {0xa4, 0x12, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 14, ._mp_d = (mp_limb_t[]) {0xffffffff,0xffffffff,0xffffffff,0xffffffff,0x65f69451,0x3d346b3e,0x193099d0,0x7ef3ecff,0xead91477,0x56ff93fa,0xc50d17a5,0xc6124673,0xf4180690,0x6a0b}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 10, ._mp_d = (mp_limb_t[]) {0x6d43ce83,0x4448853d,0x7cc4fc10,0xd5b67724,0x1d119d60,0xfe552bfc,0xa0b0e00d,0x883eeee6,0x47869a66,0x2a0bbd}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[4] = {{{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0xbfb760d1,0x3bcc0395,0x44534bc4,0x42d0134e,0xe45c628d,0x8ee98d97,0x63e5b558,0x486eb34b,0xf}}}, {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x2dace9d9,0x378}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1d}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 10, ._mp_d = (mp_limb_t[]) {0x32fb4a29,0x1e9a359f,0x8c984ce8,0xbf79f67f,0x756c8a3b,0xab7fc9fd,0xe2868bd2,0x63092339,0xfa0c0348,0x3505}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xcb}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x2}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0xbfb760d1,0x3bcc0395,0x44534bc4,0x42d0134e,0xe45c628d,0x8ee98d97,0x63e5b558,0x486eb34b,0xf}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 13, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x7f6ec1a2,0x7798072b,0x88a69788,0x85a0269c,0xc8b8c51a,0x1dd31b2f,0xc7cb6ab1,0x90dd6696,0x1e}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 10, ._mp_d = (mp_limb_t[]) {0x6d43ce83,0x4448853d,0x7cc4fc10,0xd5b67724,0x1d119d60,0xfe552bfc,0xa0b0e00d,0x883eeee6,0x47869a66,0x2a0bbd}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 10, ._mp_d = (mp_limb_t[]) {0x32fb4a29,0x1e9a359f,0x8c984ce8,0xbf79f67f,0x756c8a3b,0xab7fc9fd,0xe2868bd2,0x63092339,0xfa0c0348,0x3505}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xcb}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x2}}};
#elif 8*DIGIT_LEN == 64
const uint64_t TORSION_PLUS_EVEN_POWER = 0x81;
const uint64_t TORSION_ODD_PRIMES[4] = {0x3, 0x5, 0x7, 0x1d};
const uint64_t TORSION_ODD_POWERS[4] = {0xa4, 0x12, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[2] = {0x3, 0x5};
const size_t TORSION_PLUS_ODD_POWERS[2] = {0xa4, 0x12};
const uint64_t TORSION_MINUS_ODD_PRIMES[2] = {0x7, 0x1d};
const size_t TORSION_MINUS_ODD_POWERS[2] = {0x1, 0x1};
const size_t DEGREE_COMMITMENT_POWERS[4] = {0xa4, 0x12, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 7, ._mp_d = (mp_limb_t[]) {0xffffffffffffffff,0xffffffffffffffff,0x3d346b3e65f69451,0x7ef3ecff193099d0,0x56ff93faead91477,0xc6124673c50d17a5,0x6a0bf4180690}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x4448853d6d43ce83,0xd5b677247cc4fc10,0xfe552bfc1d119d60,0x883eeee6a0b0e00d,0x2a0bbd47869a66}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[4] = {{{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x3bcc0395bfb760d1,0x42d0134e44534bc4,0x8ee98d97e45c628d,0x486eb34b63e5b558,0xf}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3782dace9d9}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1d}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x1e9a359f32fb4a29,0xbf79f67f8c984ce8,0xab7fc9fd756c8a3b,0x63092339e2868bd2,0x3505fa0c0348}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xcb}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x2}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x3bcc0395bfb760d1,0x42d0134e44534bc4,0x8ee98d97e45c628d,0x486eb34b63e5b558,0xf}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 7, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x7798072b7f6ec1a2,0x85a0269c88a69788,0x1dd31b2fc8b8c51a,0x90dd6696c7cb6ab1,0x1e}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x4448853d6d43ce83,0xd5b677247cc4fc10,0xfe552bfc1d119d60,0x883eeee6a0b0e00d,0x2a0bbd47869a66}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0x1e9a359f32fb4a29,0xbf79f67f8c984ce8,0xab7fc9fd756c8a3b,0x63092339e2868bd2,0x3505fa0c0348}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xcb}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x2}}};
#endif
