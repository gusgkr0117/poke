#include "poke.h"
#include <time.h>
#include <stdlib.h>
#include <rng.h>
#include <bench.h>

int test_poke_kem(int bench_loops) {
    poke_sk_t sk = {0};
    poke_pk_t pk = {0};
    poke_ct_t ct = {0};
    unsigned char encap_key[32] = {0};
    unsigned char decap_key[32] = {0};
    size_t m_len = 0;
    uint64_t cycles1, cycles2;
    uint64_t cycle_runs[3] = {0};
    uint8_t s[32] = {0};

    for(int i = 0; i < bench_loops; i++) {
        cycles1 = cpucycles();
        randombytes(s, 32);
        keygen(&sk, &pk);
        cycles2 = cpucycles();
        cycle_runs[0] += cycles2 - cycles1;
        
        cycles1 = cpucycles();
        encaps(encap_key, &ct, &pk);
        cycles2 = cpucycles();
        cycle_runs[1] += cycles2 - cycles1;

        cycles1 = cpucycles();
        decaps(decap_key, &ct, &pk, &sk, s);
        cycles2 = cpucycles();
        cycle_runs[2] += cycles2 - cycles1;

        for (int j = 0; j < sizeof(encap_key); j++) {
            if (encap_key[j] != decap_key[j]) return 1;
        }
    }

    printf("test loops : %d\n", bench_loops);
    printf("  gen takes .................................... %.6f %s\n",
            (double)(cycle_runs[0])/(bench_loops), BENCH_UNITS);
    printf("  ecaps takes .................................... %.6f %s\n",
            (double)(cycle_runs[1])/(bench_loops), BENCH_UNITS);
    printf("  decaps takes .................................... %.6f %s\n",
            (double)(cycle_runs[2])/(bench_loops), BENCH_UNITS);

    return 0;
}

int test_poke(int bench_loops) {
    poke_sk_t sk = {0};
    poke_pk_t pk = {0};
    poke_ct_t ct = {0};
    unsigned char m[32] = {0};
    unsigned char dec_m[32] = {0};
    size_t m_len = 0;
    uint64_t cycles1, cycles2;
    uint64_t cycle_runs[3] = {0};

    for(int i = 0; i < bench_loops; i++) {
        randombytes(m, 32);
        cycles1 = cpucycles();
        keygen(&sk, &pk);
        cycles2 = cpucycles();
        cycle_runs[0] += cycles2 - cycles1;
        
        cycles1 = cpucycles();
        encrypt(&ct, &pk, m, 32, NULL, 0);
        cycles2 = cpucycles();
        cycle_runs[1] += cycles2 - cycles1;

        cycles1 = cpucycles();
        decrypt(dec_m, &m_len, &ct, &sk);
        cycles2 = cpucycles();
        cycle_runs[2] += cycles2 - cycles1;

        for (int j = 0; j < sizeof(m); j++) {
            if (m[j] != dec_m[j]) return 1;
        }
    }

    printf("test loops : %d\n", bench_loops);
    printf("  keygen takes .................................... %.6f %s\n",
            (double)(cycle_runs[0])/(bench_loops), BENCH_UNITS);
    printf("  encrypt takes .................................... %.6f %s\n",
            (double)(cycle_runs[1])/(bench_loops), BENCH_UNITS);
    printf("  decrypt takes .................................... %.6f %s\n",
            (double)(cycle_runs[2])/(bench_loops), BENCH_UNITS);

    return 0;
}

int main(int argc, char* argv[]) {
    int loops = 10;
    if (argc > 1) {
        loops = atoi(argv[1]);
    }
    return test_poke_kem(loops);
} 