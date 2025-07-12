#include <poke.h>
#include <time.h>
#include <stdlib.h>
#include <rng.h>

#define BENCH_LOOPS 10

static inline int64_t
cpucycles(void)
{
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
}

int test_poke() {
    poke_sk_t sk = {0};
    poke_pk_t pk = {0};
    poke_ct_t ct = {0};
    unsigned char m[32] = {0};
    unsigned char dec_m[32] = {0};
    size_t m_len = 0;
    uint64_t cycles1, cycles2;
    uint64_t cycle_runs[3] = {0};

    for(int i = 0; i < BENCH_LOOPS; i++) {
        randombytes(m, 32);
        cycles1 = cpucycles();
        keygen(&sk, &pk);
        cycles2 = cpucycles();
        cycle_runs[0] += cycles2 - cycles1;
        
        cycles1 = cpucycles();
        encrypt(&ct, &pk, m, 32);
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

    printf("  keygen takes .................................... %.6f msec\n",
            (double)(cycle_runs[0])/(1000000 * BENCH_LOOPS));
    printf("  encrypt takes .................................... %.6f msec\n",
            (double)(cycle_runs[1])/(1000000 * BENCH_LOOPS));
    printf("  decrypt takes .................................... %.6f msec\n",
            (double)(cycle_runs[2])/(1000000 * BENCH_LOOPS));

    return 0;
}

int main(int argc, char* argv[]) {

    return test_poke();
} 