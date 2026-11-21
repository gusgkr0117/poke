// SPDX-License-Identifier: Apache-2.0
/* Shared, checked PKE measurement code. BENCH_INKE_COMPLETE selects the
 * central-degree fast generation with same-q batches and exact-order PKE APIs.
 * A timeout invalidates the run; it never causes a fresh ideal to be drawn. */
#ifndef _WIN32
#define _POSIX_C_SOURCE 200809L
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <inttypes.h>
#ifndef BENCH_VARIANT
#define BENCH_STRINGIFY_(x) #x
#define BENCH_STRINGIFY(x) BENCH_STRINGIFY_(x)
#define BENCH_VARIANT BENCH_STRINGIFY(BENCH_VARIANT_TOKEN)
#endif
#include <rng.h>
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <cpuid.h>
#include <x86intrin.h>
#define BENCH_X86_TSC 1
#endif
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#if defined(BENCH_INKE)
#include <inke.h>
#if defined(BENCH_INKE_UNIFORM)
#include <inke_keygen.h>
#ifdef BENCH_INKE_COMPLETE
#include "b-complete-keygen.h"
#endif
#endif
typedef inke_sk_t secret_key;
typedef inke_pk_t public_key;
typedef inke_ct_t ciphertext;
#define SCHEME_NAME "INKE"
#elif defined(BENCH_POKE)
#include <poke.h>
typedef poke_sk_t secret_key;
typedef poke_pk_t public_key;
typedef poke_ct_t ciphertext;
#define SCHEME_NAME "POKE"
#elif defined(BENCH_PIKE)
#include <pike.h>
typedef pike_sk_t secret_key;
typedef pike_pk_t public_key;
typedef pike_ct_t ciphertext;
#define SCHEME_NAME "PIKE"
#else
#error Select exactly one scheme
#endif
#include <ec_params.h>

extern const uint64_t p[NWORDS_FIELD];

#if defined(BENCH_INKE_COMPLETE)
#define IMPLEMENTATION_NAME "central_fast_same_q"
static b_complete_public complete_public;
static double public_setup_seconds;
#elif defined(BENCH_INKE_UNIFORM)
#define IMPLEMENTATION_NAME "uniform_ideal_mixed_conversion"
#else
#define IMPLEMENTATION_NAME "upstream_heuristic"
#endif

typedef struct {
    unsigned long windows, prime_candidates, norm_candidates, unit_candidates;
    unsigned long mixed_candidates, witnesses, conversion_failures, kani_failures;
    unsigned long gaussian_attempts,general_attempts,product_fallbacks;
    unsigned long fast_candidates,fast_witnesses,fast_conversion_failures;
    unsigned long fast_guard_misses,fast_interval_misses;
    int product_steps,conversion_origin;
} keygen_statistics;

typedef struct { uint64_t ticks; unsigned int cpu; } timestamp;
static int tsc_available;

static timestamp serialized_timestamp(void)
{
    timestamp result = {0, 0};
#ifdef BENCH_X86_TSC
    if (tsc_available) {
        _mm_lfence();
        result.ticks = __rdtscp(&result.cpu);
        _mm_lfence();
    }
#endif
    return result;
}

static int invariant_rdtscp_available(void)
{
#ifdef BENCH_X86_TSC
    unsigned int a,b,c,d;
    if (__get_cpuid_max(0x80000000, NULL) < 0x80000007) return 0;
    __cpuid(0x80000001, a,b,c,d);
    if (!(d & (1U << 27))) return 0;
    __cpuid(0x80000007, a,b,c,d);
    return (d & (1U << 8)) != 0;
#else
    return 0;
#endif
}

static unsigned int prime_bits(void)
{
    for (int word=NWORDS_FIELD-1;word>=0;--word) {
        uint64_t top=p[word];
        if (top) { unsigned int bits=0; while (top) { ++bits; top>>=1; }
            return 64U*(unsigned int)word+bits; }
    }
    return 0;
}

static double monotonic_seconds(void)
{
#ifdef _WIN32
    LARGE_INTEGER counter, frequency;
    if (!QueryPerformanceFrequency(&frequency) || !QueryPerformanceCounter(&counter)) {
        fputs("Performance timer failed\n", stderr); exit(1);
    }
    return (double)counter.QuadPart / (double)frequency.QuadPart;
#else
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now)) {
        perror("clock_gettime"); exit(1);
    }
    return (double)now.tv_sec + (double)now.tv_nsec * 1e-9;
#endif
}

static int checked_keygen(secret_key *sk,public_key *pk,double deadline,
    keygen_statistics *statistics)
{
#ifdef BENCH_INKE_COMPLETE
    b_complete_keygen ctx={0};
    if(b_complete_keygen_init(&ctx,&complete_public,0,0,0,0)!=INKE_SAMPLE_OK)return 0;
    unsigned long product_before=complete_public.product.calls;int status;
    do {
        ++statistics->windows;status=b_complete_keygen_step(sk,pk,&ctx,NULL,NULL);
        if(monotonic_seconds()>deadline){
            fprintf(stderr,"Keygen time budget exceeded: phase=%d stage=%d; run invalid, original q retained.\n",ctx.key.phase,ctx.key.last_stage);
            status=INKE_SAMPLE_INVALID;break;
        }
        if(status==INKE_SAMPLE_RETRY){fprintf(stderr,"Keygen retry: window=%lu phase=%d stage=%d; retaining original q and any committed receiver map.\n",statistics->windows,ctx.key.phase,ctx.key.last_stage);fflush(stderr);}
    }while(status==INKE_SAMPLE_RETRY);
    statistics->kani_failures=ctx.key.finish_failures;
    statistics->product_fallbacks=complete_public.product.calls-product_before;
    statistics->product_steps=ctx.key.has_conversion?ctx.key.conversion.product_steps:0;
    /* Origin 4 denotes the norm-search receiver in the archived CSV format. */
    statistics->conversion_origin=ctx.key.has_conversion?4:0;
    statistics->fast_candidates=ctx.fast_counts.norm_candidates;
    statistics->fast_witnesses=ctx.fast_counts.witnesses;
    statistics->fast_conversion_failures=ctx.fast_counts.conversion_failures;
    statistics->fast_guard_misses=ctx.fast_counts.guard_misses;
    statistics->fast_interval_misses=ctx.fast_counts.interval_misses;
    if(status!=INKE_SAMPLE_OK)fprintf(stderr,"Fast Keygen failed: status=%d phase=%d stage=%d; original q and committed map retained.\n",status,ctx.key.phase,ctx.key.last_stage);
    b_complete_keygen_clear(&ctx);return status==INKE_SAMPLE_OK;
#elif defined(BENCH_INKE_UNIFORM)
    inke_keygen_context ctx;
    inke_keygen_init(&ctx);
    int status;
    do {
        ++statistics->windows;
        status=inke_keygen_uniform(sk,pk,&ctx,100000);
        if (monotonic_seconds()>deadline) {
            fprintf(stderr,"Keygen time budget exceeded: stage=%d windows=%lu witnesses=%lu; run invalid, input not replaced.\n",
                ctx.last_stage,statistics->windows,ctx.witnesses);
            status=INKE_SAMPLE_INVALID;
            break;
        }
        if (status==INKE_SAMPLE_RETRY) {
            fprintf(stderr,"Keygen retry: window=%lu stage=%d candidates=%lu witnesses=%lu; retaining context.\n",
                statistics->windows,ctx.last_stage,ctx.mixed.attempts,ctx.witnesses);
            fflush(stderr);
        }
    } while (status==INKE_SAMPLE_RETRY);
    statistics->prime_candidates=ctx.sampling.prime_candidates;
    statistics->norm_candidates=ctx.sampling.norm_candidates;
    statistics->unit_candidates=ctx.sampling.unit_candidates;
    statistics->mixed_candidates=ctx.mixed.attempts;
    statistics->witnesses=ctx.witnesses;
    statistics->conversion_failures=ctx.conversion_failures;
    statistics->kani_failures=ctx.kani_failures;
    statistics->product_steps=ctx.product_steps;
    inke_keygen_clear(&ctx);
    return status==INKE_SAMPLE_OK;
#else
    (void)statistics;
    int status=keygen(sk,pk);
    if (monotonic_seconds()>deadline) {
        fputs("Keygen time budget exceeded; run invalid.\n",stderr);
        return 0;
    }
    return status==1;
#endif
}

static void print_metadata(long loops,double budget)
{
    int affinity_cpu=-1;
    int active_five_exponent=POWER_OF_5;
#if defined(BENCH_INKE_UNIFORM) && defined(INKE_MANUSCRIPT_PARAMETERS)
    /* The inherited header retains a placeholder for unused 5-torsion. */
    active_five_exponent=0;
#endif
    unsigned int public_two_exponent=0;
#ifdef BENCH_INKE_UNIFORM
    inke_keygen_context parameter_context;
    inke_keygen_init(&parameter_context);
    const unsigned int kani_exponent=parameter_context.a;
    public_two_exponent=kani_exponent+2;
    inke_keygen_clear(&parameter_context);
#endif
#ifdef _WIN32
    DWORD_PTR allowed,system;
    if (!GetProcessAffinityMask(GetCurrentProcess(),&allowed,&system) || !allowed) {
        fputs("Cannot inspect process affinity\n",stderr); exit(1);
    }
    DWORD_PTR selected=allowed & (~allowed+1);
    if (!SetThreadAffinityMask(GetCurrentThread(),selected)) {
        fputs("Cannot pin benchmark thread\n",stderr); exit(1);
    }
    affinity_cpu=0;
    while (!(selected & 1)) { ++affinity_cpu; selected>>=1; }
#endif
    tsc_available=invariant_rdtscp_available();
    double tsc_hz=0;
    if (tsc_available) {
        double t0=monotonic_seconds(); timestamp c0=serialized_timestamp();
        double t1;
        do { t1=monotonic_seconds(); } while (t1-t0<0.1);
        timestamp c1=serialized_timestamp();
        if (c0.cpu!=c1.cpu || c1.ticks<=c0.ticks) tsc_available=0;
        else tsc_hz=(double)(c1.ticks-c0.ticks)/(t1-t0);
    }
    fprintf(stderr,"BENCH_METADATA {\"scheme\":\"%s\",\"variant\":\"%s\",\"implementation\":\"%s\",\"prime_bits\":%u,\"prime_hex\":\"",
        SCHEME_NAME,BENCH_VARIANT,IMPLEMENTATION_NAME,prime_bits());
    for (int k=NWORDS_FIELD-1;k>=0;--k) fprintf(stderr,"%016" PRIx64,p[k]);
#ifdef BENCH_INKE_UNIFORM
    fprintf(stderr,"\",\"kani_exponent\":%u",kani_exponent);
#else
    fputc('"',stderr);
#endif
    fprintf(stderr,",\"power_of_2\":%d,\"power_of_3\":%d,\"power_of_5\":%d,\"public_two_exponent\":%u,\"iterations\":%ld,\"warmup_iterations\":1,\"keygen_timeout_seconds\":%.3f,\"thread_affinity_cpu\":%d,\"invariant_rdtscp\":%s,\"tsc_hz_measured\":%.3f,\"tsc_unit\":\"reference ticks, not core cycles\",\"compiler\":\"%s\"}\n",
        POWER_OF_2,POWER_OF_3,active_five_exponent,public_two_exponent,loops,budget,affinity_cpu,
        tsc_available?"true":"false",tsc_hz,__VERSION__);
    fflush(stderr);
}

int main(int argc, char **argv)
{
    char *end = NULL;
    errno = 0;
    long loops = argc >= 2 ? strtol(argv[1], &end, 10) : 100;
    if (argc > 3 || errno || (argc >= 2 && (!*argv[1] || *end)) ||
        loops < 1 || loops > 1000000) {
        fputs("Usage: compare_<variant> [iterations <= 1000000] [keygen timeout seconds]\n", stderr);
        return 2;
    }
    errno=0;
    double budget=argc==3?strtod(argv[2],&end):3600;
    if (errno || (argc==3 && (!*argv[2] || *end)) || !(budget>0 && budget<=86400)) return 2;
    print_metadata(loops,budget);
#ifdef BENCH_INKE_COMPLETE
    double setup_start=monotonic_seconds();
    if(b_complete_public_init(&complete_public,root_system_random,NULL,10000)!=INKE_SAMPLE_OK){fputs("Public field setup failed\n",stderr);return 1;}
    public_setup_seconds=monotonic_seconds()-setup_start;
    fprintf(stderr,"BENCH_PUBLIC_SETUP {\"milliseconds\":%.9f,\"fast_candidate_batch\":%lu,\"receiver_distribution\":\"central_q_opposite_parity_fast\",\"receiver_ideal_fallback\":false,\"field_backend\":\"native_fp2\",\"scope\":\"per-parameter field and complete 2D quotient setup; excluded from per-key timer\"}\n",1000*public_setup_seconds,B_FAST_DEFAULT_TRIALS);fflush(stderr);
#endif
    puts("scheme,parameter_set,implementation,prime_bits,iteration,keygen_ms,encrypt_ms,decrypt_ms,keygen_tsc_ticks,encrypt_tsc_ticks,decrypt_tsc_ticks,cpu_migration,keygen_windows,prime_candidates,norm_candidates,unit_candidates,mixed_candidates,witnesses,conversion_failures,kani_failures,product_steps,conversion_origin,gaussian_attempts,general_attempts,keygen_product_fallbacks,decrypt_product_fallbacks,fast_candidates,fast_witnesses,fast_conversion_failures,fast_guard_misses,fast_interval_misses");
    fflush(stdout);
    /* One complete, checked warm-up; all measured iterations are also checked. */
    for (long iteration = -1; iteration < loops; ++iteration) {
        secret_key sk = {0};
        public_key pk = {0};
        ciphertext ct = {0};
        unsigned char message[32], recovered[32] = {0};
        keygen_statistics statistics={0};
        if (randombytes(message, sizeof(message))) return 1;
        double start = monotonic_seconds();
        timestamp key_start=serialized_timestamp();
        int key_ok = checked_keygen(&sk,&pk,start+budget,&statistics);
        timestamp key_end=serialized_timestamp();
        double after_keygen = monotonic_seconds();
        if (key_ok != 1) { fputs("Keygen failed\n", stderr); return 1; }
        double encrypt_start=monotonic_seconds();
        timestamp enc_start=serialized_timestamp();
#ifdef BENCH_INKE_UNIFORM
        int enc_ok=inke_encrypt_uniform(&ct,&pk,message,NULL,0);
#else
        int enc_ok = encrypt(&ct, &pk, message, sizeof(message), NULL, 0);
#endif
        timestamp enc_end=serialized_timestamp();
        double after_encrypt = monotonic_seconds();
        if (enc_ok != 1) { fputs("Encryption failed\n", stderr); return 1; }
        double decrypt_start=monotonic_seconds();
        timestamp dec_start=serialized_timestamp();
        unsigned long decrypt_fallbacks=0;
#ifdef BENCH_INKE_COMPLETE
        unsigned long calls_before=complete_public.product.calls;int decrypt_stage;
        int dec_ok=b_complete_decrypt(recovered,&ct,&sk,&complete_public,&decrypt_stage);
        decrypt_fallbacks=complete_public.product.calls-calls_before;
        size_t recovered_length=sizeof(message);
#elif defined(BENCH_INKE_UNIFORM)
        int dec_ok=inke_decrypt_uniform(recovered,&ct,&sk);
        size_t recovered_length=sizeof(message);
#else
        size_t recovered_length=0;
        int dec_ok = decrypt(recovered, &recovered_length, &ct, &sk);
#endif
        timestamp dec_end=serialized_timestamp();
        double after_decrypt = monotonic_seconds();
        if (dec_ok != 1 || recovered_length != sizeof(message) ||
            memcmp(message, recovered, sizeof(message))) {
            fprintf(stderr, "Round-trip failed at iteration %ld; no valid comparison run\n", iteration);
            return 1;
        }
        if (iteration >= 0) {
            printf("%s,%s,%s,%u,%ld,%.9f,%.9f,%.9f,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%d,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%d,%d,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu\n",
                   SCHEME_NAME, BENCH_VARIANT, IMPLEMENTATION_NAME,prime_bits(),iteration,
                   1000 * (after_keygen - start),
                   1000 * (after_encrypt - encrypt_start),
                   1000 * (after_decrypt - decrypt_start),
                   key_end.ticks-key_start.ticks,enc_end.ticks-enc_start.ticks,dec_end.ticks-dec_start.ticks,
                   key_start.cpu!=key_end.cpu || enc_start.cpu!=enc_end.cpu || dec_start.cpu!=dec_end.cpu,
                   statistics.windows,statistics.prime_candidates,statistics.norm_candidates,
                   statistics.unit_candidates,statistics.mixed_candidates,statistics.witnesses,
                   statistics.conversion_failures,statistics.kani_failures,statistics.product_steps,
                   statistics.conversion_origin,statistics.gaussian_attempts,statistics.general_attempts,
                   statistics.product_fallbacks,decrypt_fallbacks,
                   statistics.fast_candidates,statistics.fast_witnesses,statistics.fast_conversion_failures,
                   statistics.fast_guard_misses,statistics.fast_interval_misses);
            fflush(stdout);
        }
        fprintf(stderr,"Verified iteration %ld of %ld\n",iteration,loops);
        fflush(stderr);
    }
#ifdef BENCH_INKE_COMPLETE
    b_complete_public_clear(&complete_public);
#endif
    return 0;
}
