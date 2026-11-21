# POKÉ, INKE and PIKE

This repository contains the C implementations used to compare key generation, encryption and decryption for the three isogeny-based public-key encryption schemes.

| Scheme | Source files and entry points used in the comparison |
| --- | --- |
| POKÉ | [poke.c](src/poke/ref/pokex/poke.c) contains `keygen`, `encrypt` and `decrypt`. |
| INKE | [inke.c](src/poke/ref/inkex/inke.c) contains `b_complete_keygen_init`, `b_complete_keygen_step` and `b_complete_decrypt`. [inke_keygen.c](src/poke/ref/inkex/inke_keygen.c) contains `inke_encrypt_uniform` and the code that finishes constructing the public key. |
| PIKE | [pike.c](src/poke/ref/pikex/pike.c) contains `keygen`, `encrypt` and `decrypt`. |

The current INKE implementation calls the two-dimensional isogeny routines in
[complete/](src/poke/ref/inkex/complete/). The earlier INKE implementation is retained in
[inke_legacy.c](src/poke/ref/inkex/inke_legacy.c) for the original parameter sets.

The timing loops show which functions were called for [INKE](bench/comparison.c) and for [POKÉ and PIKE](bench/comparison-baseline.c). The latter two used an earlier revision of the timing driver. The original project's license and attribution are in [LICENSE](LICENSE) and [NOTICE](NOTICE).

# How-to-use

## Install Clang

The reported timings were obtained on Windows x86-64 using llvm-mingw Clang 23.1.5 (UCRT), CMake 4.5.0, Ninja 1.13.2, GMP 6.3.0 and FLINT 3.6.0. Use GMP and FLINT built for the same UCRT toolchain. Add the directories containing `clang`, `cmake`, `ninja` and the FLINT runtime DLLs to `PATH`.

## Run POKÉ, INKE and PIKE

From the repository root, configure the comparison and build the three level-I executables. The example uses GMP and FLINT installed in `C:/deps/ucrt64`. Adjust those paths to match your installation.

```sh
cmake -S . -B build -G Ninja -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DINKE_COMPLETE_BENCHMARKS=ON -DGMP=C:/deps/ucrt64/lib/libgmp.a -DGMP_INCLUDE=C:/deps/ucrt64/include -DFLINT_LIBRARY=C:/deps/ucrt64/lib/libflint.dll.a -DFLINT_INCLUDE_DIR=C:/deps/ucrt64/include
cmake --build build --target compare_poke_lvl1 compare_complete_inke_lvl1 compare_pike_lvl1
```

Each command below generates a warm-up key, then times one new key and checks that decryption recovers the 32-byte message. The arguments give the number of timed keys and the key-generation time limit in seconds. To build and run levels III and V, use `lvl3` and `lvl5` in the target and executable names.

```sh
./build/bench/compare_poke_lvl1.exe 1 86400
./build/bench/compare_complete_inke_lvl1.exe 1 86400
./build/bench/compare_pike_lvl1.exe 1 86400
```

## Output

The executables write one CSV row per timed key to standard output and write correctness checks to standard error. The paper reports the following mean times in milliseconds on a single pinned thread of an Intel Core i7-9700. Each entry is the mean for twenty keys across two blocks of ten, with a separate warm-up per block.

| NIST level | Scheme | Key generation (ms) | Encryption (ms) | Decryption (ms) |
| --- | --- | ---: | ---: | ---: |
| I | POKÉ | 313.0 | 39.0 | 63.1 |
| I | INKE | 651.0 | 88.2 | 133.3 |
| I | PIKE | 402.9 | 40.5 | 42.0 |
| III | POKÉ | 1191.7 | 181.2 | 286.6 |
| III | INKE | 2068.0 | 274.8 | 415.7 |
| III | PIKE | 1305.8 | 149.8 | 130.5 |
| V | POKÉ | 2196.8 | 348.3 | 602.5 |
| V | INKE | 4011.2 | 731.5 | 906.7 |
| V | PIKE | 1652.0 | 315.6 | 216.1 |

# Reference

- [POKÉ paper](https://eprint.iacr.org/2024/624)
- [PIKE paper](https://eprint.iacr.org/2026/473) and [reference implementation](https://github.com/Kaizhan-Lin/PIKE-C-Implementation)
- [SQIsign implementation](https://github.com/SQISign/sqisign2d-west-ac24)
- [Fiat Cryptography](https://github.com/mit-plv/fiat-crypto)
