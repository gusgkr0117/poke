# POKÉ, INKE and PIKE

This repository contains C implementations of the three isogeny-based public-key encryption schemes and benchmarks for key generation, encryption and decryption.

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

Install Clang, CMake, Ninja and GMP. On Windows, use GMP built for the same UCRT toolchain as Clang. The complete INKE quotient uses the project's own quadratic finite-field arithmetic.

## Run POKÉ, INKE and PIKE

From the repository root, configure the comparison and build the three level-I executables. The example uses GMP installed in `C:/deps/ucrt64`. Adjust the paths for your GMP installation.

```sh
cmake -S . -B build -G Ninja -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DINKE_COMPLETE_BENCHMARKS=ON -DGMP=C:/deps/ucrt64/lib/libgmp.a -DGMP_INCLUDE=C:/deps/ucrt64/include
cmake --build build --target compare_poke_lvl1 compare_complete_inke_lvl1 compare_pike_lvl1
```

Each command below generates a warm-up key, then times one new key and checks that decryption recovers the 32-byte message. The arguments give the number of timed keys and the key-generation time limit in seconds. To build and run levels III and V, use `lvl3` and `lvl5` in the target and executable names.

```sh
./build/bench/compare_poke_lvl1.exe 1 86400
./build/bench/compare_complete_inke_lvl1.exe 1 86400
./build/bench/compare_pike_lvl1.exe 1 86400
```

To check the exceptional cases, build and run the following tests. Each checks both factor labels and their torsion images on a known quotient, then forces a full key generation and decryption through the complete quotient.

```sh
cmake --build build --target test_complete_quotient_inke_lvl1 test_complete_quotient_inke_lvl3 test_complete_quotient_inke_lvl5
ctest --test-dir build -R test_complete_quotient --output-on-failure
```

## Output

The executables write one CSV row per timed key to standard output and correctness checks to standard error. The table records the paper's measurements on a single pinned thread of an Intel Core i7-9700 using llvm-mingw Clang 23.1.5, CMake 4.5.0, Ninja 1.13.2 and GMP 6.3.0. The measured revision used FLINT 3.6.0 for the complete quotient. All sixty timed INKE keys used the fast formulas, and public setup was timed separately. Each entry is the mean in milliseconds for twenty keys across two blocks of ten, with a separate warm-up per block.

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
