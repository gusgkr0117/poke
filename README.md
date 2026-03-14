# 🥗POKÉ : POint-based Key Exchange and ✒️INKE : INtermediate-curve-based Key Exchange
This implementation is based on [SQIsign](https://github.com/SQISign/sqisign2d-west-ac24) (Apache License 2.0) and has been modified.

(2026.03.13 update) The [PIKE](https://eprint.iacr.org/2026/473) code is optimized and based on [this repository](https://github.com/Kaizhan-Lin/PIKE-C-Implementation).
+ C-Implementation of the most efficient isogeny-based PKE protocols
+ GMP must be installed
+ Clang should be used

# How-to-use
## Install Clang
```bash
sudo bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"
```

## Run POKE
```bash
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release ..
make precomp
make
ctest -V -R "poke-test*" -E "hard"
```
## Output
Run in Apple M2 CPU 3.49GHz
```
test 13
    Start 13: poke-test_poke_lvl1

13: Test command: poke/build/src/poke/ref/poke_lvl1/test/poke-test_poke_lvl1 "100"
13: Working Directory: poke/build/src/poke/ref/poke_lvl1/test
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 71929159.680000 nsec
13:   encrypt takes .................................... 14066419.200000 nsec
13:   decrypt takes .................................... 23130859.520000 nsec
1/6 Test #13: poke-test_poke_lvl1 ..............   Passed   11.47 sec
test 16
    Start 16: poke-test_poke_lvl3

16: Test command: poke/build/src/poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
16: Working Directory: poke/build/src/poke/ref/poke_lvl3/test
16: Test timeout computed to be: 1500
16: test loops : 100
16:   keygen takes .................................... 245184696.320000 nsec
16:   encrypt takes .................................... 49611545.600000 nsec
16:   decrypt takes .................................... 82577244.160000 nsec
2/6 Test #16: poke-test_poke_lvl3 ..............   Passed   37.98 sec
test 19
    Start 19: poke-test_poke_lvl5

19: Test command: poke/build/src/poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
19: Working Directory: poke/build/src/poke/ref/poke_lvl5/test
19: Test timeout computed to be: 1500
19: test loops : 100
19:   keygen takes .................................... 522953984.000000 nsec
19:   encrypt takes .................................... 102804003.840000 nsec
19:   decrypt takes .................................... 168853844.480000 nsec
3/6 Test #19: poke-test_poke_lvl5 ..............   Passed   79.70 sec
test 22
    Start 22: poke-test_inke_lvl1

22: Test command: poke/build/src/poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
22: Working Directory: poke/build/src/poke/ref/inke_lvl1/test
22: Test timeout computed to be: 1500
22: test loops : 100
22:   keygen takes .................................... 67043822.080000 nsec
22:   encrypt takes .................................... 13784990.720000 nsec
22:   decrypt takes .................................... 10592750.080000 nsec
4/6 Test #22: poke-test_inke_lvl1 ..............   Passed    9.48 sec
test 25
    Start 25: poke-test_inke_lvl3

25: Test command: poke/build/src/poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
25: Working Directory: poke/build/src/poke/ref/inke_lvl3/test
25: Test timeout computed to be: 1500
25: test loops : 100
25:   keygen takes .................................... 188886853.120000 nsec
25:   encrypt takes .................................... 43448732.160000 nsec
25:   decrypt takes .................................... 33142425.600000 nsec
5/6 Test #25: poke-test_inke_lvl3 ..............   Passed   26.78 sec
test 28
    Start 28: poke-test_inke_lvl5

28: Test command: poke/build/src/poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
28: Working Directory: poke/build/src/poke/ref/inke_lvl5/test
28: Test timeout computed to be: 1500
28: test loops : 100
28:   keygen takes .................................... 410722214.400000 nsec
28:   encrypt takes .................................... 93010219.520000 nsec
28:   decrypt takes .................................... 71618705.920000 nsec
6/6 Test #28: poke-test_inke_lvl5 ..............   Passed   57.77 sec

The following tests passed:
	poke-test_poke_lvl1
	poke-test_poke_lvl3
	poke-test_poke_lvl5
	poke-test_inke_lvl1
	poke-test_inke_lvl3
	poke-test_inke_lvl5

100% tests passed, 0 tests failed out of 6

Total Test time (real) = 223.19 sec
```

# Reference
+ [INKE paper](https://eprint.iacr.org/2025/1458)
+ [PIKE paper](https://eprint.iacr.org/2026/473) and [C-implementation](https://github.com/Kaizhan-Lin/PIKE-C-Implementation)
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)