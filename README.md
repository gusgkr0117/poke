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

13: Test command: poke/ref/poke_lvl1/test/poke-test_poke_lvl1 "100"
13: Working Directory: poke/build/src/poke/ref/poke_lvl1/test
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 53593587.200000 nsec
13:   encrypt takes .................................... 12392670.720000 nsec
13:   decrypt takes .................................... 17681781.760000 nsec
1/9 Test #13: poke-test_poke_lvl1 ..............   Passed    8.85 sec
test 16
    Start 16: poke-test_poke_lvl3

16: Test command: poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
16: Working Directory: poke/build/src/poke/ref/poke_lvl3/test
16: Test timeout computed to be: 1500
16: test loops : 100
16:   keygen takes .................................... 191975690.240000 nsec
16:   encrypt takes .................................... 44050311.680000 nsec
16:   decrypt takes .................................... 64517629.440000 nsec
2/9 Test #16: poke-test_poke_lvl3 ..............   Passed   30.30 sec
test 19
    Start 19: poke-test_poke_lvl5

19: Test command: poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
19: Working Directory: poke/build/src/poke/ref/poke_lvl5/test
19: Test timeout computed to be: 1500
19: test loops : 100
19:   keygen takes .................................... 392021795.840000 nsec
19:   encrypt takes .................................... 89838517.760000 nsec
19:   decrypt takes .................................... 130401054.720000 nsec
3/9 Test #19: poke-test_poke_lvl5 ..............   Passed   61.49 sec
test 22
    Start 22: poke-test_inke_lvl1

22: Test command: poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
22: Working Directory: poke/build/src/poke/ref/inke_lvl1/test
22: Test timeout computed to be: 1500
22: test loops : 100
22:   keygen takes .................................... 44252917.760000 nsec
22:   encrypt takes .................................... 11813091.840000 nsec
22:   decrypt takes .................................... 8816678.400000 nsec
4/9 Test #22: poke-test_inke_lvl1 ..............   Passed    6.86 sec
test 25
    Start 25: poke-test_inke_lvl3

25: Test command: poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
25: Working Directory: poke/build/src/poke/ref/inke_lvl3/test
25: Test timeout computed to be: 1500
25: test loops : 100
25:   keygen takes .................................... 131215429.120000 nsec
25:   encrypt takes .................................... 36028861.440000 nsec
25:   decrypt takes .................................... 26732697.600000 nsec
5/9 Test #25: poke-test_inke_lvl3 ..............   Passed   19.63 sec
test 28
    Start 28: poke-test_inke_lvl5

28: Test command: poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
28: Working Directory: poke/build/src/poke/ref/inke_lvl5/test
28: Test timeout computed to be: 1500
28: test loops : 100
28:   keygen takes .................................... 281477163.520000 nsec
28:   encrypt takes .................................... 76106129.920000 nsec
28:   decrypt takes .................................... 56573406.720000 nsec
6/9 Test #28: poke-test_inke_lvl5 ..............   Passed   41.65 sec
test 31
    Start 31: poke-test_pike_lvl1

31: Test command: poke/ref/pike_lvl1/test/poke-test_pike_lvl1 "100"
31: Working Directory: poke/build/src/poke/ref/pike_lvl1/test
31: Test timeout computed to be: 1500
31: test loops : 100
31:   keygen takes .................................... 57891683.840000 nsec
31:   encrypt takes .................................... 11052577.280000 nsec
31:   decrypt takes .................................... 9862195.200000 nsec
7/9 Test #31: poke-test_pike_lvl1 ..............   Passed    8.25 sec
test 33
    Start 33: poke-test_pike_lvl3

33: Test command: poke/ref/pike_lvl3/test/poke-test_pike_lvl3 "100"
33: Working Directory: poke/build/src/poke/ref/pike_lvl3/test
33: Test timeout computed to be: 1500
33: test loops : 100
33:   keygen takes .................................... 177139399.680000 nsec
33:   encrypt takes .................................... 34158737.920000 nsec
33:   decrypt takes .................................... 30239623.680000 nsec
8/9 Test #33: poke-test_pike_lvl3 ..............   Passed   24.40 sec
test 35
    Start 35: poke-test_pike_lvl5

35: Test command: poke/ref/pike_lvl5/test/poke-test_pike_lvl5 "100"
35: Working Directory: poke/build/src/poke/ref/pike_lvl5/test
35: Test timeout computed to be: 1500
35: test loops : 100
35:   keygen takes .................................... 376436930.560000 nsec
35:   encrypt takes .................................... 72581629.440000 nsec
35:   decrypt takes .................................... 63933672.960000 nsec
9/9 Test #35: poke-test_pike_lvl5 ..............   Passed   51.57 sec

The following tests passed:
	poke-test_poke_lvl1
	poke-test_poke_lvl3
	poke-test_poke_lvl5
	poke-test_inke_lvl1
	poke-test_inke_lvl3
	poke-test_inke_lvl5
	poke-test_pike_lvl1
	poke-test_pike_lvl3
	poke-test_pike_lvl5

100% tests passed, 0 tests failed out of 9

Total Test time (real) = 253.03 sec
```

# Reference
+ [INKE paper](https://eprint.iacr.org/2025/1458)
+ [PIKE paper](https://eprint.iacr.org/2026/473) and [C-implementation](https://github.com/Kaizhan-Lin/PIKE-C-Implementation)
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)