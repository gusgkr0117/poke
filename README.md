# 🥗POKÉ : POint-based Key Exchange and ✒️INKE : INtermediate-curve-based Key Exchange
+ C-Implementation of the most efficient isogeny-based PKE protocols
+ Clang should be used

# How-to-use
## Install Clang
```bash
sudo bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"
```

## Run POKE
```bash
cd build
cmake -DCMAKE_C_COMPILER=clang ..
make precomp
make
ctest -V -R "poke-test*"
```

## Output
Run in Apple M2 CPU 3.49GHz
```
test 13
    Start 13: poke-test_poke_lvl1

13: Test command: /poke/build/src/poke/ref/poke_lvl1/test/poke-test_poke_lvl1 "100"
13: Working Directory: /poke/build/src/poke/ref/poke_lvl1/test
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 163999764.480000 nsec
13:   encrypt takes .................................... 23237585.920000 nsec
13:   decrypt takes .................................... 42991874.560000 nsec
1/6 Test #13: poke-test_poke_lvl1 ..............   Passed   23.03 sec
test 14
    Start 14: poke-test_poke_lvl3

14: Test command: /poke/build/src/poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
14: Working Directory: /poke/build/src/poke/ref/poke_lvl3/test
14: Test timeout computed to be: 1500
14: test loops : 100
14:   keygen takes .................................... 560542681.600000 nsec
14:   encrypt takes .................................... 78860666.880000 nsec
14:   decrypt takes .................................... 158033397.760000 nsec
2/6 Test #14: poke-test_poke_lvl3 ..............   Passed   80.11 sec
test 15
    Start 15: poke-test_poke_lvl5

15: Test command: /poke/build/src/poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
15: Working Directory: /poke/build/src/poke/ref/poke_lvl5/test
15: Test timeout computed to be: 1500
15: test loops : 100
15:   keygen takes .................................... 1149711866.880000 nsec
15:   encrypt takes .................................... 161327846.400000 nsec
15:   decrypt takes .................................... 310472757.760000 nsec
3/6 Test #15: poke-test_poke_lvl5 ..............   Passed  162.53 sec
test 16
    Start 16: poke-test_inke_lvl1

16: Test command: /poke/build/src/poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
16: Working Directory: /poke/build/src/poke/ref/inke_lvl1/test
16: Test timeout computed to be: 1500
16: test loops : 100
16:   keygen takes .................................... 124150205.440000 nsec
16:   encrypt takes .................................... 17316300.800000 nsec
16:   decrypt takes .................................... 9882741.760000 nsec
4/6 Test #16: poke-test_inke_lvl1 ..............   Passed   15.14 sec
test 17
    Start 17: poke-test_inke_lvl3

17: Test command: /poke/build/src/poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
17: Working Directory: /poke/build/src/poke/ref/inke_lvl3/test
17: Test timeout computed to be: 1500
17: test loops : 100
17:   keygen takes .................................... 378664624.640000 nsec
17:   encrypt takes .................................... 53339699.200000 nsec
17:   decrypt takes .................................... 29989557.760000 nsec
5/6 Test #17: poke-test_inke_lvl3 ..............   Passed   46.55 sec
test 18
    Start 18: poke-test_inke_lvl5

18: Test command: /poke/build/src/poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
18: Working Directory: /poke/build/src/poke/ref/inke_lvl5/test
18: Test timeout computed to be: 1500
18: test loops : 100
18:   keygen takes .................................... 808202485.760000 nsec
18:   encrypt takes .................................... 115029649.920000 nsec
18:   decrypt takes .................................... 62912957.440000 nsec
6/6 Test #18: poke-test_inke_lvl5 ..............   Passed   98.95 sec

The following tests passed:
        poke-test_poke_lvl1
        poke-test_poke_lvl3
        poke-test_poke_lvl5
        poke-test_inke_lvl1
        poke-test_inke_lvl3
        poke-test_inke_lvl5

100% tests passed, 0 tests failed out of 6

Total Test time (real) = 426.33 sec
```

# Reference
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)