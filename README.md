# 🥗POKÉ : POint-based Key Exchange and ✒️INKE : INtermediate-curve-based Key Exchange
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
13:   keygen takes .................................... 134850519.040000 nsec
13:   encrypt takes .................................... 14405957.120000 nsec
13:   decrypt takes .................................... 23145118.720000 nsec
1/6 Test #13: poke-test_poke_lvl1 ..............   Passed   17.63 sec
test 15
    Start 15: poke-test_poke_lvl3

15: Test command: /poke/build/src/poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
15: Working Directory: /poke/build/src/poke/ref/poke_lvl3/test
15: Test timeout computed to be: 1500
15: test loops : 100
15:   keygen takes .................................... 488484400.640000 nsec
15:   encrypt takes .................................... 51317742.080000 nsec
15:   decrypt takes .................................... 91717053.440000 nsec
2/6 Test #15: poke-test_poke_lvl3 ..............   Passed   63.40 sec
test 17
    Start 17: poke-test_poke_lvl5

17: Test command: /poke/build/src/poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
17: Working Directory: /poke/build/src/poke/ref/poke_lvl5/test
17: Test timeout computed to be: 1500
17: test loops : 100
17:   keygen takes .................................... 992289254.400000 nsec
17:   encrypt takes .................................... 103828646.400000 nsec
17:   decrypt takes .................................... 176275458.560000 nsec
3/6 Test #17: poke-test_poke_lvl5 ..............   Passed  127.57 sec
test 19
    Start 19: poke-test_inke_lvl1

19: Test command: /poke/build/src/poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
19: Working Directory: /poke/build/src/poke/ref/inke_lvl1/test
19: Test timeout computed to be: 1500
19: test loops : 100
19:   keygen takes .................................... 110149841.920000 nsec
19:   encrypt takes .................................... 14632714.240000 nsec
19:   decrypt takes .................................... 8628451.840000 nsec
4/6 Test #19: poke-test_inke_lvl1 ..............   Passed   13.73 sec
test 21
    Start 21: poke-test_inke_lvl3

21: Test command: /poke/build/src/poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
21: Working Directory: /poke/build/src/poke/ref/inke_lvl3/test
21: Test timeout computed to be: 1500
21: test loops : 100
21:   keygen takes .................................... 337052848.640000 nsec
21:   encrypt takes .................................... 45076684.800000 nsec
21:   decrypt takes .................................... 26370856.960000 nsec
5/6 Test #21: poke-test_inke_lvl3 ..............   Passed   41.09 sec
test 23
    Start 23: poke-test_inke_lvl5

23: Test command: /poke/build/src/poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
23: Working Directory: /poke/build/src/poke/ref/inke_lvl5/test
23: Test timeout computed to be: 1500
23: test loops : 100
23:   keygen takes .................................... 713062837.760000 nsec
23:   encrypt takes .................................... 96841850.880000 nsec
23:   decrypt takes .................................... 55803056.640000 nsec
6/6 Test #23: poke-test_inke_lvl5 ..............   Passed   86.82 sec

The following tests passed:
        poke-test_poke_lvl1
        poke-test_poke_lvl3
        poke-test_poke_lvl5
        poke-test_inke_lvl1
        poke-test_inke_lvl3
        poke-test_inke_lvl5

100% tests passed, 0 tests failed out of 6

Total Test time (real) = 350.26 sec
```

# Reference
+ [INKE paper](https://eprint.iacr.org/2025/1458)
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)