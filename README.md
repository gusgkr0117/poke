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
13:   keygen takes .................................... 146915607.040000 nsec
13:   encrypt takes .................................... 22424665.600000 nsec
13:   decrypt takes .................................... 35386342.400000 nsec
1/6 Test #13: poke-test_poke_lvl1 ..............   Passed   20.81 sec
test 14
    Start 14: poke-test_poke_lvl3

14: Test command: /poke/build/src/poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
14: Working Directory: /poke/build/src/poke/ref/poke_lvl3/test
14: Test timeout computed to be: 1500
14: test loops : 100
14:   keygen takes .................................... 520169525.760000 nsec
14:   encrypt takes .................................... 78447600.640000 nsec
14:   decrypt takes .................................... 132346401.280000 nsec
2/6 Test #14: poke-test_poke_lvl3 ..............   Passed   73.33 sec
test 15
    Start 15: poke-test_poke_lvl5

15: Test command: /poke/build/src/poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
15: Working Directory: /poke/build/src/poke/ref/poke_lvl5/test
15: Test timeout computed to be: 1500
15: test loops : 100
15:   keygen takes .................................... 1093837680.640000 nsec
15:   encrypt takes .................................... 161229027.840000 nsec
15:   decrypt takes .................................... 261533240.320000 nsec
3/6 Test #15: poke-test_poke_lvl5 ..............   Passed  152.03 sec
test 16
    Start 16: poke-test_inke_lvl1

16: Test command: /poke/build/src/poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
16: Working Directory: /poke/build/src/poke/ref/inke_lvl1/test
16: Test timeout computed to be: 1500
16: test loops : 100
16:   keygen takes .................................... 123087206.400000 nsec
16:   encrypt takes .................................... 17271086.080000 nsec
16:   decrypt takes .................................... 9744752.640000 nsec
4/6 Test #16: poke-test_inke_lvl1 ..............   Passed   15.02 sec
test 17
    Start 17: poke-test_inke_lvl3

17: Test command: /poke/build/src/poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
17: Working Directory: /poke/build/src/poke/ref/inke_lvl3/test
17: Test timeout computed to be: 1500
17: test loops : 100
17:   keygen takes .................................... 375681758.720000 nsec
17:   encrypt takes .................................... 53056145.920000 nsec
17:   decrypt takes .................................... 29667691.520000 nsec
5/6 Test #17: poke-test_inke_lvl3 ..............   Passed   45.85 sec
test 18
    Start 18: poke-test_inke_lvl5

18: Test command: /poke/build/src/poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
18: Working Directory: /poke/build/src/poke/ref/inke_lvl5/test
18: Test timeout computed to be: 1500
18: test loops : 100
18:   keygen takes .................................... 805159116.800000 nsec
18:   encrypt takes .................................... 112560581.120000 nsec
18:   decrypt takes .................................... 62083266.560000 nsec
6/6 Test #18: poke-test_inke_lvl5 ..............   Passed   97.99 sec

The following tests passed:
	poke-test_poke_lvl1
	poke-test_poke_lvl3
	poke-test_poke_lvl5
	poke-test_inke_lvl1
	poke-test_inke_lvl3
	poke-test_inke_lvl5

100% tests passed, 0 tests failed out of 6

Total Test time (real) = 405.03 sec
```

# Reference
+ [INKE paper](https://eprint.iacr.org/2025/1458)
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)