# 🥗POKÉ : POint-based Key Exchange and ✒️INKE : INtermediate-curve-based Key Exchange
This implementation is based on [SQIsign](https://github.com/SQISign/sqisign2d-west-ac24) (Apache License 2.0) and has been modified.

+ C-Implementation of the most efficient isogeny-based PKE protocols
+ GMP must be installed
+ Clang should be used

# How-to-use
## Install Clang
```bash
sudo bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"
```

## Run INKE/POKE/PIKE
### PKE test
```bash
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release ..
make
ctest -V -R "poke-test*" -E "hard"
```
### KEM test
```bash
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release ..
make
ctest -V -R "poke-kem-test*"
```
## Example output of the PKE test
Run in i7-9700 Coffe Lake CPU 3GHz
```
test 13
    Start 13: poke-test_poke_lvl1

13: Test command: poke/build/src/poke/ref/poke_lvl1/test/poke-test_poke_lvl1 "100"
13: Working Directory: poke/build/src/poke/ref/poke_lvl1/test
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 359055538.670000 cycles
13:   encrypt takes .................................... 83294860.290000 cycles
13:   decrypt takes .................................... 120975449.660000 cycles
1/9 Test #13: poke-test_poke_lvl1 ..............   Passed   18.79 sec
test 16
    Start 16: poke-test_poke_lvl3

16: Test command: poke/build/src/poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
16: Working Directory: poke/build/src/poke/ref/poke_lvl3/test
16: Test timeout computed to be: 1500
16: test loops : 100
16:   keygen takes .................................... 1381593828.220000 cycles
16:   encrypt takes .................................... 323184558.900000 cycles
16:   decrypt takes .................................... 472341775.390000 cycles
2/9 Test #16: poke-test_poke_lvl3 ..............   Passed   72.60 sec
test 19
    Start 19: poke-test_poke_lvl5

19: Test command: poke/build/src/poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
19: Working Directory: poke/build/src/poke/ref/poke_lvl5/test
19: Test timeout computed to be: 1500
19: test loops : 100
19:   keygen takes .................................... 2965943719.350000 cycles
19:   encrypt takes .................................... 688563599.450000 cycles
19:   decrypt takes .................................... 1006279158.090000 cycles
3/9 Test #19: poke-test_poke_lvl5 ..............   Passed  155.39 sec
test 22
    Start 22: poke-test_inke_lvl1

22: Test command: poke/build/src/poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
22: Working Directory: poke/build/src/poke/ref/inke_lvl1/test
22: Test timeout computed to be: 1500
22: test loops : 100
22:   keygen takes .................................... 712084209.780000 cycles
22:   encrypt takes .................................... 161225456.630000 cycles
22:   decrypt takes .................................... 187144204.560000 cycles
4/9 Test #22: poke-test_inke_lvl1 ..............   Passed   35.37 sec
test 25
    Start 25: poke-test_inke_lvl3

25: Test command: poke/build/src/poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
25: Working Directory: poke/build/src/poke/ref/inke_lvl3/test
25: Test timeout computed to be: 1500
25: test loops : 100
25:   keygen takes .................................... 2084186731.820000 cycles
25:   encrypt takes .................................... 469743586.410000 cycles
25:   decrypt takes .................................... 545816541.790000 cycles
5/9 Test #25: poke-test_inke_lvl3 ..............   Passed  103.36 sec
test 28
    Start 28: poke-test_inke_lvl5

28: Test command: poke/build/src/poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
28: Working Directory: poke/build/src/poke/ref/inke_lvl5/test
28: Test timeout computed to be: 1500
28: test loops : 100
28:   keygen takes .................................... 4947984349.820000 cycles
28:   encrypt takes .................................... 1070281654.600000 cycles
28:   decrypt takes .................................... 1244129586.980000 cycles
6/9 Test #28: poke-test_inke_lvl5 ..............   Passed  242.15 sec
test 31
    Start 31: poke-test_pike_lvl1

31: Test command: poke/build/src/poke/ref/pike_lvl1/test/poke-test_pike_lvl1 "100"
31: Working Directory: poke/build/src/poke/ref/pike_lvl1/test
31: Test timeout computed to be: 1500
31: test loops : 100
31:   keygen takes .................................... 375954946.760000 cycles
31:   encrypt takes .................................... 72628156.580000 cycles
31:   decrypt takes .................................... 65325810.010000 cycles
7/9 Test #31: poke-test_pike_lvl1 ..............   Passed   17.15 sec
test 33
    Start 33: poke-test_pike_lvl3

33: Test command: poke/build/src/poke/ref/pike_lvl3/test/poke-test_pike_lvl3 "100"
33: Working Directory: poke/build/src/poke/ref/pike_lvl3/test
33: Test timeout computed to be: 1500
33: test loops : 100
33:   keygen takes .................................... 1269569847.470000 cycles
33:   encrypt takes .................................... 248879607.410000 cycles
33:   decrypt takes .................................... 220747974.730000 cycles
8/9 Test #33: poke-test_pike_lvl3 ..............   Passed   58.00 sec
test 35
    Start 35: poke-test_pike_lvl5

35: Test command: poke/build/src/poke/ref/pike_lvl5/test/poke-test_pike_lvl5 "100"
35: Working Directory: poke/build/src/poke/ref/pike_lvl5/test
35: Test timeout computed to be: 1500
35: test loops : 100
35:   keygen takes .................................... 2789332163.190000 cycles
35:   encrypt takes .................................... 546039505.450000 cycles
35:   decrypt takes .................................... 481278694.740000 cycles
9/9 Test #35: poke-test_pike_lvl5 ..............   Passed  127.26 sec

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

Total Test time (real) = 830.07 sec
```

# Reference
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)