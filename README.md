# 🥗POKÉ : POint-based Key Exchange
+ C-Implementation of the most efficient isogeny-based PKE protocol

# How-to-use
## Run POKE
```bash
cd build
cmake ..
make precomp
make
ctest -V -R poke-test_orig
ctest -V -R poke-test_inke*
```
## Output
```
13: Test command: /poke/build/src/poke/ref/poke_lvl1/test/poke-test_poke_lvl1 "100"
13: Working Directory: /poke/build/src/poke/ref/poke_lvl1/test
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 202783772.160000 nsec
13:   encrypt takes .................................... 114953876.480000 nsec
13:   decrypt takes .................................... 43927664.640000 nsec
1/6 Test #13: poke-test_poke_lvl1 ..............   Passed   36.17 sec
test 14
    Start 14: poke-test_poke_lvl3

14: Test command: /poke/build/src/poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
14: Working Directory: /poke/build/src/poke/ref/poke_lvl3/test
14: Test timeout computed to be: 1500
14: test loops : 100
14:   keygen takes .................................... 796999966.720000 nsec
14:   encrypt takes .................................... 555472747.520000 nsec
14:   decrypt takes .................................... 164178350.080000 nsec
2/6 Test #14: poke-test_poke_lvl3 ..............   Passed  151.67 sec
test 15
    Start 15: poke-test_poke_lvl5

15: Test command: /poke/build/src/poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
15: Working Directory: /poke/build/src/poke/ref/poke_lvl5/test
15: Test timeout computed to be: 1500
15: test loops : 100
15:   keygen takes .................................... 1888557698.560000 nsec
15:   encrypt takes .................................... 1474805593.600000 nsec
15:   decrypt takes .................................... 329211745.280000 nsec
3/6 Test #15: poke-test_poke_lvl5 ..............   Passed  369.34 sec
test 16
    Start 16: poke-test_inke_lvl1

16: Test command: /poke/build/src/poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
16: Working Directory: /poke/build/src/poke/ref/inke_lvl1/test
16: Test timeout computed to be: 1500
16: test loops : 100
16:   keygen takes .................................... 169351104.000000 nsec
16:   encrypt takes .................................... 150820976.640000 nsec
16:   decrypt takes .................................... 10704919.040000 nsec
4/6 Test #16: poke-test_inke_lvl1 ..............   Passed   33.13 sec
test 17
    Start 17: poke-test_inke_lvl3

17: Test command: /poke/build/src/poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
17: Working Directory: /poke/build/src/poke/ref/inke_lvl3/test
17: Test timeout computed to be: 1500
17: test loops : 100
17:   keygen takes .................................... 587249397.760000 nsec
17:   encrypt takes .................................... 661838397.440000 nsec
17:   decrypt takes .................................... 32551234.560000 nsec
5/6 Test #17: poke-test_inke_lvl3 ..............   Passed  128.22 sec
test 18
    Start 18: poke-test_inke_lvl5

18: Test command: /poke/build/src/poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
18: Working Directory: /poke/build/src/poke/ref/inke_lvl5/test
18: Test timeout computed to be: 1500
18: test loops : 100
18:   keygen takes .................................... 1358049735.680000 nsec
18:   encrypt takes .................................... 1808848657.920000 nsec
18:   decrypt takes .................................... 68191406.080000 nsec
6/6 Test #18: poke-test_inke_lvl5 ..............   Passed  325.92 sec
```

# Reference
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)