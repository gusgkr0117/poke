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
```plaintext
==== POKE ====
4: test loops : 20
4:   keygen takes .................................... 1045711921.500000 cycles
4:   encrypt takes .................................... 547468016.050000 cycles
4:   decrypt takes .................................... 216827209.900000 cycles


==== INKE ====
...
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 785189390.650000 cycles
13:   encrypt takes .................................... 705573401.540000 cycles
13:   decrypt takes .................................... 50798305.540000 cycles
1/3 Test #13: poke-test_inke_lvl1 ..............   Passed   46.63 sec
test 14
    Start 14: poke-test_inke_lvl3

...
14: Test timeout computed to be: 1500
14: test loops : 100
14:   keygen takes .................................... 2864873504.670000 cycles
14:   encrypt takes .................................... 3306575142.910000 cycles
14:   decrypt takes .................................... 162224564.080000 cycles
2/3 Test #14: poke-test_inke_lvl3 ..............   Passed  194.20 sec
test 15
    Start 15: poke-test_inke_lvl5

...
15: Test timeout computed to be: 1500
15: test loops : 100
15:   keygen takes .................................... 7117719648.750000 cycles
15:   encrypt takes .................................... 9487551912.600000 cycles
15:   decrypt takes .................................... 357037187.360000 cycles
3/3 Test #15: poke-test_inke_lvl5 ..............   Passed  520.93 sec
```

# Reference
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)