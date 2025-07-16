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
ctest -V -R poke-test_v1
```
## Output
```plaintext
==== original ====
4: test loops : 20
4:   keygen takes .................................... 1045711921.500000 cycles
4:   encrypt takes .................................... 547468016.050000 cycles
4:   decrypt takes .................................... 216827209.900000 cycles

==== small prime ====
5: test loops : 20
5:   keygen takes .................................... 826072202.350000 cycles
5:   encrypt takes .................................... 761635015.450000 cycles
5:   decrypt takes .................................... 64544412.450000 cycles
```

# Reference
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)