# 🥗POKÉ : POint-based Key Exchange
+ C-Implementation of the most efficient isogeny-based PKE protocol

# How-to-use
## Run POKE
```bash
cd build
cmake ..
make precomp
make
ctest -V -R poke-test
```
## Output
```bash
keygen takes .................................... 208.235213 msec
encrypt takes .................................... 112.797798 msec
decrypt takes .................................... 44.286976 msec
```

# Reference
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)