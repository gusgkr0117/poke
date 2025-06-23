add_executable(poke_test poke.c)
target_include_directories(poke_test PRIVATE ${INC_GF} ${INC_COMMON} ${INC_PRECOMP} ${INC_HD} ${INC_EC} ${INC_INTBIG} ${INC_PUBLIC} ${INC_QUATERNION} ${INC_KLPT} include)
target_link_libraries(poke_test PRIVATE lib_gf poke_common_sys lib_precomp lib_hd lib_ec lib_intbig lib_quaternion lib_klpt gmp)