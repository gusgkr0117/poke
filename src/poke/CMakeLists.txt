add_executable(poke_test poke.c)
target_include_directories(poke_test PRIVATE ${INC_GF} ${INC_COMMON} ${INC_PRECOMP} ${INC_HD} ${INC_EC} ${INC_INTBIG} ${INC_PUBLIC} ${INC_QUATERNION} ${INC_KLPT} ${INC_ID2ISO} include)
target_link_libraries(poke_test PRIVATE lib_ec lib_hd lib_intbig lib_klpt lib_quaternion lib_id2iso)
target_compile_options(poke_test PRIVATE ${C_OPT_FLAGS})