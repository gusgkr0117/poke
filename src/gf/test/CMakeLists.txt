add_executable(poke_test_gf_fp test_fp.c test_utils.c sha3.c)
target_link_libraries(poke_test_gf_fp lib_gf)
target_include_directories(poke_test_gf_fp PRIVATE ${INC_GF} ${INC_COMMON} ${INC_PRECOMP})