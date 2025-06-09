add_executable(basis-gen-test basis-gen-test.c)
target_include_directories(basis-gen-test PRIVATE ${INC_GF} ${INC_COMMON} ${INC_PRECOMP} ${INC_INTBIG} ../include)
target_link_libraries(basis-gen-test lib_intbig lib_ec)

add_test(ec_basis_gen_test basis-gen-test 10)