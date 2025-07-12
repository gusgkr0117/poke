add_executable(poke-test poke-test.c)
target_include_directories(poke-test PRIVATE ${INC_POKE} ${INC_EC} ${INC_GF} ${INC_COMMON} ${INC_PRECOMP} ${INC_PUBLIC})
target_link_libraries(poke-test lib_poke)

add_test(poke-test poke-test)