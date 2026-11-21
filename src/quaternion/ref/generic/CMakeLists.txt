set(SOURCE_FILES_QUATERNION_GENERIC_REF
    algebra.c
    ideal.c
    dim4.c
    dim2.c
    integers.c
    lattice.c
    finit.c
    printer.c
    lll.c
    matkermod.c
    inke_uniform.c
    inke_mixed.c
)

add_library(${LIB_QUATERNION} ${SOURCE_FILES_QUATERNION_GENERIC_REF})
target_include_directories(${LIB_QUATERNION} PRIVATE common ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_QUATERNION})
target_link_libraries(${LIB_QUATERNION} ${LIB_INTBIG})
target_compile_options(${LIB_QUATERNION} PRIVATE ${C_OPT_FLAGS})

if (ENABLE_TESTS)
    add_executable(inke-test_uniform test/test_inke_uniform.c)
    target_include_directories(inke-test_uniform PRIVATE
        ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_QUATERNION})
    target_link_libraries(inke-test_uniform ${LIB_QUATERNION} ${GMP} poke_common_sys)
    add_test(NAME inke-test_uniform COMMAND inke-test_uniform)
    add_executable(inke-test_mixed test/test_inke_mixed.c)
    target_include_directories(inke-test_mixed PRIVATE
        ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_QUATERNION})
    target_link_libraries(inke-test_mixed ${LIB_QUATERNION} ${GMP} poke_common_sys)
    add_test(NAME inke-test_mixed COMMAND inke-test_mixed)
endif()
