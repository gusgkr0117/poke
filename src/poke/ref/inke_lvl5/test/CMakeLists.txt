add_executable(poke-test_${SVARIANT_LOWER} ${INKEX_DIR}/test/inke-test.c)
target_include_directories(poke-test_${SVARIANT_LOWER} PRIVATE ${INC_POKE} ${INC_POKE_${SVARIANT_UPPER}} ${INC_EC} ${INC_GF} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC})
target_link_libraries(poke-test_${SVARIANT_LOWER} ${LIB_POKE_${SVARIANT_UPPER}} ${LIB_ID2ISO_${SVARIANT_UPPER}} ${LIB_HD_${SVARIANT_UPPER}} ${LIB_INTBIG} ${GMP} poke_common_sys)

add_test(poke-test_${SVARIANT_LOWER} poke-test_${SVARIANT_LOWER} 100)