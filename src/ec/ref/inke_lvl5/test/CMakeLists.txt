add_executable(biextension.test_${SVARIANT_LOWER} biextension-test.c)
target_include_directories(biextension.test_${SVARIANT_LOWER} PRIVATE ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_QUATERNION} ${INC_EC} ${INC_PUBLIC})
target_link_libraries(biextension.test_${SVARIANT_LOWER} ${LIB_INTBIG} ${LIB_EC_${SVARIANT_UPPER}} ${GMP} poke_common_sys)

add_test(biextension.test_${SVARIANT_LOWER} biextension.test_${SVARIANT_LOWER} 10)