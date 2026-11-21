if(INKE_COMPLETE_BENCHMARKS)
    set(inke_dir "${PROJECT_SOURCE_DIR}/src/poke/ref/inkex")
    set(complete_dir "${inke_dir}/complete")
    set(complete_sources)
    foreach(source
            root-extension.c root-field-embedding.c a-extension-ec.c
            a-extension-quaternion.c a-extension-pairing.c a-kani-two-kernel.c
            root-theta-fast-add.c b-theta4-two.c b-theta4-odd.c root-ec-theta4.c
            root-theta-symplectic.c root-theta-label.c root-theta-factor.c
            root-kani-factor.c root-tate-projector.c root-quaternion-splitting.c
            root-ec-tate.c root-ec-halving.c a-theta-basechange.c a-theta-multifield.c
            b-theta-balanced.c a-theta-oddbasechange.c a-theta-oddfields.c
            a-theta-oddbalanced.c a-kani-fullchain.c a-kani-geometry.c
            root-poke-field-bridge.c b-gaussian-preparation.c b-poke-basis.c
            b-gaussian-convert.c a-general-kernel-lift.c root-theta-product.c
            b-product-fallback.c a-clapoti-preparation.c a-clapoti-convert.c
            b-uniform-keygen.c b-fast-keygen.c)
        list(APPEND complete_sources "${complete_dir}/${source}")
    endforeach()
endif()

foreach(VARIANT ${SVARIANT_S})
    string(TOUPPER ${VARIANT} VARIANT_UPPER)
    string(REGEX REPLACE "_lvl[135]$" "" SCHEME ${VARIANT})
    string(TOUPPER ${SCHEME} SCHEME_UPPER)
    if(NOT (INKE_MANUSCRIPT_PARAMETERS AND SCHEME STREQUAL "inke"))
        if(INKE_COMPLETE_BENCHMARKS AND NOT SCHEME STREQUAL "inke")
            # These baselines were measured with the earlier version of the driver.
            set(comparison_source comparison-baseline.c)
        else()
            set(comparison_source comparison.c)
        endif()
        add_executable(compare_${VARIANT} ${comparison_source})
        target_compile_definitions(compare_${VARIANT} PRIVATE
            BENCH_${SCHEME_UPPER} BENCH_VARIANT="${VARIANT}")
        target_include_directories(compare_${VARIANT} PRIVATE
            ${INC_POKE} ${INC_EC} ${INC_GF} ${INC_COMMON}
            ${INC_PRECOMP_${VARIANT_UPPER}} ${INC_PUBLIC})
        target_link_libraries(compare_${VARIANT}
            ${LIB_POKE_${VARIANT_UPPER}} ${LIB_ID2ISO_${VARIANT_UPPER}}
            ${LIB_HD_${VARIANT_UPPER}} ${GMP} poke_common_sys)
    endif()
    if(INKE_COMPLETE_BENCHMARKS AND SCHEME STREQUAL "inke")
        set(benchmark compare_complete_${VARIANT})
        add_executable(${benchmark} comparison.c "${inke_dir}/inke.c"
            "${inke_dir}/inke_keygen.c" ${complete_sources})
        target_compile_definitions(${benchmark} PRIVATE
            BENCH_INKE BENCH_INKE_UNIFORM BENCH_INKE_COMPLETE BENCH_VARIANT_TOKEN=${VARIANT})
        target_include_directories(${benchmark} PRIVATE
            ${GMP_INCLUDE} ${FLINT_INCLUDE_DIR} ${complete_dir} ${inke_dir}
            ${INC_POKE} ${INC_ID2ISO} ${INC_QUATERNION} ${INC_INTBIG} ${INC_EC}
            ${INC_GF} ${INC_COMMON} ${INC_PUBLIC} ${INC_HD} ${PROJECT_SOURCE_DIR}/src/hd/ref/hdx
            ${INC_KLPT} ${INC_PRECOMP_${VARIANT_UPPER}})
        if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
            target_compile_options(${benchmark} PRIVATE -Wall -Wextra)
        endif()
        target_link_libraries(${benchmark} PRIVATE
            ${LIB_POKE_${VARIANT_UPPER}} ${LIB_ID2ISO_${VARIANT_UPPER}} ${LIB_HD_${VARIANT_UPPER}}
            ${LIB_QUATERNION} poke_common_sys ${LIB_EC_${VARIANT_UPPER}} ${LIB_PRECOMP_${VARIANT_UPPER}}
            ${LIB_GF_${VARIANT_UPPER}} ${LIB_KLPT_${VARIANT_UPPER}} ${LIB_INTBIG} ${GMP} ${FLINT_LIBRARY})
        if(WIN32)
            target_link_libraries(${benchmark} PRIVATE bcrypt advapi32 kernel32 user32)
        endif()
    endif()
    if(SCHEME STREQUAL "inke")
        add_executable(compare_uniform_${VARIANT} comparison.c)
        target_compile_definitions(compare_uniform_${VARIANT} PRIVATE
            BENCH_INKE BENCH_INKE_UNIFORM BENCH_VARIANT="${VARIANT}")
        target_include_directories(compare_uniform_${VARIANT} PRIVATE
            ${INC_POKE} ${INC_ID2ISO} ${INC_QUATERNION} ${INC_INTBIG}
            ${INC_EC} ${INC_GF} ${INC_COMMON} ${INC_PRECOMP_${VARIANT_UPPER}} ${INC_PUBLIC})
        target_link_libraries(compare_uniform_${VARIANT}
            ${LIB_POKE_${VARIANT_UPPER}} ${LIB_ID2ISO_${VARIANT_UPPER}}
            ${LIB_HD_${VARIANT_UPPER}} ${LIB_QUATERNION} ${GMP} poke_common_sys)
    endif()
    if(ENABLE_TESTS AND SCHEME STREQUAL "inke")
        add_executable(test_curve_model_${VARIANT}
            ../src/poke/ref/inkex/test_curve_model.c)
        target_include_directories(test_curve_model_${VARIANT} PRIVATE
            ../src/poke/ref/inkex ${INC_POKE} ${INC_INTBIG} ${INC_EC} ${INC_GF} ${INC_COMMON}
            ${INC_PRECOMP_${VARIANT_UPPER}} ${INC_PUBLIC})
        target_link_libraries(test_curve_model_${VARIANT}
            ${LIB_POKE_${VARIANT_UPPER}} ${LIB_EC_${VARIANT_UPPER}} ${GMP} poke_common_sys)
        add_test(NAME test_curve_model_${VARIANT} COMMAND test_curve_model_${VARIANT})
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../src/hd/ref/hdx/test/test_inke_chain_strategy.c")
        add_executable(test_chain_strategy_${VARIANT}
            ../src/hd/ref/hdx/test/test_inke_chain_strategy.c)
        target_include_directories(test_chain_strategy_${VARIANT} PRIVATE
            ../src/hd/ref/hdx ${INC_HD} ${INC_ID2ISO} ${INC_QUATERNION} ${INC_INTBIG}
            ${INC_EC} ${INC_GF} ${INC_COMMON} ${INC_PUBLIC} ${INC_PRECOMP_${VARIANT_UPPER}})
        target_link_libraries(test_chain_strategy_${VARIANT}
            ${LIB_ID2ISO_${VARIANT_UPPER}} ${LIB_HD_${VARIANT_UPPER}}
            ${LIB_QUATERNION} ${GMP} poke_common_sys)
        add_test(NAME test_chain_strategy_${VARIANT} COMMAND test_chain_strategy_${VARIANT})
        endif()
        if(INKE_MANUSCRIPT_PARAMETERS)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_precomp_reference.c")
            add_executable(test_precomp_reference_${VARIANT} test_precomp_reference.c)
            target_include_directories(test_precomp_reference_${VARIANT} PRIVATE
                ${INC_HD} ${INC_ID2ISO} ${INC_QUATERNION} ${INC_INTBIG} ${INC_EC}
                ${INC_GF} ${INC_COMMON} ${INC_PRECOMP_${VARIANT_UPPER}} ${INC_PUBLIC})
            target_link_libraries(test_precomp_reference_${VARIANT}
                ${LIB_ID2ISO_${VARIANT_UPPER}} ${LIB_HD_${VARIANT_UPPER}}
                ${LIB_QUATERNION} ${GMP} poke_common_sys)
            add_test(NAME test_precomp_reference_${VARIANT} COMMAND test_precomp_reference_${VARIANT})
            endif()
            add_executable(test_field_reference_${VARIANT}
                ../src/gf/ref/gfx/test/test_field_reference.c)
            target_include_directories(test_field_reference_${VARIANT} PRIVATE
                ${INC_GF} ${INC_COMMON} ${INC_PRECOMP_${VARIANT_UPPER}} ${INC_PUBLIC})
            target_link_libraries(test_field_reference_${VARIANT}
                ${LIB_GF_${VARIANT_UPPER}} ${GMP} poke_common_sys)
            add_test(NAME test_field_reference_${VARIANT} COMMAND test_field_reference_${VARIANT})
        endif()
        add_executable(test_theta3_${VARIANT}
            ../src/hd/ref/hdx/test/test_theta3.c)
        target_include_directories(test_theta3_${VARIANT} PRIVATE
            ../src/hd/ref/hdx ${INC_HD} ${INC_EC} ${INC_GF}
            ${INC_GF_${VARIANT_UPPER}} ${INC_COMMON} ${INC_INTBIG} ${INC_QUATERNION}
            ${INC_PRECOMP_${VARIANT_UPPER}} ${INC_PUBLIC})
        target_link_libraries(test_theta3_${VARIANT}
            ${LIB_HD_${VARIANT_UPPER}} ${LIB_PRECOMP_${VARIANT_UPPER}}
            ${GMP} poke_common_sys)
        add_test(NAME test_theta3_${VARIANT} COMMAND test_theta3_${VARIANT})
        add_executable(test_conversion_${VARIANT}
            ../src/id2iso/ref/id2isox/test/test_inke_conversion.c)
        target_include_directories(test_conversion_${VARIANT} PRIVATE
            ${INC_ID2ISO} ${INC_QUATERNION} ${INC_INTBIG} ${INC_EC} ${INC_GF}
            ${INC_COMMON} ${INC_PUBLIC} ${INC_PRECOMP_${VARIANT_UPPER}})
        target_link_libraries(test_conversion_${VARIANT}
            ${LIB_ID2ISO_${VARIANT_UPPER}} ${LIB_HD_${VARIANT_UPPER}}
            ${LIB_QUATERNION} ${GMP} poke_common_sys)
        add_executable(test_uniform_keygen_${VARIANT}
            ../src/poke/ref/inkex/test_uniform_keygen.c)
        target_include_directories(test_uniform_keygen_${VARIANT} PRIVATE
            ${INC_POKE} ${INC_ID2ISO} ${INC_QUATERNION} ${INC_INTBIG} ${INC_EC} ${INC_GF}
            ${INC_COMMON} ${INC_PUBLIC} ${INC_PRECOMP_${VARIANT_UPPER}})
        target_link_libraries(test_uniform_keygen_${VARIANT}
            ${LIB_POKE_${VARIANT_UPPER}} ${LIB_ID2ISO_${VARIANT_UPPER}} ${LIB_HD_${VARIANT_UPPER}}
            ${LIB_QUATERNION} ${GMP} poke_common_sys)
        add_test(NAME test_uniform_keygen_${VARIANT} COMMAND test_uniform_keygen_${VARIANT})
    endif()
endforeach()
