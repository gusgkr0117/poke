if(INKE_MANUSCRIPT_PARAMETERS)
    set(INKE_FIELD_SOURCE ${GFX_DIR}/fp_montgomery_generic.c)
else()
    set(INKE_FIELD_SOURCE fp_p392.c)
endif()
set(SOURCE_FILES_GF_${SVARIANT_UPPER}_REF
    ${GFX_DIR}/fp.c
    ${GFX_DIR}/fp2.c
    ${GFX_DIR}/mp.c
    ${INKE_FIELD_SOURCE}
)

add_library(${LIB_GF_${SVARIANT_UPPER}} ${SOURCE_FILES_GF_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_GF_${SVARIANT_UPPER}} PRIVATE ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} include ${PROJECT_SOURCE_DIR}/include)
target_compile_options(${LIB_GF_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
