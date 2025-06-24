set(SOURCE_FILES_ID2ISO_GENERIC_REF
    id2iso.c
)

add_library(lib_id2iso ${SOURCE_FILES_ID2ISO_GENERIC_REF})
target_include_directories(lib_id2iso PRIVATE common ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_HD} ${INC_ID2ISO} ${INC_COMMON} ${INC_GF})
target_link_libraries(lib_id2iso lib_klpt lib_ec)
target_compile_options(lib_id2iso PRIVATE ${C_OPT_FLAGS})

include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})