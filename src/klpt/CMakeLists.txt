set(SOURCE_FILES_KLPT_GENERIC_REF
    tools.c
)

add_library(lib_klpt ${SOURCE_FILES_KLPT_GENERIC_REF})
target_include_directories(lib_klpt PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_PRECOMP} ${INC_QUATERNION} ${INC_KLPT} include)
target_link_libraries(lib_klpt lib_quaternion)
target_compile_options(lib_klpt PRIVATE ${C_OPT_FLAGS})


include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})
