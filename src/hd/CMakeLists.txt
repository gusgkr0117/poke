set(SOURCE_FILES_HD_GENERIC_REF
    hd.c
    theta_structure.c 
    theta_isogenies.c 
)

add_library(lib_hd ${SOURCE_FILES_HD_GENERIC_REF})
target_include_directories(lib_hd PRIVATE ${INC_INTBIG} ${INC_PRECOMP} ${INC_PUBLIC} ${INC_GF} ${INC_COMMON} ${INC_EC} include)
target_link_libraries(lib_hd lib_intbig lib_ec)
target_compile_options(lib_hd PRIVATE ${C_OPT_FLAGS})

include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})

add_subdirectory(test)