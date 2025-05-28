set(SOURCE_FILES_GF
    fp.c
    fp2.c
    mp.c
    fp_p431.c
)

add_library(poke_gf ${SOURCE_FILES_GF})
target_include_directories(poke_gf PRIVATE include ${INC_COMMON} ${INC_PRECOMP})
target_compile_options(poke_gf PRIVATE ${C_OPT_FLAGS})

include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})
add_subdirectory(test)