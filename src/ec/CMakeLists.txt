set(SOURCE_FILES_EC
    ec.c
    kps.c
    xisog.c
    xeval.c
    isog_chains.c
    basis.c
    bitextension.c
)

add_library(lib_ec ${SOURCE_FILES_EC})
target_include_directories(lib_ec PRIVATE ${INC_GF} ${INC_COMMON} ${INC_PRECOMP} ${INC_INTBIG} include)
target_link_libraries(lib_ec lib_precomp lib_gf)

include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})

add_subdirectory(test)