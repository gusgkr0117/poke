set(SOURCE_FILES_PRECOMP
    gf_constants.c
)
include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})

add_library(lib_precomp ${SOURCE_FILES_PRECOMP})
target_include_directories(lib_precomp PRIVATE ${INC_COMMON} ${PROJECT_SOURCE_DIR}/src/gf/include include)

find_program(SAGEMATH sage)
add_custom_target(precomp
    DEPENDS
        "./sqisign_parameters.txt"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_gf_constants.sage"
    WORKING_DIRECTORY
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

set_target_properties(precomp PROPERTIES EXCLUDE_FROM_ALL TRUE)