set(SOURCE_FILES_PRECOMP
    gf_constants.c
    torsion_constants.c
)
include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})

add_library(lib_precomp ${SOURCE_FILES_PRECOMP})
target_include_directories(lib_precomp PRIVATE ${INC_COMMON} ${INC_INTBIG} ${PROJECT_SOURCE_DIR}/src/gf/include include)

add_custom_command(
    OUTPUT
        "${CMAKE_CURRENT_SOURCE_DIR}/gf_constants.c"
        # "${CMAKE_CURRENT_SOURCE_DIR}/torsion_constants.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/torsion_constants.h"
        # "${CMAKE_CURRENT_SOURCE_DIR}/include/klpt_constants.h"
        # "${CMAKE_CURRENT_SOURCE_DIR}/klpt_constants.c"
        # "${CMAKE_CURRENT_SOURCE_DIR}/include/encoded_sizes.h"
        # "${CMAKE_CURRENT_SOURCE_DIR}/quaternion_data.c"
        # "${CMAKE_CURRENT_SOURCE_DIR}/include/quaternion_data.h"
        # "${CMAKE_CURRENT_SOURCE_DIR}/endomorphism_action.c"
        # "${CMAKE_CURRENT_SOURCE_DIR}/include/endomorphism_action.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/ec_params.h"
    COMMAND
        echo "Please run manually: make precomp"
)

find_program(SAGEMATH sage)
add_custom_target(precomp
    DEPENDS
        "./sqisign_parameters.txt"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_torsion_constants.sage"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/ec_params.py"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_gf_constants.sage"
    WORKING_DIRECTORY
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM true)
set_target_properties(precomp PROPERTIES EXCLUDE_FROM_ALL TRUE)