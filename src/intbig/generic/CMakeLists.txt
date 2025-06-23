set(SOURCE_FILES_INTBIG_GENERIC_REF
    intbig.c
)

add_library(lib_intbig ${SOURCE_FILES_INTBIG_GENERIC_REF})
target_include_directories(lib_intbig PRIVATE common ${INC_PUBLIC} ${INC_COMMON} include gmp)
target_link_libraries(lib_intbig PRIVATE gmp poke_common_sys lib_precomp lib_hd lib_ec lib_gf)
target_compile_options(lib_intbig PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
