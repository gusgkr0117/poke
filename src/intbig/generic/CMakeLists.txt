set(SOURCE_FILES_INTBIG_GENERIC_REF
    intbig.c
)

add_library(lib_intbig ${SOURCE_FILES_INTBIG_GENERIC_REF})
target_include_directories(lib_intbig PRIVATE common ${INC_PUBLIC} ${INC_COMMON} include gmp)
target_compile_options(lib_intbig PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
