set(SOURCE_FILES_QUATERNION_GENERIC_REF
    algebra.c
    # ideal.c
    dim4.c
    # dim2.c
    integers.c
    lattice.c
    finit.c
    printer.c
    # lll.c
    # matkermod.c
)

add_library(lib_quaternion ${SOURCE_FILES_QUATERNION_GENERIC_REF})
target_include_directories(lib_quaternion PRIVATE common ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_QUATERNION})
target_link_libraries(lib_quaternion lib_intbig)
target_compile_options(lib_quaternion PRIVATE ${C_OPT_FLAGS})
