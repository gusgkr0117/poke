add_library(lib_poke poke.c)
target_include_directories(lib_poke PRIVATE ${INC_GF} ${INC_COMMON} ${INC_PRECOMP} ${INC_HD} ${INC_EC} ${INC_INTBIG} ${INC_PUBLIC} ${INC_QUATERNION} ${INC_KLPT} ${INC_ID2ISO} include)
target_link_libraries(lib_poke PRIVATE lib_ec lib_hd lib_intbig lib_klpt lib_quaternion lib_id2iso)
target_compile_options(lib_poke PRIVATE ${C_OPT_FLAGS})

include(${SELECT_IMPL_TYPE})
include(${SELECT_POKE_VARIANT})

add_subdirectory(test)