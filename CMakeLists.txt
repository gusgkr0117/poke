cmake_minimum_required(VERSION 3.5)
project(POKE VERSION 1.0 LANGUAGES C)

set(POKE_VERSION "0")
set(CMAKE_C_STANDARD 99)

include(CTest)

option(ENABLE_STRICT "Build with strict compile options." OFF)
option(ENABLE_TESTS  "Enable compilation of tests." ON)
option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)
option(ENABLE_GMP_BUILD "Download and build external version of GMP" OFF)

SET(SVARIANT_S "orig;v1;")

set(SELECT_POKE_VARIANT ${PROJECT_SOURCE_DIR}/.cmake/poke_variant.cmake)
set(SELECT_IMPL_TYPE ${PROJECT_SOURCE_DIR}/.cmake/impl_type.cmake)

set(INC_PUBLIC ${PROJECT_SOURCE_DIR}/include)

include(.cmake/flags.cmake)
include(.cmake/target.cmake)
include(.cmake/gmpconfig.cmake)

add_subdirectory(src)