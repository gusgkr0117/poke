cmake_minimum_required(VERSION 3.5)
project(POKE VERSION 1.0 LANGUAGES C)

set(POKE_VERSION "0")
set(CMAKE_C_STANDARD 99)

include(CTest)

option(ENABLE_STRICT "Build with strict compile options." OFF)
option(ENABLE_TESTS  "Enable compilation of tests." ON)
option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)
option(ENABLE_GMP_BUILD "Download and build external version of GMP" OFF)
option(INKE_MANUSCRIPT_PARAMETERS "Use the manuscript's primes and bases for revised INKE" OFF)
option(INKE_COMPLETE_BENCHMARKS "Build the INKE, POKE and PIKE benchmarks used for the paper" OFF)
if (INKE_COMPLETE_BENCHMARKS)
    set(INKE_MANUSCRIPT_PARAMETERS ON)
endif()

if (NOT DEFINED POKE_BUILD_TYPE)
	SET(POKE_BUILD_TYPE "ref")
endif()

SET(SVARIANT_S "poke_lvl1;poke_lvl3;poke_lvl5;inke_lvl1;inke_lvl3;inke_lvl5;pike_lvl1;pike_lvl3;pike_lvl5")

set(SELECT_POKE_VARIANT ${PROJECT_SOURCE_DIR}/.cmake/poke_variant.cmake)
set(SELECT_IMPL_TYPE ${PROJECT_SOURCE_DIR}/.cmake/impl_type.cmake)

set(INC_PUBLIC ${PROJECT_SOURCE_DIR}/include)

include(.cmake/flags.cmake)
include(.cmake/sanitizers.cmake)
include(.cmake/target.cmake)
include(.cmake/gmpconfig.cmake)

# The high security levels use large automatic arrays. Give every Windows
# executable the same stack reservation (the MinGW default is only 1 MiB).
if (MINGW)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,16777216")
endif()

add_subdirectory(src)

if (ENABLE_TESTS OR INKE_COMPLETE_BENCHMARKS)
    add_subdirectory(bench)
endif()
