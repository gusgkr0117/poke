cmake_minimum_required(VERSION 3.5)
project(POKE VERSION 1.0 LANGUAGES C)

set(POKE_VERSION "0")
set(CMAKE_C_STANDARD 99)

set(SELECT_POKE_VARIANT ${PROJECT_SOURCE_DIR}/.cmake/poke_variant.cmake)
set(SELECT_IMPL_TYPE ${PROJECT_SOURCE_DIR}/.cmake/impl_type.cmake)

include(.cmake/target.cmake)

add_subdirectory(src)