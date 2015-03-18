project(chebyshev_algebra)
cmake_minimum_required(VERSION 2.8)

set(PROJECT_PATH                               "${CMAKE_CURRENT_SOURCE_DIR}")
set(INCLUDE_PATH                               "${PROJECT_PATH}")
include_directories(${INCLUDE_PATH})

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

