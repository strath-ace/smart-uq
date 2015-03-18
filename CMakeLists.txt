project(chebyshev_algebra)
cmake_minimum_required(VERSION 2.8)

set(PROJECT_PATH                               "${CMAKE_CURRENT_SOURCE_DIR}")
set(INCLUDE_PATH                               "${PROJECT_PATH}")
SET(TARGET_LIB                                 "${PROJECT_NAME}")
include_directories(${INCLUDE_PATH})

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${TARGET_LIB} lapack blas)
