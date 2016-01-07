cmake_minimum_required(VERSION 2.8)

# Change this section to correspond with your project info.
project(smart-uq)

set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ")

set(PROJECT_SUMMARY                            "Strathclyde Mechanical and Aerospace Researt Toolbox for Uncertainity Quantification (SMART-UQ)")
set(${CMAKE_PROJECT_NAME}_VENDOR_NAME          "spaceART")
set(${CMAKE_PROJECT_NAME}_VENDOR_CONTACT       "spaceART@strathclyde.ac.uk")
set(${CMAKE_PROJECT_NAME}_VERSION_MAJOR        0.1)
set(${CMAKE_PROJECT_NAME}_VERSION_MINOR        0.1)
set(${CMAKE_PROJECT_NAME}_VERSION_PATCH        0.1)
set(PROJECT_VERSION
  "${${CMAKE_PROJECT_NAME}_VERSION_MAJOR}.${${CMAKE_PROJECT_NAME}_VERSION_MINOR}.${${CMAKE_PROJECT_NAME}_VERSION_PATCH}")

set(PROJECT_PATH                               "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}     "${PROJECT_PATH}/cmake/Modules")
set(INCLUDE_PATH                               "${PROJECT_PATH}/include")
set(EIGEN_INCLUDE_PATH                         "${PROJECT_PATH}/include/Math")
set(SRC_PATH                                   "${PROJECT_PATH}/src")
set(TEST_SRC_PATH                              "${PROJECT_PATH}/test")
if(NOT EXTERNAL_PATH)
  set(EXTERNAL_PATH                            "${PROJECT_PATH}/external")
endif(NOT EXTERNAL_PATH)
if(NOT DOCS_PATH)
  set(DOCS_PATH                                "${PROJECT_PATH}/docs")
endif(NOT DOCS_PATH)
set(LIB_PATH                                   "${PROJECT_BINARY_DIR}/lib")
set(LIB_NAME                                   "smart-uq")
set(BIN_PATH                                   "${PROJECT_BINARY_DIR}/bin")
set(MAIN_NAME                                  "main")
set(TEST_PATH                                  "${PROJECT_BINARY_DIR}/test")
set(TEST_NAME                                  "test")

option(BUILD_MAIN                              "Build main function"            ON)
option(BUILD_DOCS                              "Build docs"                     OFF)
option(BUILD_TESTS                             "Build tests"                    OFF)

find_library(FFTW_LIBRARY fftw3)
find_library(FFTWF_LIBRARY fftw3f)
find_library(FFTWL_LIBRARY fftw3l)
if (FFTW_LIBRARY AND (FFTWF_LIBRARY AND FFTWL_LIBRARY))
  set(FFTW_LIBRARIES "${FFTW_LIBRARY}" "${FFTWF_LIBRARY}" "${FFTWL_LIBRARY}")
  set(CHEBYSHEV_DCT_MULTIPLICATION ON)
else (FFTW_LIBRARY AND (FFTWF_LIBRARY AND FFTWL_LIBRARY))
  set(CHEBYSHEV_DCT_MULTIPLICATION OFF)
  message(STATUS "WARNING: fftw3 libraries not found; always using direct multiplication for Chebyshev Polynomials.")
endif(FFTW_LIBRARY AND (FFTWF_LIBRARY AND FFTWL_LIBRARY))

file(GLOB SRC
  "${SRC_PATH}/Sampling/*.cpp"
  "${SRC_PATH}/Non-Intrusive/*.cpp"
  "${SRC_PATH}/Polynomial/*.cpp"
  "${SRC_PATH}/main_accuracy.cpp"
  "${SRC_PATH}/main_AKP_ni.cpp"
  "${SRC_PATH}/main_AKP_ni_sparse.cpp"
  "${SRC_PATH}/main_AKP_i_chebyshev.cpp"
  "${SRC_PATH}/main_AKP_i_canonical.cpp"
  "${SRC_PATH}/main_collision.cpp"
  "${SRC_PATH}/main_multiphase.cpp"
  "${SRC_PATH}/main_vanderpol_i_chebyshev.cpp"
  "${SRC_PATH}/main_vanderpol_i_canonical.cpp"
  "${SRC_PATH}/main_vanderpol_ni.cpp"
  "${SRC_PATH}/main_vanderpol_ni_sparse.cpp"  
  "${SRC_PATH}/main_lotka_volterra_i_chebyshev.cpp"
  "${SRC_PATH}/main_lotka_volterra_i_canonical.cpp"
  "${SRC_PATH}/main_lotka_volterra_ni.cpp"
  "${SRC_PATH}/main_lotka_volterra_ni_sparse.cpp"

  "${SRC_PATH}/main_spring_mass_ni.cpp"
  "${SRC_PATH}/main_spring_mass_ni_sparse.cpp"
  "${SRC_PATH}/main_spring_mass_i_canonical.cpp"
  "${SRC_PATH}/main_spring_mass_i_chebyshev.cpp"

  "${SRC_PATH}/main_intersect.cpp"
  "${SRC_PATH}/main_test_dct.cpp"
  "${SRC_PATH}/main_test_canonical_multiplication.cpp"
  "${SRC_PATH}/main_test_sampling.cpp"
  "${SRC_PATH}/main_tests.cpp"
)

set(MAIN_SRC
  "${SRC_PATH}/main.cpp"
)

file(GLOB TEST_SRC
  "${TEST_SRC_PATH}/Sampling/*.cpp"
  "${TEST_SRC_PATH}/Non-Intrusive/*.cpp"
  "${TEST_SRC_PATH}/Polynomial/*.cpp"
)

# Set CMake build-type. If it not supplied by the user, the default built type is "Release".
if(((NOT CMAKE_BUILD_TYPE)
  AND (NOT BUILD_COVERAGE_ANALYSIS))
  OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif((CMAKE_BUILD_TYPE STREQUAL "Debug") OR (BUILD_COVERAGE_ANALYSIS))
    message(STATUS "WARNING: building debug version!")
endif(((NOT CMAKE_BUILD_TYPE)
  AND (NOT BUILD_COVERAGE_ANALYSIS))
  OR (CMAKE_BUILD_TYPE STREQUAL "Release"))

# Set platform-specific compiler flags.
if(WIN32)
  if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # We need to force MSVC to use static runtime:
    # (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
    foreach(flag_var
      CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
  else(MSVC)
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor")
  endif(MSVC)
else(WIN32)
  set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor")
endif(WIN32)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif(CMAKE_COMPILER_IS_GNUCXX)

configure_file (
  "${SRC_PATH}/config.h.in"
  "${INCLUDE_PATH}/config.h"
  )

include_directories(${INCLUDE_PATH} ${EIGEN_INCLUDE_PATH})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_PATH})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH})
add_library(${LIB_NAME} ${SRC})

if(BUILD_MAIN)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_PATH})
  add_executable(${MAIN_NAME} ${MAIN_SRC})
  target_link_libraries(${MAIN_NAME} ${LIB_NAME} ${FFTW_LIBRARIES})
endif(BUILD_MAIN)

if(BUILD_DOCS)
  find_package(Doxygen)

  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
            "Doxygen is needed to build the documentation. Please make sure it is installed.")
  endif(NOT DOXYGEN_FOUND)

  configure_file(${PROJECT_PATH}/Doxyfile.in
                 ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # Add a custom target to run Doxygen whenever the project is built.
  # If you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the command below.
  add_custom_target(docs ALL COMMAND ${DOXYGEN_EXECUTABLE}
                    ${PROJECT_BINARY_DIR}/Doxyfile
                    SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif(BUILD_DOCS)

if(BUILD_TESTS)
  enable_testing()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_PATH})

  add_executable(${TEST_NAME} ${TEST_SRC})
  target_link_libraries(${TEST_NAME} ${LIB_NAME})
  add_test(NAME ${TEST_NAME} COMMAND "${TEST_PATH}/${TEST_NAME}")

  if(BUILD_COVERAGE_ANALYSIS)
    include(CodeCoverage)
    set(COVERAGE_EXTRACT '${PROJECT_PATH}/include/*' '${PROJECT_PATH}/src/*')
    setup_target_for_coverage(coverage "${TEST_NAME}" coverage)
  endif(BUILD_COVERAGE_ANALYSIS)
endif(BUILD_TESTS)

# Install header files and library.
# Destination is set by CMAKE_INSTALL_PREFIX and defaults to usual locations, unless overridden by
# user.
install(DIRECTORY ${INCLUDE_PATH}/${CMAKE_PROJECT_NAME}
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
install(TARGETS ${LIB_NAME} DESTINATION lib)
install(TARGETS ${BIN_NAME} DESTINATION bin)


# Set up packager.
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CMAKE_PROJECT_NAME}")
set(CPACK_PACKAGE_VENDOR              "${CMAKE_PROJECT_NAME}_VENDOR_NAME")
set(CPACK_PACKAGE_DESCRIPTION_FILE    "${PROJECT_PATH}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE       "${PROJECT_PATH}/LICENSE.md")
set(CPACK_PACKAGE_VERSION_MAJOR       "${${CMAKE_PROJECT_NAME}_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR       "${${CMAKE_PROJECT_NAME}_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH       "${${CMAKE_PROJECT_NAME}_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY   "${CMAKE_PROJECT_NAME}")
if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  set(CPACK_NSIS_DISPLAY_NAME         "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
  set(CPACK_NSIS_CONTACT              "${CMAKE_PROJECT_NAME}_VENDOR_CONTACT")
  set(CPACK_NSIS_MODIFY_PATH ON)
else(WIN32 AND NOT UNIX)
endif(WIN32 AND NOT UNIX)
include(CPack)
