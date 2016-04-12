
add_executable(example_intrusive example_intrusive.cpp)
target_link_libraries(example_intrusive ${LIB_NAME} ${FFTW_LIBRARIES})

add_executable(example_taylor example_taylor.cpp)
target_link_libraries(example_taylor ${LIB_NAME} ${FFTW_LIBRARIES})

add_executable(example_nonintrusive example_nonintrusive.cpp)
target_link_libraries(example_nonintrusive ${LIB_NAME} ${FFTW_LIBRARIES})

