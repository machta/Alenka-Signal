# Build Google Test library.
set(SRC_GTEST googletest/googletest/src/gtest-all.cc googletest/googletest/src/gtest_main.cc)
include_directories(googletest/googletest/include googletest/googletest)

# The tests.
set(TEST_OPENCL_PLATFORM 0 CACHE STRING "OpenCL platform id to use for the test.")
set(TEST_OPENCL_DEVICE 0 CACHE STRING "OpenCL device id to use for the test.")
add_definitions(-DOPENCL_PLATFORM=${TEST_OPENCL_PLATFORM} -DOPENCL_DEVICE=${TEST_OPENCL_DEVICE})

file(GLOB SRC *.cpp)
add_executable(unit-test ${SRC} ${SRC_GTEST})

find_package (Threads)
target_link_libraries(unit-test ${LIBS_TO_LINK_ALENKA_SIGNAL} ${LIBS_TO_LINK_ALENKA_FILE} ${CMAKE_THREAD_LIBS_INIT})

if(MSVC)
	set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS "-W4")
else()
	set_source_files_properties(${SRC1} PROPERTIES COMPILE_FLAGS "-Wall -pedantic")# -Wextra -Wconversion")
endif()
