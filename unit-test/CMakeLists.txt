# Build Google Test library.
set(BUILD_GMOCK OFF)
set(BUILD_GTEST ON)
add_subdirectory(googletest)
include_directories(../include googletest/googletest/include)

# alenka-file
set(Boost_INCLUDE_DIR "")
set(ALENKA_FILE_DIR ../../Alenka-File CACHE PATH "Path to Alenka-File source.")
add_subdirectory(${ALENKA_FILE_DIR} alenka-file)
include_directories(${ALENKA_FILE_DIR}/include)

# The tests.
set(TEST_OPENCL_PLATFORM 0 CACHE STRING "OpenCL platform id to use for the test.")
set(TEST_OPENCL_DEVICE 0 CACHE STRING "OpenCL device id to use for the test.")
add_definitions(-DOPENCL_PLATFORM=${TEST_OPENCL_PLATFORM} -DOPENCL_DEVICE=${TEST_OPENCL_DEVICE})

file(GLOB SRC *.cpp)
add_executable(unit-test ${SRC})
target_link_libraries(unit-test ${LIBS_TO_LINK_ALENKA_SIGNAL} ${LIBS_TO_LINK_ALENKA_FILE} gtest gtest_main)

if(MSVC)
	set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS "-W4")
else()
	set_source_files_properties(${SRC1} PROPERTIES COMPILE_FLAGS "-Wall -pedantic")# -Wextra -Wconversion")
endif()

# Download test files if needed. This is a workaround for GitHub's screwed up LFS system.
add_custom_command(TARGET PRE_BUILD unit-test COMMAND bash download-data.sh WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/unit-test/data)

