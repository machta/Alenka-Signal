cmake_minimum_required(VERSION 3.1)
project(Alenka-Signal)

# Options.
option(BUILD_TESTS_ALENKA_SIGNAL "Builds unit tests for Alenka-Signal." OFF)
option(USE_CL_1_2 "Use OpenCL 1.2" ON)
option(BUILD_CLFFT "Should a local build of clFFT library be used." OFF)
set(CMAKE_CXX_STANDARD 11)
set(BUILD_SHARED_LIBS off CACHE BOOL "")

add_definitions(-DCL_USE_DEPRECATED_OPENCL_1_2_APIS)

if(MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

if(USE_CL_1_2)
	add_definitions(-DCL_1_2)
endif()

find_package(OpenCL)
include_directories(${OpenCL_INCLUDE_DIRS})

# First build all the library dependencies.
set(BUILD_LOADLIBRARIES off CACHE BOOL "")
set(BUILD_EXAMPLES off CACHE BOOL "")
add_subdirectory(clFFT/src)
include_directories(clFFT/src/include ${CMAKE_BINARY_DIR}/clFFT/src/include)

# If you want to use this library, you need to link to these libraries.
set(LIBS_TO_LINK_ALENKA_SIGNAL alenka-signal clFFT ${OpenCL_LIBRARIES})
set(LIBS_TO_LINK_ALENKA_SIGNAL ${LIBS_TO_LINK_ALENKA_SIGNAL} PARENT_SCOPE)

# The library.
file(GLOB SRC1 src/*)
file(GLOB SRC2 alglib/src/*.cpp)
add_library(alenka-signal STATIC ${SRC1} ${SRC2})

if(MSVC)
	set_source_files_properties(${SRC1} PROPERTIES COMPILE_FLAGS "-W4")
else()
	set_source_files_properties(${SRC1} PROPERTIES COMPILE_FLAGS "-Wall -pedantic")# -Wextra -Wconversion")
endif()

include_directories(include alglib/src eigen)

if(BUILD_TESTS_ALENKA_SIGNAL)
	add_subdirectory(unit-test)
endif()

# Install rule.
if(BUILD_CLFFT)
	add_custom_target(install-alenka-signal
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:clFFT> $<TARGET_FILE:alenka-signal> ${CMAKE_INSTALL_PREFIX}/lib
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/include ${CMAKE_INSTALL_PREFIX}/include
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/clFFT/src/include/clFFT.h ${CMAKE_BINARY_DIR}/clFFT/src/include/clFFT.version.h ${CMAKE_INSTALL_PREFIX}/include)
else()
	add_custom_target(install-alenka-signal
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:alenka-signal> ${CMAKE_INSTALL_PREFIX}/lib
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/include ${CMAKE_INSTALL_PREFIX}/include)
endif()
