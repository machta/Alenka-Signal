cmake_minimum_required(VERSION 3.1)
project(Alenka-Signal)

# Options.
option(BUILD_TESTS_ALENKA_SIGNAL "Builds unit tests for Alenka-Signal." OFF)
option(USE_CL_1_2 "Use OpenCL 1.2" ON)
set(CMAKE_CXX_STANDARD 11)

add_definitions(-DCL_USE_DEPRECATED_OPENCL_1_2_APIS)

if(MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
endif()

if(USE_CL_1_2)
	add_definitions(-DCL_1_2)
endif()

find_package(OpenCL REQUIRED)
set(INC ${INC} ${OpenCL_INCLUDE_DIRS})

# First build all the library dependencies.
set(BUILD_LOADLIBRARIES false)
set(BUILD_EXAMPLES false)
set(BUILD_SHARED_LIBS false)
add_subdirectory(clFFT/src)
set(INC ${INC} clFFT/src/include ${CMAKE_BINARY_DIR}/clFFT/src/include)

# If you want to use this library, you need to link to these libraries.
set(LIBS_TO_LINK_ALENKA_SIGNAL alenka-signal clFFT ${OpenCL_LIBRARIES})

# The library.
if(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wconversion -pedantic")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

include_directories(include ${INC} alglib/src)

if(BUILD_TESTS_ALENKA_SIGNAL)
	add_subdirectory (unit-test)
endif()

file(GLOB SRC src/* alglib/src/fasttransforms.cpp alglib/src/interpolation.cpp)
add_library(alenka-signal ${SRC})
