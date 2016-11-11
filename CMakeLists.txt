cmake_minimum_required(VERSION 3.1)
project(Alenka-Signal)

# Options.
option(BUILD_TESTS_ALENKA_SIGNAL "Builds unit tests for Alenka-Signal." OFF)
option(USE_CL_1_2 "Use OpenCL 1.2" ON)
set(CMAKE_CXX_STANDARD 11)
set(<variable> <value>... CACHE <type> <docstring> [FORCE])

add_definitions(-DCL_USE_DEPRECATED_OPENCL_1_2_APIS)

if(MSVC)
	find_path(OPENCL_DIR include/CL/cl_gl.h PATHS "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v8.0" "C:/Program Files (x86)/AMD APP SDK/3.0")

	if(CMAKE_CL_64)
		find_library(OPENCL_LIB OpenCL PATHS ${OPENCL_DIR} PATH_SUFFIXES lib/x64 lib/x86_64)
		
		set(CLFFT_DIR "C:/Program Files/clFFT" CACHE PATH "Path to clFFT library")
		find_library(CLFFT_LIB clFFT PATHS ${CLFFT_DIR} PATH_SUFFIXES lib64/import)
	else()
		find_library(OPENCL_LIB OpenCL PATHS ${OPENCL_DIR} PATH_SUFFIXES lib/Win32 lib/x86)
		
		set(CLFFT_DIR "C:/Program Files (x86)/clFFT" CACHE PATH "Path to clFFT library")
		find_library(CLFFT_LIB clFFT PATHS ${CLFFT_DIR} PATH_SUFFIXES lib/import)
	endif()

	set(INC ${INC} ${OPENCL_DIR}/include)
	set(INC ${INC} ${CLFFT_DIR}/include)
	
	add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
endif()

if(USE_CL_1_2)
	add_definitions(-DCL_1_2)
endif()

# If you include this library into your project, you will want to link to these libraries.
if(MSVC)
	set(LIBS_TO_LINK alenka-signal ${OPENCL_LIB} ${CLFFT_LIB})
else()
	set(LIBS_TO_LINK alenka-signal clFFT OpenCL)
endif()

# The library.
if(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wconversion -pedantic")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

include_directories(include ${INC})

if(BUILD_TESTS_ALENKA_SIGNAL)
	add_subdirectory (unit-test)
endif()

file(GLOB SRC src/*)
add_library(alenka-signal ${SRC})
