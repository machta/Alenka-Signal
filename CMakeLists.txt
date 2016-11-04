cmake_minimum_required(VERSION 3.1)
project(Alenka-Signal)

# Options.
option(BUILD_TESTS_ALENKA_SIGNAL "Builds unit tests for Alenka-Signal." OFF)
option(USE_CL_1_2 "Use OpenCL 1.2" ON)
set(CMAKE_CXX_STANDARD 11)

if (MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

if (USE_CL_1_2)
	add_definitions(-DCL_1_2)
endif()

# If you include this library into your project, you will want to link to these libraries.
set(LIBS_TO_LINK alenka-signal clFFT OpenCL)

# The library.
if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wconversion -pedantic")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

include_directories(include)

if(BUILD_TESTS_ALENKA_SIGNAL)
	add_subdirectory (unit-test)
endif (BUILD_TESTS_ALENKA_SIGNAL)

file(GLOB SRC src/*)
add_library(alenka-signal ${SRC})
