# Build Google Test library.
set(BUILD_GMOCK OFF)
set(BUILD_GTEST ON)
add_subdirectory(googletest)
include_directories(../include googletest/googletest/include)

# alenka-file
set(ALENKA_FILE_DIR ../../Alenka-File CACHE PATH "Path to Alenka-File source.")
add_subdirectory(${ALENKA_FILE_DIR} alenka-file)
include_directories(${ALENKA_FILE_DIR}/include)

# The tests.
file(GLOB SRC *.cpp)
add_executable(unit-test ${SRC})
target_link_libraries(unit-test ${LIBS_TO_LINK_ALENKA_SIGNAL} ${LIBS_TO_LINK_ALENKA_FILE} gtest gtest_main)

message("link: ${LIBS_TO_LINK_ALENKA_FILE}")

if(MSVC)
	find_file(CLFFT_DLL clFFT.dll PATHS ${CLFFT_DIR} PATH_SUFFIXES bin)
	if(CLFFT_DLL)
		add_custom_command(TARGET unit-test POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CLFFT_DLL} $<TARGET_FILE_DIR:unit-test>)
	endif()
endif()

# Download test files if needed. This is a workaround for GitHub's screwed up LFS system.
add_custom_command(TARGET PRE_BUILD unit-test COMMAND bash download-data.sh WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/unit-test/data)

