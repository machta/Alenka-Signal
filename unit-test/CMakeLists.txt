# Build Google Test library.
set(BUILD_GMOCK OFF)
set(BUILD_GTEST ON)
add_subdirectory(googletest)
include_directories(../include googletest/googletest/include)

# The tests.
file(GLOB SRC *)
add_executable (unit-test ${SRC})
target_link_libraries (unit-test ${LIBS_TO_LINK_ALENKA_SIGNAL} gtest gtest_main)

if(MSVC)
	find_file(CLFFT_DLL clFFT.dll PATHS ${CLFFT_DIR} PATH_SUFFIXES bin)
	if(CLFFT_DLL)
		add_custom_command(TARGET unit-test POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CLFFT_DLL} $<TARGET_FILE_DIR:unit-test>)
	endif()
endif()
