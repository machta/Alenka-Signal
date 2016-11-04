# Build Google Test library.
set(BUILD_GMOCK OFF)
set(BUILD_GTEST ON)
add_subdirectory(googletest)
include_directories(../include googletest/googletest/include)

# The tests.
file(GLOB SRC *)
add_executable (unit-test ${SRC})
target_link_libraries (unit-test ${LIBS_TO_LINK} gtest gtest_main)
