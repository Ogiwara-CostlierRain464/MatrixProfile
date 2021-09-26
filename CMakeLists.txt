cmake_minimum_required(VERSION 3.13)

#set(ENV{CC} gcc)
#set(ENV{CXX} g++)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -march=native")

project(MatrixProfile)

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE PROJECT_SOURCES src/*.cpp)
file(GLOB_RECURSE PROJECT_HEADERS src/*.h)

add_executable(${PROJECT_NAME}
    ${PROJECT_SOURCES}
    ${PROJECT_HEADERS}
)


target_link_libraries(
    ${PROJECT_NAME}
)


enable_testing()

include(FetchContent)
FetchContent_Declare(
    googletest
    # Specify the commit you depend on and update it regularly.
    URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

file(GLOB_RECURSE TEST_SOURCES test/*.cpp)
file(GLOB_RECURSE TEST_HEADERS test/*.h)
# main.cppのみ外す
list(REMOVE_ITEM PROJECT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

add_executable(tests
    ${TEST_SOURCES}
    ${TEST_HEADERS}
    ${PROJECT_SOURCES}
    ${PROJECT_HEADERS}
)

target_link_libraries(tests gtest_main)
add_test(NAME example_test COMMAND tests)