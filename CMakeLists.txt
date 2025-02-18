cmake_minimum_required(VERSION 3.16)
project(dnegri.jiarray C CXX Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(GTest)


set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -DJIARRAY_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")

include_directories(/opt/homebrew/include)

# Add test source files
file(GLOB SOURCES ./test/simple_test.cpp ./test/slice_test.cpp)

add_executable(${PROJECT_NAME} ./test/iter_test.cpp  ./test/range_test.cpp ./test/main.cpp)

if(GTest_FOUND)
  enable_testing()

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -DJIARRAY_DEBUG")

  add_executable(${PROJECT_NAME}.test ${SOURCES})
  target_link_libraries(
    ${PROJECT_NAME}.test
    PRIVATE GTest::gtest_main
            GTest::gtest
            GTest::gmock)
  gtest_discover_tests(${PROJECT_NAME}.test)
endif()
