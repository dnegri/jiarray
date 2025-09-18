cmake_minimum_required(VERSION 3.20)
project(dnegri.jiarray C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(GTest)
option(BUILD_SUBPROJECT_TESTS "Build subproject tests" ON)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -DJIARRAY_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")

include_directories(/opt/homebrew/include)

# CRITICAL FIX: Find and link Blitz++ library
find_library(BLITZ_LIBRARY 
    NAMES blitz
    PATHS /opt/homebrew/lib
    REQUIRED
)

add_executable(${PROJECT_NAME}.column_major ./test/column_major_iter_test.cpp  ./test/main.cpp)

target_link_libraries(${PROJECT_NAME}.column_major PRIVATE ${BLITZ_LIBRARY})

# Optional: Suppress warnings
target_compile_options(${PROJECT_NAME}.column_major PRIVATE -Wno-deprecated-declarations)

add_executable(${PROJECT_NAME}.row_major ./test/row_major_iter_test.cpp  ./test/main.cpp)
target_compile_options(${PROJECT_NAME}.row_major PRIVATE -Wno-deprecated-declarations -DJIARRAY_COLUMN_MAJOR=0 )
target_link_libraries(${PROJECT_NAME}.row_major PRIVATE ${BLITZ_LIBRARY})

if(GTest_FOUND)
  enable_testing()

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -DJIARRAY_DEBUG")

  add_executable(${PROJECT_NAME}.column_major_test ./test/column_major_test.cpp )

  target_compile_options(${PROJECT_NAME}.column_major_test PRIVATE -Wno-deprecated-declarations)

  target_link_libraries(
    ${PROJECT_NAME}.column_major_test 
    PRIVATE ${BLITZ_LIBRARY}
    GTest::gtest_main
    GTest::gtest
    GTest::gmock
  )
  
  gtest_discover_tests(${PROJECT_NAME}.column_major_test)

  add_executable(${PROJECT_NAME}.row_major_test ./test/row_major_test.cpp )
  target_compile_options(${PROJECT_NAME}.row_major_test PRIVATE -Wno-deprecated-declarations -DJIARRAY_COLUMN_MAJOR=0 )

  target_link_libraries(
    ${PROJECT_NAME}.row_major_test
    PRIVATE ${BLITZ_LIBRARY}
    GTest::gtest_main
    GTest::gtest
    GTest::gmock
  )

  gtest_discover_tests(${PROJECT_NAME}.row_major_test)

endif()
