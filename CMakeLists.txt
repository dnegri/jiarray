cmake_minimum_required(VERSION 3.21)
project(dnegri.jiarray C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

add_library(jiarray INTERFACE)
target_include_directories(jiarray INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

option(JIARRAY_BUILD_TESTS "Build jiarray tests" ${PROJECT_IS_TOP_LEVEL})

if(JIARRAY_BUILD_TESTS)
  find_package(GTest)


  # CRITICAL FIX: Find and link Blitz++ library
  find_library(BLITZ_LIBRARY
    NAMES blitz
    PATHS /opt/homebrew/lib
    REQUIRED
  )
  include_directories(/opt/homebrew/include)


  option(BUILD_SUBPROJECT_TESTS "Build subproject tests" ON)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -DJIARRAY_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")


  add_executable(${PROJECT_NAME}.column_major ./test/column_major_iter_test.cpp ./test/main.cpp)

  target_link_libraries(${PROJECT_NAME}.column_major PRIVATE jiarray ${BLITZ_LIBRARY})

  # Optional: Suppress warnings
  target_compile_options(${PROJECT_NAME}.column_major PRIVATE -Wno-deprecated-declarations)

  add_executable(${PROJECT_NAME}.row_major ./test/row_major_iter_test.cpp ./test/main.cpp)
  target_compile_options(${PROJECT_NAME}.row_major PRIVATE -Wno-deprecated-declarations -DJIARRAY_COLUMN_MAJOR=0)
  target_link_libraries(${PROJECT_NAME}.row_major PRIVATE jiarray ${BLITZ_LIBRARY})

  if(GTest_FOUND)
    enable_testing()

    add_executable(${PROJECT_NAME}.column_major_test ./test/column_major_test.cpp)
    target_compile_definitions(${PROJECT_NAME}.column_major_test PRIVATE JIARRAY_DEBUG)
    target_compile_options(${PROJECT_NAME}.column_major_test PRIVATE -Wno-deprecated-declarations)

    target_link_libraries(
      ${PROJECT_NAME}.column_major_test
      PRIVATE jiarray ${BLITZ_LIBRARY}
      GTest::gtest_main
      GTest::gtest
      GTest::gmock
    )

    gtest_discover_tests(${PROJECT_NAME}.column_major_test)

    add_executable(${PROJECT_NAME}.row_major_test ./test/row_major_test.cpp)
    target_compile_definitions(${PROJECT_NAME}.row_major_test PRIVATE JIARRAY_DEBUG)
    target_compile_options(${PROJECT_NAME}.row_major_test PRIVATE -Wno-deprecated-declarations -DJIARRAY_COLUMN_MAJOR=0)

    target_link_libraries(
      ${PROJECT_NAME}.row_major_test
      PRIVATE jiarray ${BLITZ_LIBRARY}
      GTest::gtest_main
      GTest::gtest
      GTest::gmock
    )

    gtest_discover_tests(${PROJECT_NAME}.row_major_test)

    add_executable(${PROJECT_NAME}.fast_array_test ./test/fast_array_test.cpp)
    target_compile_definitions(${PROJECT_NAME}.fast_array_test PRIVATE JIARRAY_DEBUG)
    target_compile_options(${PROJECT_NAME}.fast_array_test PRIVATE -Wno-deprecated-declarations)
    target_link_libraries(
      ${PROJECT_NAME}.fast_array_test
      PRIVATE jiarray
      GTest::gtest_main
      GTest::gtest
    )
    gtest_discover_tests(${PROJECT_NAME}.fast_array_test)

    add_executable(${PROJECT_NAME}.jivector_test ./test/jivector_test.cpp)
    target_compile_definitions(${PROJECT_NAME}.jivector_test PRIVATE JIARRAY_DEBUG)
    target_compile_options(${PROJECT_NAME}.jivector_test PRIVATE -Wno-deprecated-declarations)
    target_link_libraries(
      ${PROJECT_NAME}.jivector_test
      PRIVATE jiarray
      GTest::gtest_main
      GTest::gtest
    )
    gtest_discover_tests(${PROJECT_NAME}.jivector_test)

    add_executable(${PROJECT_NAME}.jiarray_coverage_test ./test/jiarray_coverage_test.cpp)
    target_compile_definitions(${PROJECT_NAME}.jiarray_coverage_test PRIVATE JIARRAY_DEBUG)
    target_compile_options(${PROJECT_NAME}.jiarray_coverage_test PRIVATE -Wno-deprecated-declarations)
    target_link_libraries(
      ${PROJECT_NAME}.jiarray_coverage_test
      PRIVATE jiarray
      GTest::gtest_main
      GTest::gtest
    )
    gtest_discover_tests(${PROJECT_NAME}.jiarray_coverage_test)

  endif()
endif()
