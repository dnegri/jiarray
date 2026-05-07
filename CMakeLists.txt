cmake_minimum_required(VERSION 3.21)
project(dnegri.jiarray C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

add_library(dnegri.jiarray INTERFACE)
target_include_directories(dnegri.jiarray INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
add_library(jiarray ALIAS dnegri.jiarray)

# `ffor(i, b, e, step)` in JIArray.h uses C++20 __VA_OPT__, which on
# MSVC requires the conforming preprocessor.  Propagate the flag to
# every consumer that links jiarray; it is a no-op on GCC / Clang.
if (MSVC)
    target_compile_options(dnegri.jiarray INTERFACE /Zc:preprocessor)
endif()

option(JIARRAY_BUILD_TESTS "Build jiarray tests" ${PROJECT_IS_TOP_LEVEL})

if(JIARRAY_BUILD_TESTS)
  find_package(GTest)

  # Blitz++ is used only by the legacy iterator benchmark test files
  # (column_major_iter_test / row_major_iter_test).  Linux/Windows CI
  # typically doesn't have it installed; Mac Homebrew has it under
  # /opt/homebrew/lib.  Make it opt-in so the core test suite builds
  # unconditionally.
  option(JIARRAY_WITH_BLITZ "Link Blitz++ benchmark comparison (Mac Homebrew)" OFF)
  if(JIARRAY_WITH_BLITZ)
    find_library(BLITZ_LIBRARY NAMES blitz PATHS /opt/homebrew/lib /usr/local/lib)
    if(NOT BLITZ_LIBRARY)
      message(WARNING "Blitz++ not found; disabling JIARRAY_WITH_BLITZ")
      set(JIARRAY_WITH_BLITZ OFF)
    else()
      include_directories(/opt/homebrew/include /usr/local/include)
    endif()
  endif()

  option(BUILD_SUBPROJECT_TESTS "Build subproject tests" ON)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -DJIARRAY_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")

  # Blitz-dependent benchmark executables are gated by JIARRAY_WITH_BLITZ
  if(JIARRAY_WITH_BLITZ)
    add_executable(${PROJECT_NAME}.column_major ./test/column_major_iter_test.cpp ./test/main.cpp)
    target_link_libraries(${PROJECT_NAME}.column_major PRIVATE jiarray ${BLITZ_LIBRARY})
    target_compile_options(${PROJECT_NAME}.column_major PRIVATE -Wno-deprecated-declarations)

    add_executable(${PROJECT_NAME}.row_major ./test/row_major_iter_test.cpp ./test/main.cpp)
    target_compile_options(${PROJECT_NAME}.row_major PRIVATE -Wno-deprecated-declarations -DJIARRAY_COLUMN_MAJOR=0)
    target_link_libraries(${PROJECT_NAME}.row_major PRIVATE jiarray ${BLITZ_LIBRARY})
  endif()

  if(GTest_FOUND)
    enable_testing()

    # Core gtest suites — no Blitz++ required
    function(_ji_add_gtest name src)
      add_executable(${PROJECT_NAME}.${name} ${src})
      target_compile_definitions(${PROJECT_NAME}.${name} PRIVATE JIARRAY_DEBUG)
      target_compile_options   (${PROJECT_NAME}.${name} PRIVATE -Wno-deprecated-declarations)
      target_link_libraries    (${PROJECT_NAME}.${name}
        PRIVATE jiarray GTest::gtest_main GTest::gtest GTest::gmock)
      gtest_discover_tests(${PROJECT_NAME}.${name})
    endfunction()

    _ji_add_gtest(fast_array_test          ./test/fast_array_test.cpp)
    _ji_add_gtest(jivector_test            ./test/jivector_test.cpp)
    _ji_add_gtest(jiarray_coverage_test    ./test/jiarray_coverage_test.cpp)
    _ji_add_gtest(move_lifecycle_test      ./test/move_lifecycle_test.cpp)
    _ji_add_gtest(ffor_test                ./test/ffor_test.cpp)

    # Column-major (default) + row-major variants of the main matrix suite
    _ji_add_gtest(column_major_test        ./test/column_major_test.cpp)
    _ji_add_gtest(row_major_test           ./test/row_major_test.cpp)
    target_compile_options(${PROJECT_NAME}.row_major_test
      PRIVATE -DJIARRAY_COLUMN_MAJOR=0)

  endif()
endif()
