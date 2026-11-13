cmake_minimum_required(VERSION 3.21)
project(dnegri.jiarray VERSION 0.10.0 LANGUAGES C CXX)

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
# A CUDA consumer gets it through -Xcompiler: nvcc reads a bare /Zc:... as an
# input file, while the cl compiling a .cu file's host side needs the flag.
if (MSVC)
    target_compile_options(dnegri.jiarray INTERFACE
        $<$<COMPILE_LANGUAGE:C,CXX>:/Zc:preprocessor>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>)
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
    _ji_add_gtest(fast_array_rowmajor_test ./test/fast_array_test.cpp)
    target_compile_options(${PROJECT_NAME}.fast_array_rowmajor_test
      PRIVATE -DJIARRAY_COLUMN_MAJOR=0)
    _ji_add_gtest(jivector_test            ./test/jivector_test.cpp)
    _ji_add_gtest(jiarray_coverage_test    ./test/jiarray_coverage_test.cpp)
    _ji_add_gtest(move_lifecycle_test      ./test/move_lifecycle_test.cpp)
    _ji_add_gtest(ffor_test                ./test/ffor_test.cpp)

    # FastArray-indexed access under both index bases and both storage orders
    _ji_add_gtest(fast_array_index_test           ./test/fast_array_index_test.cpp)
    _ji_add_gtest(fast_array_index_rowmajor_test  ./test/fast_array_index_test.cpp)
    target_compile_options(${PROJECT_NAME}.fast_array_index_rowmajor_test PRIVATE -DJIARRAY_COLUMN_MAJOR=0)
    _ji_add_gtest(fast_array_index_base0_test     ./test/fast_array_index_test.cpp)
    target_compile_options(${PROJECT_NAME}.fast_array_index_base0_test PRIVATE -DJIARRAY_OFFSET=0)
    _ji_add_gtest(fast_array_index_base0_rowmajor_test ./test/fast_array_index_test.cpp)
    target_compile_options(${PROJECT_NAME}.fast_array_index_base0_rowmajor_test
      PRIVATE -DJIARRAY_OFFSET=0 -DJIARRAY_COLUMN_MAJOR=0)

    # init0 range and large-offset addressing: both index bases x both storage orders,
    # each with and without JIARRAY_DEBUG (-U undoes the Debug-config -DJIARRAY_DEBUG).
    foreach(_ji_base 0 1)
      foreach(_ji_order col row)
        foreach(_ji_debug debug nodebug)
          set(_ji_name init0_base${_ji_base}_${_ji_order}_${_ji_debug}_test)
          _ji_add_gtest(${_ji_name} ./test/init0_test.cpp)
          target_compile_definitions(${PROJECT_NAME}.${_ji_name} PRIVATE JIARRAY_OFFSET=${_ji_base})
          if(_ji_order STREQUAL "row")
            target_compile_definitions(${PROJECT_NAME}.${_ji_name} PRIVATE JIARRAY_COLUMN_MAJOR=0)
          endif()
          if(_ji_debug STREQUAL "nodebug")
            get_target_property(_ji_defs ${PROJECT_NAME}.${_ji_name} COMPILE_DEFINITIONS)
            list(REMOVE_ITEM _ji_defs JIARRAY_DEBUG)
            set_target_properties(${PROJECT_NAME}.${_ji_name} PROPERTIES COMPILE_DEFINITIONS "${_ji_defs}")
            target_compile_options(${PROJECT_NAME}.${_ji_name} PRIVATE -UJIARRAY_DEBUG)
          endif()
        endforeach()
      endforeach()
    endforeach()

    # Column-major (default) + row-major variants of the main matrix suite
    _ji_add_gtest(column_major_test        ./test/column_major_test.cpp)
    _ji_add_gtest(row_major_test           ./test/row_major_test.cpp)
    target_compile_options(${PROJECT_NAME}.row_major_test
      PRIVATE -DJIARRAY_COLUMN_MAJOR=0)

    # JIDeviceArray suites need nvcc and a GPU, so they are opt-in.
    option(JIARRAY_WITH_CUDA_TESTS "Build the JIDeviceArray CUDA tests and benchmark" OFF)
    if(JIARRAY_WITH_CUDA_TESTS)
      if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        if(CMAKE_VERSION VERSION_LESS 3.24)
          message(FATAL_ERROR "Set CMAKE_CUDA_ARCHITECTURES (e.g. 89) for JIARRAY_WITH_CUDA_TESTS")
        endif()
        set(CMAKE_CUDA_ARCHITECTURES native)
      endif()
      enable_language(CUDA)
      find_package(CUDAToolkit REQUIRED)

      function(_ji_add_cuda_gtest name src)
        add_executable(${PROJECT_NAME}.${name} ${src})
        target_compile_definitions(${PROJECT_NAME}.${name} PRIVATE JIARRAY_DEBUG)
        set_target_properties(${PROJECT_NAME}.${name} PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
        target_link_libraries(${PROJECT_NAME}.${name}
          PRIVATE jiarray CUDA::cudart GTest::gtest_main GTest::gtest)
        # The leak check reads device-wide free memory, so device tests must not run concurrently.
        gtest_discover_tests(${PROJECT_NAME}.${name} PROPERTIES RESOURCE_LOCK cuda_device)
      endfunction()

      _ji_add_cuda_gtest(device_array_test          ./test/device_array_test.cu)
      _ji_add_cuda_gtest(device_array_rowmajor_test ./test/device_array_test.cu)
      target_compile_definitions(${PROJECT_NAME}.device_array_rowmajor_test
        PRIVATE JIARRAY_COLUMN_MAJOR=0)
      # Same header through the host compiler only (no nvcc).
      _ji_add_cuda_gtest(device_array_host_test     ./test/device_array_host_test.cpp)

      # init0 and device addressing on JICudaArray / JIDeviceArray: both index bases x both storage
      # orders, each with and without JIARRAY_DEBUG.
      foreach(_ji_base 0 1)
        foreach(_ji_order col row)
          foreach(_ji_debug debug nodebug)
            set(_ji_name cuda_array_init0_base${_ji_base}_${_ji_order}_${_ji_debug}_test)
            _ji_add_cuda_gtest(${_ji_name} ./test/cuda_array_init0_test.cu)
            target_compile_definitions(${PROJECT_NAME}.${_ji_name} PRIVATE JIARRAY_OFFSET=${_ji_base})
            if(_ji_order STREQUAL "row")
              target_compile_definitions(${PROJECT_NAME}.${_ji_name} PRIVATE JIARRAY_COLUMN_MAJOR=0)
            endif()
            if(_ji_debug STREQUAL "nodebug")
              get_target_property(_ji_defs ${PROJECT_NAME}.${_ji_name} COMPILE_DEFINITIONS)
              list(REMOVE_ITEM _ji_defs JIARRAY_DEBUG)
              set_target_properties(${PROJECT_NAME}.${_ji_name} PROPERTIES COMPILE_DEFINITIONS "${_ji_defs}")
              target_compile_options(${PROJECT_NAME}.${_ji_name} PRIVATE -UJIARRAY_DEBUG)
            endif()
          endforeach()
        endforeach()
      endforeach()

      # Raw-pointer vs view timing; not a ctest (no JIARRAY_DEBUG, run by hand).
      add_executable(${PROJECT_NAME}.device_array_bench ./test/device_array_bench.cu)
      set_target_properties(${PROJECT_NAME}.device_array_bench PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
      target_compile_options(${PROJECT_NAME}.device_array_bench PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 -lineinfo>)
      target_link_libraries(${PROJECT_NAME}.device_array_bench PRIVATE jiarray CUDA::cudart)
    endif()

  endif()
endif()
