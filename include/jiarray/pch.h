#pragma once

#include <assert.h>

#include <algorithm>
#include <cstddef> // For std::ptrdiff_t
#include <iterator>
#include <iterator> // For std::forward_iterator_tag
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace dnegri::jiarray {

// CUDA host/device annotation macro
#ifdef __CUDACC__
  #define JIARRAY_HD __host__ __device__
#else
  #define JIARRAY_HD
#endif

#ifndef JIARRAY_OFFSET
#define JIARRAY_OFFSET 1
#endif

#ifndef JIARRAY_COLUMN_MAJOR
#define JIARRAY_COLUMN_MAJOR 1
#endif

// Portable loop-unroll pragma. Issued inside functions to hint the
// compiler. Silent no-op on compilers that don't understand the
// specific pragma (avoids -Wunknown-pragmas on gcc/clang/msvc).
#if defined(__CUDACC__) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  #define JIARRAY_UNROLL _Pragma("unroll")
#elif defined(__clang__)
  #define JIARRAY_UNROLL _Pragma("clang loop unroll(enable)")
#elif defined(__GNUC__) && __GNUC__ >= 8
  #define JIARRAY_UNROLL _Pragma("GCC unroll 16")
#else
  #define JIARRAY_UNROLL
#endif

// Portable SIMD loop pragma. Requires -fopenmp/-qopenmp to actually
// vectorize on gcc/clang/icc; MSVC <= 2022 lacks `omp simd` so we fall
// back to the native `loop(ivdep)` hint. Without OpenMP flag the
// directive is silently omitted (no -Wunknown-pragmas noise).
#if defined(_MSC_VER) && !defined(__clang__) && !defined(__INTEL_LLVM_COMPILER)
  #define JIARRAY_PRAGMA_STR(...) #__VA_ARGS__
  #define JIARRAY_SIMD_LOOP __pragma(loop(ivdep))
  #define JIARRAY_SIMD_REDUCTION(op, var) __pragma(loop(ivdep))
#elif defined(_OPENMP)
  #define JIARRAY_PRAGMA_STR(...) #__VA_ARGS__
  #define JIARRAY_SIMD_LOOP _Pragma("omp simd")
  #define JIARRAY_SIMD_REDUCTION(op, var) \
      _Pragma(JIARRAY_PRAGMA_STR(omp simd reduction(op : var)))
#else
  #define JIARRAY_PRAGMA_STR(...) #__VA_ARGS__
  #define JIARRAY_SIMD_LOOP
  #define JIARRAY_SIMD_REDUCTION(op, var)
#endif

#if defined(JIARRAY_DEBUG) && !defined(__CUDA_ARCH__)
#define JIARRAY_CHECK_NOT_ALLOCATED()                                                              \
    if (allocated)                                                                                 \
    throw std::invalid_argument("Already allocated")
#define JIARRAY_CHECK_BOUND(i, beg, end)                                                           \
    if (i < beg || i > end)                                                                        \
    throw std::out_of_range("Index out of bounds")
#define JIARRAY_CHECK_RANK(r1, r2)                                                                 \
    if (r1 != r2)                                                                                  \
    throw std::invalid_argument("Rank mismatch")
#define JIARRAY_CHECK_SIZE(r1, r2)                                                                 \
    if (r1 != r2)                                                                                  \
    throw std::invalid_argument("Size mismatch")
#else
#define JIARRAY_CHECK_NOT_ALLOCATED()
#define JIARRAY_CHECK_BOUND(i, beg, end)
#define JIARRAY_CHECK_RANK(r1, r2)
#define JIARRAY_CHECK_SIZE(r1, r2)
#endif

#define JIARRAY_ALLOCATED_NONE 0
#define JIARRAY_ALLOCATED_MEMORY 1
#define JIARRAY_ALLOCATED_RANKSIZE 2
#define JIARRAY_ALLOCATED_OFFSET 4
#define JIARRAY_ALLOCATED_ALL 7
#define JIARRAY_ALLOCATED_RANKSIZE_OFFSET 6
}; // namespace dnegri::jiarray
