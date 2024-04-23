#pragma once

#include <assert.h>

#include <algorithm>
#include <cstddef> // For std::ptrdiff_t
#include <iterator>
#include <iterator> // For std::forward_iterator_tag
#include <string>
#include <tuple>
#include <utility>

namespace dnegri::jiarray {

#ifndef JIARRAY_OFFSET
    #define JIARRAY_OFFSET 1
#endif

#ifdef JIARRAY_DEBUG
    #define JIARRAY_CHECK_NOT_ALLOCATED()    assert(allocated == false)
    #define JIARRAY_CHECK_BOUND(i, beg, end) assert(i >= beg && i <= end)
    #define JIARRAY_CHECK_RANK(r1, r2)       assert(r1 == r2)
    #define JIARRAY_CHECK_SIZE(r1, r2)       assert(r1 == r2)
#else
    #define JIARRAY_CHECK_NOT_ALLOCATED()
    #define JIARRAY_CHECK_BOUND(i, beg, end)
    #define JIARRAY_CHECK_RANK(r1, r2)
    #define JIARRAY_CHECK_SIZE(r1, r2)
#endif

#define JIARRAY_ALLOCATED_NONE            0
#define JIARRAY_ALLOCATED_MEMORY          1
#define JIARRAY_ALLOCATED_RANKSIZE        2
#define JIARRAY_ALLOCATED_OFFSET          4
#define JIARRAY_ALLOCATED_ALL             7
#define JIARRAY_ALLOCATED_RANKSIZE_OFFSET 6
}; // namespace dnegri::jiarray