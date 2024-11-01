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

#ifndef JIARRAY_OFFSET
    #define JIARRAY_OFFSET 1
#endif

#ifdef JIARRAY_DEBUG
    #define JIARRAY_CHECK_NOT_ALLOCATED()    if(allocated) throw std::invalid_argument("Already allocated")
    #define JIARRAY_CHECK_BOUND(i, beg, end) if(i < beg || i > end) throw std::out_of_range("Index out of bounds")
    #define JIARRAY_CHECK_RANK(r1, r2)       if(r1 != r2) throw std::invalid_argument("Rank mismatch")
    #define JIARRAY_CHECK_SIZE(r1, r2)       if(r1 != r2) throw std::invalid_argument("Size mismatch")
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
