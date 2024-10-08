#pragma once

#include "pch.h"

namespace dnegri::jiarray {

template <class T, std::size_t SIZE, std::size_t OFFSET = JIARRAY_OFFSET>
class FastArray {
public:
    T memory[SIZE]{};

public:
    FastArray() {
    }

    FastArray(std::initializer_list<T> arrays) {
        int idx = -1;
        for (const auto& value : arrays) {
            memory[++idx] = value;
        }
    }

    FastArray(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = val;
        }
    }

    template <size_t array_size>
    FastArray(T (&a)[array_size]) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = a[i];
        }
    }

    inline int findFirst(const T& val) {
        int pos = OFFSET - 1;
        for (size_t i = 0; i < SIZE; i++) {
            if (memory[i] == val) {
                pos = i + OFFSET;
                break;
            }
        }

        return pos;
    }

    inline T min() const {
        T mx = memory[0];
        for (int i = 1; i < SIZE; ++i) {
            if (memory[i] < mx) mx = memory[i];
        }
        return mx;
    }

    inline T max() const {
        T mx = memory[0];
        for (int i = 1; i < SIZE; ++i) {
            if (memory[i] > mx) mx = memory[i];
        }
        return mx;
    }


    inline T& operator()(int i) {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return memory[i - OFFSET];
    }
    inline const T& operator()(int i) const {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return memory[i - OFFSET];
    }
    inline FastArray<T, SIZE, OFFSET>& operator=(const FastArray<T, SIZE, OFFSET>& array) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = array.memory[i];
        }
        return *this;
    }
    inline FastArray<T, SIZE, OFFSET>& operator=(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = val;
        }
        return *this;
    }
};

template <class T, std::size_t SIZE1, std::size_t SIZE2, std::size_t OFFSET = JIARRAY_OFFSET>
class FastArray2D {
public:
    const int SIZE = SIZE1 * SIZE2;
    T         memory[SIZE1 * SIZE2]{};

public:
    FastArray2D() {
    }

    FastArray2D(std::initializer_list<T> arrays) {
        int idx = -1;
        for (const auto& value : arrays) {
            memory[++idx] = value;
        }
    }

    FastArray2D(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = val;
        }
    }

    template <size_t array_size>
    FastArray2D(T (&a)[array_size]) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = a[i];
        }
    }

    inline T& operator()(int i, int j) {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE1 - 1);
        JIARRAY_CHECK_BOUND(j, OFFSET, OFFSET + SIZE2 - 1);
        return memory[(j - OFFSET) * SIZE1 + i - OFFSET];
    }
    inline const T& operator()(int i, int j) const {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE1 - 1);
        JIARRAY_CHECK_BOUND(j, OFFSET, OFFSET + SIZE2 - 1);
        return memory[(j - OFFSET) * SIZE1 + i - OFFSET];
    }
    inline FastArray2D<T, SIZE1, SIZE2, OFFSET>& operator=(const FastArray2D<T, SIZE1, SIZE2, OFFSET>& array) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = array.memory[i];
        }
        return *this;
    }
    inline FastArray2D<T, SIZE1, SIZE2, OFFSET>& operator=(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = val;
        }
        return *this;
    }
};

template <std::size_t SIZE, std::size_t OFFSET = JIARRAY_OFFSET>
class StringFastArray {
public:
    std::string memory[SIZE];

public:
    StringFastArray() {
    }
    // template<size_t array_size>
    // template<typename... Ts>
    StringFastArray(std::initializer_list<const char*> strings) {
        int idx = -1;
        for (const auto& string : strings) {
            memory[++idx] = string;
        }
    }

    inline std::string& operator()(int i) {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return memory[i - OFFSET];
    }
    inline const std::string& operator()(int i) const {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return memory[i - OFFSET];
    }
    inline StringFastArray<SIZE, OFFSET>& operator=(const StringFastArray<SIZE, OFFSET>& array) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = array.memory[i];
        }
        return *this;
    }
    inline StringFastArray<SIZE, OFFSET>& operator=(const std::string& val) {
        for (size_t i = 0; i < SIZE; i++) {
            memory[i] = val;
        }
        return *this;
    }
};

template <int I>
using bool1d = dnegri::jiarray::FastArray<bool, I>;

template <int I>
using int1d = dnegri::jiarray::FastArray<int, I>;

template <int I>
using float1d = dnegri::jiarray::FastArray<float, I>;

template <int I>
using double1d = dnegri::jiarray::FastArray<double, I>;

template <int I>
using string1d = dnegri::jiarray::StringFastArray<I>;

template <int I, int J>
using bool2d = dnegri::jiarray::FastArray<bool, I, J>;

template <int I, int J>
using int2d = dnegri::jiarray::FastArray<int, I, J>;

template <int I, int J>
using float2d = dnegri::jiarray::FastArray<float, I, J>;

template <int I, int J>
using double2d = dnegri::jiarray::FastArray<double, I, J>;

template <class T, int I>
using farray = dnegri::jiarray::FastArray<T, I>;

} // namespace dnegri::jiarray