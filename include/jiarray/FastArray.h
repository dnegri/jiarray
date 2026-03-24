#pragma once

#include "pch.h"

namespace dnegri::jiarray {

template <class T, std::size_t SIZE, std::size_t OFFSET = JIARRAY_OFFSET>
class FastArray {
public:
    T mm[SIZE]{};

public:
    FastArray() {
    }

    FastArray(const FastArray<T, SIZE, OFFSET>& array) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = array.mm[i];
        }
    };

    explicit FastArray(std::initializer_list<T> arrays) {
        int idx = -1;
        for (const auto& value : arrays) {
            mm[++idx] = value;
        }
    }

    explicit FastArray(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = val;
        }
    }

    template <size_t array_size>
    explicit FastArray(T (&a)[array_size]) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = a[i];
        }
    }

    inline int findFirst(const T& val) {
        int pos = OFFSET - 1;
        for (size_t i = 0; i < SIZE; i++) {
            if (mm[i] == val) {
                pos = i + OFFSET;
                break;
            }
        }

        return pos;
    }

    inline T min() const {
        T mx = mm[0];
        for (int i = 1; i < SIZE; ++i) {
            if (mm[i] < mx)
                mx = mm[i];
        }
        return mx;
    }

    inline T max() const {
        T mx = mm[0];
        for (int i = 1; i < SIZE; ++i) {
            if (mm[i] > mx)
                mx = mm[i];
        }
        return mx;
    }

    inline T& operator()(int i) {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return mm[i - OFFSET];
    }
    inline const T& operator()(int i) const {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return mm[i - OFFSET];
    }

    inline FastArray<T, SIZE, OFFSET>& operator=(const FastArray<T, SIZE, OFFSET>& array) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = array.mm[i];
        }
        return *this;
    }

    inline FastArray<T, SIZE, OFFSET>& operator=(const T* val) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = val[i];
        }
        return *this;
    }

    inline FastArray<T, SIZE, OFFSET>& operator=(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = val;
        }
        return *this;
    }

    inline FastArray<T, SIZE, OFFSET>& operator+=(const T& val) {
        for (size_t i = 0; i < SIZE; ++i) {
            this->mm[i] += val;
        }
        return *this;
    }

    inline FastArray<T, SIZE, OFFSET>& operator*=(const T& val) {
        for (size_t i = 0; i < SIZE; ++i) {
            this->mm[i] *= val;
        }
        return *this;
    }

    inline FastArray<T, SIZE, OFFSET>& operator/=(const T& val) {
        for (size_t i = 0; i < SIZE; ++i) {
            this->mm[i] /= val;
        }
        return *this;
    }

    inline FastArray<T, SIZE, OFFSET> operator*(const T& val) {
        FastArray<T, SIZE, OFFSET> result = *this;
        for (size_t i = 0; i < SIZE; ++i) {
            result.mm[i] *= val;
        }
        return result;
    }

    // ========== Iterator support ==========
    using value_type             = T;
    using size_type              = std::size_t;
    using difference_type        = std::ptrdiff_t;
    using reference              = T&;
    using const_reference        = const T&;
    using pointer                = T*;
    using const_pointer          = const T*;
    using iterator               = T*;
    using const_iterator         = const T*;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    // Iterator functions
    iterator begin() noexcept {
        return mm;
    }
    const_iterator begin() const noexcept {
        return mm;
    }
    const_iterator cbegin() const noexcept {
        return mm;
    }

    iterator end() noexcept {
        return mm + SIZE;
    }
    const_iterator end() const noexcept {
        return mm + SIZE;
    }
    const_iterator cend() const noexcept {
        return mm + SIZE;
    }

    reverse_iterator rbegin() noexcept {
        return reverse_iterator(end());
    }
    const_reverse_iterator rbegin() const noexcept {
        return const_reverse_iterator(end());
    }
    const_reverse_iterator crbegin() const noexcept {
        return const_reverse_iterator(end());
    }

    reverse_iterator rend() noexcept {
        return reverse_iterator(begin());
    }
    const_reverse_iterator rend() const noexcept {
        return const_reverse_iterator(begin());
    }
    const_reverse_iterator crend() const noexcept {
        return const_reverse_iterator(begin());
    }

    // Size functions
    constexpr size_type size() const noexcept {
        return SIZE;
    }
    constexpr size_type max_size() const noexcept {
        return SIZE;
    }
    constexpr bool empty() const noexcept {
        return SIZE == 0;
    }

    // Element access (0-based, for STL compatibility)
    reference operator[](size_type i) noexcept {
        return mm[i];
    }
    const_reference operator[](size_type i) const noexcept {
        return mm[i];
    }

    reference front() noexcept {
        return mm[0];
    }
    const_reference front() const noexcept {
        return mm[0];
    }

    reference back() noexcept {
        return mm[SIZE - 1];
    }
    const_reference back() const noexcept {
        return mm[SIZE - 1];
    }

    pointer data() noexcept {
        return mm;
    }
    const_pointer data() const noexcept {
        return mm;
    }
};

template <class T, std::size_t SIZE1, std::size_t SIZE2, std::size_t OFFSET = JIARRAY_OFFSET>
class FastArray2D {
public:
    const int SIZE = SIZE1 * SIZE2;
    T         mm[SIZE1 * SIZE2]{};

public:
    FastArray2D() {
    }

    FastArray2D(std::initializer_list<T> arrays) {
        int idx = -1;
        for (const auto& value : arrays) {
            mm[++idx] = value;
        }
    }

    FastArray2D(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = val;
        }
    }

    template <size_t array_size>
    FastArray2D(T (&a)[array_size]) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = a[i];
        }
    }

    inline T& operator()(int i, int j) {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE1 - 1);
        JIARRAY_CHECK_BOUND(j, OFFSET, OFFSET + SIZE2 - 1);
        return mm[(j - OFFSET) * SIZE1 + i - OFFSET];
    }
    inline const T& operator()(int i, int j) const {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE1 - 1);
        JIARRAY_CHECK_BOUND(j, OFFSET, OFFSET + SIZE2 - 1);
        return mm[(j - OFFSET) * SIZE1 + i - OFFSET];
    }
    inline FastArray2D<T, SIZE1, SIZE2, OFFSET>&
    operator=(const FastArray2D<T, SIZE1, SIZE2, OFFSET>& array) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = array.mm[i];
        }
        return *this;
    }
    inline FastArray2D<T, SIZE1, SIZE2, OFFSET>& operator=(const T& val) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = val;
        }
        return *this;
    }
};

template <std::size_t SIZE, std::size_t OFFSET = JIARRAY_OFFSET>
class StringFastArray {
public:
    std::string mm[SIZE];

public:
    StringFastArray() {
    }
    // template<size_t array_size>
    // template<typename... Ts>
    StringFastArray(std::initializer_list<const char*> strings) {
        int idx = -1;
        for (const auto& string : strings) {
            mm[++idx] = string;
        }
    }

    inline std::string& operator()(int i) {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return mm[i - OFFSET];
    }
    inline const std::string& operator()(int i) const {
        JIARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return mm[i - OFFSET];
    }
    inline StringFastArray<SIZE, OFFSET>& operator=(const StringFastArray<SIZE, OFFSET>& array) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = array.mm[i];
        }
        return *this;
    }
    inline StringFastArray<SIZE, OFFSET>& operator=(const std::string& val) {
        for (size_t i = 0; i < SIZE; i++) {
            mm[i] = val;
        }
        return *this;
    }
};

template <int I, int OFFSET = JIARRAY_OFFSET>
using fbool1d = dnegri::jiarray::FastArray<bool, I>;

template <int I, int OFFSET = JIARRAY_OFFSET>
using fint1d = dnegri::jiarray::FastArray<int, I>;

template <int I, int OFFSET = JIARRAY_OFFSET>
using ffloat1d = dnegri::jiarray::FastArray<float, I>;

// template <int I> using fdouble1d = dnegri::jiarray::FastArray<double, I>;
template <int I, int OFFSET = JIARRAY_OFFSET>
using fdouble1d = dnegri::jiarray::FastArray<double, I, OFFSET>;

template <int I, int OFFSET = JIARRAY_OFFSET>
using fstring1d = dnegri::jiarray::StringFastArray<I, OFFSET>;

template <int I, int J, int OFFSET = JIARRAY_OFFSET>
using fbool2d = dnegri::jiarray::FastArray2D<bool, I, J, OFFSET>;

template <int I, int J, int OFFSET = JIARRAY_OFFSET>
using fint2d = dnegri::jiarray::FastArray2D<int, I, J, OFFSET>;

template <int I, int J, int OFFSET = JIARRAY_OFFSET>
using ffloat2d = dnegri::jiarray::FastArray2D<float, I, J, OFFSET>;

template <int I, int J, int OFFSET = JIARRAY_OFFSET>
using fdouble2d = dnegri::jiarray::FastArray2D<double, I, J, OFFSET>;

template <class T, int I, int OFFSET = JIARRAY_OFFSET>
using farray = dnegri::jiarray::FastArray<T, I, OFFSET>;

} // namespace dnegri::jiarray