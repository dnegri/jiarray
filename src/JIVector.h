#pragma once

#include "JIArray.h"
#include <algorithm>
#include <initializer_list>
#include <vector>

namespace dnegri::jiarray {

template <typename _Tp, typename _Alloc = std::allocator<_Tp>>
class JIVector : public std::vector<_Tp, _Alloc> {
public:
    using Base = std::vector<_Tp, _Alloc>;
    using typename Base::const_iterator;
    using typename Base::const_reference;
    using typename Base::iterator;
    using typename Base::reference;

    // Constructors
    JIVector() = default;

    JIVector(std::initializer_list<_Tp> il)
        : Base(il) {
    }

    JIVector(const std::vector<_Tp>& vec)
        : Base(vec) {
    }

    JIVector(size_t size, const _Tp& value = _Tp())
        : Base(size, value) {
    }

    // Accessors with custom offset
    reference operator[](int index) {
        return Base::operator[](index - JIARRAY_OFFSET);
    }

    const_reference operator[](int index) const {
        return Base::operator[](index - JIARRAY_OFFSET);
    }

    reference at(int index) {
        return Base::at(index - JIARRAY_OFFSET);
    }

    const_reference at(int index) const {
        return Base::at(index - JIARRAY_OFFSET);
    }

    reference operator()(int index) {
        return this->at(index);
    }

    const_reference operator()(int index) const {
        return this->at(index);
    }

    _Tp sum() const {
        return std::accumulate(Base::begin(), Base::end(), 0.0);
    }

    iterator begin() {
        return Base::begin();
    }

    const_iterator begin() const {
        return Base::begin();
    }

    iterator end() {
        return Base::end();
    }

    const_iterator end() const {
        return Base::end();
    }

    // Memory access
    _Tp* getMemory(int index) {
        return &Base::at(index);
    }

    _Tp* get_pointer() {
        return Base::data();
    }

    // Assignment operators
    void operator=(const std::vector<_Tp>& other) {
        Base::operator=(other);
    }

    void operator=(const JIArray<_Tp, 1>& other) {
        ffor(i, 1, other.size()) {
            this->push_back(other(i));
        }
    }

    void operator*=(const _Tp& value) {
        for (auto& val : *this) {
            val *= value;
        }
    }

    void operator/=(const _Tp& value) {
        for (auto& val : *this) {
            val /= value;
        }
    }

    int count(const _Tp& value) const {
        return std::count(Base::begin(), Base::end(), value);
    }

    // Contains method
    bool contains(const _Tp& value) const {
        return std::find(Base::begin(), Base::end(), value) != Base::end();
    }

    int index(const _Tp& value) const {
        auto it = std::find(Base::begin(), Base::end(), value);
        if (it != Base::end()) {
            return std::distance(Base::begin(), it) + JIARRAY_OFFSET;
        }
        return -1 + JIARRAY_OFFSET; // Not found
    }

    // Insert methods
    void insert(const_iterator iter, std::initializer_list<_Tp> il) {
        Base::insert(iter, il);
    }

    void insert(const_iterator iter, const JIVector<_Tp>& other) {
        Base::insert(iter, other.begin(), other.end());
    }

    void insert(const_iterator iter, const std::vector<_Tp>& other) {
        Base::insert(iter, other.begin(), other.end());
    }

    void insert(const_iterator iter, const JIArray<_Tp>& other) {
        Base::insert(iter, other.begin(), other.end());
    }

    template <typename _InputIterator, typename = std::_RequireInputIter<_InputIterator>>
    void insert(const_iterator iter, _InputIterator first, _InputIterator last) {
        Base::insert(iter, first, last);
    }
};

template <typename T>
using zvector = JIVector<T>;

using zstrings = JIVector<std::string>;
using zdoubles = JIVector<double>;
using zfloats  = JIVector<float>;
using zints    = JIVector<int>;

}; // namespace dnegri::jiarray