#pragma once

#include "JIArray.h"
#include <vector>

namespace dnegri::jiarray {

template <typename _Tp, typename _Alloc = std::allocator<_Tp>>
class JIVector {
private:
    typedef typename std::vector<_Tp>::reference              reference;
    typedef typename std::vector<_Tp>::const_reference        const_reference;
    typedef typename std::vector<_Tp, _Alloc>::const_iterator const_iterator;

    std::vector<_Tp, _Alloc> vec_impl;

public:
    JIVector() {
    }

    JIVector(std::initializer_list<_Tp> il) {
        vec_impl = il;
    }

    reference front() {
        return vec_impl.front();
    }

    const_reference front() const {
        return vec_impl.front();
    }

    reference operator[](int index) {
        return vec_impl[index - JIARRAY_OFFSET];
    }

    void push_back(const _Tp& value) {
        vec_impl.push_back(value);
    }

    reference at(int index) {
        return vec_impl.at(index - JIARRAY_OFFSET);
    }

    const_reference at(int index) const {
        return vec_impl.at(index - JIARRAY_OFFSET);
    }

    reference operator()(int index) {
        return this->at(index);
    }

    const_reference operator()(int index) const {
        return this->at(index);
    }

    size_t size() const {
        return vec_impl.size();
    }

    void resize(size_t count) {
        vec_impl.resize(count);
    }

    const_iterator begin() const {
        return vec_impl.begin();
    }

    const_iterator end() const {
        return vec_impl.end();
    }

    const std::vector<_Tp>& std() const {
        return vec_impl;
    }

    const _Tp* data() const {
        return vec_impl.data();
    }

    void operator=(const JIArray<_Tp, 1>& other) {
        ffor(i, 1, other.size()) {
            this->push_back(other(i));
        }
    }

    const bool contains(const _Tp& value) const {
        return std::find(vec_impl.begin(), vec_impl.end(), value) != vec_impl.end();
    }
};
template <typename T>
using zvector = JIVector<T>;

using zstrings = JIVector<std::string>;
using zdoubles = JIVector<double>;
using zfloats = JIVector<float>;
}; // namespace dnegri::jiarray