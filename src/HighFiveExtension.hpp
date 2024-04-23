#pragma once
#include <highfive/H5File.hpp>
#include <JIArray.h>

using namespace HighFive;
using namespace HighFive::details;

template <typename T, size_t N>
struct inspector<dnegri::jiarray::JIArray<T, N>> {
    using type = dnegri::jiarray::JIArray<T, N>;
    using value_type = unqualified_t<T>;
    using base_type = typename inspector<value_type>::base_type;
    using hdf5_type = typename inspector<value_type>::hdf5_type;

    static constexpr size_t ndim = 1;
    static constexpr size_t recursive_ndim = ndim;
    static constexpr bool is_trivially_copyable = std::is_trivially_copyable<value_type>::value &&
                                                  inspector<value_type>::is_trivially_copyable;

    static std::vector<size_t> getDimensions(const type& val) {
        std::vector<size_t> sizes;
        sizes.push_back(val.getSize());        
        return sizes;
    }

    static size_t getSizeVal(const type& val) {
        return val.getSize();
    }

    static size_t getSize(const std::vector<size_t>& dims) {
        return compute_total_size(dims);
    }

    static void prepare(type& val, const std::vector<size_t>& dims) {
        // val.setSize(dims);
        // std::vector<size_t> next_dims(dims.begin() + 1, dims.end());
        // for (auto&& e: val) {
        //     inspector<value_type>::prepare(e, next_dims);
        // }
    }

    static hdf5_type* data(type& val) {
        return inspector<value_type>::data(*val.data());
    }

    static const hdf5_type* data(const type& val) {
        return inspector<value_type>::data(*val.data());
    }

    static void serialize(const type& val, hdf5_type* m) {
        size_t subsize = inspector<value_type>::getSizeVal(*val.data());
        ffor(i, 1, val.getSize()) {
            inspector<value_type>::serialize(val(i), m);
            m += subsize;
        }
    }

    static void unserialize(const hdf5_type* vec_align,
                            const std::vector<size_t>& dims,
                            type& val) {
        std::vector<size_t> next_dims(dims.begin() + 1, dims.end());
        size_t next_size = compute_total_size(next_dims);
        for (size_t i = 0; i < dims[0]; ++i) {
            inspector<value_type>::unserialize(vec_align + i * next_size, next_dims, *val.getMemory(i));
        }
    }
};

template <typename T>
struct inspector<dnegri::jiarray::JIVector<T>> {
    using type = dnegri::jiarray::JIVector<T>;
    using value_type = unqualified_t<T>;
    using base_type = typename inspector<value_type>::base_type;
    using hdf5_type = typename inspector<value_type>::hdf5_type;

    static constexpr size_t ndim = 1;
    static constexpr size_t recursive_ndim = ndim;
    static constexpr bool is_trivially_copyable = std::is_trivially_copyable<value_type>::value &&
                                                  inspector<value_type>::is_trivially_copyable;

    static std::vector<size_t> getDimensions(const type& val) {
        std::vector<size_t> sizes;
        sizes.push_back(val.size());        
        return sizes;
    }

    static size_t getSizeVal(const type& val) {
        return val.size();
    }

    static size_t getSize(const std::vector<size_t>& dims) {
        return compute_total_size(dims);
    }

    static void prepare(type& val, const std::vector<size_t>& dims) {
        // val.setSize(dims);
        // std::vector<size_t> next_dims(dims.begin() + 1, dims.end());
        // for (auto&& e: val) {
        //     inspector<value_type>::prepare(e, next_dims);
        // }
    }

    static hdf5_type* data(type& val) {
        return inspector<value_type>::data(*val.data());
    }

    static const hdf5_type* data(const type& val) {
        return inspector<value_type>::data(*val.data());
    }

    static void serialize(const type& val, hdf5_type* m) {
        size_t subsize = inspector<value_type>::getSizeVal(*val.data());
        ffor(i, 1, val.size()) {
            inspector<value_type>::serialize(val(i), m);
            m += subsize;
        }
    }

    static void unserialize(const hdf5_type* vec_align,
                            const std::vector<size_t>& dims,
                            type& val) {
        std::vector<size_t> next_dims(dims.begin() + 1, dims.end());
        size_t next_size = compute_total_size(next_dims);
        for (size_t i = 0; i < dims[0]; ++i) {
            inspector<value_type>::unserialize(vec_align + i * next_size, next_dims, *val.getMemory(i));
        }
    }
};
