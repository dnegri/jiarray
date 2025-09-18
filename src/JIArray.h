#pragma once

#include "FastArray.h"
#include "pch.h"
#include <type_traits>
#include <vector>

namespace dnegri::jiarray {

template <class T, int RANK = 1, class = std::make_index_sequence<RANK>> class JIArray;

template <class T, int RANK, int... INTS> class JIArray<T, RANK, std::index_sequence<INTS...>> {
  private:
    constexpr static bool is_row_major = (JIARRAY_COLUMN_MAJOR == 0);

    int nn = 0;
    T *mm = nullptr;
    int allocated = JIARRAY_ALLOCATED_NONE;
    int rankSize[RANK]{};
    int offset[RANK]{};
    int sumOfOffset = 0;
    int sizes[RANK]{};

    // Index calculation helpers - optimized for each storage order
    template <typename... Args> inline int calculateIndex(Args... indices) const {
        static_assert(sizeof...(indices) == RANK, "Number of indices must match array rank");

        int idx[] = {static_cast<int>(indices)...};
        int pos = -sumOfOffset;

        for (int i = 0; i < RANK; i++) {
            JIARRAY_CHECK_BOUND(idx[i], offset[i], offset[i] + sizes[i] - 1);
            pos += rankSize[i] * idx[i];
        }

        return pos;
    }

    // Calculate rank sizes based on storage order for optimal memory layout
    void calculateRankSize(const int dimensions[RANK]) {
        if constexpr (is_row_major) {
            // Row-major: rightmost dimension has stride 1
            rankSize[RANK - 1] = 1;
            for (int i = RANK - 2; i >= 0; i--) {
                rankSize[i] = rankSize[i + 1] * dimensions[i + 1];
            }
        } else {
            // Column-major: leftmost dimension has stride 1 (original behavior)
            rankSize[0] = 1;
            for (int i = 1; i < RANK; i++) {
                rankSize[i] = rankSize[i - 1] * dimensions[i - 1];
            }
        }
    }

  public:
    JIArray() = default;

    JIArray(std::initializer_list<T> il) : JIArray(static_cast<int>(il.size())) {
        std::copy(il.begin(), il.end(), mm);
    }

    JIArray(const JIArray<T, RANK> &array) {
        nn = array.nn;
        mm = array.mm;
        std::copy(array.rankSize, array.rankSize + RANK, rankSize);
        std::copy(array.offset, array.offset + RANK, offset);
        allocated = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
        std::copy(array.sizes, array.sizes + RANK, sizes);
    }

    virtual ~JIArray() { destroy(); }

    template <typename... Args, typename = std::enable_if_t<(std::is_integral_v<Args> && ...)>>
    JIArray(Args... args) {
        init(args...);
    }

    JIArray(T *memory_, decltype(INTS)... args) { init(args..., memory_); }

    // Initialize with range specifications (min, max pairs)
    template <typename... INTS2, typename = std::enable_if_t<(sizeof...(INTS2) == 2 * RANK)>>
    void init0(INTS2... array_sizes) {
        int temp_size[] = {static_cast<int>(array_sizes)...};

        int dimensions[RANK];
        int offsets[RANK];

        for (int i = 0; i < RANK; i++) {
            dimensions[i] = temp_size[2 * i + 1] - temp_size[2 * i] + 1;
            offsets[i] = temp_size[2 * i];
        }

        init(dimensions);
        setOffsets(offsets);
    }

    template <typename... Args, typename = std::enable_if_t<(std::is_integral_v<Args> && ...)>,
              typename = std::enable_if_t<(sizeof...(Args) == RANK)>>
    void init(Args... array_sizes) {
        int dimensions[] = {static_cast<int>(array_sizes)...};
        init(dimensions);
    }

    void init(int dimensions[RANK]) {
        destroy();

        for (int i = 0; i < RANK; i++) {
            if (dimensions[i] == 0)
                return;
        }

        // Calculate rank sizes based on storage order for optimal performance
        calculateRankSize(dimensions);

        // Initialize offsets
        for (int i = 0; i < RANK; i++) {
            offset[i] = JIARRAY_OFFSET;
        }

        // Calculate total size and sum of offsets
        nn = 1;
        for (int i = 0; i < RANK; i++) {
            nn *= dimensions[i];
        }

        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++) {
            sumOfOffset += offset[i] * rankSize[i];
        }

        mm = new T[nn]{};
        allocated = JIARRAY_ALLOCATED_ALL;
        std::copy(dimensions, dimensions + RANK, sizes);
    }

    void initByRankSize(int size, const int *rankSizes, const int *offsets, T *memory) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        nn = size;
        sumOfOffset = 0;

        for (int i = 0; i < RANK; i++) {
            rankSize[i] = rankSizes[i];
            offset[i] = offsets[i];
            sumOfOffset += offset[i] * rankSize[i];
        }

        mm = memory;
        allocated = JIARRAY_ALLOCATED_NONE;

        // Reconstruct sizes from rank information based on storage order
        if constexpr (is_row_major) {
            sizes[0] = (nn != 0 && rankSize[0] != 0) ? nn / rankSize[0] : 0;

            for (int i = 1; i < RANK; i++) {
                sizes[i] = (rankSize[i] != 0) ? rankSize[i - 1] / rankSize[i] : 0;
            }
        } else {
            // Column-major (original behavior)
            for (int i = 0; i < RANK - 1; i++) {
                sizes[i] = (rankSize[i] != 0) ? rankSize[i + 1] / rankSize[i] : 0;
            }
            sizes[RANK - 1] = (nn != 0 && rankSize[RANK - 1] != 0) ? nn / rankSize[RANK - 1] : 0;
        }
    }

    void initByRankSize(int size, const int *rankSizes, const int *offsets) {
        destroy();
        if (size == 0)
            return;

        mm = new T[size];
        initByRankSize(size, rankSizes, offsets, mm);
        allocated = JIARRAY_ALLOCATED_ALL;
    }

    void init(decltype(INTS)... array_sizes, T *memory_) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        int dimensions[] = {static_cast<int>(array_sizes)...};
        nn = (array_sizes * ...);
        mm = memory_;

        calculateRankSize(dimensions);

        for (int i = 0; i < RANK; i++) {
            offset[i] = JIARRAY_OFFSET;
        }

        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++) {
            sumOfOffset += offset[i] * rankSize[i];
        }

        allocated = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
        std::copy(dimensions, dimensions + RANK, sizes);
    }

    void destroy() {
        if (allocated != JIARRAY_ALLOCATED_NONE) {
            if ((allocated & JIARRAY_ALLOCATED_MEMORY) != 0 && mm != nullptr) {
                delete[] mm;
            }
            mm = nullptr;
            allocated = JIARRAY_ALLOCATED_NONE;
            nn = 0;
        }
    }

    void erase() {
        destroy();
        std::fill(rankSize, rankSize + RANK, 0);
        std::fill(offset, offset + RANK, 0);
        std::fill(sizes, sizes + RANK, 0);
        sumOfOffset = 0;
    }

    inline void setOffsets(decltype(INTS)... offsets) {
        int tempOffsets[] = {static_cast<int>(offsets)...};
        setOffsets(tempOffsets);
    }

    inline void setOffsets(int offsets[RANK]) {
        std::copy(offsets, offsets + RANK, offset);
        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++) {
            sumOfOffset += rankSize[i] * offset[i];
        }
    }

    void setSize(decltype(INTS)... array_sizes) {
        int dimensions[] = {static_cast<int>(array_sizes)...};

        // Check if sizes are the same
        bool same = true;
        for (int i = 0; i < RANK; i++) {
            if (sizes[i] != dimensions[i]) {
                same = false;
                break;
            }
        }

        if (same)
            return;

        destroy();
        init(array_sizes...);
    }

    // Storage-order-aware slicing
    template <typename... INDEX>
    inline JIArray<T, RANK - sizeof...(INDEX)> slice(INDEX... index) const {
        constexpr int RANK2 = RANK - sizeof...(INDEX);
        constexpr int num_idx = sizeof...(INDEX);

        int idx[] = {static_cast<int>(index)...};
        int p_mm = 0;

        if constexpr (is_row_major) {
            // Row-major slicing: indices from left to right
            for (int i = 0; i < num_idx; i++) {
                JIARRAY_CHECK_BOUND(idx[i], offset[i], offset[i] + sizes[i] - 1);
                p_mm += rankSize[i] * (idx[i] - offset[i]);
            }

            auto array = JIArray<T, RANK2>();
            array.initByRankSize(rankSize[num_idx - 1], rankSize + num_idx, offset + num_idx,
                                 mm + p_mm);
            return array;
        } else {
            // Column-major slicing (original behavior)
            int rank = RANK;
            for (int i = num_idx - 1; i >= 0; i--) {
                rank--;
                JIARRAY_CHECK_BOUND(idx[i], offset[rank], offset[rank] + sizes[rank] - 1);
                p_mm += rankSize[rank] * (idx[i] - offset[rank]);
            }

            auto array = JIArray<T, RANK2>();
            array.initByRankSize(rankSize[RANK2], rankSize, offset, mm + p_mm);
            return array;
        }
    }

    // Accessors and properties
    inline bool isAllocated() const { return allocated != JIARRAY_ALLOCATED_NONE; }
    inline T *getMemory(int idx = 0) { return mm + idx; }
    inline const T *getMemory(int idx = 0) const { return mm + idx; }
    inline T *data() { return mm; }
    inline const T *data() const { return mm; }
    inline const int &size() const { return nn; }
    inline int getSize() const { return nn; }
    inline const int *getRankSize() const { return rankSize; }
    inline const int *getSizeOfRank() const { return sizes; }
    inline const int &getSize(int rank) const { return sizes[rank - JIARRAY_OFFSET]; }
    inline const int *getOffset() const { return offset; }
    inline const int &getOffset(int rank) const { return offset[rank - JIARRAY_OFFSET]; }

    // Optimized element access operators
    inline T &at(decltype(INTS)... index) {
        JIARRAY_CHECK_RANK(RANK, sizeof...(index));
        int pos = calculateIndex(index...);
        return mm[pos];
    }

    inline const T &at(decltype(INTS)... index) const {
        return const_cast<JIArray<T, RANK> &>(*this).at(index...);
    }

    inline const T &operator()(decltype(INTS)... index) const { return at(index...); }

    inline T &operator()(decltype(INTS)... index) { return at(index...); }

    inline T &operator()(const FastArray<int, RANK> &idx) {
        int pos = -sumOfOffset;

        for (int i = 0; i < RANK; i++) {
            JIARRAY_CHECK_BOUND(idx(i + 1), offset[i], offset[i] + sizes[i] - 1);
            pos += rankSize[i] * idx(i + 1);
        }

        return mm[pos];
    }

    // Storage-order-aware data access
    template <typename... INDEX> inline T *data(INDEX... index) {
        int num_idx = sizeof...(index);
        JIARRAY_CHECK_BOUND(num_idx, 1, RANK);

        int idx[] = {static_cast<int>(index)...};
        int pos = 0;

        if constexpr (is_row_major) {
            for (int i = 0; i < num_idx; i++) {
                pos += rankSize[i] * (idx[i] - offset[i]);
            }
        } else {
            // Column-major (original behavior)
            int rank = RANK;
            for (int i = num_idx - 1; i >= 0; i--) {
                --rank;
                pos += rankSize[rank] * (idx[i] - offset[rank]);
            }
        }

        return mm + pos;
    }

    template <typename... INDEX> inline const T *data(INDEX... index) const {
        return const_cast<JIArray<T, RANK> &>(*this).data(index...);
    }

    // Reshape with same storage order
    template <typename... INTS2> inline JIArray<T, sizeof...(INTS2)> reshape(INTS2... sizes) {
        constexpr int RANK2 = sizeof...(INTS2);
        JIArray<T, RANK2> array;
        array.init(sizes..., getMemory());
        return array;
    }

    // Assignment operators (remain mostly the same)
    inline JIArray<T, RANK> &operator=(const std::initializer_list<T> &list) {
        JIARRAY_CHECK_SIZE(nn, list.size());
        int idx = 0;
        for (const T &val : list) {
            mm[idx++] = val;
        }
        return *this;
    }

    inline JIArray<T, RANK> &operator=(const T &val) {
        assert(nn > 0);
        for (int i = 0; i < nn; ++i) {
            mm[i] = val;
        }
        return *this;
    }

    inline JIArray<T, RANK> &operator=(const std::vector<T> &val) {
        if (nn == 0) {
            int dimensions[RANK];
            std::fill(dimensions, dimensions + RANK, 1);
            dimensions[RANK - 1] = static_cast<int>(val.size());
            init(dimensions);
        } else {
            JIARRAY_CHECK_SIZE(nn, val.size());
        }

        int i = -1;
        for (auto &v : val) {
            mm[++i] = v;
        }
        return *this;
    }

    inline JIArray<T, RANK> &operator=(const T *array) {
        assert(nn > 0);
        for (int i = 0; i < nn; ++i) {
            mm[i] = array[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> &operator=(const JIArray<T, RANK> &array) {
        if (allocated == JIARRAY_ALLOCATED_NONE && mm == nullptr) {
            initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        } else {
            JIARRAY_CHECK_SIZE(nn, array.getSize());
            auto arraySizeOfRank = array.getSizeOfRank();
            auto arrayOffset = array.getOffset();
            for (int rank = 0; rank < RANK; ++rank) {
                assert(sizes[rank] == arraySizeOfRank[rank]);
                assert(offset[rank] == arrayOffset[rank]);
            }
        }

        auto arrayMemory = array.data();
        std::copy(arrayMemory, arrayMemory + nn, mm);
        return *this;
    }

    // Statistical operations (unchanged)
    inline T average() const {
        T result = T{};
        for (int i = 0; i < nn; ++i) {
            result += mm[i];
        }
        return result / nn;
    }

    inline T sum() const {
        T result = T{};
        for (int i = 0; i < nn; ++i) {
            result += mm[i];
        }
        return result;
    }

    inline T max() const {
        T mx = mm[0];
        for (int i = 1; i < nn; ++i) {
            if (mm[i] > mx)
                mx = mm[i];
        }
        return mx;
    }

    inline T min() const {
        T mn = mm[0];
        for (int i = 1; i < nn; ++i) {
            if (mm[i] < mn)
                mn = mm[i];
        }
        return mn;
    }

    // Arithmetic operators (updated to preserve ORDER)
    inline bool operator==(const JIArray<T, RANK> &array) const {
        if (nn != array.nn)
            return false;
        for (int i = 0; i < nn; ++i) {
            if (array.mm[i] != mm[i])
                return false;
        }
        return true;
    }

    inline JIArray<T, RANK> &operator+=(const T &val) {
        for (int i = 0; i < nn; ++i) {
            mm[i] += val;
        }
        return *this;
    }

    inline JIArray<T, RANK> &operator+=(const JIArray<T, RANK> &array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        for (int i = 0; i < nn; ++i) {
            mm[i] += array.mm[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> operator+(const JIArray<T, RANK> &array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIArray<T, RANK> result;
        result.initByRankSize(getSize(), getRankSize(), getOffset());
        for (int i = 0; i < nn; ++i) {
            result.mm[i] = mm[i] + array.mm[i];
        }
        return result;
    }

    inline JIArray<T, RANK> &operator-=(const JIArray<T, RANK> &array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        for (int i = 0; i < nn; ++i) {
            mm[i] -= array.mm[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> operator-(const JIArray<T, RANK> &array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] - array.mm[i];
        }
        return out;
    }

    inline JIArray<T, RANK> friend operator-(const JIArray<T, RANK> &array) {
        JIArray<T, RANK> out;
        out = array;
        for (int i = 0; i < out.nn; ++i) {
            out.mm[i] = -out.mm[i];
        }
        return out;
    }

    inline JIArray<T, RANK> operator-(const JIArray<T, RANK> &array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] - array.mm[i];
        }
        return out;
    }

    inline JIArray<T, RANK> operator*(const JIArray<T, RANK> &array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] * array.mm[i];
        }
        return out;
    }

    inline JIArray<T, RANK> &operator*=(const T &val) {
        for (int i = 0; i < nn; ++i) {
            this->mm[i] *= val;
        }
        return *this;
    }

    inline JIArray<T, RANK> &operator*=(const JIArray<T, RANK> &array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);

        for (int i = 0; i < nn; ++i) {
            this->mm[i] *= array.mm[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> operator*(const JIArray<T, RANK> &array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIArray<T, RANK> result;

        result.nn = array.nn;
        result.mm = new T[result.nn];
        result.allocated = JIARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, result.rankSize);
        std::copy(array.offset, array.offset + RANK, result.offset);
        for (int i = 0; i < nn; ++i) {
            result.mm[i] = this->mm[i] * array.mm[i];
        }
        // #ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            result.sizes[i] = array.sizes[i];
        // #endif

        return result;
    }

    inline void operator/=(const T &val) {
        auto rval = 1.0 / val;
        for (int i = 0; i < nn; ++i) {
            mm[i] *= rval;
        }
    }

    inline JIArray<T, RANK> operator/=(const JIArray<T, RANK> &array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);

        for (int i = 0; i < nn; ++i) {
            this->mm[i] /= array.mm[i];
        }
        return *this;
    }

    friend inline JIArray<T, RANK> operator/(const JIArray<T, RANK> &lhs,
                                             const JIArray<T, RANK> &rhs) {
        JIARRAY_CHECK_SIZE(lhs.nn, rhs.nn);
        JIArray<T, RANK> out;
        out.initByRankSize(lhs.getSize(), lhs.getRankSize(), lhs.getOffset());
        for (int i = 0; i < lhs.nn; ++i) {
            out.mm[i] = lhs.mm[i] / rhs.mm[i];
        }
        return out;
    }

    template <typename Scalar, typename = std::enable_if_t<
                                   std::is_arithmetic<Scalar>::value &&
                                   !std::is_same<std::decay_t<Scalar>, JIArray<T, RANK>>::value>>
    friend inline JIArray<T, RANK> operator/(const Scalar &val, const JIArray<T, RANK> &array) {
        JIArray<T, RANK> result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val / array.mm[i];
        }
        return result;
    }

    template <typename Scalar, typename = std::enable_if_t<
                                   std::is_arithmetic<Scalar>::value &&
                                   !std::is_same<std::decay_t<Scalar>, JIArray<T, RANK>>::value>>
    friend inline JIArray<T, RANK> operator/(const JIArray<T, RANK> &array, const Scalar &val) {
        JIArray<T, RANK> result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = array.mm[i] / val;
        }
        return result;
    }

    inline double sqsum() const {
        double result = 0;
        for (int i = 0; i < nn; ++i) {
            result += mm[i] * mm[i];
        }

        return result;
    }

    inline JIArray<T, RANK> copy() const {
        JIArray<T, RANK> array;
        array = *this;
        return array;
    }

    inline JIArray<T, RANK> friend operator+(const T &val, const JIArray<T, RANK> &array) {
        JIArray<T, RANK> result;

        result.nn = array.nn;
        result.mm = new T[result.nn];
        result.allocated = JIARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, result.rankSize);
        std::copy(array.offset, array.offset + RANK, result.offset);
        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val + array.mm[i];
        }
        // #ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            result.sizes[i] = array.sizes[i];
        // #endif

        return result;
    }

    template <typename Scalar, typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
    friend inline JIArray<T, RANK> operator*(const Scalar &val, const JIArray<T, RANK> &array) {
        JIArray<T, RANK> result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());

        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val * array.mm[i];
        }

        return result;
    }

    template <typename Scalar, typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
    friend inline JIArray<T, RANK> operator*(const JIArray<T, RANK> &array, const Scalar &val) {
        JIArray<T, RANK> result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());

        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val * array.mm[i];
        }

        return result;
    }

    // friend inline JIArray<T, RANK> operator/(const T &val, const JIArray<T, RANK> &array) {
    //     JIArray<T, RANK> result;

    //     result.nn = array.nn;
    //     result.mm = new T[result.nn];
    //     result.isAllocated() = JIARRAY_ALLOCATED_ALL;

    //     std::copy(array.rankSize, array.rankSize + RANK, result.rankSize);
    //     std::copy(array.offset, array.offset + RANK, result.offset);
    //     for (int i = 0; i < result.nn; ++i) {
    //         result.mm[i] = val / array.mm[i];
    //     }
    //     // #ifdef JIARRAY_DEBUG
    //     for (int i = 0; i < RANK; i++)
    //         result.sizes[i] = array.sizes[i];
    //     // #endif

    //     return result;
    // }

    friend inline T dot(const JIArray<T, RANK> &array1, const JIArray<T, RANK> &array2) {
        // #ifdef JIARRAY_DEBUG
        for (int rank = 0; rank < RANK; ++rank) {
            assert(array1.sizes[rank] == array2.sizes[rank]);
            assert(array1.offset[rank] == array2.offset[rank]);
            JIARRAY_CHECK_SIZE(array1.nn, array2.nn);
        }
        // #endif

        T result = 0.0;
        for (int i = 0; i < array1.nn; ++i) {
            result += array1.mm[i] * array2.mm[i];
        }
        return result;
    }

    inline bool contains(const T &item) const {
        for (int i = 0; i < nn; ++i) {
            if (mm[i] == item)
                return true;
        }
        return false;
    }

    inline auto findFirst(const T &item) const {
        if constexpr (RANK == 1) {
            for (int i = 0; i < this->nn; ++i) {
                if (mm[i] == item) {
                    return i + offset[0];
                }
            }
            return (-1 + JIARRAY_OFFSET);
        } else {
            int loc = -1;
            for (int i = 0; i < this->nn; ++i) {
                if (mm[i] == item) {
                    loc = i;
                    break;
                }
            }

            FastArray<int, RANK> location;
            location = (-1 + JIARRAY_OFFSET);

            if (loc == -1)
                return location;

            for (int rank = RANK - 1; rank > 0; --rank) {
                location.mm[rank] = loc / rankSize[rank];
                loc -= (location.mm[rank]) * rankSize[rank];
            }
            location.mm[0] = loc;

            for (int rank = 0; rank < RANK; ++rank) {
                location.mm[rank] += offset[rank];
            }

            return location;
        }
    }

    inline FastArray<int, RANK> maxloc() {
        int maxloc = 0;
        T maxval = mm[0];
        for (int i = 1; i < this->nn; ++i) {
            if (mm[i] > maxval) {
                maxloc = i;
                maxval = mm[i];
            }
        }

        FastArray<int, RANK> maxLocation;

        for (int rank = RANK - 1; rank > 0; --rank) {
            maxLocation.mm[rank] = maxloc / rankSize[rank];
            maxloc -= (maxLocation.mm[rank]) * rankSize[rank];
        }
        maxLocation.mm[0] = maxloc;

        for (int rank = 0; rank < RANK; ++rank) {
            maxLocation.mm[rank] += offset[rank];
        }

        return maxLocation;
    }

    inline FastArray<int, RANK> maxloc(int from, int to) {
        int maxloc = 0;
        T maxval = mm[from - JIARRAY_OFFSET];
        for (int i = from - JIARRAY_OFFSET + 1; i < to - JIARRAY_OFFSET; ++i) {
            if (mm[i] > maxval) {
                maxloc = i;
                maxval = mm[i];
            }
        }

        FastArray<int, RANK> maxLocation;

        for (int rank = RANK - 1; rank > 0; ++rank) {
            maxLocation.mm[rank] = maxloc / rankSize[rank] + JIARRAY_OFFSET;
            maxloc -= maxLocation.mm[rank] * rankSize[rank];
        }
        maxLocation.mm[0] = maxloc + JIARRAY_OFFSET;
        return maxLocation;
    }

    inline void shareWith(const JIArray<T, RANK> &array) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        nn = array.nn;
        mm = array.mm;
        std::copy((int *)array.rankSize, ((int *)array.rankSize) + RANK, rankSize);
        std::copy((int *)array.offset, ((int *)array.offset) + RANK, offset);
        allocated = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
        // #ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            this->sizes[i] = array.sizes[i];
        // #endif
    }

    // Iterator implementation (unchanged - operates on linear memory)
    struct Iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;
        using pointer = T *;
        using reference = T &;

        Iterator(pointer ptr) : m_ptr(ptr) {}

        Iterator operator+(difference_type n) const { return Iterator(m_ptr + n); }
        Iterator operator-(difference_type n) const { return Iterator(m_ptr - n); }
        Iterator &operator+=(difference_type n) {
            m_ptr += n;
            return *this;
        }
        Iterator &operator-=(difference_type n) {
            m_ptr -= n;
            return *this;
        }

        T &operator[](difference_type n) { return *(m_ptr + n); }
        const T &operator[](difference_type n) const { return *(m_ptr + n); }

        difference_type operator-(const Iterator &other) const { return m_ptr - other.m_ptr; }

        bool operator<(const Iterator &other) const { return m_ptr < other.m_ptr; }
        bool operator<=(const Iterator &other) const { return m_ptr <= other.m_ptr; }
        bool operator>(const Iterator &other) const { return m_ptr > other.m_ptr; }
        bool operator>=(const Iterator &other) const { return m_ptr >= other.m_ptr; }

        reference operator*() { return *m_ptr; }
        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }
        const pointer operator->() const { return m_ptr; }

        Iterator &operator++() {
            m_ptr++;
            return *this;
        }
        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }
        Iterator &operator--() {
            --m_ptr;
            return *this;
        }
        Iterator operator--(int) {
            Iterator temp = *this;
            --(*this);
            return temp;
        }

        friend bool operator==(const Iterator &a, const Iterator &b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!=(const Iterator &a, const Iterator &b) { return a.m_ptr != b.m_ptr; }

      private:
        pointer m_ptr;
    };

  public:
    Iterator begin() { return Iterator(mm); }
    Iterator end() { return Iterator(mm + nn); }
    const Iterator begin() const { return Iterator(mm); }
    const Iterator end() const { return Iterator(mm + nn); }

    std::vector<T> convertToVector() {
        std::vector<T> vec(mm, mm + nn);
        return vec;
    }
};

#define zbool1 JIArray<bool, 1>
#define zbool2 JIArray<bool, 2>
#define zbool3 JIArray<bool, 3>
#define zbool4 JIArray<bool, 4>
#define zbool5 JIArray<bool, 5>
#define zint1 JIArray<int, 1>
#define zint2 JIArray<int, 2>
#define zint3 JIArray<int, 3>
#define zint4 JIArray<int, 4>
#define zint5 JIArray<int, 5>
#define zdouble1 JIArray<double, 1>
#define zdouble2 JIArray<double, 2>
#define zdouble3 JIArray<double, 3>
#define zdouble4 JIArray<double, 4>
#define zdouble5 JIArray<double, 5>
#define zdouble6 JIArray<double, 6>

#define zfloat1 JIArray<float, 1>
#define zfloat2 JIArray<float, 2>
#define zfloat3 JIArray<float, 3>
#define zfloat4 JIArray<float, 4>
#define zfloat5 JIArray<float, 5>

#define zstring1 JIArray<string, 1>
#define zstring2 JIArray<string, 2>
#define zstring3 JIArray<string, 3>
#define zstring4 JIArray<string, 4>
#define zstring5 JIArray<string, 5>

#if JIARRAY_OFFSET == 0
#define ffor(i, begin, end) for (int i = begin; i < end; ++i)
#define ffor_back(i, begin, end) for (int i = begin; i >= end; --i)
#else
#define ffor(i, begin, end) for (int i = begin; i <= end; ++i)
#define ffor_back(i, begin, end) for (int i = begin; i >= end; --i)
#endif

#define zfor(i, end) ffor(i, JIARRAY_OFFSET, end)

template <typename Type, int N = 1> using zarray = JIArray<Type, N>;

} // namespace dnegri::jiarray
