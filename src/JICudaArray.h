#pragma once

#include "FastArray.h"
#include "JIArray.h"
#include "pch.h"
#include <numeric>
#include <vector>

namespace dnegri::jiarray {

template <class T, std::size_t RANK = 1, class = std::make_index_sequence<RANK>>
class JICudaArray;

template <class T, std::size_t RANK, std::size_t... INTS>
class JICudaArray<T, RANK, std::index_sequence<INTS...>> {
private:
    int nn        = 0;
    T*  mm        = nullptr;
    int allocated = JIARRAY_ALLOCATED_NONE;
    int rankSize[RANK]{};
    int offset[RANK]{};
    int sumOfOffset = 0;

    // private:
    // #ifdef JIARRAY_DEBUG
    int sizes[RANK]{};
    // #endif

public:
    __host__ __device__ JICudaArray() {
    }

    __host__ __device__ JICudaArray(std::initializer_list<T> il)
        : JICudaArray(il.size()) {
        std::copy(il.begin(), il.end(), mm);
    }

    __host__ __device__ JICudaArray(const JICudaArray<T, RANK>& array) {
        nn = array.nn;
        mm = array.mm;
        std::copy((int*)array.rankSize, ((int*)array.rankSize) + RANK, rankSize);
        std::copy((int*)array.offset, ((int*)array.offset) + RANK, offset);
        allocated   = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
        // #ifdef JIARRAY_DEBUG
        for (auto i = 0; i < RANK; i++)
            this->sizes[i] = array.sizes[i];
        // #endif
    }

    __host__ __device__ virtual ~JICudaArray() {
        destroy();
    }

    template <typename... Args,
              typename = std::enable_if_t<(std::is_integral_v<Args> && ...)>>
    JICudaArray(Args... args) {
        init(args...);
    }

    JICudaArray(T* memory_, decltype(INTS)... args) {
        init(args..., memory_);
    }

    template <typename... INTS2, typename = std::enable_if_t<(sizeof...(INTS2) == 2 * RANK)>>
    void init0(INTS2... array_sizes) {
        int temp_size[] = {static_cast<int>(array_sizes)...};

        int sizes[RANK];
        int offsets[RANK];

        for (int i = 0; i < RANK; i++) {
            sizes[i]   = temp_size[2 * i + 1] - temp_size[2 * i] + JIARRAY_OFFSET;
            offsets[i] = temp_size[2 * i];
        }

        init(sizes);
        setOffsets(offsets);
    }

    template <typename... Args, typename = std::enable_if_t<(std::is_integral_v<Args> && ...)>, typename = std::enable_if_t<(sizeof...(Args) == RANK)>>
    void init(Args... array_sizes) {
        int sizes[] = {static_cast<int>(array_sizes)...};
        init(sizes);
    }

    void init(int sizes[RANK]) {
        destroy();

        for (int i = 0; i < RANK; i++) {
            if (sizes[i] == 0) return;
        }

        int ioffset = 0;

        rankSize[0] = 1;
        offset[0]   = JIARRAY_OFFSET;
        nn          = sizes[0];
        sumOfOffset = rankSize[0] * offset[0];

        for (int i = 1; i < RANK; i++) {
            rankSize[i] = rankSize[i - 1] * sizes[i - 1];
            offset[i]   = JIARRAY_OFFSET;
            nn *= sizes[i];
            sumOfOffset += offset[i] * rankSize[i];
        }

        // mm = new T[nn]{};
        cudaMallocManaged(&mm, nn * sizeof(T));

        allocated = JIARRAY_ALLOCATED_ALL;
        // #ifdef JIARRAY_DEBUG
        std::copy(sizes, sizes + RANK, this->sizes);
        // #endif
    }

    __host__ __device__ void initByRankSize(int size, const int* rankSizes, const int* offsets, T* memory) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        nn = size;

        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++) {
            rankSize[i] = rankSizes[i];
            offset[i]   = offsets[i];
            sumOfOffset += offset[i] * rankSize[i];
        }

        mm = memory;

        allocated = JIARRAY_ALLOCATED_NONE;
        // #ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK - 1; i++) {
            this->sizes[i] = 0;
            if (rankSize[i] != 0) this->sizes[i] = rankSize[i + 1] / rankSize[i];
        }
        this->sizes[RANK - 1] = 0;
        if (rankSize[RANK - 1] != 0) this->sizes[RANK - 1] = size / rankSize[RANK - 1];
        // #endif
    }

    __host__ __device__ void initByRankSize(int size, const int* rankSizes, const int* offsets) {
        destroy();
        if (size == 0) return;

        mm = new T[size];
        initByRankSize(size, rankSizes, offsets, mm);
        allocated = JIARRAY_ALLOCATED_ALL;
    }

    __host__ __device__ void init(decltype(INTS)... array_sizes, T* memory_) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        nn = (array_sizes * ...);
        mm = memory_;

        int ioffset = 0;

        int sizes[] = {static_cast<int>(array_sizes)...};

        rankSize[0] = 1;
        offset[0]   = JIARRAY_OFFSET;

        sumOfOffset = rankSize[0] * offset[0];
        for (int i = 1; i < RANK; i++) {
            offset[i]   = JIARRAY_OFFSET;
            rankSize[i] = rankSize[i - 1] * sizes[i - 1];
            sumOfOffset += offset[i] * rankSize[i];
        }

        allocated = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
        // #ifdef JIARRAY_DEBUG
        std::copy(sizes, sizes + RANK, this->sizes);
        // #endif
    }

    __host__ __device__ void destroy() {
        if (allocated != JIARRAY_ALLOCATED_NONE) {
            if ((allocated & JIARRAY_ALLOCATED_MEMORY) != 0 && mm != nullptr) {
                cudaFree(mm);
            }
            mm        = nullptr;
            allocated = JIARRAY_ALLOCATED_NONE;
            nn        = 0;
        }
    }

    void erase() {
        destroy();
        std::fill(rankSize, rankSize + RANK, 0);
        std::fill(offset, offset + RANK, 0);
        sumOfOffset = 0;
        // #ifdef JIARRAY_DEBUG
        std::fill(this->sizes, this->sizes + RANK, 0);
        // #endif
    }

    template <typename... Args,
              typename = std::enable_if_t<(std::is_integral_v<Args> && ...)>>
    inline void setOffsets(Args... offsets) {
        static_assert(sizeof...(offsets) > 0, "Offsets cannot be empty");

        int tempOffsets[] = {static_cast<int>(offsets)...}; // Expand to an array
        setOffsets(tempOffsets);                            // Call overloaded function
    }

    inline void setOffsets(int offsets[RANK]) {
        std::copy(offsets, offsets + RANK, offset);

        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++)
            sumOfOffset += rankSize[i] * offset[i];
    }

    void setSize(decltype(INTS)... array_sizes) {
        size_t sizes[] = {static_cast<size_t>(array_sizes)...};

        bool same = true;

        for (int i = 0; i < RANK; i++) {
            if (this->sizes[i] != sizes[i]) {
                same = false;
                break;
            }
        }

        if (same) return;

        auto new_array = JICudaArray<T, RANK>();
        new_array.init(array_sizes...);

        // complete copy this to new_array with copy_size
        for (int new_pos = 0; new_pos < new_array.nn; new_pos++) {
            int idx[RANK];
            int temp = new_pos;
            for (int j = RANK - 1; j >= 0; j--) {
                idx[j] = temp % sizes[j];
                temp /= sizes[j];
            }

            int old_pos = 0;
            for (int j = 0; j < RANK; j++) {
                old_pos += idx[j] * rankSize[j];
            }

            new_array.mm[new_pos] = mm[old_pos];
        }

        destroy();
        initByRankSize(new_array.getSize(), new_array.getRankSize(), new_array.getOffset(), new_array.getMemory());
        new_array.allocated = JIARRAY_ALLOCATED_NONE;
        allocated           = JIARRAY_ALLOCATED_ALL;
    }

    template <typename... INDEX>
    __host__ __device__ inline JICudaArray<T, RANK - sizeof...(INDEX)> slice(INDEX... index) const {
        constexpr int RANK2   = RANK - sizeof...(INDEX);
        constexpr int num_idx = sizeof...(INDEX);

        int idx[] = {static_cast<int>(index)...};

        int p_mm = 0;
        int rank = RANK;
        for (int i = num_idx - 1; i >= 0; i--) {
            rank--;
            JIARRAY_CHECK_BOUND(idx[i], offset[rank], offset[rank] + this->sizes[rank] - 1);
            p_mm += rankSize[rank] * (idx[i] - offset[rank]);
        }

        auto array = JICudaArray<T, RANK2>();
        array.initByRankSize(rankSize[RANK2], rankSize, offset, mm + p_mm);
        return array;
    }

    inline bool isAllocated() const {
        return allocated != JIARRAY_ALLOCATED_NONE;
    }

    __host__ __device__ inline T* getMemory(int idx = 0) {
        return mm + idx;
    }

    __host__ __device__ inline const T* getMemory(int idx = 0) const {
        return mm + idx;
    }

    inline T* data() {
        return mm;
    }

    __host__ __device__ inline const T* data() const {
        return mm;
    }

    __host__ __device__ inline const int& size() const {
        return nn;
    }

    inline T average() const {
        T result = 0.0;
        for (int i = 0; i < nn; ++i) {
            result += mm[i];
        }
        return result / nn;
    }

    inline T sum() const {
        T result = std::accumulate(mm, mm + nn, 0.0);
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
        T mx = mm[0];
        for (int i = 1; i < nn; ++i) {
            if (mm[i] < mx)
                mx = mm[i];
        }
        return mx;
    }

    inline int getSize() const {
        return nn;
    }

    __host__ __device__ inline const int* getRankSize() const {
        return rankSize;
    }

    __host__ __device__ inline const int* getSizeOfRank() const {
        return sizes;
    }

    // #ifdef JIARRAY_DEBUG
    __host__ __device__ inline const int& getSize(int rank) const {
        return this->sizes[rank - JIARRAY_OFFSET];
    }
    // #endif

    __host__ __device__ inline const int* getOffset() const {
        return offset;
    }

    __host__ __device__ inline const int& getOffset(int rank) const {
        return offset[rank - JIARRAY_OFFSET];
    }

    __host__ __device__ inline T& at(decltype(INTS)... index) {
        JIARRAY_CHECK_RANK(RANK, sizeof...(index));

        size_t idx[] = {index...};

        int pos = -sumOfOffset;
        for (int i = 0; i < sizeof...(index); i++) {
            JIARRAY_CHECK_BOUND(idx[i], offset[i], offset[i] + this->sizes[i] - 1);
            pos += rankSize[i] * idx[i];
        }

        return mm[pos];
    }

    __host__ __device__ inline const T& at(decltype(INTS)... index) const {
        return const_cast<std::remove_const_t<JICudaArray<T, RANK>>&>(*this).at(index...);
    }

    __host__ __device__ inline const T& operator[](size_t index) const {
        return at(index);
    }

    __host__ __device__ inline T& operator[](size_t index) {
        return at(index);
    }

    __host__ __device__ inline const T& operator()(decltype(INTS)... index) const {
        return at(index...);
    }

    __host__ __device__ inline T& operator()(decltype(INTS)... index) {
        return at(index...);
    }

    __host__ __device__ inline T& operator()(const FastArray<int, RANK>& idx) {
        int pos = idx(1) - sumOfOffset;
        JIARRAY_CHECK_BOUND(idx(1), offset[0], offset[0] + this->sizes[0] - 1);
        for (int i = 1; i < RANK; i++) {
            JIARRAY_CHECK_BOUND(idx(i + 1), offset[i], offset[i] + this->sizes[i] - 1);
            pos += rankSize[i] * idx(i + 1);
        }

        return mm[pos];
    }

    template <typename... INDEX>
    inline T* data(INDEX... index) {
        int num_idx = sizeof...(index);
        JIARRAY_CHECK_BOUND(num_idx, 1, RANK);

        int idx[] = {static_cast<int>(index)...};

        int rank = RANK;
        int pos  = 0;
        for (int i = num_idx - 1; i >= 0; i--) {
            --rank;
            pos += rankSize[rank] * (idx[i] - offset[rank]);
        }

        return mm + pos;
    }

    template <typename... INDEX>
    __host__ __device__ inline const T* data(INDEX... index) const {
        return const_cast<std::remove_const_t<JICudaArray<T, RANK>>&>(*this).data(index...);
    }

    template <typename... INTS2>
    inline JICudaArray<T, sizeof...(INTS2)> reshape(INTS2... sizes) {
        constexpr int RANK2 = sizeof...(INTS2);

        JICudaArray<T, RANK2> array;
        array.init(sizes..., getMemory());

        return array;
    }

    inline JICudaArray<T, RANK>& operator=(const std::initializer_list<T>& list) {
        JIARRAY_CHECK_SIZE(this->nn, list.size());
        int idx = 0;
        for (const T& val : list) {
            mm[idx++] = val;
        }
        return *this;
    }

    inline JICudaArray<T, RANK>& operator=(const T& val) {
        assert(nn > 0);
        for (int i = 0; i < nn; ++i) {
            mm[i] = val;
        }

        return *this;
    }

    template <typename T2>
    inline JICudaArray<T, RANK>& operator=(const std::vector<T2>& val) {
        if (nn == 0) {
            int size[RANK];
            std::fill(size, size + RANK, 1);
            size[RANK - 1] = val.size();
            init(size);
        } else {
            JIARRAY_CHECK_SIZE(nn, val.size());
        }

        for (int i = 0; i < nn; ++i) {
            mm[i] = val[i];
        }

        return *this;
    }

    inline JICudaArray<T, RANK>& operator=(const T* array) {
        assert(nn > 0);

        for (int i = 0; i < nn; ++i) {
            mm[i] = array[i];
        }

        return *this;
    }

    inline JICudaArray<T, RANK>& operator=(const JICudaArray<T, RANK>& array) {
        if (allocated == JIARRAY_ALLOCATED_NONE && mm == nullptr) {
            initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());

        } else {
            // #ifdef JIARRAY_DEBUG
            JIARRAY_CHECK_SIZE(nn, array.getSize());

            auto arraySizeOfRank = array.getSizeOfRank();
            auto arrayOffset     = array.getOffset();
            for (int rank = 0; rank < RANK; ++rank) {
                assert(this->sizes[rank] == arraySizeOfRank[rank]);
                assert(offset[rank] == arrayOffset[rank]);
            }
            // #endif
        }

        auto arrayMemory = array.data();
        std::copy(arrayMemory, arrayMemory + nn, mm);

        return *this;
    }

    template <typename T2>
    inline JICudaArray<T, RANK>& operator=(const JIArray<T2, RANK>& array) {
        if (allocated == JIARRAY_ALLOCATED_NONE && mm == nullptr) {
            initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());

        } else {
            // #ifdef JIARRAY_DEBUG
            JIARRAY_CHECK_SIZE(nn, array.getSize());

            auto arraySizeOfRank = array.getSizeOfRank();
            auto arrayOffset     = array.getOffset();
            for (int rank = 0; rank < RANK; ++rank) {
                assert(this->sizes[rank] == arraySizeOfRank[rank]);
                assert(offset[rank] == arrayOffset[rank]);
            }
            // #endif
        }

        auto arrayMemory = array.data();
        std::copy(arrayMemory, arrayMemory + nn, mm);

        return *this;
    }

    inline bool operator==(const JICudaArray<T, RANK>& array) const {
        if (nn != array.nn)
            return false;

        for (int i = 0; i < nn; ++i) {
            if (array.mm[i] != this->mm[i])
                return false;
        }

        return true;
    }

    inline JICudaArray<T, RANK>& operator+=(const T& val) {
        for (int i = 0; i < nn; ++i) {
            this->mm[i] += val;
        }
        return *this;
    }

    inline JICudaArray<T, RANK> operator+=(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);

        for (int i = 0; i < nn; ++i) {
            this->mm[i] += array.mm[i];
        }
        return *this;
    }

    inline JICudaArray<T, RANK> operator+(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JICudaArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] + array.mm[i];
        }
        return out;
    }

    inline JICudaArray<T, RANK> operator+(const JICudaArray<T, RANK>& array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JICudaArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] + array.mm[i];
        }
        return out;
    }

    inline JICudaArray<T, RANK>& operator-=(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        for (int i = 0; i < nn; ++i) {
            mm[i] -= array.mm[i];
        }
        return *this;
    }

    inline JICudaArray<T, RANK> operator-(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JICudaArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] - array.mm[i];
        }
        return out;
    }

    inline JICudaArray<T, RANK> friend operator-(const JICudaArray<T, RANK>& array) {
        JICudaArray<T, RANK> out;
        out = array;
        for (int i = 0; i < out.nn; ++i) {
            out.mm[i] = -out.mm[i];
        }
        return out;
    }

    inline JICudaArray<T, RANK> operator-(const JICudaArray<T, RANK>& array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JICudaArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] - array.mm[i];
        }
        return out;
    }

    inline JICudaArray<T, RANK> operator*(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JICudaArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] * array.mm[i];
        }
        return out;
    }

    inline JICudaArray<T, RANK>& operator*=(const T& val) {
        for (int i = 0; i < nn; ++i) {
            this->mm[i] *= val;
        }
        return *this;
    }

    inline JICudaArray<T, RANK>& operator*=(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);

        for (int i = 0; i < nn; ++i) {
            this->mm[i] *= array.mm[i];
        }
        return *this;
    }

    inline JICudaArray<T, RANK> operator*(const JICudaArray<T, RANK>& array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JICudaArray<T, RANK> result;

        result.nn        = array.nn;
        result.mm        = new T[result.nn];
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

    inline void operator/=(const T& val) {
        auto rval = 1.0 / val;
        for (int i = 0; i < nn; ++i) {
            mm[i] *= rval;
        }
    }

    inline JICudaArray<T, RANK> operator/=(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);

        for (int i = 0; i < nn; ++i) {
            this->mm[i] /= array.mm[i];
        }
        return *this;
    }

    inline JICudaArray<T, RANK> operator/(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JICudaArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] / array.mm[i];
        }
        return out;
    }

    inline JICudaArray<T, RANK> operator/(const T& val) {
        JICudaArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < nn; ++i) {
            out.mm[i] = this->mm[i] / val;
        }
        return out;
    }

    inline double sqsum() const {
        double result = 0;
        for (int i = 0; i < nn; ++i) {
            result += mm[i] * mm[i];
        }

        return result;
    }

    inline JICudaArray<T, RANK> copy() const {
        JICudaArray<T, RANK> array;
        array = *this;
        return array;
    }

    inline JICudaArray<T, RANK> friend operator+(const T& val, const JICudaArray<T, RANK>& array) {
        JICudaArray<T, RANK> result;

        result.nn        = array.nn;
        result.mm        = new T[result.nn];
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

    template <typename Scalar,
              typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
    friend inline JICudaArray<T, RANK> operator*(const Scalar& val, const JICudaArray<T, RANK>& array) {
        JICudaArray<T, RANK> result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());

        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val * array.mm[i];
        }

        return result;
    }

    template <typename Scalar,
              typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
    friend inline JICudaArray<T, RANK> operator*(const JICudaArray<T, RANK>& array, const Scalar& val) {
        JICudaArray<T, RANK> result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());

        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val * array.mm[i];
        }

        return result;
    }

    friend inline JICudaArray<T, RANK> operator/(const T& val, const JICudaArray<T, RANK>& array) {
        JICudaArray<T, RANK> result;

        result.nn            = array.nn;
        result.mm            = new T[result.nn];
        result.isAllocated() = JIARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, result.rankSize);
        std::copy(array.offset, array.offset + RANK, result.offset);
        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val / array.mm[i];
        }
        // #ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            result.sizes[i] = array.sizes[i];
        // #endif

        return result;
    }

    friend inline T dot(const JICudaArray<T, RANK>& array1, const JICudaArray<T, RANK>& array2) {
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

    inline bool contains(const T& item) const {
        for (int i = 0; i < nn; ++i) {
            if (mm[i] == item)
                return true;
        }
        return false;
    }

    inline FastArray<int, RANK> findFirst(const T& item) const {
        int loc = -1;
        T   value;
        for (int i = 0; i < this->nn; ++i) {
            if (mm[i] == item) {
                loc   = i;
                value = mm[i];
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

    inline FastArray<int, RANK> maxloc() {
        int maxloc = 0;
        T   maxval = mm[0];
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
        T   maxval = mm[from - JIARRAY_OFFSET];
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

    inline void shareWith(const JICudaArray<T, RANK>& array) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        nn = array.nn;
        mm = array.mm;
        std::copy((int*)array.rankSize, ((int*)array.rankSize) + RANK, rankSize);
        std::copy((int*)array.offset, ((int*)array.offset) + RANK, offset);
        allocated   = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
        // #ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            this->sizes[i] = array.sizes[i];
        // #endif
    }

    struct Iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*; // or also value_type*
        using reference         = T&; // or also value_type&

        Iterator(pointer ptr)
            : m_ptr(ptr) {};

        Iterator operator+(difference_type n) const {
            return Iterator(m_ptr + n);
        }

        Iterator operator-(difference_type n) const {
            return Iterator(m_ptr - n);
        }

        Iterator& operator+=(difference_type n) {
            m_ptr += n;
            return *this;
        }

        Iterator& operator-=(difference_type n) {
            m_ptr -= n;
            return *this;
        }

        T& operator[](difference_type n) {
            return *(m_ptr + n);
        }

        const T& operator[](difference_type n) const {
            return *(m_ptr + n);
        }

        difference_type operator-(const Iterator& other) const {
            return m_ptr - other.m_ptr;
        }

        bool operator<(const Iterator& other) const {
            return m_ptr < other.m_ptr;
        }
        bool operator<=(const Iterator& other) const {
            return m_ptr <= other.m_ptr;
        }
        bool operator>(const Iterator& other) const {
            return m_ptr > other.m_ptr;
        }
        bool operator>=(const Iterator& other) const {
            return m_ptr >= other.m_ptr;
        }

        reference operator*() {
            return *m_ptr;
        }

        reference operator*() const {
            return *m_ptr;
        }

        pointer operator->() {
            return m_ptr;
        }

        const pointer operator->() const {
            return m_ptr;
        }

        // Prefix increment
        Iterator& operator++() {
            m_ptr++;
            return *this;
        }

        // Postfix increment
        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        const Iterator operator++(int) const {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        Iterator& operator--() {
            --m_ptr;
            return *this;
        }

        Iterator operator--(int) {
            Iterator temp = *this;
            --(*this);
            return temp;
        }

        const Iterator operator--(int) const {
            Iterator temp = *this;
            --(*this);
            return temp;
        }

        friend bool operator==(const Iterator& a, const Iterator& b) {
            return a.m_ptr == b.m_ptr;
        }
        friend bool operator!=(const Iterator& a, const Iterator& b) {
            return a.m_ptr != b.m_ptr;
        }

    private:
        pointer m_ptr;
    };

public:
    Iterator begin() {
        return Iterator(mm);
    }
    Iterator end() {
        return Iterator(mm + nn);
    } // 200 is out of bounds

    const Iterator begin() const {
        return Iterator(mm);
    }
    const Iterator end() const {
        return Iterator(mm + nn);
    } // 200 is out of bounds

public:
    std::vector<T> convertToVector() {
        std::vector<T> vec(mm, mm + nn);
        return vec;
    }
};

#define cbool1   JICudaArray<bool, 1>
#define cbool2   JICudaArray<bool, 2>
#define cbool3   JICudaArray<bool, 3>
#define cbool4   JICudaArray<bool, 4>
#define cbool5   JICudaArray<bool, 5>
#define cint1    JICudaArray<int, 1>
#define cint2    JICudaArray<int, 2>
#define cint3    JICudaArray<int, 3>
#define cint4    JICudaArray<int, 4>
#define cint5    JICudaArray<int, 5>
#define cdouble1 JICudaArray<double, 1>
#define cdouble2 JICudaArray<double, 2>
#define cdouble3 JICudaArray<double, 3>
#define cdouble4 JICudaArray<double, 4>
#define cdouble5 JICudaArray<double, 5>
#define cdouble6 JICudaArray<double, 6>

#define cfloat1 JICudaArray<float, 1>
#define cfloat2 JICudaArray<float, 2>
#define cfloat3 JICudaArray<float, 3>
#define cfloat4 JICudaArray<float, 4>
#define cfloat5 JICudaArray<float, 5>

#define cstring1 JICudaArray<string, 1>
#define cstring2 JICudaArray<string, 2>
#define cstring3 JICudaArray<string, 3>
#define cstring4 JICudaArray<string, 4>
#define cstring5 JICudaArray<string, 5>

#if JIARRAY_OFFSET == 0
    #define ffor(i, begin, end)      for (int i = begin; i < end; ++i)
    #define ffor_back(i, begin, end) for (int i = begin; i >= end; --i)
#else
    #define ffor(i, begin, end)      for (int i = begin; i <= end; ++i)
    #define ffor_back(i, begin, end) for (int i = begin; i >= end; --i)
#endif

#define zfor(i, end) ffor(i, JIARRAY_OFFSET, end)

template <typename Type, int N = 1>
using carray = JICudaArray<Type, N>;

} // namespace dnegri::jiarray
