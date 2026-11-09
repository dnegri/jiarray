// Non-owning, trivially copyable N-D view of device memory with JIArray indexing (passed to kernels by value).
#pragma once

/**
 * @file JIDeviceView.h
 * @brief `JIDeviceView<T, RANK, Index>` — the kernel-side half of the device array.
 *
 * A view carries a device pointer plus the shape, offsets and strides of the
 * array it looks at.  Indexing follows JIArray exactly (`JIARRAY_OFFSET`,
 * `JIARRAY_COLUMN_MAJOR`, per-dimension offsets from `init0`/`setOffsets`,
 * `slice` peeling the slowest dimension), while the cost of `a(i, j, k)` is
 * that of a hand-written 0-based raw index:
 *
 * - the offsets are folded into a shifted base pointer (`origin`) when the
 *   view is built, so no `- offset` term survives into the access;
 * - the fastest dimension always has unit stride (every view is contiguous:
 *   allocation, `slice` and `reshape` all preserve that), so its multiply is
 *   dropped at compile time;
 * - the linear index is computed in `Index` (32-bit `int` by default).
 *
 * The whole header is empty when neither nvcc nor the CUDA runtime header is
 * available, so host-only translation units may include it unconditionally.
 */

#if defined(__CUDACC__) || __has_include(<cuda_runtime.h>)

#include "FastArray.h"
#include "JIArray.h"
#include "pch.h"

#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <type_traits>

// Bounds/size checks of the device views.  Host side throws the same
// exceptions as JIArray (`std::out_of_range`, `std::invalid_argument`);
// device side cannot throw, so it prints the offending values and traps,
// which aborts the kernel and leaves cudaErrorLaunchFailure on the context.
#if defined(JIARRAY_DEBUG)
#if defined(__CUDA_ARCH__)
#define JIARRAY_DEVICE_CHECK_BOUND(i, beg, end)                                                                        \
    do {                                                                                                               \
        if ((i) < (beg) || (i) > (end)) {                                                                              \
            printf("JIDeviceView: index %lld outside [%lld, %lld] (%s:%d)\n", static_cast<long long>(i),               \
                   static_cast<long long>(beg), static_cast<long long>(end), __FILE__, __LINE__);                      \
            __trap();                                                                                                  \
        }                                                                                                              \
    } while (0)
#define JIARRAY_DEVICE_CHECK_SIZE(a, b)                                                                                \
    do {                                                                                                               \
        if ((a) != (b)) {                                                                                              \
            printf("JIDeviceView: size %lld != %lld (%s:%d)\n", static_cast<long long>(a), static_cast<long long>(b),  \
                   __FILE__, __LINE__);                                                                                \
            __trap();                                                                                                  \
        }                                                                                                              \
    } while (0)
#else
#define JIARRAY_DEVICE_CHECK_BOUND(i, beg, end)                                                                        \
    do {                                                                                                               \
        if ((i) < (beg) || (i) > (end)) {                                                                              \
            throw std::out_of_range("Index out of bounds");                                                            \
        }                                                                                                              \
    } while (0)
#define JIARRAY_DEVICE_CHECK_SIZE(a, b)                                                                                \
    do {                                                                                                               \
        if ((a) != (b)) {                                                                                              \
            throw std::invalid_argument("Size mismatch");                                                              \
        }                                                                                                              \
    } while (0)
#endif
#else
#define JIARRAY_DEVICE_CHECK_BOUND(i, beg, end) ((void)0)
#define JIARRAY_DEVICE_CHECK_SIZE(a, b) ((void)0)
#endif

namespace dnegri::jiarray {

template <class T, size_t RANK = 1, class Index = int>
class JIDeviceView;

template <class T, size_t RANK, class Index>
class JIDeviceArray;

/**
 * @class JIDeviceView
 * @brief Non-owning view of a contiguous device array, indexed like JIArray.
 *
 * Trivially copyable: pass it to a kernel by value.  Like `std::span`, the
 * constness of the view object does not propagate to the elements — a
 * `const` view still writes; use `JIDeviceView<const T, ...>` (the owner's
 * `view() const`, alias `ConstView`) for read-only access.  A
 * `JIDeviceView<T>` converts implicitly to `JIDeviceView<const T>`.
 *
 * The view never frees memory and must not outlive its owner
 * (`JIDeviceArray`) or the external buffer it wraps.
 *
 * @code
 * __global__ void scale(ddouble3::View a, double factor) {
 *     const int i = blockIdx.x * blockDim.x + threadIdx.x + JIARRAY_OFFSET;
 *     if (i > a.getSize(1)) return;
 *     zfor (k, a.getSize(3)) {
 *         zfor (j, a.getSize(2)) {
 *             a(i, j, k) *= factor;
 *         }
 *     }
 * }
 * @endcode
 *
 * @tparam T     Element type (may be const-qualified).
 * @tparam RANK  Number of dimensions (≥ 1).
 * @tparam Index Signed integer type of sizes, strides and the linear index:
 *               `int` (default, 32-bit address arithmetic) or
 *               `long long`/`std::int64_t` for arrays beyond 2^31 elements.
 */
template <class T, size_t RANK, class Index>
class JIDeviceView {
    static_assert(RANK > 0, "Array rank must be at least 1");
    static_assert(std::is_integral_v<Index> && std::is_signed_v<Index>, "Index must be a signed integer type");

    template <class, size_t, class>
    friend class JIDeviceView;
    template <class, size_t, class>
    friend class JIDeviceArray;

    static constexpr bool is_row_major = (JIARRAY_COLUMN_MAJOR == 0);

    /// Dimension with unit stride: the first (column-major) or the last (row-major).
    static constexpr size_t kFastDim = is_row_major ? RANK - 1 : 0;

public:
    using value_type = T;
    using index_type = Index;
    using size_type = Index;

    static constexpr size_t rank = RANK;
    static constexpr bool row_major = is_row_major;

    /// @brief Empty view (null pointer, zero size).
    JIDeviceView() = default;

    /**
     * @brief Wrap an external contiguous device buffer with default offsets (`JIARRAY_OFFSET`).
     * @param device_data Device pointer to the first element; not owned.
     * @param dims        Size of each dimension.
     */
    template <typename... Dims, typename = std::enable_if_t<all_integral_v<Dims...> && (sizeof...(Dims) == RANK)>>
    JIARRAY_HD JIDeviceView(T* device_data, Dims... dims) {
        const Index dimensions[RANK] = {static_cast<Index>(dims)...};
        initLayout(device_data, dimensions);
    }

    /// @brief Read-only view of a mutable one (`View` → `ConstView`).
    template <class U, typename = std::enable_if_t<std::is_same_v<const U, T> && !std::is_same_v<U, T>>>
    JIARRAY_HD JIDeviceView(const JIDeviceView<U, RANK, Index>& other)
        : mm(other.mm), origin(other.origin), nn(other.nn) {
        for (size_t d = 0; d < RANK; ++d) {
            sizes[d] = other.sizes[d];
            offset[d] = other.offset[d];
            rankSize[d] = other.rankSize[d];
        }
    }

    // ------------------------------------------------------------------------
    // Element access
    // ------------------------------------------------------------------------

    /**
     * @brief Element at offset-relative indices (one per dimension), as `JIArray::operator()`.
     * @note Costs the same instructions as a 0-based raw index: offsets are
     *       pre-folded into `origin` and the fastest dimension has unit stride.
     */
    template <typename... Args, typename = std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK)>>
    JIARRAY_HD inline T& operator()(Args... index) const {
        const Index idx[RANK] = {static_cast<Index>(index)...};
        return origin[computeLinearIndex(idx)];
    }

    /// @brief Element at a FastArray of indices, as `JIArray::operator()(const FastArray<int, RANK>&)`.
    JIARRAY_HD inline T& operator()(const FastArray<int, RANK>& index) const {
        Index idx[RANK];
        for (size_t d = 0; d < RANK; ++d) {
            idx[d] = static_cast<Index>(index[d]);
        }
        return origin[computeLinearIndex(idx)];
    }

    /// @brief Same as operator().
    template <typename... Args, typename = std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK)>>
    JIARRAY_HD inline T& at(Args... index) const {
        return (*this)(index...);
    }

    /// @brief 0-based flat access over contiguous storage, as `JIArray::operator[]`.
    JIARRAY_HD inline T& operator[](Index flat_index) const {
        return mm[flat_index];
    }

    /// @brief Device pointer to the first stored element.
    JIARRAY_HD inline T* data() const noexcept {
        return mm;
    }

    /**
     * @brief Device pointer at partial indices, as `JIArray::data(INDEX...)`.
     *
     * With RANK indices it points at that element; with fewer it points at
     * the start of `slice(index...)` (column-major fixes the trailing
     * dimensions, row-major the leading ones).
     */
    template <typename... Args>
    JIARRAY_HD inline std::enable_if_t<all_integral_v<Args...>, T*> data(Args... index) const {
        constexpr size_t kNumIndex = sizeof...(Args);
        static_assert(kNumIndex >= 1 && kNumIndex <= RANK, "Number of indices must be between 1 and RANK");
        if constexpr (kNumIndex == RANK) {
            return &(*this)(index...);
        } else {
            return slice(index...).data();
        }
    }

    // ------------------------------------------------------------------------
    // Size and offset queries
    // ------------------------------------------------------------------------

    /// @brief Total number of elements.
    JIARRAY_HD inline Index size() const noexcept {
        return nn;
    }

    /// @overload
    JIARRAY_HD inline Index getSize() const noexcept {
        return nn;
    }

    /// @brief Size of dimension `dim` (numbered from JIARRAY_OFFSET), as `JIArray::getSize(int)`.
    JIARRAY_HD inline Index getSize(int dim) const {
        JIARRAY_DEVICE_CHECK_BOUND(dim, JIARRAY_OFFSET, JIARRAY_OFFSET + static_cast<int>(RANK) - 1);
        return sizes[dim - JIARRAY_OFFSET];
    }

    /// @brief Sizes of all dimensions (RANK entries).
    JIARRAY_HD inline const Index* getSizeOfRank() const noexcept {
        return sizes;
    }

    /// @brief Strides of all dimensions (RANK entries), as `JIArray::getRankSize()`.
    JIARRAY_HD inline const Index* getRankSize() const noexcept {
        return rankSize;
    }

    /// @brief First index of every dimension (RANK entries).
    JIARRAY_HD inline const Index* getOffset() const noexcept {
        return offset;
    }

    /// @brief First index of dimension `dim` (numbered from JIARRAY_OFFSET).
    JIARRAY_HD inline Index getOffset(int dim) const {
        JIARRAY_DEVICE_CHECK_BOUND(dim, JIARRAY_OFFSET, JIARRAY_OFFSET + static_cast<int>(RANK) - 1);
        return offset[dim - JIARRAY_OFFSET];
    }

    /// @brief Whether the view points at memory.
    JIARRAY_HD inline bool isAllocated() const noexcept {
        return mm != nullptr;
    }

    // ------------------------------------------------------------------------
    // Layout changes (views only — no data moves)
    // ------------------------------------------------------------------------

    /// @brief Set the first index of every dimension, as `JIArray::setOffsets`.
    template <typename... Args>
    JIARRAY_HD inline std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK), void>
    setOffsets(Args... offsets) noexcept {
        const Index new_offsets[RANK] = {static_cast<Index>(offsets)...};
        setOffsets(new_offsets);
    }

    /// @overload
    JIARRAY_HD inline void setOffsets(const Index new_offsets[RANK]) noexcept {
        for (size_t d = 0; d < RANK; ++d) {
            offset[d] = new_offsets[d];
        }
        foldOrigin();
    }

    /**
     * @brief Lower-rank view with the given indices fixed, as `JIArray::slice`.
     *
     * Column-major fixes the trailing dimensions (`a.slice(k)` is `a(:, :, k)`),
     * row-major the leading ones (`a.slice(i)` is `a(i, :, :)`).  The kept
     * dimensions keep their sizes and offsets.  Callable in kernels.
     */
    template <typename... Args>
    JIARRAY_HD inline JIDeviceView<T, RANK - sizeof...(Args), Index> slice(Args... index) const {
        constexpr size_t kNumFixed = sizeof...(Args);
        static_assert(all_integral_v<Args...>, "Slice indices must be integral");
        static_assert(kNumFixed >= 1 && kNumFixed < RANK, "Number of slice indices must be in [1, RANK)");
        constexpr size_t kRank2 = RANK - kNumFixed;
        constexpr size_t kFirstFix = is_row_major ? 0 : kRank2;
        constexpr size_t kFirstKeep = is_row_major ? kNumFixed : 0;

        const Index fixed[kNumFixed] = {static_cast<Index>(index)...};
        Index shift = 0;
        for (size_t f = 0; f < kNumFixed; ++f) {
            const size_t d = kFirstFix + f;
            JIARRAY_DEVICE_CHECK_BOUND(fixed[f], offset[d], offset[d] + sizes[d] - 1);
            shift += fixed[f] * rankSize[d];
        }

        JIDeviceView<T, kRank2, Index> out;
        out.origin = origin + shift;
        out.nn = 1;
        for (size_t r = 0; r < kRank2; ++r) {
            out.sizes[r] = sizes[kFirstKeep + r];
            out.offset[r] = offset[kFirstKeep + r];
            out.rankSize[r] = rankSize[kFirstKeep + r];
            out.nn *= out.sizes[r];
        }
        out.mm = out.origin + out.computeOriginShift();
        return out;
    }

    /**
     * @brief View of the same elements with another shape and default offsets, as `JIArray::reshape`.
     * @note The total size must not change (checked under JIARRAY_DEBUG).
     */
    template <typename... Dims>
    JIARRAY_HD inline std::enable_if_t<all_integral_v<Dims...>, JIDeviceView<T, sizeof...(Dims), Index>>
    reshape(Dims... dims) const {
        const Index dimensions[sizeof...(Dims)] = {static_cast<Index>(dims)...};
        JIDeviceView<T, sizeof...(Dims), Index> out;
        out.initLayout(mm, dimensions);
        JIARRAY_DEVICE_CHECK_SIZE(nn, out.nn);
        return out;
    }

private:
    /// @brief Set pointer, sizes, default offsets and contiguous strides.
    JIARRAY_HD inline void initLayout(T* device_data, const Index dimensions[RANK]) {
        mm = device_data;
        nn = 1;
        for (size_t d = 0; d < RANK; ++d) {
            sizes[d] = dimensions[d];
            offset[d] = JIARRAY_OFFSET;
            nn *= dimensions[d];
        }
        if constexpr (is_row_major) {
            rankSize[RANK - 1] = 1;
            for (size_t d = RANK - 1; d > 0; --d) {
                rankSize[d - 1] = rankSize[d] * sizes[d];
            }
        } else {
            rankSize[0] = 1;
            for (size_t d = 1; d < RANK; ++d) {
                rankSize[d] = rankSize[d - 1] * sizes[d - 1];
            }
        }
        foldOrigin();
    }

    /// @brief Linear position of the first element relative to `origin` (sum of offset × stride).
    JIARRAY_HD inline Index computeOriginShift() const {
        Index shift = 0;
        for (size_t d = 0; d < RANK; ++d) {
            shift += offset[d] * rankSize[d];
        }
        return shift;
    }

    // The shifted pointer may lie outside the allocation (it is never
    // dereferenced there); this is the classic Fortran-style base shift and
    // is what makes an offset-relative access a single multiply-add chain.
    JIARRAY_HD inline void foldOrigin() {
        origin = mm - computeOriginShift();
    }

    /// @brief Offset-relative indices → position relative to `origin`; fully unrolled for constant RANK.
    JIARRAY_HD inline Index computeLinearIndex(const Index idx[RANK]) const {
        Index pos = idx[kFastDim];
        JIARRAY_DEVICE_CHECK_BOUND(idx[kFastDim], offset[kFastDim], offset[kFastDim] + sizes[kFastDim] - 1);
        for (size_t d = 0; d < RANK; ++d) {
            if (d != kFastDim) {
                JIARRAY_DEVICE_CHECK_BOUND(idx[d], offset[d], offset[d] + sizes[d] - 1);
                pos += idx[d] * rankSize[d];
            }
        }
        return pos;
    }

    T* mm = nullptr;         ///< Device pointer to the first stored element.
    T* origin = nullptr;     ///< mm shifted back by Σ offset·stride: element (i, j, …) is origin[Σ index·stride].
    Index nn = 0;            ///< Total number of elements.
    Index sizes[RANK]{};     ///< Size of each dimension.
    Index offset[RANK]{};    ///< First index of each dimension.
    Index rankSize[RANK]{};  ///< Stride of each dimension; rankSize[kFastDim] == 1.
};

}  // namespace dnegri::jiarray

#endif  // __CUDACC__ || cuda_runtime.h available
