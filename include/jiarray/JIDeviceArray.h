// Owning device-memory array (cudaMalloc / cudaMallocAsync) with JIArray shape, offsets and upload/download.
#pragma once

/**
 * @file JIDeviceArray.h
 * @brief `JIDeviceArray<T, RANK, Index>` — host-side owner of a device array.
 *
 * The owner allocates device-only memory (`cudaMalloc`, or stream-ordered
 * `cudaMallocAsync`), is move-only (RAII, no virtual), and moves data to and
 * from host `JIArray`s.  Kernels never see the owner: they take its
 * `view()` (`JIDeviceView`, trivially copyable) by value and index it with
 * JIArray syntax.
 *
 * Unlike `JICudaArray` (managed memory, element access from the host), this
 * array is never dereferenced on the host.
 *
 * Usable from nvcc and from host-only compilers that see `cuda_runtime.h`
 * (link `CUDA::cudart`); `fillAsync` needs nvcc because it launches a kernel.
 * Without the CUDA runtime header the file is empty.
 */

#include "JIDeviceView.h"

#if defined(__CUDACC__) || __has_include(<cuda_runtime.h>)

#include "JIArray.h"
#include "pch.h"

#include <cuda_runtime.h>

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace dnegri::jiarray {

namespace detail {

/// @brief Throw std::runtime_error on a CUDA runtime failure (host-side setup/transfer paths only).
inline void checkCuda(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("JIDeviceArray: ") + what + ": " + cudaGetErrorString(status));
    }
}

#ifdef __CUDACC__
/// @brief Grid-stride fill of `count` elements.
template <class T, class Index>
__global__ void fillDeviceKernel(T* data, Index count, T value) {
    const long long stride = static_cast<long long>(blockDim.x) * gridDim.x;
    for (long long i = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x; i < count; i += stride) {
        data[i] = value;
    }
}
#endif

}  // namespace detail

/// @brief How the current device buffer of a JIDeviceArray was allocated (selects its release call).
enum class DeviceAllocation {
    None,          ///< No memory owned.
    Blocking,      ///< cudaMalloc → cudaFree.
    StreamOrdered  ///< cudaMallocAsync(stream) → cudaFreeAsync(stream).
};

/**
 * @class JIDeviceArray
 * @brief Move-only owner of a contiguous device array with JIArray shape and offsets.
 *
 * @code
 * zdouble3 flux_host(ng, nxy, nz);           // host JIArray
 * ddouble3 flux;                             // empty device array
 * flux = flux_host;                          // allocate + upload (shape and offsets copied)
 * dfloat3 flux32;
 * flux32.uploadAs(flux_host);                // FP64 → FP32 on upload
 * scale<<<grid, block>>>(flux.view(), 2.0);  // kernel indexes flux(ig, ixy, iz)
 * flux.copyTo(flux_host);                    // download
 * @endcode
 *
 * Shape and offsets follow JIArray: `init(n1, n2, ...)` starts every
 * dimension at JIARRAY_OFFSET, `init0(lo1, hi1, ...)` and `setOffsets(...)`
 * choose per-dimension first indices.  New memory is zero-filled (the
 * JIArray `init` contract).
 *
 * @tparam T     Element type (trivially copyable).
 * @tparam RANK  Number of dimensions (≥ 1).
 * @tparam Index `int` (default) or a 64-bit signed integer; see JIDeviceView.
 */
template <class T, size_t RANK = 1, class Index = int>
class JIDeviceArray {
    static_assert(std::is_trivially_copyable_v<T>, "Device array elements must be trivially copyable");

public:
    using value_type = T;
    using index_type = Index;
    using View = JIDeviceView<T, RANK, Index>;             ///< Kernel argument, mutable elements.
    using ConstView = JIDeviceView<const T, RANK, Index>;  ///< Kernel argument, read-only elements.

    static_assert(std::is_trivially_copyable_v<View>, "JIDeviceView must stay trivially copyable");

    static constexpr size_t rank = RANK;

    /// @brief Empty array (no allocation).
    JIDeviceArray() = default;

    /// @brief Allocate (cudaMalloc, zero-filled) with the given sizes.
    template <typename... Dims, typename = std::enable_if_t<all_integral_v<Dims...> && (sizeof...(Dims) == RANK)>>
    explicit JIDeviceArray(Dims... dims) {
        init(dims...);
    }

    ~JIDeviceArray() {
        destroy();
    }

    JIDeviceArray(const JIDeviceArray&) = delete;
    JIDeviceArray& operator=(const JIDeviceArray&) = delete;

    JIDeviceArray(JIDeviceArray&& other) noexcept
        : layout(other.layout), stream(other.stream), allocation(other.allocation) {
        other.forget();
    }

    JIDeviceArray& operator=(JIDeviceArray&& other) noexcept {
        if (this != &other) {
            destroy();
            layout = other.layout;
            stream = other.stream;
            allocation = other.allocation;
            other.forget();
        }
        return *this;
    }

    // ------------------------------------------------------------------------
    // Allocation
    // ------------------------------------------------------------------------

    /// @brief Allocate with cudaMalloc; every dimension starts at JIARRAY_OFFSET; zero-filled.
    template <typename... Dims>
    std::enable_if_t<all_integral_v<Dims...> && (sizeof...(Dims) == RANK), void> init(Dims... dims) {
        const long long dimensions[RANK] = {static_cast<long long>(dims)...};
        allocate(dimensions, nullptr, DeviceAllocation::Blocking);
    }

    /**
     * @brief Allocate with cudaMallocAsync on `alloc_stream`; zero-filled on that stream.
     * @note The memory is released with cudaFreeAsync on the same stream, so the stream must outlive the
     *       allocation: destroy the array before destroying the stream.
     */
    template <typename... Dims>
    std::enable_if_t<all_integral_v<Dims...> && (sizeof...(Dims) == RANK), void> initAsync(cudaStream_t alloc_stream,
                                                                                           Dims... dims) {
        const long long dimensions[RANK] = {static_cast<long long>(dims)...};
        allocate(dimensions, alloc_stream, DeviceAllocation::StreamOrdered);
    }

    /**
     * @brief Allocate from (first, end) index pairs, as `JIArray::init0`: each range is
     *        [first, end) when JIARRAY_OFFSET is 0 and [first, end] when it is 1.
     * @code
     * dint2 a;
     * a.init0(0, 2, -1, 3);  // 1-based: a(0..2, -1..3); 0-based: a(0..1, -1..2)
     * @endcode
     */
    template <typename... Bounds>
    std::enable_if_t<all_integral_v<Bounds...> && (sizeof...(Bounds) == 2 * RANK), void> init0(Bounds... bounds) {
        const long long pairs[2 * RANK] = {static_cast<long long>(bounds)...};
        long long dimensions[RANK];
        long long first[RANK];
        for (size_t d = 0; d < RANK; ++d) {
            first[d] = pairs[2 * d];
            dimensions[d] = pairs[2 * d + 1] - pairs[2 * d] + JIARRAY_OFFSET;
        }
        allocate(dimensions, nullptr, DeviceAllocation::Blocking);
        applyOffsets(first);
    }

    /**
     * @brief Release the memory (cudaFree / cudaFreeAsync) and reset to empty.
     * @note Never throws; a release failure is left on the CUDA error state.
     */
    void destroy() noexcept {
        if (allocation == DeviceAllocation::Blocking) {
            cudaFree(layout.mm);
        } else if (allocation == DeviceAllocation::StreamOrdered) {
            cudaFreeAsync(layout.mm, stream);
        }
        forget();
    }

    // ------------------------------------------------------------------------
    // Offsets
    // ------------------------------------------------------------------------

    /// @brief Set the first index of every dimension, as `JIArray::setOffsets`.
    template <typename... Args>
    std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK), void> setOffsets(Args... offsets) {
        const long long first[RANK] = {static_cast<long long>(offsets)...};
        applyOffsets(first);
    }

    // ------------------------------------------------------------------------
    // Host ↔ device transfer
    // ------------------------------------------------------------------------

    /**
     * @brief Upload a host JIArray (blocking).  An empty array is allocated to its shape and offsets.
     * @throws std::invalid_argument if this array is allocated with another shape or other offsets.
     * @note An empty host array empties this one, as JIArray copy assignment does.
     */
    JIDeviceArray& operator=(const JIArray<T, RANK>& host) {
        uploadAs(host);
        return *this;
    }

    /**
     * @brief Upload a host JIArray of another element type, converting each element to T (e.g. FP64 → FP32).
     * @throws std::invalid_argument on a shape/offset mismatch with an allocated array.
     */
    template <class U>
    void uploadAs(const JIArray<U, RANK>& host) {
        if (host.size() == 0) {
            destroy();
            return;
        }
        if (!isAllocated()) {
            allocateLike(host);
        } else {
            checkSameShape(host);
        }
        if constexpr (std::is_same_v<U, T>) {
            detail::checkCuda(cudaMemcpy(layout.mm, host.data(), getBytes(), cudaMemcpyHostToDevice), "upload");
        } else {
            std::vector<T> staging(static_cast<size_t>(layout.nn));
            for (size_t i = 0; i < staging.size(); ++i) {
                staging[i] = static_cast<T>(host.data()[i]);
            }
            detail::checkCuda(cudaMemcpy(layout.mm, staging.data(), getBytes(), cudaMemcpyHostToDevice), "upload");
        }
    }

    /**
     * @brief Download into a host JIArray (blocking), converting to its element type if it differs.
     * @param[out] host Allocated to this shape and offsets if empty, else must match them.
     * @throws std::invalid_argument on a shape/offset mismatch.
     */
    template <class U>
    void copyTo(JIArray<U, RANK>& host) const {
        if (!isAllocated()) {
            host.destroy();
            return;
        }
        if (host.size() == 0) {
            int dimensions[RANK];
            int first[RANK];
            for (size_t d = 0; d < RANK; ++d) {
                dimensions[d] = static_cast<int>(layout.sizes[d]);
                first[d] = static_cast<int>(layout.offset[d]);
            }
            host.init(dimensions);
            host.setOffsets(first);
        } else {
            checkSameShape(host);
        }
        if constexpr (std::is_same_v<U, T>) {
            detail::checkCuda(cudaMemcpy(host.data(), layout.mm, getBytes(), cudaMemcpyDeviceToHost), "download");
        } else {
            std::vector<T> staging(static_cast<size_t>(layout.nn));
            detail::checkCuda(cudaMemcpy(staging.data(), layout.mm, getBytes(), cudaMemcpyDeviceToHost), "download");
            for (size_t i = 0; i < staging.size(); ++i) {
                host.data()[i] = static_cast<U>(staging[i]);
            }
        }
    }

    /**
     * @brief Enqueue an upload on `copy_stream` (cudaMemcpyAsync); both arrays must already have the same shape.
     * @note Overlaps with host work only when `host` lives in pinned memory; `host` must stay alive and
     *       unchanged until the stream reaches the copy.
     * @throws std::invalid_argument if this array is empty or the shapes differ.
     */
    void copyFromAsync(const JIArray<T, RANK>& host, cudaStream_t copy_stream) {
        checkAllocated();
        checkSameShape(host);
        detail::checkCuda(cudaMemcpyAsync(layout.mm, host.data(), getBytes(), cudaMemcpyHostToDevice, copy_stream),
                          "copyFromAsync");
    }

    /**
     * @brief Enqueue a download on `copy_stream`; `host` must already have this shape.
     * @note `host` may be read only after the stream is synchronized.
     * @throws std::invalid_argument if this array is empty or the shapes differ.
     */
    void copyToAsync(JIArray<T, RANK>& host, cudaStream_t copy_stream) const {
        checkAllocated();
        checkSameShape(host);
        detail::checkCuda(cudaMemcpyAsync(host.data(), layout.mm, getBytes(), cudaMemcpyDeviceToHost, copy_stream),
                          "copyToAsync");
    }

#ifdef __CUDACC__
    /// @brief Enqueue `value` into every element on `fill_stream` (a grid-stride kernel; nvcc only).
    void fillAsync(const T& value, cudaStream_t fill_stream) {
        if (!isAllocated()) {
            return;
        }
        constexpr int kThreads = 256;
        constexpr long long kMaxBlocks = 4096;
        const long long blocks = std::min<long long>((layout.nn + kThreads - 1) / kThreads, kMaxBlocks);
        detail::fillDeviceKernel<<<static_cast<unsigned>(blocks), kThreads, 0, fill_stream>>>(layout.mm, layout.nn,
                                                                                              value);
        detail::checkCuda(cudaGetLastError(), "fillAsync");
    }
#endif

    // ------------------------------------------------------------------------
    // Views and queries (host side; the pointers are device pointers)
    // ------------------------------------------------------------------------

    /// @brief Kernel argument with mutable elements.
    View view() noexcept {
        return layout;
    }

    /// @brief Kernel argument with read-only elements.
    ConstView view() const noexcept {
        return ConstView(layout);
    }

    /// @brief View with the given trailing (column-major) or leading (row-major) indices fixed.
    template <typename... Args>
    JIDeviceView<T, RANK - sizeof...(Args), Index> slice(Args... index) {
        return layout.slice(index...);
    }

    /// @overload
    template <typename... Args>
    JIDeviceView<const T, RANK - sizeof...(Args), Index> slice(Args... index) const {
        return view().slice(index...);
    }

    /// @brief View of the same elements with another shape and default offsets.
    template <typename... Dims>
    JIDeviceView<T, sizeof...(Dims), Index> reshape(Dims... dims) {
        return layout.reshape(dims...);
    }

    /// @brief Device pointer to the first element.
    T* data() noexcept {
        return layout.mm;
    }

    /// @overload
    const T* data() const noexcept {
        return layout.mm;
    }

    /// @brief Device pointer at partial indices, as `JIArray::data(INDEX...)`.
    template <typename... Args>
    std::enable_if_t<all_integral_v<Args...>, T*> data(Args... index) {
        return layout.data(index...);
    }

    Index size() const noexcept {
        return layout.nn;
    }
    Index getSize() const noexcept {
        return layout.nn;
    }
    Index getSize(int dim) const {
        return layout.getSize(dim);
    }
    const Index* getSizeOfRank() const noexcept {
        return layout.sizes;
    }
    const Index* getRankSize() const noexcept {
        return layout.rankSize;
    }
    const Index* getOffset() const noexcept {
        return layout.offset;
    }
    Index getOffset(int dim) const {
        return layout.getOffset(dim);
    }
    bool isAllocated() const noexcept {
        return allocation != DeviceAllocation::None;
    }
    /// @brief How the current buffer was allocated.
    DeviceAllocation getAllocation() const noexcept {
        return allocation;
    }

private:
    /// @brief Bytes of the whole buffer.
    size_t getBytes() const noexcept {
        return static_cast<size_t>(layout.nn) * sizeof(T);
    }

    /// @brief Reset to empty without releasing (after destroy or a move).
    void forget() noexcept {
        layout = View();
        stream = nullptr;
        allocation = DeviceAllocation::None;
    }

    /// @brief Allocate a zero-filled buffer of the given sizes with default offsets.
    void allocate(const long long dimensions[RANK], cudaStream_t alloc_stream, DeviceAllocation kind) {
        destroy();
        long long count = 1;
        for (size_t d = 0; d < RANK; ++d) {
            if (dimensions[d] < 0) {
                throw std::invalid_argument("JIDeviceArray: negative dimension");
            }
            count *= dimensions[d];
        }
        if (count == 0) {
            return;
        }
        checkFitsIndex(count);

        T* device_data = nullptr;
        if (kind == DeviceAllocation::StreamOrdered) {
            detail::checkCuda(cudaMallocAsync(reinterpret_cast<void**>(&device_data), count * sizeof(T), alloc_stream),
                              "cudaMallocAsync");
            detail::checkCuda(cudaMemsetAsync(device_data, 0, count * sizeof(T), alloc_stream), "cudaMemsetAsync");
        } else {
            detail::checkCuda(cudaMalloc(reinterpret_cast<void**>(&device_data), count * sizeof(T)), "cudaMalloc");
            detail::checkCuda(cudaMemset(device_data, 0, count * sizeof(T)), "cudaMemset");
        }

        Index sizes[RANK];
        for (size_t d = 0; d < RANK; ++d) {
            sizes[d] = static_cast<Index>(dimensions[d]);
        }
        layout.initLayout(device_data, sizes);
        stream = alloc_stream;
        allocation = kind;
    }

    /// @brief Allocate to the shape and offsets of a host JIArray.
    template <class U>
    void allocateLike(const JIArray<U, RANK>& host) {
        long long dimensions[RANK];
        long long first[RANK];
        for (size_t d = 0; d < RANK; ++d) {
            dimensions[d] = host.getSizeOfRank()[d];
            first[d] = host.getOffset()[d];
        }
        allocate(dimensions, nullptr, DeviceAllocation::Blocking);
        applyOffsets(first);
    }

    /// @brief Set offsets after checking that the folded index range fits in Index.
    void applyOffsets(const long long first[RANK]) {
        long long lowest = 0;
        long long highest = 0;
        for (size_t d = 0; d < RANK; ++d) {
            lowest += first[d] * layout.rankSize[d];
            highest += (first[d] + layout.sizes[d] - 1) * layout.rankSize[d];
        }
        checkFitsIndex(lowest);
        checkFitsIndex(highest);
        Index offsets[RANK];
        for (size_t d = 0; d < RANK; ++d) {
            offsets[d] = static_cast<Index>(first[d]);
        }
        layout.setOffsets(offsets);
    }

    /// @brief Reject sizes and folded indices that overflow Index (choose a 64-bit Index for those).
    static void checkFitsIndex(long long value) {
        if (value > static_cast<long long>(std::numeric_limits<Index>::max()) ||
            value < static_cast<long long>(std::numeric_limits<Index>::min())) {
            throw std::length_error("JIDeviceArray: index range exceeds the Index type; use a 64-bit Index");
        }
    }

    void checkAllocated() const {
        if (!isAllocated()) {
            throw std::invalid_argument("JIDeviceArray: array is not allocated");
        }
    }

    /// @brief Same sizes and offsets as a host JIArray.
    template <class U>
    void checkSameShape(const JIArray<U, RANK>& host) const {
        for (size_t d = 0; d < RANK; ++d) {
            if (host.getSizeOfRank()[d] != layout.sizes[d]) {
                throw std::invalid_argument("JIDeviceArray: size mismatch with host array");
            }
            if (host.getOffset()[d] != layout.offset[d]) {
                throw std::invalid_argument("JIDeviceArray: offset mismatch with host array");
            }
        }
    }

    View layout;                                           ///< Pointer, sizes, offsets and strides.
    cudaStream_t stream = nullptr;                         ///< Stream of a stream-ordered allocation.
    DeviceAllocation allocation = DeviceAllocation::None;  ///< Selects the release call.
};

// ============================================================================
// Type aliases: `d` (device) + element type + rank.  `c*` names are taken by
// the JICudaArray macros (managed memory) and `f*` by FastArray.
// ============================================================================

/// @name Integer device array aliases
/// @{
using dint1 = JIDeviceArray<int, 1>;
using dint2 = JIDeviceArray<int, 2>;
using dint3 = JIDeviceArray<int, 3>;
using dint4 = JIDeviceArray<int, 4>;
using dint5 = JIDeviceArray<int, 5>;
/// @}

/// @name Single-precision device array aliases
/// @{
using dfloat1 = JIDeviceArray<float, 1>;
using dfloat2 = JIDeviceArray<float, 2>;
using dfloat3 = JIDeviceArray<float, 3>;
using dfloat4 = JIDeviceArray<float, 4>;
using dfloat5 = JIDeviceArray<float, 5>;
/// @}

/// @name Double-precision device array aliases
/// @{
using ddouble1 = JIDeviceArray<double, 1>;
using ddouble2 = JIDeviceArray<double, 2>;
using ddouble3 = JIDeviceArray<double, 3>;
using ddouble4 = JIDeviceArray<double, 4>;
using ddouble5 = JIDeviceArray<double, 5>;
using ddouble6 = JIDeviceArray<double, 6>;
/// @}

/// @brief Generic device array alias, as `zarray`.
template <typename Type, size_t N = 1, class Index = int>
using darray = JIDeviceArray<Type, N, Index>;

}  // namespace dnegri::jiarray

#endif  // __CUDACC__ || cuda_runtime.h available
