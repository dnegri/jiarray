#pragma once

#include "pch.h"

#include <array>
#include <cstddef>
#include <initializer_list>
#include <type_traits>

namespace dnegri::jiarray {

/**
 * @brief Compile-time fixed-size N-dimensional array.
 *
 * A single variadic class template that subsumes the former
 * `FastArray<T,N>`, `FastArray2D<T,I,J>`, and `StringFastArray<N>`.
 *
 * Storage order follows the global `JIARRAY_COLUMN_MAJOR` macro
 * (same convention as `JIArray`):
 *   - `JIARRAY_COLUMN_MAJOR == 0` → row-major (last dim stride 1)
 *   - `JIARRAY_COLUMN_MAJOR != 0` → column-major (first dim stride 1)
 *
 * Index base follows `JIARRAY_OFFSET` (0 or 1). `operator()` takes
 * exactly `RANK` offset-relative indices; `operator[]` provides
 * 0-based flat access.
 *
 * @tparam T     Element type.
 * @tparam Dims  Per-dimension extents (compile-time). Rank must be ≥ 1.
 *
 * Example:
 * @code
 *   fint<6>       v;          // 1D
 *   fint<3, 4>   m;          // 2D
 *   fint<3, 4, 5> c;          // 3D
 *   fstring<10>  names;       // FastArray<std::string, 10>
 * @endcode
 */
template <typename T, std::size_t... Dims>
class FastArray {
    static_assert(sizeof...(Dims) > 0,
                  "FastArray requires at least one dimension");

public:
    // Public type / size info
    static constexpr std::size_t RANK = sizeof...(Dims);
    static constexpr std::size_t SIZE = (std::size_t{1} * ... * Dims);
    static constexpr bool is_row_major = (JIARRAY_COLUMN_MAJOR == 0);

    T mm[SIZE]{};

private:
    // Compile-time dimension table.  Host-only constexpr member (CUDA device
    // code cannot link static class members), so strides are recomputed in
    // each `linearize` call — the compiler inlines this into a constant.
    static constexpr std::size_t dims_arr_[RANK] = {Dims...};

    template <typename... Idx>
    JIARRAY_HD inline std::size_t linearize(Idx... idx) const {
        static_assert(sizeof...(Idx) == RANK,
                      "operator() requires exactly RANK indices");
        constexpr std::size_t dims_local[RANK] = {Dims...};

        // Compute strides locally — `constexpr` so the compiler folds them.
        std::size_t strides[RANK]{};
        if constexpr (is_row_major) {
            strides[RANK - 1] = 1;
            for (std::size_t i = RANK - 1; i > 0; --i) {
                strides[i - 1] = strides[i] * dims_local[i];
            }
        } else {
            strides[0] = 1;
            for (std::size_t i = 1; i < RANK; ++i) {
                strides[i] = strides[i - 1] * dims_local[i - 1];
            }
        }

        const std::size_t arr[RANK] = {
            static_cast<std::size_t>(idx - static_cast<int>(JIARRAY_OFFSET))...
        };
#if defined(JIARRAY_DEBUG) && !defined(__CUDA_ARCH__)
        {
            const int raw[RANK] = { static_cast<int>(idx)... };
            for (std::size_t d = 0; d < RANK; ++d) {
                JIARRAY_CHECK_BOUND(raw[d],
                                    static_cast<int>(JIARRAY_OFFSET),
                                    static_cast<int>(JIARRAY_OFFSET + dims_local[d] - 1));
            }
        }
#endif
        std::size_t pos = 0;
        for (std::size_t d = 0; d < RANK; ++d) {
            pos += arr[d] * strides[d];
        }
        return pos;
    }

public:
    // ---- constructors ----
    FastArray() = default;

    FastArray(std::initializer_list<T> il) {
        std::size_t i = 0;
        for (auto const& v : il) {
            if (i >= SIZE) break;
            mm[i++] = v;
        }
    }

    explicit FastArray(const T& val) {
        for (std::size_t i = 0; i < SIZE; ++i) mm[i] = val;
    }

    template <std::size_t N>
    explicit FastArray(const T (&a)[N]) {
        static_assert(N == SIZE, "C-array size must match FastArray SIZE");
        for (std::size_t i = 0; i < SIZE; ++i) mm[i] = a[i];
    }

    FastArray(const FastArray&)            = default;
    FastArray(FastArray&&)                 = default;
    FastArray& operator=(const FastArray&) = default;
    FastArray& operator=(FastArray&&)      = default;

    // ---- element access ----
    template <typename... Idx>
    JIARRAY_HD inline T& operator()(Idx... idx) {
        return mm[linearize(idx...)];
    }
    template <typename... Idx>
    JIARRAY_HD inline const T& operator()(Idx... idx) const {
        return mm[linearize(idx...)];
    }

    // 0-based flat access (rank-agnostic)
    JIARRAY_HD inline T&       operator[](std::size_t i)       { return mm[i]; }
    JIARRAY_HD inline const T& operator[](std::size_t i) const { return mm[i]; }

    JIARRAY_HD inline T*       data() noexcept       { return mm; }
    JIARRAY_HD inline const T* data() const noexcept { return mm; }

    // ---- scalar / pointer assignment ----
    inline FastArray& operator=(const T& val) {
        for (std::size_t i = 0; i < SIZE; ++i) mm[i] = val;
        return *this;
    }
    inline FastArray& operator=(const T* p) {
        for (std::size_t i = 0; i < SIZE; ++i) mm[i] = p[i];
        return *this;
    }

    // ---- arithmetic (only valid when T supports the operator) ----
    inline FastArray& operator+=(const T& val) {
        for (std::size_t i = 0; i < SIZE; ++i) mm[i] += val;
        return *this;
    }
    inline FastArray& operator*=(const T& val) {
        for (std::size_t i = 0; i < SIZE; ++i) mm[i] *= val;
        return *this;
    }
    inline FastArray& operator/=(const T& val) {
        for (std::size_t i = 0; i < SIZE; ++i) mm[i] /= val;
        return *this;
    }
    inline FastArray operator*(const T& val) const {
        FastArray r = *this;
        r *= val;
        return r;
    }

    // ---- equality ----
    inline bool operator==(const FastArray& other) const {
        for (std::size_t i = 0; i < SIZE; ++i) {
            if (!(mm[i] == other.mm[i])) return false;
        }
        return true;
    }
    inline bool operator!=(const FastArray& other) const {
        return !(*this == other);
    }

    // ---- searches / reductions ----
    inline int findFirst(const T& val) const {
        for (std::size_t i = 0; i < SIZE; ++i) {
            if (mm[i] == val) return static_cast<int>(i + JIARRAY_OFFSET);
        }
        return static_cast<int>(JIARRAY_OFFSET) - 1;
    }

    inline T min() const {
        T m = mm[0];
        for (std::size_t i = 1; i < SIZE; ++i) {
            if (mm[i] < m) m = mm[i];
        }
        return m;
    }
    inline T max() const {
        T m = mm[0];
        for (std::size_t i = 1; i < SIZE; ++i) {
            if (m < mm[i]) m = mm[i];
        }
        return m;
    }
    JIARRAY_HD inline T sum() const {
        T s = T{};
        for (std::size_t i = 0; i < SIZE; ++i) s += mm[i];
        return s;
    }

    // ---- STL-compatible interface (flat) ----
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

    iterator       begin()        noexcept { return mm; }
    const_iterator begin()  const noexcept { return mm; }
    const_iterator cbegin() const noexcept { return mm; }
    iterator       end()          noexcept { return mm + SIZE; }
    const_iterator end()    const noexcept { return mm + SIZE; }
    const_iterator cend()   const noexcept { return mm + SIZE; }

    reverse_iterator       rbegin()        noexcept { return reverse_iterator(end()); }
    const_reverse_iterator rbegin()  const noexcept { return const_reverse_iterator(end()); }
    const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(end()); }
    reverse_iterator       rend()          noexcept { return reverse_iterator(begin()); }
    const_reverse_iterator rend()    const noexcept { return const_reverse_iterator(begin()); }
    const_reverse_iterator crend()   const noexcept { return const_reverse_iterator(begin()); }

    constexpr size_type size()     const noexcept { return SIZE; }
    constexpr size_type max_size() const noexcept { return SIZE; }
    constexpr bool      empty()    const noexcept { return SIZE == 0; }

    reference       front()       noexcept { return mm[0]; }
    const_reference front() const noexcept { return mm[0]; }
    reference       back()        noexcept { return mm[SIZE - 1]; }
    const_reference back()  const noexcept { return mm[SIZE - 1]; }
};

// Out-of-class definition for static constexpr C-array member.  In C++17 the
// in-class declaration is itself an inline definition, but supplying this
// keeps older toolchains (and CUDA host/device linker) happy.
template <typename T, std::size_t... Dims>
constexpr std::size_t FastArray<T, Dims...>::dims_arr_[FastArray<T, Dims...>::RANK];

// ============================================================================
// Preferred type aliases (variadic — match dimension list directly)
// ============================================================================
template <std::size_t... Dims> using fbool   = FastArray<bool,        Dims...>;
template <std::size_t... Dims> using fchar   = FastArray<char,        Dims...>;
template <std::size_t... Dims> using fshort  = FastArray<short,       Dims...>;
template <std::size_t... Dims> using fint    = FastArray<int,         Dims...>;
template <std::size_t... Dims> using flong   = FastArray<long,        Dims...>;
template <std::size_t... Dims> using ffloat  = FastArray<float,       Dims...>;
template <std::size_t... Dims> using fdouble = FastArray<double,      Dims...>;
template <std::size_t... Dims> using fstring = FastArray<std::string, Dims...>;

template <typename T, std::size_t... Dims>
using farray = FastArray<T, Dims...>;

// ============================================================================
// Class-name compatibility (pre-0.8.0): FastArray2D and StringFastArray were
// separate types.  Now thin aliases over the unified FastArray.
// ============================================================================
template <typename T, std::size_t I, std::size_t J, std::size_t = JIARRAY_OFFSET>
using FastArray2D = FastArray<T, I, J>;

template <std::size_t N, std::size_t = JIARRAY_OFFSET>
using StringFastArray = FastArray<std::string, N>;

// ============================================================================
// Backward-compatibility aliases (rank-suffixed; deprecated, kept for
// downstream code still using the pre-0.8.0 names).
// ============================================================================
template <int I, int = JIARRAY_OFFSET>
using fbool1d   = FastArray<bool,        static_cast<std::size_t>(I)>;
template <int I, int = JIARRAY_OFFSET>
using fint1d    = FastArray<int,         static_cast<std::size_t>(I)>;
template <int I, int = JIARRAY_OFFSET>
using ffloat1d  = FastArray<float,       static_cast<std::size_t>(I)>;
template <int I, int = JIARRAY_OFFSET>
using fdouble1d = FastArray<double,      static_cast<std::size_t>(I)>;
template <int I, int = JIARRAY_OFFSET>
using fstring1d = FastArray<std::string, static_cast<std::size_t>(I)>;

template <int I, int J, int = JIARRAY_OFFSET>
using fbool2d   = FastArray<bool,   static_cast<std::size_t>(I), static_cast<std::size_t>(J)>;
template <int I, int J, int = JIARRAY_OFFSET>
using fint2d    = FastArray<int,    static_cast<std::size_t>(I), static_cast<std::size_t>(J)>;
template <int I, int J, int = JIARRAY_OFFSET>
using ffloat2d  = FastArray<float,  static_cast<std::size_t>(I), static_cast<std::size_t>(J)>;
template <int I, int J, int = JIARRAY_OFFSET>
using fdouble2d = FastArray<double, static_cast<std::size_t>(I), static_cast<std::size_t>(J)>;

} // namespace dnegri::jiarray
