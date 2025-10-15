#pragma once

/**
 * @file JIArray.h
 * @brief High-performance multidimensional array container with flexible indexing
 *
 * This file contains the JIArray template class, a feature-rich N-dimensional array
 * implementation optimized for scientific computing and numerical analysis. It supports
 * both row-major and column-major storage orders, customizable indexing offsets,
 * array slicing, and vectorized operations.
 *
 * Key Features:
 * - Configurable storage order (row-major/column-major)
 * - Customizable index offsets (0-based or 1-based)
 * - Efficient slicing with view semantics
 * - Conditional SIMD vectorization
 * - STL-compatible iterators
 * - Compile-time type safety
 *
 * @author dnegri
 * @version 0.6
 * @date 2025-10-06
 */

#include "FastArray.h"
#include "pch.h"
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

namespace dnegri::jiarray {

// ============================================================================
// Compile-time configuration for OpenMP usage
// ============================================================================

/**
 * @def JIARRAY_USE_SIMD
 * @brief Control SIMD vectorization in JIArray operations
 *
 * Set to 1 to enable SIMD directives, 0 to disable.
 * Should be disabled (0) when JIArray methods are called from within
 * OpenMP parallel regions to avoid overhead.
 *
 * Usage:
 * @code
 * // For use inside OpenMP parallel regions
 * #define JIARRAY_USE_SIMD 0
 *
 * // For standalone use with SIMD
 * #define JIARRAY_USE_SIMD 1
 * @endcode
 */
#ifndef JIARRAY_USE_SIMD
    #define JIARRAY_USE_SIMD 1
#endif

// Helper macros for conditional SIMD - uses _Pragma for portability
#if JIARRAY_USE_SIMD
    #define JIARRAY_SIMD_LOOP _Pragma("omp simd")
// #define JIARRAY_SIMD_REDUCTION(op, var) _Pragma("omp simd reduction(" #op ":" #var ")")
#else
    #define JIARRAY_SIMD_LOOP
// #define JIARRAY_SIMD_REDUCTION(op, var)
#endif

// ============================================================================
// Compile-time type traits and helpers
// ============================================================================

/**
 * @brief Type trait to detect JIArray types at compile time
 * @tparam T Type to check
 *
 * Base case: returns false_type for non-JIArray types.
 */
template <typename T>
struct is_jiarray : std::false_type {};

/// Forward declaration of JIArray template class
template <class T, size_t RANK = 1, class SEQ = std::make_index_sequence<RANK>>
class JIArray;

/**
 * @brief Specialization of is_jiarray for JIArray types
 * @tparam T Element type
 * @tparam RANK Array rank (dimensionality)
 * @tparam SEQ Index sequence
 */
template <class T, size_t RANK, class SEQ>
struct is_jiarray<JIArray<T, RANK, SEQ>> : std::true_type {};

/**
 * @brief Helper variable template for is_jiarray
 * @tparam T Type to check
 */
template <typename T>
inline constexpr bool is_jiarray_v = is_jiarray<T>::value;

/**
 * @brief Type trait to check if all types are integral
 * @tparam Args Variadic types to check
 */
template <typename... Args>
struct all_integral : std::conjunction<std::is_integral<Args>...> {};

/**
 * @brief Helper variable template for all_integral
 * @tparam Args Variadic types to check
 */
template <typename... Args>
inline constexpr bool all_integral_v = all_integral<Args...>::value;

/**
 * @brief Type trait to check if type is arithmetic (for scalar operations)
 * @tparam T Type to check
 */
template <typename T>
inline constexpr bool is_scalar_v = std::is_arithmetic_v<std::decay_t<T>>;

/**
 * @brief Compile-time index sequence validation helper
 * @tparam N Number of indices provided
 * @tparam RANK Expected array rank
 */
template <size_t N, size_t RANK>
struct validate_index_count {
    static_assert(N == RANK, "Number of indices must match array rank");
    static constexpr bool value = (N == RANK);
};

// ============================================================================
// Main JIArray class
// ============================================================================

/**
 * @class JIArray
 * @brief High-performance N-dimensional array container with flexible memory layout
 *
 * JIArray is a template class providing efficient multidimensional array operations
 * with the following features:
 * - Configurable storage order (row-major or column-major)
 * - Customizable index offsets (0-based or 1-based indexing)
 * - Efficient slicing operations with view semantics
 * - Conditional vectorized arithmetic operations
 * - STL-compatible iterators
 * - Compile-time type safety
 *
 * @tparam T Element type (must be default constructible)
 * @tparam RANK Array dimensionality (must be > 0)
 * @tparam INTS Parameter pack for index sequence
 *
 * @section example_usage Example Usage
 * @code
 * // Create a 3x4 matrix with 1-based indexing
 * zint2 matrix;
 * matrix.init(3, 4);
 *
 * // Access elements
 * matrix(1, 1) = 42;
 *
 * // Slice to get a row
 * zint1 row = matrix.slice(2);
 *
 * // Arithmetic operations
 * matrix *= 2;
 * auto sum = matrix + matrix;
 * @endcode
 *
 * @section openmp_usage OpenMP Usage
 * @code
 * // For use inside OpenMP parallel regions, disable SIMD
 * #define JIARRAY_USE_SIMD 0
 * #include "JIArray.h"
 *
 * #pragma omp parallel for
 * for (int i = 0; i < n; i++) {
 *     array(i) = compute(i);  // Safe: no nested parallelism
 * }
 * @endcode
 *
 * @note The storage order is determined by the JIARRAY_COLUMN_MAJOR macro:
 *       - JIARRAY_COLUMN_MAJOR=0: Row-major (C-style)
 *       - JIARRAY_COLUMN_MAJOR=1: Column-major (Fortran-style)
 *
 * @warning When using JIArray inside OpenMP parallel regions, set JIARRAY_USE_SIMD=0
 *          to avoid nested parallelism overhead.
 */
template <class T, size_t RANK, size_t... INTS>
class JIArray<T, RANK, std::index_sequence<INTS...>> {
private:
    // Compile-time checks
    static_assert(RANK > 0, "Array rank must be at least 1");
    static_assert(std::is_default_constructible_v<T>, "Element type must be default constructible");

    /// Storage order flag: true for row-major, false for column-major
    constexpr static bool is_row_major = (JIARRAY_COLUMN_MAJOR == 0);

    // Member variables
    int nn        = 0;                      ///< Total number of elements
    T*  mm        = nullptr;                ///< Pointer to data storage
    int allocated = JIARRAY_ALLOCATED_NONE; ///< Allocation status flag
    int rankSize[RANK]{};                   ///< Stride for each dimension
    int offset[RANK]{};                     ///< Index offset for each dimension
    int sumOfOffset = 0;                    ///< Precomputed sum of offsets
    int sizes[RANK]{};                      ///< Size of each dimension

    // ========================================================================
    // Private helper methods
    // ========================================================================

    /**
     * @brief Calculate linear index from multidimensional indices
     * @tparam Args Index types (must be integral)
     * @param indices Multidimensional indices
     * @return Linear index in memory
     * @note Uses manual loop unrolling for small ranks (≤ 4) for better performance
     */
    template <typename... Args>
    inline int calculateIndex(Args... indices) const {
        static_assert(sizeof...(indices) == RANK, "Number of indices must match array rank");
        static_assert(all_integral_v<Args...>, "All indices must be integral types");

        const int idx[] = {static_cast<int>(indices)...};
        int       pos   = -sumOfOffset;

        // Manual unrolling for small ranks
        if constexpr (RANK <= 4) {
#pragma unroll
            for (int i = 0; i < RANK; i++) {
                JIARRAY_CHECK_BOUND(idx[i], offset[i], offset[i] + sizes[i] - 1);
                pos += rankSize[i] * idx[i];
            }
        } else {
            for (int i = 0; i < RANK; i++) {
                JIARRAY_CHECK_BOUND(idx[i], offset[i], offset[i] + sizes[i] - 1);
                pos += rankSize[i] * idx[i];
            }
        }

        return pos;
    }

    /**
     * @brief Calculate stride sizes for each dimension based on storage order
     * @param dimensions Array of dimension sizes
     * @note For row-major: rightmost dimension has stride 1
     *       For column-major: leftmost dimension has stride 1
     */
    constexpr void calculateRankSize(const int dimensions[RANK]) {
        if constexpr (is_row_major) {
            // Row-major: rightmost dimension has stride 1
            rankSize[RANK - 1] = 1;
            if constexpr (RANK > 1) {
#pragma unroll
                for (int i = static_cast<int>(RANK) - 2; i >= 0; i--) {
                    rankSize[i] = rankSize[i + 1] * dimensions[i + 1];
                }
            }
        } else {
            // Column-major: leftmost dimension has stride 1
            rankSize[0] = 1;
            if constexpr (RANK > 1) {
#pragma unroll
                for (int i = 1; i < RANK; i++) {
                    rankSize[i] = rankSize[i - 1] * dimensions[i - 1];
                }
            }
        }
    }

    /**
     * @brief Calculate memory offset for slice operation
     * @tparam INDEX Index types (must be integral)
     * @param index Indices specifying the slice position
     * @return Memory offset for the sliced subarray
     */
    template <typename... INDEX>
    inline int calculateSliceOffset(INDEX... index) const {
        constexpr int num_idx = sizeof...(INDEX);
        static_assert(num_idx < RANK, "Slice indices must be less than array rank");

        const int idx[] = {static_cast<int>(index)...};
        int       p_mm  = 0;

        if constexpr (is_row_major) {
#pragma unroll
            for (int i = 0; i < num_idx; i++) {
                JIARRAY_CHECK_BOUND(idx[i], offset[i], offset[i] + sizes[i] - 1);
                p_mm += rankSize[i] * (idx[i] - offset[i]);
            }
        } else {
            int rank = RANK;
            for (int i = num_idx - 1; i >= 0; i--) {
                rank--;
                JIARRAY_CHECK_BOUND(idx[i], offset[rank], offset[rank] + sizes[rank] - 1);
                p_mm += rankSize[rank] * (idx[i] - offset[rank]);
            }
        }

        return p_mm;
    }

public:
    // ========================================================================
    // Type aliases for better interface
    // ========================================================================

    using value_type      = T;              ///< Element type
    using size_type       = int;            ///< Size type
    using difference_type = std::ptrdiff_t; ///< Difference type for iterators
    using reference       = T&;             ///< Reference type
    using const_reference = const T&;       ///< Const reference type
    using pointer         = T*;             ///< Pointer type
    using const_pointer   = const T*;       ///< Const pointer type

    /// Type alias for this instantiation
    using this_type = JIArray<T, RANK, std::index_sequence<INTS...>>;

    static constexpr size_t rank      = RANK;         ///< Array dimensionality
    static constexpr bool   row_major = is_row_major; ///< Storage order flag

    // ========================================================================
    // Constructors and destructor
    // ========================================================================

    /**
     * @brief Default constructor creates an empty array
     */
    JIArray() = default;

    /**
     * @brief Construct 1D array from initializer list
     * @param il Initializer list with elements
     * @note Only available for RANK == 1
     */
    JIArray(std::initializer_list<T> il) : JIArray(static_cast<int>(il.size())) {
        static_assert(RANK == 1, "Initializer list constructor only available for 1D arrays");
        std::copy(il.begin(), il.end(), mm);
    }

    /**
     * @brief Copy constructor (creates a view, not a deep copy)
     * @param array Source array
     * @note This creates a shallow copy that shares memory with the source
     */
    JIArray(const this_type& array) noexcept {
        nn = array.nn;
        mm = array.mm;
        std::copy(array.rankSize, array.rankSize + RANK, rankSize);
        std::copy(array.offset, array.offset + RANK, offset);
        allocated   = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
        std::copy(array.sizes, array.sizes + RANK, sizes);
    }

    /**
     * @brief Virtual destructor
     * @note Automatically deallocates memory if owned by this array
     */
    virtual ~JIArray() {
        destroy();
    }

    /**
     * @brief Construct array with specified dimensions
     * @tparam Args Dimension types (must be integral)
     * @param args Dimension sizes
     * @note Number of arguments must match RANK
     */
    template <typename... Args,
              typename = std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK)>>
    explicit JIArray(Args... args) {
        init(args...);
    }

    /**
     * @brief Construct array wrapping external memory
     * @param memory_ Pointer to external memory
     * @param args Dimension sizes
     */
    JIArray(T* memory_, decltype(INTS)... args) {
        static_assert(sizeof...(INTS) == RANK, "Number of dimensions must match rank");
        init(args..., memory_);
    }

    // ========================================================================
    // Initialization methods
    // ========================================================================

    /**
     * @brief Initialize with range specifications (min, max pairs)
     * @tparam INTS2 Integer types for range specifications
     * @param array_sizes Alternating min and max values for each dimension
     * @note Expects 2*RANK arguments: min1, max1, min2, max2, ...
     *
     * @code
     * zint2 arr;
     * arr.init0(1, 3, 1, 4);  // Creates 3x4 array with indices [1:3, 1:4]
     * @endcode
     */
    template <typename... INTS2>
    std::enable_if_t<(sizeof...(INTS2) == 2 * RANK), void>
    init0(INTS2... array_sizes) {
        static_assert(all_integral_v<INTS2...>, "All size parameters must be integral");

        const int temp_size[] = {static_cast<int>(array_sizes)...};
        int       dimensions[RANK];
        int       offsets[RANK];

#pragma unroll
        for (int i = 0; i < RANK; i++) {
            dimensions[i] = temp_size[2 * i + 1] - temp_size[2 * i] + 1;
            offsets[i]    = temp_size[2 * i];
        }

        init(dimensions);
        setOffsets(offsets);
    }

    /**
     * @brief Initialize array with specified dimensions
     * @tparam Args Dimension types (must be integral)
     * @param array_sizes Size of each dimension
     * @note Allocates new memory and uses default offset (JIARRAY_OFFSET)
     */
    template <typename... Args>
    std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK), void>
    init(Args... array_sizes) {
        const int dimensions[] = {static_cast<int>(array_sizes)...};
        init(dimensions);
    }

    /**
     * @brief Initialize array with dimension array
     * @param dimensions Array of dimension sizes
     * @note Destroys existing data if any, allocates new memory
     */
    void init(const int dimensions[RANK]) {
        destroy();

        // Check for zero dimensions
        for (int i = 0; i < RANK; i++) {
            if (dimensions[i] == 0)
                return;
        }

        calculateRankSize(dimensions);

        // Initialize offsets
        std::fill(offset, offset + RANK, JIARRAY_OFFSET);

        // Calculate total size and sum of offsets
        nn = 1;
        for (int i = 0; i < RANK; i++) {
            nn *= dimensions[i];
        }

        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++) {
            sumOfOffset += offset[i] * rankSize[i];
        }

        mm        = new T[nn]{};
        allocated = JIARRAY_ALLOCATED_ALL;
        std::copy(dimensions, dimensions + RANK, sizes);
    }

    /**
     * @brief Initialize array to wrap external memory
     * @param array_sizes Dimension sizes
     * @param memory_ Pointer to external memory buffer
     * @note Array will not own the memory and won't deallocate it
     */
    void init(decltype(INTS)... array_sizes, T* memory_) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        const int dimensions[] = {static_cast<int>(array_sizes)...};
        nn                     = (array_sizes * ...);
        mm                     = memory_;

        calculateRankSize(dimensions);

        std::fill(offset, offset + RANK, JIARRAY_OFFSET);

        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++) {
            sumOfOffset += offset[i] * rankSize[i];
        }

        allocated = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
        std::copy(dimensions, dimensions + RANK, sizes);
    }

    /**
     * @brief Initialize by directly specifying rank sizes and offsets (advanced)
     * @param size Total number of elements
     * @param rankSizes Array of stride values for each dimension
     * @param offsets Array of index offsets for each dimension
     * @param memory Pointer to memory buffer
     * @note This is an advanced method for creating array views
     */
    void initByRankSize(int size, const int* rankSizes, const int* offsets, T* memory) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        nn          = size;
        sumOfOffset = 0;

        for (int i = 0; i < RANK; i++) {
            rankSize[i] = rankSizes[i];
            offset[i]   = offsets[i];
            sumOfOffset += offset[i] * rankSize[i];
        }

        mm        = memory;
        allocated = JIARRAY_ALLOCATED_NONE;

        // Reconstruct sizes from rank information
        if constexpr (is_row_major) {
            sizes[0] = (nn != 0 && rankSize[0] != 0) ? nn / rankSize[0] : 0;
            for (int i = 1; i < RANK; i++) {
                sizes[i] = (rankSize[i] != 0) ? rankSize[i - 1] / rankSize[i] : 0;
            }
        } else {
            for (int i = 0; i < RANK - 1; i++) {
                sizes[i] = (rankSize[i] != 0) ? rankSize[i + 1] / rankSize[i] : 0;
            }
            sizes[RANK - 1] = (nn != 0 && rankSize[RANK - 1] != 0) ? nn / rankSize[RANK - 1] : 0;
        }
    }

    /**
     * @brief Initialize by rank sizes with new memory allocation
     * @param size Total number of elements
     * @param rankSizes Array of stride values
     * @param offsets Array of index offsets
     */
    void initByRankSize(int size, const int* rankSizes, const int* offsets) {
        destroy();
        if (size == 0)
            return;

        mm = new T[size]{};
        initByRankSize(size, rankSizes, offsets, mm);
        allocated = JIARRAY_ALLOCATED_ALL;
    }

    // ========================================================================
    // Memory management
    // ========================================================================

    /**
     * @brief Destroy array and deallocate memory if owned
     * @note Only deallocates if this array owns the memory
     */
    void destroy() {
        if (allocated != JIARRAY_ALLOCATED_NONE) {
            if ((allocated & JIARRAY_ALLOCATED_MEMORY) != 0 && mm != nullptr) {
                delete[] mm;
            }
            mm        = nullptr;
            allocated = JIARRAY_ALLOCATED_NONE;
            nn        = 0;
        }
    }

    /**
     * @brief Erase array completely, resetting all metadata
     * @note Deallocates memory and resets all dimension information
     */
    void erase() {
        destroy();
        std::fill(rankSize, rankSize + RANK, 0);
        std::fill(offset, offset + RANK, 0);
        std::fill(sizes, sizes + RANK, 0);
        sumOfOffset = 0;
    }

    // ========================================================================
    // Offset and size management
    // ========================================================================

    /**
     * @brief Set index offsets for all dimensions
     * @tparam Args Offset types (must be integral)
     * @param offsets Offset value for each dimension
     * @note Updates internal stride calculations
     */
    template <typename... Args>
    std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK), void> inline setOffsets(Args... offsets) noexcept {
        const int tempOffsets[] = {static_cast<int>(offsets)...};
        setOffsets(tempOffsets);
    }

    /**
     * @brief Set index offsets from array
     * @param offsets Array of offset values
     */
    inline void setOffsets(const int offsets[RANK]) noexcept {
        std::copy(offsets, offsets + RANK, offset);
        sumOfOffset = 0;
        for (int i = 0; i < RANK; i++) {
            sumOfOffset += rankSize[i] * offset[i];
        }
    }

    /**
     * @brief Resize array if dimensions differ
     * @tparam Args Size types (must be integral)
     * @param array_sizes New dimension sizes
     * @note Only reallocates if sizes actually changed
     */
    template <typename... Args>
    std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK), void>
    setSize(Args... array_sizes) {
        const int dimensions[] = {static_cast<int>(array_sizes)...};

        // Check if sizes are the same
        bool same = std::equal(sizes, sizes + RANK, dimensions);

        if (same)
            return;

        destroy();
        init(array_sizes...);
    }

    // ========================================================================
    // Slicing operations with proper const-correctness
    // ========================================================================

    /**
     * @brief Create a mutable view (slice) of the array
     * @tparam INDEX Index types (must be integral)
     * @param index Indices specifying which slice to extract
     * @return Lower-dimensional array view
     * @note Returns a view that shares memory with the original array
     *
     * @code
     * zint3 arr(10, 20, 30);
     * zint2 slice1 = arr.slice(5);      // Extract arr(5, :, :)
     * zint1 slice2 = arr.slice(5, 10);  // Extract arr(5, 10, :)
     * @endcode
     */
    template <typename... INDEX>
    inline JIArray<T, RANK - sizeof...(INDEX)> slice(INDEX... index) {
        constexpr int RANK2 = RANK - sizeof...(INDEX);
        static_assert(RANK2 > 0, "Slice rank must be at least 1");
        static_assert(sizeof...(INDEX) < RANK, "Number of slice indices must be less than rank");
        constexpr int num_idx = sizeof...(INDEX);

        int p_mm = calculateSliceOffset(index...);

        auto array = JIArray<T, RANK2>();
        if constexpr (is_row_major) {
            array.initByRankSize(rankSize[num_idx - 1], rankSize + num_idx, offset + num_idx,
                                 mm + p_mm);
        } else {
            array.initByRankSize(rankSize[RANK2], rankSize, offset, mm + p_mm);
        }
        return array;
    }

    /**
     * @brief Create a const view (slice) of the array
     * @tparam INDEX Index types (must be integral)
     * @param index Indices specifying which slice to extract
     * @return Lower-dimensional const array view
     * @note Returns a const view that prevents modification of the original data
     */
    template <typename... INDEX>
    inline const JIArray<T, RANK - sizeof...(INDEX)> slice(INDEX... index) const {
        constexpr int RANK2 = RANK - sizeof...(INDEX);
        static_assert(RANK2 > 0, "Slice rank must be at least 1");
        static_assert(sizeof...(INDEX) < RANK, "Number of slice indices must be less than rank");
        constexpr int num_idx = sizeof...(INDEX);

        int p_mm = calculateSliceOffset(index...);

        auto array = JIArray<T, RANK2>();
        if constexpr (is_row_major) {
            array.initByRankSize(rankSize[num_idx - 1], rankSize + num_idx, offset + num_idx,
                                 const_cast<T*>(mm) + p_mm);
        } else {
            array.initByRankSize(rankSize[RANK2], rankSize, offset, const_cast<T*>(mm) + p_mm);
        }
        return array;
    }

    // ========================================================================
    // Accessors and properties
    // ========================================================================

    /**
     * @brief Check if array has been allocated
     * @return true if array owns or references memory
     */
    inline bool isAllocated() const noexcept {
        return allocated != JIARRAY_ALLOCATED_NONE;
    }

    /**
     * @brief Get pointer to memory at specified offset
     * @param idx Linear offset from array start
     * @return Pointer to memory location
     */
    inline T* getMemory(int idx = 0) noexcept {
        return mm + idx;
    }

    /// @overload
    inline const T* getMemory(int idx = 0) const noexcept {
        return mm + idx;
    }

    /**
     * @brief Get raw pointer to data
     * @return Pointer to first element
     */
    inline T* data() noexcept {
        return mm;
    }

    /// @overload
    inline const T* data() const noexcept {
        return mm;
    }

    /**
     * @brief Get total number of elements
     * @return Total element count
     */
    inline const int& size() const noexcept {
        return nn;
    }

    /// @overload
    inline int getSize() const noexcept {
        return nn;
    }

    /**
     * @brief Get array of stride values
     * @return Pointer to rankSize array
     */
    inline const int* getRankSize() const noexcept {
        return rankSize;
    }

    /**
     * @brief Get array of dimension sizes
     * @return Pointer to sizes array
     */
    inline const int* getSizeOfRank() const noexcept {
        return sizes;
    }

    /**
     * @brief Get size of specific dimension
     * @param rank Dimension number (1-based or 0-based depending on JIARRAY_OFFSET)
     * @return Size of specified dimension
     */
    inline const int& getSize(int rank) const {
        JIARRAY_CHECK_BOUND(rank, JIARRAY_OFFSET, JIARRAY_OFFSET + RANK - 1);
        return sizes[rank - JIARRAY_OFFSET];
    }

    /**
     * @brief Get array of index offsets
     * @return Pointer to offset array
     */
    inline const int* getOffset() const noexcept {
        return offset;
    }

    /**
     * @brief Get offset for specific dimension
     * @param rank Dimension number
     * @return Offset value for specified dimension
     */
    inline const int& getOffset(int rank) const noexcept {
        return offset[rank - JIARRAY_OFFSET];
    }

    // ========================================================================
    // Element access operators
    // ========================================================================

    /**
     * @brief Access element with bounds checking
     * @tparam Args Index types (must be integral)
     * @param index Multidimensional indices
     * @return Reference to element
     * @note Performs bounds checking in debug mode
     */
    template <typename... Args,
              typename = std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK)>>
    inline T& at(Args... index) {
        int pos = calculateIndex(index...);
        return mm[pos];
    }

    /// @overload
    template <typename... Args,
              typename = std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK)>>
    inline const T&
    at(Args... index) const {
        int pos = calculateIndex(index...);
        return mm[pos];
    }

    /**
     * @brief Function call operator for element access
     * @tparam Args Index types (must be integral)
     * @param index Multidimensional indices
     * @return Reference to element
     */
    template <typename... Args,
              typename = std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK)>>
    inline const T& operator()(Args... index) const {
        return at(index...);
    }

    /// @overload
    template <typename... Args,
              typename = std::enable_if_t<all_integral_v<Args...> && (sizeof...(Args) == RANK)>>
    inline T& operator()(Args... index) {
        return at(index...);
    }

    /**
     * @brief Access element using FastArray index
     * @param idx FastArray containing indices
     * @return Reference to element
     */
    inline T& operator()(const FastArray<int, RANK>& idx) {
        int pos = -sumOfOffset;

        for (int i = 0; i < RANK; i++) {
            JIARRAY_CHECK_BOUND(idx(i + 1), offset[i], offset[i] + sizes[i] - 1);
            pos += rankSize[i] * idx(i + 1);
        }

        return mm[pos];
    }

    /// @overload
    inline const T& operator()(const FastArray<int, RANK>& idx) const {
        return const_cast<this_type&>(*this)(idx);
    }

    /**
     * @brief Get pointer to subarray at specified indices
     * @tparam INDEX Index types (must be integral)
     * @param index Partial indices
     * @return Pointer to memory location
     * @note Number of indices must be between 1 and RANK
     */
    template <typename... INDEX>
    inline std::enable_if_t<all_integral_v<INDEX...>, T*>
    data(INDEX... index) {
        constexpr int num_idx = sizeof...(index);
        static_assert(num_idx >= 1 && num_idx <= RANK, "Number of indices must be between 1 and RANK");

        const int idx[] = {static_cast<int>(index)...};
        int       pos   = 0;

        if constexpr (is_row_major) {
            for (int i = 0; i < num_idx; i++) {
                pos += rankSize[i] * (idx[i] - offset[i]);
            }
        } else {
            int rank = RANK;
            for (int i = num_idx - 1; i >= 0; i--) {
                --rank;
                pos += rankSize[rank] * (idx[i] - offset[rank]);
            }
        }

        return mm + pos;
    }

    /// @overload
    template <typename... INDEX>
    inline std::enable_if_t<all_integral_v<INDEX...>, const T*>
    data(INDEX... index) const {
        return const_cast<this_type&>(*this).data(index...);
    }

    // ========================================================================
    // Reshape
    // ========================================================================

    /**
     * @brief Reshape array to different dimensionality
     * @tparam INTS2 New dimension types
     * @param sizes New dimension sizes
     * @return Reshaped array view
     * @note Total element count must remain the same
     */
    template <typename... INTS2>
    inline std::enable_if_t<all_integral_v<INTS2...>, JIArray<T, sizeof...(INTS2)>>
    reshape(INTS2... sizes) {

        int total_size = (1 * ... * sizes);

        JIARRAY_CHECK_SIZE(nn, total_size);

        constexpr int     RANK2 = sizeof...(INTS2);
        JIArray<T, RANK2> array;
        array.init(sizes..., getMemory());
        return array;
    }

    // ========================================================================
    // Assignment operators
    // ========================================================================

    /**
     * @brief Assign from initializer list
     * @param list Initializer list of values
     * @return Reference to this array
     */
    inline this_type& operator=(const std::initializer_list<T>& list) {
        JIARRAY_CHECK_SIZE(nn, list.size());
        std::copy(list.begin(), list.end(), mm);
        return *this;
    }

    /**
     * @brief Assign scalar value to all elements
     * @param val Value to assign
     * @return Reference to this array
     */
    inline this_type& operator=(const T& val) {
        assert(nn > 0);
        std::fill(mm, mm + nn, val);
        return *this;
    }

    /**
     * @brief Assign from std::vector
     * @param val Vector of values
     * @return Reference to this array
     * @note If array is empty, initializes to match vector size
     */
    inline this_type& operator=(const std::vector<T>& val) {
        if (nn == 0) {
            int dimensions[RANK];
            std::fill(dimensions, dimensions + RANK, 1);
            dimensions[RANK - 1] = static_cast<int>(val.size());
            init(dimensions);
        } else {
            JIARRAY_CHECK_SIZE(nn, val.size());
        }

        std::copy(val.begin(), val.end(), mm);
        return *this;
    }

    /**
     * @brief Assign from C-style array
     * @param array Pointer to source array
     * @return Reference to this array
     */
    inline this_type& operator=(const T* array) {
        assert(nn > 0);
        std::copy(array, array + nn, mm);
        return *this;
    }

    /**
     * @brief Assign from another JIArray (deep copy)
     * @param array Source array
     * @return Reference to this array
     * @note Allocates memory if this array is uninitialized
     */
    inline this_type& operator=(const this_type& array) {
        if (allocated == JIARRAY_ALLOCATED_NONE && mm == nullptr) {
            initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        } else {
            JIARRAY_CHECK_SIZE(nn, array.getSize());
            auto arraySizeOfRank = array.getSizeOfRank();
            auto arrayOffset     = array.getOffset();
            for (int rank = 0; rank < RANK; ++rank) {
                assert(sizes[rank] == arraySizeOfRank[rank]);
                assert(offset[rank] == arrayOffset[rank]);
            }
        }

        auto arrayMemory = array.data();
        std::copy(arrayMemory, arrayMemory + nn, mm);
        return *this;
    }

    // ========================================================================
    // Statistical operations - with conditional SIMD
    // ========================================================================

    /**
     * @brief Calculate arithmetic mean of all elements
     * @return Average value
     * @note Requires T to be arithmetic type
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline T average() const {
        static_assert(std::is_arithmetic_v<T>, "Average requires arithmetic type");
        T result = T{};

        // Conditional SIMD based on JIARRAY_USE_SIMD
        // JIARRAY_SIMD_REDUCTION("+", result)
        for (int i = 0; i < nn; ++i) {
            result += mm[i];
        }
        return result / nn;
    }

    /**
     * @brief Calculate sum of all elements
     * @return Sum of elements
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline T sum() const {
        T result = T{};
        // JIARRAY_SIMD_REDUCTION("+", result)
        for (int i = 0; i < nn; ++i) {
            result += mm[i];
        }
        return result;
    }

    /**
     * @brief Find maximum element value
     * @return Maximum value
     */
    inline T max() const {
        assert(nn > 0);
        return *std::max_element(mm, mm + nn);
    }

    /**
     * @brief Find minimum element value
     * @return Minimum value
     */
    inline T min() const {
        assert(nn > 0);
        return *std::min_element(mm, mm + nn);
    }

    // ========================================================================
    // Comparison operators
    // ========================================================================

    /**
     * @brief Element-wise equality comparison
     * @param array Array to compare with
     * @return true if all elements are equal
     */
    inline bool operator==(const this_type& array) const {
        if (nn != array.nn)
            return false;
        return std::equal(mm, mm + nn, array.mm);
    }

    /**
     * @brief Element-wise inequality comparison
     * @param array Array to compare with
     * @return true if any element differs
     */
    inline bool operator!=(const this_type& array) const {
        return !(*this == array);
    }

    // ========================================================================
    // Arithmetic operators - with conditional SIMD
    // ========================================================================

    /**
     * @brief Add scalar to all elements in-place
     * @param val Scalar value to add
     * @return Reference to this array
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type& operator+=(const T& val) {
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            mm[i] += val;
        }
        return *this;
    }

    /**
     * @brief Add another array element-wise in-place
     * @param array Array to add
     * @return Reference to this array
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type& operator+=(const this_type& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            mm[i] += array.mm[i];
        }
        return *this;
    }

    /**
     * @brief Element-wise addition of two arrays
     * @param array Array to add
     * @return New array containing sum
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type operator+(const this_type& array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        this_type result;
        result.initByRankSize(getSize(), getRankSize(), getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            result.mm[i] = mm[i] + array.mm[i];
        }
        return result;
    }

    /**
     * @brief Subtract another array element-wise in-place
     * @param array Array to subtract
     * @return Reference to this array
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type& operator-=(const this_type& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            mm[i] -= array.mm[i];
        }
        return *this;
    }

    /**
     * @brief Element-wise subtraction of two arrays
     * @param array Array to subtract
     * @return New array containing difference
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type operator-(const this_type& array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        this_type result;
        result.initByRankSize(getSize(), getRankSize(), getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            result.mm[i] = mm[i] - array.mm[i];
        }
        return result;
    }

    /**
     * @brief Unary negation operator
     * @param array Array to negate
     * @return New array with negated elements
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    friend inline this_type operator-(const this_type& array) {
        this_type result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < array.nn; ++i) {
            result.mm[i] = -array.mm[i];
        }
        return result;
    }

    /**
     * @brief Multiply all elements by scalar in-place
     * @param val Scalar multiplier
     * @return Reference to this array
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type& operator*=(const T& val) {
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            mm[i] *= val;
        }
        return *this;
    }

    /**
     * @brief Element-wise multiplication in-place
     * @param array Array to multiply with
     * @return Reference to this array
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type& operator*=(const this_type& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            mm[i] *= array.mm[i];
        }
        return *this;
    }

    /**
     * @brief Element-wise multiplication of two arrays
     * @param array Array to multiply with
     * @return New array containing product
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type operator*(const this_type& array) const {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        this_type result;
        result.initByRankSize(getSize(), getRankSize(), getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            result.mm[i] = mm[i] * array.mm[i];
        }
        return result;
    }

    /**
     * @brief Divide all elements by scalar in-place
     * @param val Scalar divisor
     * @note For floating-point types, uses reciprocal multiplication for performance
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline void operator/=(const T& val) {
        static_assert(std::is_floating_point_v<T> || std::is_integral_v<T>,
                      "Division requires arithmetic type");
        if constexpr (std::is_floating_point_v<T>) {
            const T rval = T{1} / val;
            JIARRAY_SIMD_LOOP
            for (int i = 0; i < nn; ++i) {
                mm[i] *= rval;
            }
        } else {
            JIARRAY_SIMD_LOOP
            for (int i = 0; i < nn; ++i) {
                mm[i] /= val;
            }
        }
    }

    /**
     * @brief Element-wise division in-place
     * @param array Array to divide by
     * @return Reference to this array
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline this_type& operator/=(const this_type& array) {
        JIARRAY_CHECK_SIZE(nn, array.nn);
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < nn; ++i) {
            mm[i] /= array.mm[i];
        }
        return *this;
    }

    /**
     * @brief Element-wise division of two arrays
     * @param lhs Numerator array
     * @param rhs Denominator array
     * @return New array containing quotient
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    friend inline this_type operator/(const this_type& lhs, const this_type& rhs) {
        JIARRAY_CHECK_SIZE(lhs.nn, rhs.nn);
        this_type result;
        result.initByRankSize(lhs.getSize(), lhs.getRankSize(), lhs.getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < lhs.nn; ++i) {
            result.mm[i] = lhs.mm[i] / rhs.mm[i];
        }
        return result;
    }

    // ========================================================================
    // Scalar arithmetic operations - with conditional SIMD
    // ========================================================================

    /**
     * @brief Add scalar to array (commutative)
     * @tparam Scalar Arithmetic type
     * @param val Scalar value
     * @param array Array operand
     * @return New array with scalar added to each element
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    template <typename Scalar>
    friend inline std::enable_if_t<is_scalar_v<Scalar>, this_type>
    operator+(const Scalar& val, const this_type& array) {
        this_type result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val + array.mm[i];
        }
        return result;
    }

    /// @overload
    template <typename Scalar>
    friend inline std::enable_if_t<is_scalar_v<Scalar>, this_type>
    operator+(const this_type& array, const Scalar& val) {
        return val + array;
    }

    /**
     * @brief Multiply array by scalar (commutative)
     * @tparam Scalar Arithmetic type
     * @param val Scalar multiplier
     * @param array Array operand
     * @return New array with each element multiplied by scalar
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    template <typename Scalar>
    friend inline std::enable_if_t<is_scalar_v<Scalar>, this_type>
    operator*(const Scalar& val, const this_type& array) {
        this_type result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val * array.mm[i];
        }
        return result;
    }

    /// @overload
    template <typename Scalar>
    friend inline std::enable_if_t<is_scalar_v<Scalar>, this_type>
    operator*(const this_type& array, const Scalar& val) {
        return val * array;
    }

    /**
     * @brief Divide scalar by array element-wise
     * @tparam Scalar Arithmetic type
     * @param val Scalar numerator
     * @param array Array of denominators
     * @return New array where each element is val / array[i]
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    template <typename Scalar>
    friend inline std::enable_if_t<is_scalar_v<Scalar>, this_type>
    operator/(const Scalar& val, const this_type& array) {
        this_type result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());
        JIARRAY_SIMD_LOOP
        for (int i = 0; i < result.nn; ++i) {
            result.mm[i] = val / array.mm[i];
        }
        return result;
    }

    /**
     * @brief Divide array by scalar element-wise
     * @tparam Scalar Arithmetic type
     * @param array Array numerator
     * @param val Scalar denominator
     * @return New array with each element divided by scalar
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    template <typename Scalar>
    friend inline std::enable_if_t<is_scalar_v<Scalar>, this_type>
    operator/(const this_type& array, const Scalar& val) {
        this_type result;
        result.initByRankSize(array.getSize(), array.getRankSize(), array.getOffset());

        if constexpr (std::is_floating_point_v<Scalar>) {
            const Scalar rval = Scalar{1} / val;
            JIARRAY_SIMD_LOOP
            for (int i = 0; i < result.nn; ++i) {
                result.mm[i] = array.mm[i] * rval;
            }
        } else {
            JIARRAY_SIMD_LOOP
            for (int i = 0; i < result.nn; ++i) {
                result.mm[i] = array.mm[i] / val;
            }
        }
        return result;
    }

    // ========================================================================
    // Additional mathematical operations - with conditional SIMD
    // ========================================================================

    /**
     * @brief Calculate sum of squared elements
     * @return Sum of squares
     * @note Returns double for numerical stability
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    inline double sqsum() const {
        double result = 0;
        // JIARRAY_SIMD_REDUCTION("+", result)
        for (int i = 0; i < nn; ++i) {
            result += static_cast<double>(mm[i]) * static_cast<double>(mm[i]);
        }
        return result;
    }

    /**
     * @brief Compute dot product of two arrays
     * @param array1 First array
     * @param array2 Second array
     * @return Dot product value
     * @note Arrays must have same dimensions and offsets
     * @note SIMD can be controlled via JIARRAY_USE_SIMD macro
     */
    friend inline T dot(const this_type& array1, const this_type& array2) {
        JIARRAY_CHECK_SIZE(array1.nn, array2.nn);
        for (int rank = 0; rank < RANK; ++rank) {
            assert(array1.sizes[rank] == array2.sizes[rank]);
            assert(array1.offset[rank] == array2.offset[rank]);
        }

        T result = T{};
        // JIARRAY_SIMD_REDUCTION("+", result)
        for (int i = 0; i < array1.nn; ++i) {
            result += array1.mm[i] * array2.mm[i];
        }
        return result;
    }

    // ========================================================================
    // Search operations
    // ========================================================================

    /**
     * @brief Check if array contains specified value
     * @param item Value to search for
     * @return true if item is found
     */
    inline bool contains(const T& item) const {
        return std::find(mm, mm + nn, item) != (mm + nn);
    }

    /**
     * @brief Find first occurrence of value
     * @param item Value to search for
     * @return For 1D: linear index; For ND: FastArray with multidimensional indices
     * @note Returns offset-adjusted indices, or -1+OFFSET if not found
     */
    inline auto findFirst(const T& item) const {
        if constexpr (RANK == 1) {
            auto it = std::find(mm, mm + nn, item);
            if (it == mm + nn)
                return (-1 + JIARRAY_OFFSET);
            return static_cast<int>(std::distance(mm, it)) + offset[0];
        } else {
            auto it = std::find(mm, mm + nn, item);

            FastArray<int, RANK> location;
            location = (-1 + JIARRAY_OFFSET);

            if (it == mm + nn)
                return location;

            int loc = static_cast<int>(std::distance(mm, it));

            // Convert linear index to multidimensional indices
            for (int rank = RANK - 1; rank > 0; --rank) {
                location.mm[rank] = loc / rankSize[rank];
                loc -= (location.mm[rank]) * rankSize[rank];
            }
            location.mm[0] = loc;

            // Apply offsets
            for (int rank = 0; rank < RANK; ++rank) {
                location.mm[rank] += offset[rank];
            }

            return location;
        }
    }

    /**
     * @brief Find location of maximum element
     * @return FastArray containing multidimensional indices of maximum element
     */
    inline FastArray<int, RANK> maxloc() const {
        assert(nn > 0);
        auto it     = std::max_element(mm, mm + nn);
        int  maxloc = static_cast<int>(std::distance(mm, it));

        FastArray<int, RANK> maxLocation;

        // Convert linear index to multidimensional indices
        for (int rank = RANK - 1; rank > 0; --rank) {
            maxLocation.mm[rank] = maxloc / rankSize[rank];
            maxloc -= (maxLocation.mm[rank]) * rankSize[rank];
        }
        maxLocation.mm[0] = maxloc;

        // Apply offsets
        for (int rank = 0; rank < RANK; ++rank) {
            maxLocation.mm[rank] += offset[rank];
        }

        return maxLocation;
    }

    /**
     * @brief Find location of maximum element within range
     * @param from Start index (inclusive)
     * @param to End index (exclusive or inclusive depending on JIARRAY_OFFSET)
     * @return FastArray containing multidimensional indices of maximum element
     */
    inline FastArray<int, RANK> maxloc(int from, int to) const {
        assert(from >= JIARRAY_OFFSET && to <= nn + JIARRAY_OFFSET);
        auto it     = std::max_element(mm + from - JIARRAY_OFFSET, mm + to - JIARRAY_OFFSET);
        int  maxloc = static_cast<int>(std::distance(mm, it));

        FastArray<int, RANK> maxLocation;

        for (int rank = RANK - 1; rank > 0; --rank) {
            maxLocation.mm[rank] = maxloc / rankSize[rank];
            maxloc -= maxLocation.mm[rank] * rankSize[rank];
        }
        maxLocation.mm[0] = maxloc;

        for (int rank = 0; rank < RANK; ++rank) {
            maxLocation.mm[rank] += offset[rank];
        }

        return maxLocation;
    }

    // ========================================================================
    // Utility methods
    // ========================================================================

    /**
     * @brief Create a deep copy of the array
     * @return New array with copied data
     */
    inline this_type copy() const {
        this_type array;
        array = *this;
        return array;
    }

    /**
     * @brief Share memory with another array (create view)
     * @param array Source array to share memory with
     * @note This array must not be allocated before calling shareWith
     */
    inline void shareWith(const this_type& array) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        nn = array.nn;
        mm = array.mm;
        std::copy(array.rankSize, array.rankSize + RANK, rankSize);
        std::copy(array.offset, array.offset + RANK, offset);
        allocated   = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
        std::copy(array.sizes, array.sizes + RANK, sizes);
    }

    /**
     * @brief Convert array to std::vector
     * @return Vector containing copy of all elements
     */
    inline std::vector<T> convertToVector() const {
        return std::vector<T>(mm, mm + nn);
    }

    // ========================================================================
    // Iterator implementation
    // ========================================================================

    // ========================================================================
    // Iterator implementation
    // ========================================================================

    /**
     * @class Iterator
     * @brief Random-access iterator for JIArray (non-const version)
     */
    class Iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        constexpr explicit Iterator(pointer ptr) noexcept : m_ptr(ptr) {
        }

        // Arithmetic operators
        constexpr Iterator operator+(difference_type n) const noexcept {
            return Iterator(m_ptr + n);
        }
        constexpr Iterator operator-(difference_type n) const noexcept {
            return Iterator(m_ptr - n);
        }

        constexpr Iterator& operator+=(difference_type n) noexcept {
            m_ptr += n;
            return *this;
        }
        constexpr Iterator& operator-=(difference_type n) noexcept {
            m_ptr -= n;
            return *this;
        }

        constexpr reference operator[](difference_type n) noexcept {
            return *(m_ptr + n);
        }

        constexpr difference_type operator-(const Iterator& other) const noexcept {
            return m_ptr - other.m_ptr;
        }

        // Comparison operators
        constexpr bool operator<(const Iterator& other) const noexcept {
            return m_ptr < other.m_ptr;
        }
        constexpr bool operator<=(const Iterator& other) const noexcept {
            return m_ptr <= other.m_ptr;
        }
        constexpr bool operator>(const Iterator& other) const noexcept {
            return m_ptr > other.m_ptr;
        }
        constexpr bool operator>=(const Iterator& other) const noexcept {
            return m_ptr >= other.m_ptr;
        }

        // Dereference operators
        constexpr reference operator*() const noexcept {
            return *m_ptr;
        }
        constexpr pointer operator->() const noexcept {
            return m_ptr;
        }

        // Increment/decrement operators
        constexpr Iterator& operator++() noexcept {
            ++m_ptr;
            return *this;
        }
        constexpr Iterator operator++(int) noexcept {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }
        constexpr Iterator& operator--() noexcept {
            --m_ptr;
            return *this;
        }
        constexpr Iterator operator--(int) noexcept {
            Iterator tmp = *this;
            --(*this);
            return tmp;
        }

        // Equality operators
        friend constexpr bool operator==(const Iterator& a, const Iterator& b) noexcept {
            return a.m_ptr == b.m_ptr;
        }
        friend constexpr bool operator!=(const Iterator& a, const Iterator& b) noexcept {
            return a.m_ptr != b.m_ptr;
        }

    private:
        pointer m_ptr;
    };

    /**
     * @class ConstIterator
     * @brief Random-access const iterator for JIArray
     */
    class ConstIterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = const T*;
        using reference         = const T&;

        constexpr explicit ConstIterator(pointer ptr) noexcept : m_ptr(ptr) {
        }

        // Allow implicit conversion from Iterator to ConstIterator
        constexpr ConstIterator(const Iterator& it) noexcept : m_ptr(&(*it)) {
        }

        // Arithmetic operators
        constexpr ConstIterator operator+(difference_type n) const noexcept {
            return ConstIterator(m_ptr + n);
        }
        constexpr ConstIterator operator-(difference_type n) const noexcept {
            return ConstIterator(m_ptr - n);
        }

        constexpr ConstIterator& operator+=(difference_type n) noexcept {
            m_ptr += n;
            return *this;
        }
        constexpr ConstIterator& operator-=(difference_type n) noexcept {
            m_ptr -= n;
            return *this;
        }

        constexpr reference operator[](difference_type n) const noexcept {
            return *(m_ptr + n);
        }

        constexpr difference_type operator-(const ConstIterator& other) const noexcept {
            return m_ptr - other.m_ptr;
        }

        // Comparison operators
        constexpr bool operator<(const ConstIterator& other) const noexcept {
            return m_ptr < other.m_ptr;
        }
        constexpr bool operator<=(const ConstIterator& other) const noexcept {
            return m_ptr <= other.m_ptr;
        }
        constexpr bool operator>(const ConstIterator& other) const noexcept {
            return m_ptr > other.m_ptr;
        }
        constexpr bool operator>=(const ConstIterator& other) const noexcept {
            return m_ptr >= other.m_ptr;
        }

        // Dereference operators - const만!
        constexpr reference operator*() const noexcept {
            return *m_ptr;
        }
        constexpr pointer operator->() const noexcept {
            return m_ptr;
        }

        // Increment/decrement operators
        constexpr ConstIterator& operator++() noexcept {
            ++m_ptr;
            return *this;
        }
        constexpr ConstIterator operator++(int) noexcept {
            ConstIterator tmp = *this;
            ++(*this);
            return tmp;
        }
        constexpr ConstIterator& operator--() noexcept {
            --m_ptr;
            return *this;
        }
        constexpr ConstIterator operator--(int) noexcept {
            ConstIterator tmp = *this;
            --(*this);
            return tmp;
        }

        // Equality operators
        friend constexpr bool operator==(const ConstIterator& a, const ConstIterator& b) noexcept {
            return a.m_ptr == b.m_ptr;
        }
        friend constexpr bool operator!=(const ConstIterator& a, const ConstIterator& b) noexcept {
            return a.m_ptr != b.m_ptr;
        }

    private:
        pointer m_ptr;
    };

    /// @name Iterator access
    /// @{
    Iterator begin() noexcept {
        return Iterator(mm);
    }
    Iterator end() noexcept {
        return Iterator(mm + nn);
    }
    ConstIterator begin() const noexcept {
        return ConstIterator(mm);
    }
    ConstIterator end() const noexcept {
        return ConstIterator(mm + nn);
    }
    ConstIterator cbegin() const noexcept {
        return ConstIterator(mm);
    }
    ConstIterator cend() const noexcept {
        return ConstIterator(mm + nn);
    }
    /// @}
};

// ============================================================================
// Type aliases for convenience
// ============================================================================

/// @name Boolean array type aliases
/// @{
#define zbool1 JIArray<bool, 1>
#define zbool2 JIArray<bool, 2>
#define zbool3 JIArray<bool, 3>
#define zbool4 JIArray<bool, 4>
#define zbool5 JIArray<bool, 5>
/// @}

/// @name Integer array type aliases
/// @{
#define zint1 JIArray<int, 1>
#define zint2 JIArray<int, 2>
#define zint3 JIArray<int, 3>
#define zint4 JIArray<int, 4>
#define zint5 JIArray<int, 5>
/// @}

/// @name Double-precision floating-point array type aliases
/// @{
#define zdouble1 JIArray<double, 1>
#define zdouble2 JIArray<double, 2>
#define zdouble3 JIArray<double, 3>
#define zdouble4 JIArray<double, 4>
#define zdouble5 JIArray<double, 5>
#define zdouble6 JIArray<double, 6>
/// @}

/// @name Single-precision floating-point array type aliases
/// @{
#define zfloat1 JIArray<float, 1>
#define zfloat2 JIArray<float, 2>
#define zfloat3 JIArray<float, 3>
#define zfloat4 JIArray<float, 4>
#define zfloat5 JIArray<float, 5>
/// @}

/// @name String array type aliases
/// @{
#define zstring1 JIArray<std::string, 1>
#define zstring2 JIArray<std::string, 2>
#define zstring3 JIArray<std::string, 3>
#define zstring4 JIArray<std::string, 4>
#define zstring5 JIArray<std::string, 5>
/// @}

// ============================================================================
// Loop macros
// ============================================================================

/**
 * @def ffor
 * @brief Flexible for-loop macro that adapts to indexing convention
 * @param i Loop variable
 * @param begin Start index
 * @param end End index (inclusive for 1-based, exclusive for 0-based)
 */
#if JIARRAY_OFFSET == 0
    #define ffor(i, begin, end)      for (int i = begin; i < end; ++i)
    #define ffor_back(i, begin, end) for (int i = begin; i >= end; --i)
#else
    #define ffor(i, begin, end)      for (int i = begin; i <= end; ++i)
    #define ffor_back(i, begin, end) for (int i = begin; i >= end; --i)
#endif

/**
 * @def zfor
 * @brief Loop from offset to end (offset-aware)
 * @param i Loop variable
 * @param end End index
 */
#define zfor(i, end) ffor(i, JIARRAY_OFFSET, end)

/**
 * @brief Generic multidimensional array type alias
 * @tparam Type Element type
 * @tparam N Rank (default: 1)
 */
template <typename Type, int N = 1>
using zarray = JIArray<Type, N, std::make_index_sequence<N>>;

} // namespace dnegri::jiarray
