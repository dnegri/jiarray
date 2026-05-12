#include <jiarray/JIArray.h>
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <string>

using namespace dnegri::jiarray;

// ============================================================================
// FastArray tests
// ============================================================================

TEST(FastArrayTest, DefaultConstructorZeroInitialized) {
    FastArray<int, 5> arr;
    for (int i = 1; i <= 5; ++i) {
        EXPECT_EQ(arr(i), 0);
    }
}

TEST(FastArrayTest, ScalarFillConstructor) {
    FastArray<double, 4> arr(3.14);
    for (int i = 1; i <= 4; ++i) {
        EXPECT_DOUBLE_EQ(arr(i), 3.14);
    }
}

TEST(FastArrayTest, InitializerListConstructor) {
    FastArray<int, 4> arr({10, 20, 30, 40});
    EXPECT_EQ(arr(1), 10);
    EXPECT_EQ(arr(2), 20);
    EXPECT_EQ(arr(3), 30);
    EXPECT_EQ(arr(4), 40);
}

TEST(FastArrayTest, CopyConstructor) {
    FastArray<int, 3> src({5, 10, 15});
    FastArray<int, 3> copy(src);
    EXPECT_EQ(copy(1), 5);
    EXPECT_EQ(copy(2), 10);
    EXPECT_EQ(copy(3), 15);
    // Verify independent - modify source shouldn't affect copy
    src(1) = 999;
    EXPECT_EQ(copy(1), 5);
}

TEST(FastArrayTest, CArrayConstructor) {
    int raw[] = {7, 8, 9};
    FastArray<int, 3> arr(raw);
    EXPECT_EQ(arr(1), 7);
    EXPECT_EQ(arr(2), 8);
    EXPECT_EQ(arr(3), 9);
}

TEST(FastArrayTest, OneBasedAccessOperator) {
    FastArray<int, 3> arr({100, 200, 300});
    // 1-based indexing via operator()
    EXPECT_EQ(arr(1), 100);
    EXPECT_EQ(arr(3), 300);
    arr(2) = 999;
    EXPECT_EQ(arr(2), 999);
}

TEST(FastArrayTest, ConstAccessOperator) {
    const FastArray<int, 3> arr({100, 200, 300});
    EXPECT_EQ(arr(1), 100);
    EXPECT_EQ(arr(3), 300);
}

TEST(FastArrayTest, ZeroBasedBracketOperator) {
    FastArray<int, 3> arr({10, 20, 30});
    // 0-based indexing via operator[]
    EXPECT_EQ(arr[0], 10);
    EXPECT_EQ(arr[1], 20);
    EXPECT_EQ(arr[2], 30);
}

TEST(FastArrayTest, FindFirstFound) {
    FastArray<int, 5> arr({10, 20, 30, 20, 50});
    // findFirst returns 1-based index of first match
    EXPECT_EQ(arr.findFirst(30), 3);
    EXPECT_EQ(arr.findFirst(20), 2); // first occurrence
    EXPECT_EQ(arr.findFirst(10), 1);
    EXPECT_EQ(arr.findFirst(50), 5);
}

TEST(FastArrayTest, FindFirstNotFound) {
    FastArray<int, 3> arr({1, 2, 3});
    // not found returns OFFSET - 1 = 0
    EXPECT_EQ(arr.findFirst(99), 0);
}

TEST(FastArrayTest, MinMax) {
    FastArray<double, 5> arr({3.0, 1.0, 4.0, 1.5, 2.0});
    EXPECT_DOUBLE_EQ(arr.min(), 1.0);
    EXPECT_DOUBLE_EQ(arr.max(), 4.0);
}

TEST(FastArrayTest, MinMaxSingleElement) {
    FastArray<int, 1> arr({42});
    EXPECT_EQ(arr.min(), 42);
    EXPECT_EQ(arr.max(), 42);
}

TEST(FastArrayTest, AssignmentFromArray) {
    FastArray<int, 3> src({1, 2, 3});
    FastArray<int, 3> dst;
    dst = src;
    EXPECT_EQ(dst(1), 1);
    EXPECT_EQ(dst(2), 2);
    EXPECT_EQ(dst(3), 3);
}

TEST(FastArrayTest, AssignmentFromPointer) {
    int vals[] = {7, 8, 9};
    FastArray<int, 3> arr;
    arr = static_cast<const int*>(vals);
    EXPECT_EQ(arr(1), 7);
    EXPECT_EQ(arr(2), 8);
    EXPECT_EQ(arr(3), 9);
}

TEST(FastArrayTest, AssignmentFromScalar) {
    FastArray<int, 4> arr;
    arr = 42;
    for (int i = 1; i <= 4; ++i) {
        EXPECT_EQ(arr(i), 42);
    }
}

TEST(FastArrayTest, PlusEqualsScalar) {
    FastArray<int, 3> arr({10, 20, 30});
    arr += 5;
    EXPECT_EQ(arr(1), 15);
    EXPECT_EQ(arr(2), 25);
    EXPECT_EQ(arr(3), 35);
}

TEST(FastArrayTest, MultiplyEqualsScalar) {
    FastArray<double, 3> arr({2.0, 3.0, 4.0});
    arr *= 10.0;
    EXPECT_DOUBLE_EQ(arr(1), 20.0);
    EXPECT_DOUBLE_EQ(arr(2), 30.0);
    EXPECT_DOUBLE_EQ(arr(3), 40.0);
}

TEST(FastArrayTest, DivideEqualsScalar) {
    FastArray<double, 3> arr({10.0, 20.0, 30.0});
    arr /= 5.0;
    EXPECT_DOUBLE_EQ(arr(1), 2.0);
    EXPECT_DOUBLE_EQ(arr(2), 4.0);
    EXPECT_DOUBLE_EQ(arr(3), 6.0);
}

TEST(FastArrayTest, MultiplyOperator) {
    FastArray<int, 3> arr({2, 3, 4});
    auto result = arr * 10;
    EXPECT_EQ(result(1), 20);
    EXPECT_EQ(result(2), 30);
    EXPECT_EQ(result(3), 40);
    // Original unchanged
    EXPECT_EQ(arr(1), 2);
}

// ============================================================================
// FastArray Iterator tests
// ============================================================================

TEST(FastArrayTest, ForwardIterator) {
    FastArray<int, 4> arr({10, 20, 30, 40});
    std::vector<int> collected(arr.begin(), arr.end());
    EXPECT_EQ(collected, (std::vector<int>{10, 20, 30, 40}));
}

TEST(FastArrayTest, ConstIterator) {
    const FastArray<int, 3> arr({5, 10, 15});
    int sum = 0;
    for (auto it = arr.cbegin(); it != arr.cend(); ++it) {
        sum += *it;
    }
    EXPECT_EQ(sum, 30);
}

TEST(FastArrayTest, ReverseIterator) {
    FastArray<int, 4> arr({1, 2, 3, 4});
    std::vector<int> reversed(arr.rbegin(), arr.rend());
    EXPECT_EQ(reversed, (std::vector<int>{4, 3, 2, 1}));
}

TEST(FastArrayTest, ConstReverseIterator) {
    const FastArray<int, 3> arr({10, 20, 30});
    std::vector<int> reversed(arr.crbegin(), arr.crend());
    EXPECT_EQ(reversed, (std::vector<int>{30, 20, 10}));
}

TEST(FastArrayTest, RangeBasedFor) {
    FastArray<int, 4> arr({1, 2, 3, 4});
    int sum = 0;
    for (auto& v : arr) {
        sum += v;
    }
    EXPECT_EQ(sum, 10);
}

TEST(FastArrayTest, STLAlgorithmSort) {
    FastArray<int, 5> arr({5, 3, 1, 4, 2});
    std::sort(arr.begin(), arr.end());
    EXPECT_EQ(arr[0], 1);
    EXPECT_EQ(arr[1], 2);
    EXPECT_EQ(arr[2], 3);
    EXPECT_EQ(arr[3], 4);
    EXPECT_EQ(arr[4], 5);
}

TEST(FastArrayTest, STLAccumulate) {
    FastArray<double, 4> arr({1.5, 2.5, 3.5, 4.5});
    double sum = std::accumulate(arr.begin(), arr.end(), 0.0);
    EXPECT_DOUBLE_EQ(sum, 12.0);
}

// ============================================================================
// FastArray STL-like interface
// ============================================================================

TEST(FastArrayTest, SizeAndEmpty) {
    FastArray<int, 5> arr;
    EXPECT_EQ(arr.size(), 5u);
    EXPECT_EQ(arr.max_size(), 5u);
    EXPECT_FALSE(arr.empty());
}

TEST(FastArrayTest, FrontAndBack) {
    FastArray<int, 4> arr({10, 20, 30, 40});
    EXPECT_EQ(arr.front(), 10);
    EXPECT_EQ(arr.back(), 40);
    arr.front() = 99;
    EXPECT_EQ(arr(1), 99);
}

TEST(FastArrayTest, ConstFrontAndBack) {
    const FastArray<int, 3> arr({5, 10, 15});
    EXPECT_EQ(arr.front(), 5);
    EXPECT_EQ(arr.back(), 15);
}

TEST(FastArrayTest, DataPointer) {
    FastArray<int, 3> arr({1, 2, 3});
    int* p = arr.data();
    EXPECT_EQ(p[0], 1);
    EXPECT_EQ(p[2], 3);
}

TEST(FastArrayTest, ConstDataPointer) {
    const FastArray<int, 3> arr({1, 2, 3});
    const int* p = arr.data();
    EXPECT_EQ(p[0], 1);
    EXPECT_EQ(p[2], 3);
}

// ============================================================================
// FastArray2D tests
// ============================================================================

TEST(FastArray2DTest, DefaultConstructorZeroInit) {
    FastArray2D<int, 3, 4> arr;
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 4; ++j)
            EXPECT_EQ(arr(i, j), 0);
}

TEST(FastArray2DTest, ScalarFillConstructor) {
    FastArray2D<double, 2, 3> arr(5.5);
    for (int i = 1; i <= 2; ++i)
        for (int j = 1; j <= 3; ++j)
            EXPECT_DOUBLE_EQ(arr(i, j), 5.5);
}

TEST(FastArray2DTest, InitializerListConstructor) {
    // Initializer list fills mm[] in declared order.  How (i,j) maps to mm
    // depends on storage order, so check whichever layout was compiled in.
    FastArray2D<int, 2, 3> arr({1, 2, 3, 4, 5, 6});
    if (arr.is_row_major) {
        // row-major: (i,j) → mm[(i-1)*3 + (j-1)]
        EXPECT_EQ(arr(1, 1), 1); EXPECT_EQ(arr(1, 2), 2); EXPECT_EQ(arr(1, 3), 3);
        EXPECT_EQ(arr(2, 1), 4); EXPECT_EQ(arr(2, 2), 5); EXPECT_EQ(arr(2, 3), 6);
    } else {
        // column-major: (i,j) → mm[(j-1)*2 + (i-1)]
        EXPECT_EQ(arr(1, 1), 1); EXPECT_EQ(arr(2, 1), 2);
        EXPECT_EQ(arr(1, 2), 3); EXPECT_EQ(arr(2, 2), 4);
        EXPECT_EQ(arr(1, 3), 5); EXPECT_EQ(arr(2, 3), 6);
    }
}

TEST(FastArray2DTest, CArrayConstructor) {
    int raw[] = {10, 20, 30, 40};
    FastArray2D<int, 2, 2> arr(raw);
    if (arr.is_row_major) {
        EXPECT_EQ(arr(1, 1), 10); EXPECT_EQ(arr(1, 2), 20);
        EXPECT_EQ(arr(2, 1), 30); EXPECT_EQ(arr(2, 2), 40);
    } else {
        EXPECT_EQ(arr(1, 1), 10); EXPECT_EQ(arr(2, 1), 20);
        EXPECT_EQ(arr(1, 2), 30); EXPECT_EQ(arr(2, 2), 40);
    }
}

TEST(FastArray2DTest, ElementAccessAndModification) {
    FastArray2D<int, 3, 3> arr(0);
    arr(2, 3) = 42;
    EXPECT_EQ(arr(2, 3), 42);
    EXPECT_EQ(arr(1, 1), 0);
}

TEST(FastArray2DTest, ConstAccess) {
    const FastArray2D<int, 2, 2> arr({1, 2, 3, 4});
    EXPECT_EQ(arr(1, 1), 1);
    EXPECT_EQ(arr(2, 2), 4);
}

TEST(FastArray2DTest, AssignmentFromScalar) {
    FastArray2D<int, 2, 3> arr;
    arr = 7;
    for (int i = 1; i <= 2; ++i)
        for (int j = 1; j <= 3; ++j)
            EXPECT_EQ(arr(i, j), 7);
}

TEST(FastArray2DTest, CopyAssignment) {
    FastArray2D<int, 2, 2> src({1, 2, 3, 4});
    FastArray2D<int, 2, 2> dst;
    dst = src;
    EXPECT_EQ(dst(1, 1), 1);
    EXPECT_EQ(dst(2, 2), 4);
    // Verify independent copy
    src(1, 1) = 999;
    EXPECT_EQ(dst(1, 1), 1);
}

// ============================================================================
// StringFastArray tests
// ============================================================================

TEST(StringFastArrayTest, DefaultConstructorEmpty) {
    StringFastArray<3> arr;
    for (int i = 1; i <= 3; ++i) {
        EXPECT_TRUE(arr(i).empty());
    }
}

TEST(StringFastArrayTest, InitializerListConstructor) {
    StringFastArray<3> arr({"hello", "world", "test"});
    EXPECT_EQ(arr(1), "hello");
    EXPECT_EQ(arr(2), "world");
    EXPECT_EQ(arr(3), "test");
}

TEST(StringFastArrayTest, ConstAccess) {
    const StringFastArray<2> arr({"foo", "bar"});
    EXPECT_EQ(arr(1), "foo");
    EXPECT_EQ(arr(2), "bar");
}

TEST(StringFastArrayTest, Modification) {
    StringFastArray<2> arr({"a", "b"});
    arr(1) = "modified";
    EXPECT_EQ(arr(1), "modified");
    EXPECT_EQ(arr(2), "b");
}

TEST(StringFastArrayTest, CopyAssignment) {
    StringFastArray<2> src({"x", "y"});
    StringFastArray<2> dst;
    dst = src;
    EXPECT_EQ(dst(1), "x");
    EXPECT_EQ(dst(2), "y");
}

TEST(StringFastArrayTest, ScalarAssignment) {
    StringFastArray<3> arr;
    arr = std::string("same");
    for (int i = 1; i <= 3; ++i) {
        EXPECT_EQ(arr(i), "same");
    }
}

// ============================================================================
// Type alias tests
// ============================================================================

TEST(TypeAliasTest, IntAliases) {
    fint1d<3> iarr({1, 2, 3});
    EXPECT_EQ(iarr(1), 1);
    EXPECT_EQ(iarr(3), 3);
}

TEST(TypeAliasTest, DoubleAliases) {
    fdouble1d<3> darr({1.1, 2.2, 3.3});
    EXPECT_DOUBLE_EQ(darr(1), 1.1);
    EXPECT_DOUBLE_EQ(darr(3), 3.3);
}

TEST(TypeAliasTest, BoolAliases) {
    fbool1d<3> barr({true, false, true});
    EXPECT_TRUE(barr(1));
    EXPECT_FALSE(barr(2));
    EXPECT_TRUE(barr(3));
}

TEST(TypeAliasTest, FloatAliases) {
    ffloat1d<3> farr({1.0f, 2.0f, 3.0f});
    EXPECT_FLOAT_EQ(farr(1), 1.0f);
    EXPECT_FLOAT_EQ(farr(3), 3.0f);
}

TEST(TypeAliasTest, StringAliases) {
    fstring1d<2> sarr({"ab", "cd"});
    EXPECT_EQ(sarr(1), "ab");
    EXPECT_EQ(sarr(2), "cd");
}

// ============================================================================
// Bounds checking (JIARRAY_DEBUG active in test builds)
// ============================================================================

TEST(FastArrayTest, OutOfBoundsThrows) {
    FastArray<int, 3> arr({1, 2, 3});
    EXPECT_THROW(arr(0), std::out_of_range);  // below lower bound
    EXPECT_THROW(arr(4), std::out_of_range);  // above upper bound
}

TEST(FastArrayTest, ConstOutOfBoundsThrows) {
    const FastArray<int, 3> arr({1, 2, 3});
    EXPECT_THROW(arr(0), std::out_of_range);
    EXPECT_THROW(arr(4), std::out_of_range);
}

TEST(FastArray2DTest, OutOfBoundsThrows) {
    FastArray2D<int, 2, 3> arr(0);
    EXPECT_THROW(arr(0, 1), std::out_of_range);
    EXPECT_THROW(arr(1, 0), std::out_of_range);
    EXPECT_THROW(arr(3, 1), std::out_of_range);
    EXPECT_THROW(arr(1, 4), std::out_of_range);
}

TEST(FastArray2DTest, ConstOutOfBoundsThrows) {
    const FastArray2D<int, 2, 3> arr(0);
    EXPECT_THROW(arr(0, 1), std::out_of_range);
    EXPECT_THROW(arr(1, 4), std::out_of_range);
}

TEST(StringFastArrayTest, OutOfBoundsThrows) {
    StringFastArray<3> arr({"a", "b", "c"});
    EXPECT_THROW(arr(0), std::out_of_range);
    EXPECT_THROW(arr(4), std::out_of_range);
}

TEST(StringFastArrayTest, ConstOutOfBoundsThrows) {
    const StringFastArray<3> arr({"a", "b", "c"});
    EXPECT_THROW(arr(0), std::out_of_range);
    EXPECT_THROW(arr(4), std::out_of_range);
}

// ============================================================================
// Type alias tests
// ============================================================================

TEST(TypeAliasTest, Int2DAliases) {
    fint2d<2, 3> arr(0);
    arr(1, 2) = 42;
    EXPECT_EQ(arr(1, 2), 42);
}

TEST(TypeAliasTest, Double2DAliases) {
    fdouble2d<2, 2> arr({1.0, 2.0, 3.0, 4.0});
    EXPECT_DOUBLE_EQ(arr(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(arr(2, 2), 4.0);
}

// ============================================================================
// 0.8.0+ unified API: variadic aliases (fint<...>, fdouble<...>, ...),
// equality operators, 3D rank, and storage-layout verification.
// ============================================================================

TEST(VariadicAliasTest, OneDimensionalFint) {
    fint<5> v({10, 20, 30, 40, 50});
    EXPECT_EQ(v.RANK, 1u);
    EXPECT_EQ(v.SIZE, 5u);
    for (int i = 1; i <= 5; ++i) EXPECT_EQ(v(i), i * 10);
}

TEST(VariadicAliasTest, TwoDimensionalFint) {
    fint<3, 2> m(0);
    EXPECT_EQ(m.RANK, 2u);
    EXPECT_EQ(m.SIZE, 6u);
    m(1, 1) = 11;
    m(3, 2) = 32;
    EXPECT_EQ(m(1, 1), 11);
    EXPECT_EQ(m(3, 2), 32);
}

TEST(VariadicAliasTest, ThreeDimensionalFdouble) {
    fdouble<2, 3, 4> cube(0.0);
    EXPECT_EQ(cube.RANK, 3u);
    EXPECT_EQ(cube.SIZE, 24u);
    cube(2, 3, 4) = 99.5;
    EXPECT_DOUBLE_EQ(cube(2, 3, 4), 99.5);
}

TEST(VariadicAliasTest, FstringAndFarray) {
    fstring<3> names({"a", "b", "c"});
    EXPECT_EQ(names(1), "a");
    EXPECT_EQ(names(3), "c");

    farray<unsigned, 4> u({1u, 2u, 3u, 4u});
    EXPECT_EQ(u(2), 2u);
}

// ----------------- Equality -----------------

TEST(EqualityTest, FastArray1DEqual) {
    fint<4> a({1, 2, 3, 4});
    fint<4> b({1, 2, 3, 4});
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);
}

TEST(EqualityTest, FastArray1DDifferentMiddle) {
    fint<4> a({1, 2, 3, 4});
    fint<4> b({1, 9, 3, 4});
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a != b);
}

TEST(EqualityTest, FastArray2DEqual) {
    fint<2, 3> a(0);
    fint<2, 3> b(0);
    a(1, 1) = 7; b(1, 1) = 7;
    a(2, 3) = 8; b(2, 3) = 8;
    EXPECT_TRUE(a == b);

    b(2, 3) = 0;
    EXPECT_FALSE(a == b);
}

TEST(EqualityTest, StringFastArrayEqual) {
    fstring<3> a({"x", "y", "z"});
    fstring<3> b({"x", "y", "z"});
    fstring<3> c({"x", "y", "Z"});
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

// ----------------- 3D access -----------------

TEST(VariadicAliasTest, ThreeDimensionalLinearLayout) {
    // Write every cell with a unique value and read back to confirm
    // operator() and operator[] are consistent.
    fint<2, 3, 4> a(0);
    int n = 0;
    for (int k = 1; k <= 4; ++k)
        for (int j = 1; j <= 3; ++j)
            for (int i = 1; i <= 2; ++i)
                a(i, j, k) = ++n;

    EXPECT_EQ(a.SIZE, 24u);
    // Flat-iterate via operator[] and confirm we hit every value once.
    long sum = 0;
    for (std::size_t flat = 0; flat < a.SIZE; ++flat) sum += a[flat];
    EXPECT_EQ(sum, 24L * 25L / 2L);  // 1+2+...+24
}

// ----------------- Storage order (compile-time check) -----------------

TEST(StorageOrderTest, MatchesJiarrayMacro) {
    // Confirm the static flag is consistent with the global macro.
    constexpr bool expected_row_major = (JIARRAY_COLUMN_MAJOR == 0);
    constexpr bool actual             = fint<2, 2>::is_row_major;
    EXPECT_EQ(actual, expected_row_major);
}

TEST(StorageOrderTest, ColumnMajorStrideOrder) {
    // For a (3 rows × 2 cols) array under column-major, incrementing the
    // *first* index moves by 1 in flat memory; incrementing the second
    // moves by 3.  Under row-major the roles swap.
    fint<3, 2> a(0);
    a(1, 1) = 100;
    a(2, 1) = 200;
    a(1, 2) = 300;

    if (a.is_row_major) {
        // row-major: a(i,j) at flat = (i-OFFSET) * 2 + (j-OFFSET)
        EXPECT_EQ(a[0], 100);  // (1,1)
        EXPECT_EQ(a[1], 300);  // (1,2)
        EXPECT_EQ(a[2], 200);  // (2,1)
    } else {
        // column-major: a(i,j) at flat = (i-OFFSET) + (j-OFFSET) * 3
        EXPECT_EQ(a[0], 100);  // (1,1)
        EXPECT_EQ(a[1], 200);  // (2,1)
        EXPECT_EQ(a[3], 300);  // (1,2)
    }
}

TEST(StorageOrderTest, ThreeDStrides) {
    // Spot-check that 3D access agrees with the layout-implied strides.
    fint<2, 3, 4> a(0);
    a(1, 1, 1) = 1;
    a(2, 1, 1) = 2;
    a(1, 2, 1) = 3;
    a(1, 1, 2) = 4;

    if (a.is_row_major) {
        // strides (row-major): [3*4, 4, 1]  →  i*12 + j*4 + k
        EXPECT_EQ(a[0], 1);   // (1,1,1)
        EXPECT_EQ(a[12], 2);  // (2,1,1) → stride 12
        EXPECT_EQ(a[4], 3);   // (1,2,1) → stride 4
        EXPECT_EQ(a[1], 4);   // (1,1,2) → stride 1
    } else {
        // strides (column-major): [1, 2, 6]  →  i + j*2 + k*6
        EXPECT_EQ(a[0], 1);   // (1,1,1)
        EXPECT_EQ(a[1], 2);   // (2,1,1) → stride 1
        EXPECT_EQ(a[2], 3);   // (1,2,1) → stride 2
        EXPECT_EQ(a[6], 4);   // (1,1,2) → stride 6
    }
}
