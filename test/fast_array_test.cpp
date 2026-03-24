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
    // 2x3 column-major: mm stores [col1_row1, col1_row2, col2_row1, col2_row2, col3_row1, col3_row2]
    FastArray2D<int, 2, 3> arr({1, 2, 3, 4, 5, 6});
    // Column-major layout: (i,j) maps to mm[(j-1)*SIZE1 + (i-1)]
    EXPECT_EQ(arr(1, 1), 1); // mm[0]
    EXPECT_EQ(arr(2, 1), 2); // mm[1]
    EXPECT_EQ(arr(1, 2), 3); // mm[2]
    EXPECT_EQ(arr(2, 2), 4); // mm[3]
    EXPECT_EQ(arr(1, 3), 5); // mm[4]
    EXPECT_EQ(arr(2, 3), 6); // mm[5]
}

TEST(FastArray2DTest, CArrayConstructor) {
    int raw[] = {10, 20, 30, 40};
    FastArray2D<int, 2, 2> arr(raw);
    EXPECT_EQ(arr(1, 1), 10);
    EXPECT_EQ(arr(2, 1), 20);
    EXPECT_EQ(arr(1, 2), 30);
    EXPECT_EQ(arr(2, 2), 40);
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
