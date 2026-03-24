#include <jiarray/JIArray.h>
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace dnegri::jiarray;

// ============================================================================
// erase() - completely untested
// ============================================================================

TEST(JIArrayCoverageTest, EraseResetsCompletely) {
    zint2 arr(3, 4);
    arr = 5;
    EXPECT_EQ(arr.size(), 12);

    arr.erase();
    EXPECT_EQ(arr.size(), 0);
    EXPECT_EQ(arr.data(), nullptr);

    // Can re-init after erase
    arr.init(2, 3);
    EXPECT_EQ(arr.size(), 6);
    arr(1, 1) = 42;
    EXPECT_EQ(arr(1, 1), 42);
}

TEST(JIArrayCoverageTest, EraseOnAlreadyEmpty) {
    zint1 arr;
    arr.erase(); // should not crash
    EXPECT_EQ(arr.size(), 0);
}

// ============================================================================
// data(INDEX...) - partial index pointer access
// ============================================================================

TEST(JIArrayCoverageTest, DataPartialIndex2D) {
    zint2 arr(3, 4);
    // Fill with known pattern: arr(i,j) = i*10 + j
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 4; ++j)
            arr(i, j) = i * 10 + j;

    // data(1) should point to column-major column where last dim index = 1
    // In column-major 2D: data(i) gives pointer to the start of "slice at last dim = i"
    const int* p = arr.data(1);
    // For column-major: data(1) points to arr(:, 1)
    // p[0] = arr(1,1) = 11, p[1] = arr(2,1) = 21, p[2] = arr(3,1) = 31
    EXPECT_EQ(p[0], 11);
    EXPECT_EQ(p[1], 21);
    EXPECT_EQ(p[2], 31);
}

TEST(JIArrayCoverageTest, DataPartialIndex3D) {
    zint3 arr(2, 3, 4);
    // Fill with pattern
    for (int i = 1; i <= 2; ++i)
        for (int j = 1; j <= 3; ++j)
            for (int k = 1; k <= 4; ++k)
                arr(i, j, k) = i * 100 + j * 10 + k;

    // data(2) in column-major: gets a pointer for slice where last dim = 2
    const int* p = arr.data(2);
    // arr(:,:,2) starts at some offset - first element should be arr(1,1,2) = 112
    EXPECT_EQ(p[0], 112);
}

TEST(JIArrayCoverageTest, ConstDataPartialIndex) {
    zint2 arr(3, 3);
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
            arr(i, j) = i + j;

    const zint2& cref = arr;
    const int* p = cref.data(2);
    // Column 2: arr(1,2)=3, arr(2,2)=4, arr(3,2)=5
    EXPECT_EQ(p[0], 3);
    EXPECT_EQ(p[1], 4);
    EXPECT_EQ(p[2], 5);
}

// ============================================================================
// initByRankSize() - advanced initialization
// ============================================================================

TEST(JIArrayCoverageTest, InitByRankSizeWithAllocation) {
    zint1 arr;
    int rankSizes[] = {1};
    int offsets[] = {1};
    arr.initByRankSize(5, rankSizes, offsets);
    EXPECT_EQ(arr.size(), 5);
    arr(1) = 10;
    arr(5) = 50;
    EXPECT_EQ(arr(1), 10);
    EXPECT_EQ(arr(5), 50);
}

TEST(JIArrayCoverageTest, InitByRankSizeWithExternalMemory) {
    int buffer[6] = {10, 20, 30, 40, 50, 60};
    zint2 arr;
    // Column-major 2x3: rankSize = {1, 2}, offsets = {1, 1}
    int rankSizes[] = {1, 2};
    int offsets[] = {1, 1};
    arr.initByRankSize(6, rankSizes, offsets, buffer);
    // arr(1,1) = buffer[0] = 10, arr(2,1) = buffer[1] = 20
    // arr(1,2) = buffer[2] = 30, arr(2,2) = buffer[3] = 40
    EXPECT_EQ(arr(1, 1), 10);
    EXPECT_EQ(arr(2, 1), 20);
    EXPECT_EQ(arr(1, 2), 30);
    EXPECT_EQ(arr(2, 3), 60);
}

TEST(JIArrayCoverageTest, InitByRankSizeZeroSize) {
    zint1 arr;
    int rankSizes[] = {1};
    int offsets[] = {1};
    arr.initByRankSize(0, rankSizes, offsets);
    EXPECT_EQ(arr.size(), 0);
}

// ============================================================================
// Assignment from C-style array (operator=(const T*))
// ============================================================================

TEST(JIArrayCoverageTest, AssignFromCArray) {
    zint1 arr(5);
    int source[] = {10, 20, 30, 40, 50};
    arr = source;
    EXPECT_EQ(arr(1), 10);
    EXPECT_EQ(arr(3), 30);
    EXPECT_EQ(arr(5), 50);
}

TEST(JIArrayCoverageTest, AssignFromCArray2D) {
    zdouble2 arr(2, 3);
    double source[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    arr = source;
    // Column-major: first 2 elements form column 1
    EXPECT_DOUBLE_EQ(arr(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(arr(2, 1), 2.0);
    EXPECT_DOUBLE_EQ(arr(1, 2), 3.0);
}

// ============================================================================
// Assignment from FastArray
// ============================================================================

TEST(JIArrayCoverageTest, AssignFromFastArray) {
    FastArray<int, 4> fa({10, 20, 30, 40});
    zint1 arr;  // empty, should auto-init
    arr = fa;
    EXPECT_EQ(arr.size(), 4);
    EXPECT_EQ(arr(1), 10);
    EXPECT_EQ(arr(4), 40);
}

TEST(JIArrayCoverageTest, AssignFromFastArrayToExisting) {
    FastArray<int, 3> fa({7, 8, 9});
    zint1 arr(3);
    arr = fa;
    EXPECT_EQ(arr(1), 7);
    EXPECT_EQ(arr(2), 8);
    EXPECT_EQ(arr(3), 9);
}

// ============================================================================
// maxloc(from, to) - range-based maxloc
// ============================================================================

TEST(JIArrayCoverageTest, MaxlocRange1D) {
    zint1 arr{5, 30, 10, 25, 15};
    // Search in range [2, 4] (1-based) -> elements: 30, 10, 25
    auto loc = arr.maxloc(2, 4);
    EXPECT_EQ(loc(1), 2); // 30 is at index 2
}

TEST(JIArrayCoverageTest, MaxlocRangeFull) {
    zdouble1 arr(4);
    arr(1) = 1.0;
    arr(2) = 4.0;
    arr(3) = 2.0;
    arr(4) = 3.0;
    auto loc = arr.maxloc(1, 4);
    EXPECT_EQ(loc(1), 2); // max value 4.0 at position 2
}

// ============================================================================
// ConstIterator explicit usage
// ============================================================================

TEST(JIArrayCoverageTest, ConstIteratorFromConst) {
    const zint1 arr{10, 20, 30, 40};
    int sum = 0;
    for (auto it = arr.cbegin(); it != arr.cend(); ++it) {
        sum += *it;
    }
    EXPECT_EQ(sum, 100);
}

TEST(JIArrayCoverageTest, ConstIteratorArithmetic) {
    const zint1 arr{1, 2, 3, 4, 5};
    auto it = arr.cbegin();
    EXPECT_EQ(*(it + 2), 3);
    EXPECT_EQ(it[3], 4);

    auto it2 = arr.cend();
    EXPECT_EQ(it2 - it, 5);
}

TEST(JIArrayCoverageTest, ConstIteratorComparison) {
    const zint1 arr{1, 2, 3};
    auto first = arr.cbegin();
    auto last = arr.cend();
    EXPECT_TRUE(first < last);
    EXPECT_TRUE(first <= last);
    EXPECT_TRUE(last > first);
    EXPECT_TRUE(last >= first);
    EXPECT_TRUE(first <= first);
    EXPECT_TRUE(first >= first);
}

TEST(JIArrayCoverageTest, ConstIteratorDecrement) {
    const zint1 arr{10, 20, 30};
    auto it = arr.cend();
    --it;
    EXPECT_EQ(*it, 30);
    it--;
    EXPECT_EQ(*it, 20);
}

TEST(JIArrayCoverageTest, ConstIteratorCompoundAssign) {
    const zint1 arr{1, 2, 3, 4, 5};
    auto it = arr.cbegin();
    it += 3;
    EXPECT_EQ(*it, 4);
    it -= 2;
    EXPECT_EQ(*it, 2);
}

TEST(JIArrayCoverageTest, IteratorToConstIteratorConversion) {
    zint1 arr{10, 20, 30};
    auto it = arr.begin();
    // Implicit conversion from Iterator to ConstIterator
    typename zint1::ConstIterator cit = it;
    EXPECT_EQ(*cit, 10);
}

// ============================================================================
// Iterator mutable operations
// ============================================================================

TEST(JIArrayCoverageTest, IteratorArithmetic) {
    zint1 arr{1, 2, 3, 4, 5};
    auto it = arr.begin();
    EXPECT_EQ(*(it + 2), 3);
    EXPECT_EQ(it[4], 5);
    auto it2 = arr.end();
    EXPECT_EQ(it2 - it, 5);
}

TEST(JIArrayCoverageTest, IteratorCompoundAssign) {
    zint1 arr{10, 20, 30, 40};
    auto it = arr.begin();
    it += 2;
    EXPECT_EQ(*it, 30);
    it -= 1;
    EXPECT_EQ(*it, 20);
}

TEST(JIArrayCoverageTest, IteratorComparison) {
    zint1 arr{1, 2, 3};
    auto a = arr.begin();
    auto b = arr.end();
    EXPECT_TRUE(a < b);
    EXPECT_TRUE(b > a);
    EXPECT_TRUE(a <= a);
    EXPECT_TRUE(a >= a);
}

TEST(JIArrayCoverageTest, IteratorDecrement) {
    zint1 arr{10, 20, 30};
    auto it = arr.end();
    --it;
    EXPECT_EQ(*it, 30);
    it--;
    EXPECT_EQ(*it, 20);
}

TEST(JIArrayCoverageTest, IteratorPostIncrement) {
    zint1 arr{10, 20};
    auto it = arr.begin();
    auto old = it++;
    EXPECT_EQ(*old, 10);
    EXPECT_EQ(*it, 20);
}

TEST(JIArrayCoverageTest, IteratorSubtract) {
    zint1 arr{1, 2, 3, 4};
    auto it = arr.end();
    auto it2 = it - 2;
    EXPECT_EQ(*it2, 3);
}

// ============================================================================
// ffor_back macro
// ============================================================================

TEST(JIArrayCoverageTest, FforBackMacro) {
    zint1 arr{10, 20, 30, 40, 50};
    std::vector<int> reversed;
    ffor_back(i, 5, 1) {
        reversed.push_back(arr(i));
    }
    EXPECT_EQ(reversed, (std::vector<int>{50, 40, 30, 20, 10}));
}

TEST(JIArrayCoverageTest, FforBackEveryElement) {
    zint1 arr{10, 20, 30, 40, 50};
    std::vector<int> collected;
    ffor_back(i, 5, 1) {
        collected.push_back(arr(i));
    }
    EXPECT_EQ(collected, (std::vector<int>{50, 40, 30, 20, 10}));
}

// ============================================================================
// init0() range-based init
// ============================================================================

TEST(JIArrayCoverageTest, Init0CustomRanges) {
    zint2 arr;
    // init0(min1, max1, min2, max2) -> dim1=[3..5], dim2=[10..12]
    arr.init0(3, 5, 10, 12);
    EXPECT_EQ(arr.size(), 9); // 3 * 3
    arr(3, 10) = 42;
    arr(5, 12) = 99;
    EXPECT_EQ(arr(3, 10), 42);
    EXPECT_EQ(arr(5, 12), 99);
}

TEST(JIArrayCoverageTest, Init0_1D) {
    zint1 arr;
    arr.init0(5, 10); // indices [5..10], size = 6
    EXPECT_EQ(arr.size(), 6);
    arr(5) = 100;
    arr(10) = 200;
    EXPECT_EQ(arr(5), 100);
    EXPECT_EQ(arr(10), 200);
}

// ============================================================================
// setSize() - resize only if dimensions differ
// ============================================================================

TEST(JIArrayCoverageTest, SetSizeSameDimNoOp) {
    zint1 arr(5);
    arr(1) = 42;
    arr.setSize(5); // same size, should not reallocate
    EXPECT_EQ(arr(1), 42); // data preserved
}

TEST(JIArrayCoverageTest, SetSizeDifferentDimReallocates) {
    zint1 arr(3);
    arr(1) = 99;
    arr.setSize(5);
    EXPECT_EQ(arr.size(), 5);
    // After reallocation, data is zero-initialized
    EXPECT_EQ(arr(1), 0);
}

TEST(JIArrayCoverageTest, SetSize2D) {
    zint2 arr(2, 3);
    arr.setSize(2, 3); // same, no-op
    EXPECT_EQ(arr.size(), 6);

    arr.setSize(4, 5);
    EXPECT_EQ(arr.size(), 20);
}

// ============================================================================
// shareWith()
// ============================================================================

TEST(JIArrayCoverageTest, ShareWithSharesMemory) {
    zint1 original(4);
    original(1) = 10;
    original(2) = 20;
    original(3) = 30;
    original(4) = 40;

    zint1 view;
    view.shareWith(original);

    EXPECT_EQ(view(1), 10);
    EXPECT_EQ(view(4), 40);

    // Modifying via view affects original
    view(2) = 999;
    EXPECT_EQ(original(2), 999);
}

// ============================================================================
// copy() - deep copy
// ============================================================================

TEST(JIArrayCoverageTest, CopyCreatesIndependent) {
    zint1 arr{10, 20, 30};
    auto copied = arr.copy();
    EXPECT_EQ(copied(1), 10);
    EXPECT_EQ(copied(3), 30);

    // Modifying copy doesn't affect original
    copied(1) = 999;
    EXPECT_EQ(arr(1), 10);
}

// ============================================================================
// operator== / operator!=
// ============================================================================

TEST(JIArrayCoverageTest, EqualityDifferentContent) {
    zint1 a{1, 2, 3};
    zint1 b{1, 2, 4};
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a != b);
}

TEST(JIArrayCoverageTest, EqualitySameContent) {
    zint1 a{1, 2, 3};
    zint1 b{1, 2, 3};
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);
}

// ============================================================================
// Unary negation
// ============================================================================

TEST(JIArrayCoverageTest, UnaryNegation) {
    zdouble1 arr{1.0, -2.0, 3.0};
    auto neg = -arr;
    EXPECT_DOUBLE_EQ(neg(1), -1.0);
    EXPECT_DOUBLE_EQ(neg(2), 2.0);
    EXPECT_DOUBLE_EQ(neg(3), -3.0);
}

// ============================================================================
// Subtraction operators
// ============================================================================

TEST(JIArrayCoverageTest, SubtractionOperator) {
    zint1 a{10, 20, 30};
    zint1 b{3, 5, 7};
    auto result = a - b;
    EXPECT_EQ(result(1), 7);
    EXPECT_EQ(result(2), 15);
    EXPECT_EQ(result(3), 23);
}

TEST(JIArrayCoverageTest, SubtractionEqualsOperator) {
    zint1 arr{10, 20, 30};
    zint1 sub{1, 2, 3};
    arr -= sub;
    EXPECT_EQ(arr(1), 9);
    EXPECT_EQ(arr(2), 18);
    EXPECT_EQ(arr(3), 27);
}

// ============================================================================
// Division operators
// ============================================================================

TEST(JIArrayCoverageTest, DivisionByScalarInt) {
    zint1 arr{10, 20, 30};
    arr /= 5;
    EXPECT_EQ(arr(1), 2);
    EXPECT_EQ(arr(2), 4);
    EXPECT_EQ(arr(3), 6);
}

TEST(JIArrayCoverageTest, DivisionByScalarDouble) {
    zdouble1 arr{10.0, 20.0, 30.0};
    arr /= 4.0;
    EXPECT_DOUBLE_EQ(arr(1), 2.5);
    EXPECT_DOUBLE_EQ(arr(2), 5.0);
    EXPECT_DOUBLE_EQ(arr(3), 7.5);
}

TEST(JIArrayCoverageTest, ArrayDivisionEquals) {
    zdouble1 a{10.0, 20.0, 30.0};
    zdouble1 b{2.0, 4.0, 5.0};
    a /= b;
    EXPECT_DOUBLE_EQ(a(1), 5.0);
    EXPECT_DOUBLE_EQ(a(2), 5.0);
    EXPECT_DOUBLE_EQ(a(3), 6.0);
}

TEST(JIArrayCoverageTest, ArrayDivisionOperator) {
    zdouble1 a{10.0, 20.0, 30.0};
    zdouble1 b{2.0, 5.0, 6.0};
    auto result = a / b;
    EXPECT_DOUBLE_EQ(result(1), 5.0);
    EXPECT_DOUBLE_EQ(result(2), 4.0);
    EXPECT_DOUBLE_EQ(result(3), 5.0);
}

TEST(JIArrayCoverageTest, ScalarDivideByArray) {
    zdouble1 arr{2.0, 4.0, 5.0};
    auto result = 20.0 / arr;
    EXPECT_DOUBLE_EQ(result(1), 10.0);
    EXPECT_DOUBLE_EQ(result(2), 5.0);
    EXPECT_DOUBLE_EQ(result(3), 4.0);
}

TEST(JIArrayCoverageTest, ArrayDivideByScalar) {
    zdouble1 arr{10.0, 20.0, 30.0};
    auto result = arr / 5.0;
    EXPECT_DOUBLE_EQ(result(1), 2.0);
    EXPECT_DOUBLE_EQ(result(2), 4.0);
    EXPECT_DOUBLE_EQ(result(3), 6.0);
}

TEST(JIArrayCoverageTest, ArrayDivideByScalarInt) {
    zint1 arr{10, 20, 30};
    auto result = arr / 5;
    EXPECT_EQ(result(1), 2);
    EXPECT_EQ(result(2), 4);
    EXPECT_EQ(result(3), 6);
}

// ============================================================================
// Scalar + array, array + scalar
// ============================================================================

TEST(JIArrayCoverageTest, ScalarPlusArray) {
    zint1 arr{1, 2, 3};
    auto result = 10 + arr;
    EXPECT_EQ(result(1), 11);
    EXPECT_EQ(result(2), 12);
    EXPECT_EQ(result(3), 13);
}

TEST(JIArrayCoverageTest, ArrayPlusScalar) {
    zint1 arr{1, 2, 3};
    auto result = arr + 10;
    EXPECT_EQ(result(1), 11);
    EXPECT_EQ(result(2), 12);
    EXPECT_EQ(result(3), 13);
}

// ============================================================================
// Scalar * array, array * scalar
// ============================================================================

TEST(JIArrayCoverageTest, ScalarTimesArray) {
    zdouble1 arr{1.0, 2.0, 3.0};
    auto result = 3.0 * arr;
    EXPECT_DOUBLE_EQ(result(1), 3.0);
    EXPECT_DOUBLE_EQ(result(2), 6.0);
    EXPECT_DOUBLE_EQ(result(3), 9.0);
}

TEST(JIArrayCoverageTest, ArrayTimesScalar) {
    zdouble1 arr{1.0, 2.0, 3.0};
    auto result = arr * 3.0;
    EXPECT_DOUBLE_EQ(result(1), 3.0);
    EXPECT_DOUBLE_EQ(result(3), 9.0);
}

// ============================================================================
// Array multiply (element-wise) operators
// ============================================================================

TEST(JIArrayCoverageTest, ArrayMultiplyEquals) {
    zint1 a{2, 3, 4};
    zint1 b{5, 6, 7};
    a *= b;
    EXPECT_EQ(a(1), 10);
    EXPECT_EQ(a(2), 18);
    EXPECT_EQ(a(3), 28);
}

TEST(JIArrayCoverageTest, ArrayMultiplyOperator) {
    zint1 a{2, 3, 4};
    zint1 b{5, 6, 7};
    auto result = a * b;
    EXPECT_EQ(result(1), 10);
    EXPECT_EQ(result(2), 18);
    EXPECT_EQ(result(3), 28);
}

// ============================================================================
// sqsum() and dot()
// ============================================================================

TEST(JIArrayCoverageTest, SqSum) {
    zdouble1 arr{3.0, 4.0};
    // 3^2 + 4^2 = 9 + 16 = 25
    EXPECT_DOUBLE_EQ(arr.sqsum(), 25.0);
}

TEST(JIArrayCoverageTest, DotProduct) {
    zdouble1 a{1.0, 2.0, 3.0};
    zdouble1 b{4.0, 5.0, 6.0};
    // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    EXPECT_DOUBLE_EQ(dot(a, b), 32.0);
}

// ============================================================================
// findFirst() multi-dimensional
// ============================================================================

TEST(JIArrayCoverageTest, FindFirst2DFound) {
    zint2 arr(3, 3);
    arr = 0;
    arr(2, 3) = 42;
    auto loc = arr.findFirst(42);
    // Should return position (2, 3)
    EXPECT_EQ(loc(1), 2);
    EXPECT_EQ(loc(2), 3);
}

TEST(JIArrayCoverageTest, FindFirst2DNotFound) {
    zint2 arr(2, 2);
    arr = 0;
    auto loc = arr.findFirst(99);
    // Not found: all indices = JIARRAY_OFFSET - 1 = 0
    EXPECT_EQ(loc(1), 0);
    EXPECT_EQ(loc(2), 0);
}

TEST(JIArrayCoverageTest, FindFirst1DFound) {
    zint1 arr{10, 20, 30};
    EXPECT_EQ(arr.findFirst(20), 2);
    EXPECT_EQ(arr.findFirst(10), 1);
}

TEST(JIArrayCoverageTest, FindFirst1DNotFound) {
    zint1 arr{1, 2, 3};
    EXPECT_EQ(arr.findFirst(99), 0);
}

// ============================================================================
// FastArray index access
// ============================================================================

TEST(JIArrayCoverageTest, AccessWithFastArrayIndex) {
    zint2 arr(3, 4);
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 4; ++j)
            arr(i, j) = i * 10 + j;

    FastArray<int, 2> idx({2, 3});
    EXPECT_EQ(arr(idx), 23);
}

TEST(JIArrayCoverageTest, ConstAccessWithFastArrayIndex) {
    zint2 arr(3, 3);
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
            arr(i, j) = i + j;

    const zint2& cref = arr;
    FastArray<int, 2> idx({2, 3});
    EXPECT_EQ(cref(idx), 5); // 2 + 3
}

// ============================================================================
// Assignment from std::vector to empty array
// ============================================================================

TEST(JIArrayCoverageTest, AssignVectorToEmptyAutoInit) {
    zint1 arr;
    std::vector<int> v{10, 20, 30, 40};
    arr = v;
    EXPECT_EQ(arr.size(), 4);
    EXPECT_EQ(arr(1), 10);
    EXPECT_EQ(arr(4), 40);
}

// ============================================================================
// Assignment from initializer list
// ============================================================================

TEST(JIArrayCoverageTest, AssignInitializerList) {
    zint1 arr(3);
    arr = {7, 8, 9};
    EXPECT_EQ(arr(1), 7);
    EXPECT_EQ(arr(2), 8);
    EXPECT_EQ(arr(3), 9);
}

// ============================================================================
// getRankSize(), getSizeOfRank(), getSize(rank), getOffset()
// ============================================================================

TEST(JIArrayCoverageTest, GetRankSizeColumnMajor2D) {
    zint2 arr(3, 4);
    const int* rs = arr.getRankSize();
    // Column-major: rankSize[0] = 1, rankSize[1] = 3
    EXPECT_EQ(rs[0], 1);
    EXPECT_EQ(rs[1], 3);
}

TEST(JIArrayCoverageTest, GetSizeOfRank) {
    zint2 arr(3, 4);
    const int* s = arr.getSizeOfRank();
    EXPECT_EQ(s[0], 3);
    EXPECT_EQ(s[1], 4);
}

TEST(JIArrayCoverageTest, GetSizeByRank) {
    zint3 arr(2, 3, 4);
    EXPECT_EQ(arr.getSize(1), 2);
    EXPECT_EQ(arr.getSize(2), 3);
    EXPECT_EQ(arr.getSize(3), 4);
}

TEST(JIArrayCoverageTest, GetOffset) {
    zint2 arr(3, 4);
    const int* off = arr.getOffset();
    EXPECT_EQ(off[0], 1); // JIARRAY_OFFSET
    EXPECT_EQ(off[1], 1);
    EXPECT_EQ(arr.getOffset(1), 1);
    EXPECT_EQ(arr.getOffset(2), 1);
}

// ============================================================================
// setOffsets() variadic
// ============================================================================

TEST(JIArrayCoverageTest, SetOffsetsVariadic) {
    zint2 arr(3, 4);
    arr = 0;
    arr.setOffsets(0, 0);
    arr(0, 0) = 42;
    arr(2, 3) = 99;
    EXPECT_EQ(arr(0, 0), 42);
    EXPECT_EQ(arr(2, 3), 99);
}

// ============================================================================
// Const slice
// ============================================================================

TEST(JIArrayCoverageTest, ConstSlice) {
    zint2 arr(3, 4);
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 4; ++j)
            arr(i, j) = i * 10 + j;

    const zint2& cref = arr;
    auto s = cref.slice(2); // column 2 in column-major
    EXPECT_EQ(s(1), 12);
    EXPECT_EQ(s(2), 22);
    EXPECT_EQ(s(3), 32);
}

// ============================================================================
// Destroy and re-init
// ============================================================================

TEST(JIArrayCoverageTest, DestroyThenReinit) {
    zint1 arr(5);
    arr(3) = 42;
    arr.destroy();
    EXPECT_EQ(arr.size(), 0);

    arr.init(3);
    EXPECT_EQ(arr.size(), 3);
    arr(1) = 10;
    EXPECT_EQ(arr(1), 10);
}

// ============================================================================
// isAllocated()
// ============================================================================

TEST(JIArrayCoverageTest, IsAllocated) {
    zint1 arr;
    EXPECT_FALSE(arr.isAllocated());

    arr.init(5);
    EXPECT_TRUE(arr.isAllocated());

    arr.destroy();
    EXPECT_FALSE(arr.isAllocated());
}

// ============================================================================
// External memory init
// ============================================================================

TEST(JIArrayCoverageTest, ExternalMemoryInit2D) {
    double buffer[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    zdouble2 arr(buffer, 2, 3);
    // Column-major: arr(1,1)=1, arr(2,1)=2, arr(1,2)=3, arr(2,2)=4
    EXPECT_DOUBLE_EQ(arr(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(arr(2, 1), 2.0);
    EXPECT_DOUBLE_EQ(arr(1, 2), 3.0);
    EXPECT_DOUBLE_EQ(arr(2, 3), 6.0);
}

// ============================================================================
// init() with zero dimensions
// ============================================================================

TEST(JIArrayCoverageTest, InitZeroDimensionNoAlloc) {
    zint2 arr;
    arr.init(0, 5);
    EXPECT_EQ(arr.size(), 0);
}

// ============================================================================
// Reshape
// ============================================================================

TEST(JIArrayCoverageTest, Reshape1Dto2D) {
    zdouble1 arr(6);
    for (int i = 1; i <= 6; ++i)
        arr(i) = static_cast<double>(i);

    auto reshaped = arr.reshape(2, 3);
    EXPECT_EQ(reshaped.size(), 6);
    // Data is shared, column-major: (1,1)=1, (2,1)=2, (1,2)=3 ...
    EXPECT_DOUBLE_EQ(reshaped(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(reshaped(2, 1), 2.0);
    EXPECT_DOUBLE_EQ(reshaped(1, 2), 3.0);
}

// ============================================================================
// convertToVector()
// ============================================================================

TEST(JIArrayCoverageTest, ConvertToVector) {
    zint1 arr{10, 20, 30};
    auto v = arr.convertToVector();
    EXPECT_EQ(v, (std::vector<int>{10, 20, 30}));
}

// ============================================================================
// Deep copy via operator=
// ============================================================================

TEST(JIArrayCoverageTest, DeepCopyAssignment) {
    zint1 src{1, 2, 3};
    zint1 dst;
    dst = src;
    EXPECT_EQ(dst(1), 1);
    EXPECT_EQ(dst(3), 3);
    // Independent copy
    dst(1) = 999;
    EXPECT_EQ(src(1), 1);
}

TEST(JIArrayCoverageTest, DeepCopyAssignmentExisting) {
    zint1 src{10, 20, 30};
    zint1 dst(3);
    dst = src;
    EXPECT_EQ(dst(1), 10);
    EXPECT_EQ(dst(3), 30);
}

// ============================================================================
// String array
// ============================================================================

TEST(JIArrayCoverageTest, StringArray1D) {
    zstring1 arr(3);
    arr(1) = "hello";
    arr(2) = "world";
    arr(3) = "test";
    EXPECT_EQ(arr(1), "hello");
    EXPECT_EQ(arr(3), "test");
    EXPECT_TRUE(arr.contains("world"));
    EXPECT_FALSE(arr.contains("missing"));
}

// ============================================================================
// 4D and 5D arrays
// ============================================================================

TEST(JIArrayCoverageTest, Array4DAccessPattern) {
    zint4 arr(2, 3, 4, 5);
    EXPECT_EQ(arr.size(), 120);
    arr(1, 2, 3, 4) = 1234;
    arr(2, 3, 4, 5) = 2345;
    EXPECT_EQ(arr(1, 2, 3, 4), 1234);
    EXPECT_EQ(arr(2, 3, 4, 5), 2345);
}

TEST(JIArrayCoverageTest, Array4DSliceTo3D) {
    zint4 arr(2, 3, 4, 5);
    arr = 0;
    arr(1, 2, 3, 4) = 99;
    auto s = arr.slice(4);
    EXPECT_EQ(s(1, 2, 3), 99);
}

// ============================================================================
// zfor macro
// ============================================================================

TEST(JIArrayCoverageTest, ZforMacro) {
    zint1 arr(5);
    zfor(i, 5) {
        arr(i) = i * i;
    }
    EXPECT_EQ(arr(1), 1);
    EXPECT_EQ(arr(2), 4);
    EXPECT_EQ(arr(3), 9);
    EXPECT_EQ(arr(4), 16);
    EXPECT_EQ(arr(5), 25);
}

// ============================================================================
// ffor with step
// ============================================================================

TEST(JIArrayCoverageTest, FforBasic) {
    std::vector<int> collected;
    ffor(i, 3, 7) {
        collected.push_back(i);
    }
    EXPECT_EQ(collected, (std::vector<int>{3, 4, 5, 6, 7}));
}

// ============================================================================
// Type traits
// ============================================================================

TEST(JIArrayCoverageTest, IsJIArrayTrait) {
    using int1_t = JIArray<int, 1>;
    using double3_t = JIArray<double, 3>;
    using vec_t = std::vector<int>;
    EXPECT_TRUE(is_jiarray_v<int1_t>);
    EXPECT_TRUE(is_jiarray_v<double3_t>);
    EXPECT_FALSE(is_jiarray_v<int>);
    EXPECT_FALSE(is_jiarray_v<vec_t>);
}

// ============================================================================
// getMemory() with offset
// ============================================================================

TEST(JIArrayCoverageTest, GetMemoryWithOffset) {
    zint1 arr{10, 20, 30, 40};
    int* p0 = arr.getMemory(0);
    EXPECT_EQ(*p0, 10);
    int* p2 = arr.getMemory(2);
    EXPECT_EQ(*p2, 30);
}

TEST(JIArrayCoverageTest, ConstGetMemory) {
    const zint1 arr{10, 20, 30};
    const int* p = arr.getMemory(1);
    EXPECT_EQ(*p, 20);
}

// ============================================================================
// Scalar assignment to all elements
// ============================================================================

TEST(JIArrayCoverageTest, ScalarAssignment2D) {
    zdouble2 arr(3, 4);
    arr = 7.5;
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 4; ++j)
            EXPECT_DOUBLE_EQ(arr(i, j), 7.5);
}

// ============================================================================
// += with scalar and array
// ============================================================================

TEST(JIArrayCoverageTest, PlusEqualsScalar) {
    zint1 arr{1, 2, 3};
    arr += 10;
    EXPECT_EQ(arr(1), 11);
    EXPECT_EQ(arr(2), 12);
    EXPECT_EQ(arr(3), 13);
}

TEST(JIArrayCoverageTest, PlusEqualsArray) {
    zint1 a{1, 2, 3};
    zint1 b{10, 20, 30};
    a += b;
    EXPECT_EQ(a(1), 11);
    EXPECT_EQ(a(2), 22);
    EXPECT_EQ(a(3), 33);
}

// ============================================================================
// *= with scalar
// ============================================================================

TEST(JIArrayCoverageTest, MultiplyEqualsScalar) {
    zdouble1 arr{2.0, 3.0, 4.0};
    arr *= 5.0;
    EXPECT_DOUBLE_EQ(arr(1), 10.0);
    EXPECT_DOUBLE_EQ(arr(2), 15.0);
    EXPECT_DOUBLE_EQ(arr(3), 20.0);
}

// ============================================================================
// Coverage gap: operator== with different sizes -> return false (line 1086)
// ============================================================================

TEST(JIArrayCoverageTest, EqualityDifferentSize) {
    zint1 a(3);
    a = 1;
    zint1 b(5);
    b = 1;
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a != b);
}

// ============================================================================
// Coverage gap: maxloc(from, to) on multi-dimensional array (lines 1540-1542)
// ============================================================================

TEST(JIArrayCoverageTest, MaxlocRange2D) {
    zint2 arr(3, 4);
    arr = 0;
    // Place max value at (2, 3)
    arr(2, 3) = 99;
    // Search full range
    auto loc = arr.maxloc(1, 12);
    EXPECT_EQ(loc(1), 2);
    EXPECT_EQ(loc(2), 3);
}

TEST(JIArrayCoverageTest, MaxlocRange2DPartial) {
    zint2 arr(3, 3);
    arr = 0;
    arr(1, 1) = 5;
    arr(2, 2) = 10;
    arr(3, 3) = 3;
    // Search range [4, 9] (linear indices 4..9 -> covers (1,2) through (3,3))
    auto loc = arr.maxloc(4, 9);
    EXPECT_EQ(loc(1), 2);
    EXPECT_EQ(loc(2), 2);
}
