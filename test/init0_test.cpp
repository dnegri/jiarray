// init0's range contract and large-offset addressing on the host type.
// Built for both index bases and both storage orders, with and without JIARRAY_DEBUG.
//
// With JIARRAY_OFFSET 0 an init0 range is [first, end), as ffor's; with 1 it is [first, end].
// A block of rows taken from a large global row space must address exactly: the product of a
// stride and the first row passes INT_MAX long before the element count does.
#include <jiarray/JIArray.h>
#include <gtest/gtest.h>

#include <climits>
#include <cstddef>
#include <stdexcept>

using namespace dnegri::jiarray;

namespace {
constexpr int kBase       = JIARRAY_OFFSET;
constexpr int kInnerSize  = 1300;      ///< elements per row
constexpr int kLargeFirst = 2'100'000; ///< a first row whose offset x kInnerSize passes INT_MAX

/// init0's end argument for a range whose last index is `last`.
constexpr int endOf(int last) {
    return last + 1 - kBase;
}

/// Allocate `num_row` rows [first_row, first_row + num_row) of `inner` elements, the row index
/// being the slow dimension in either storage order.
void initRows(zdouble2& array, int inner, int first_row, int num_row) {
#if JIARRAY_COLUMN_MAJOR
    array.init0(0, endOf(inner - 1), first_row, endOf(first_row + num_row - 1));
#else
    array.init0(first_row, endOf(first_row + num_row - 1), 0, endOf(inner - 1));
#endif
}

/// Element `i` of row `row` of an array built by initRows.
double& element(zdouble2& array, int i, int row) {
#if JIARRAY_COLUMN_MAJOR
    return array(i, row);
#else
    return array(row, i);
#endif
}

/// The member layout before the offset sum became 64-bit, for the size check.
template <size_t RANK>
struct OldLayout {
    int     nn;
    double* mm;
    int     allocated;
    int     rank_size[RANK];
    int     offset[RANK];
    int     sum_of_offset;
    int     sizes[RANK];
};
} // namespace

TEST(Init0Test, Rank1RangeFollowsTheIndexBase) {
    zint1 a;
    a.init0(5, 9);
#if JIARRAY_OFFSET == 0
    EXPECT_EQ(a.getSize(), 4); // [5, 9)
#else
    EXPECT_EQ(a.getSize(), 5); // [5, 9]
#endif
    EXPECT_EQ(a.getOffset(kBase), 5);
    EXPECT_EQ(a.data(), &a(5));
    EXPECT_EQ(&a(8) - &a(5), 3);
}

TEST(Init0Test, Rank2RowStrideIsTheInnerSize) {
    zdouble2 x;
    initRows(x, 3, 10, 4);
    EXPECT_EQ(x.getSize(), 12);
    EXPECT_EQ(x.data(), &element(x, 0, 10));
    EXPECT_EQ(&element(x, 0, 11) - &element(x, 0, 10), 3);
    EXPECT_EQ(&element(x, 2, 13) - x.data(), 11);
}

TEST(Init0Test, Rank2ExplicitBoundsGiveTheExpectedShape) {
    zdouble2 x;
    x.init0(0, 3, 10, 14);
#if JIARRAY_OFFSET == 0
    EXPECT_EQ(x.getSize(kBase), 3);
    EXPECT_EQ(x.getSize(kBase + 1), 4);
#else
    EXPECT_EQ(x.getSize(kBase), 4);
    EXPECT_EQ(x.getSize(kBase + 1), 5);
#endif
    EXPECT_EQ(x.getOffset(kBase), 0);
    EXPECT_EQ(x.getOffset(kBase + 1), 10);
}

TEST(Init0Test, AnEqualFirstAndEndHoldsNoRowOnlyWhenZeroBased) {
    zdouble2 x;
    x.init0(0, 3, 7, 7);
#if JIARRAY_OFFSET == 0
    EXPECT_FALSE(x.isAllocated());
    EXPECT_EQ(x.getSize(), 0);
#else
    EXPECT_TRUE(x.isAllocated());
    EXPECT_EQ(x.getSize(), 4); // [0, 3] x [7, 7]
#endif
}

// Under UBSan (-fsanitize=signed-integer-overflow) an int offset sum reports the overflow here.
TEST(Init0Test, ALargeFirstRowAddressesExactly) {
    zdouble2 x;
    initRows(x, kInnerSize, kLargeFirst, 3);
    ASSERT_EQ(x.getSize(), 3 * kInnerSize);

    element(x, kInnerSize - 1, kLargeFirst + 2) = 7.0;
    EXPECT_EQ(&element(x, kInnerSize - 1, kLargeFirst + 2) - x.data(), 3 * kInnerSize - 1);
    EXPECT_EQ(&element(x, 0, kLargeFirst) - x.data(), 0);
    EXPECT_EQ(x.slice(kLargeFirst + 2)(kInnerSize - 1), 7.0);

    FastArray<int, 2> idx;
#if JIARRAY_COLUMN_MAJOR
    idx[0] = kInnerSize - 1;
    idx[1] = kLargeFirst + 2;
#else
    idx[0] = kLargeFirst + 2;
    idx[1] = kInnerSize - 1;
#endif
    EXPECT_EQ(x(idx), 7.0);
}

TEST(Init0Test, TheElementCountIsCheckedAgainstAnIntPosition) {
    const int at_limit[]   = {INT_MAX, 1};
    const int past_limit[] = {46341, 46341}; // 46341^2 = INT_MAX + 4,634
    const int three[]      = {65536, 65536, 65536};
    EXPECT_EQ(checkElementCount(at_limit, 2), INT_MAX);
    EXPECT_THROW(checkElementCount(past_limit, 2), std::length_error);
    EXPECT_THROW(checkElementCount(three, 3), std::length_error);
}

// 65536 x 65536 wraps an int to 0 and 65536 x 65537 to 65536: unchecked, both would allocate.
TEST(Init0Test, AnArrayPastAnIntPositionIsNotAllocated) {
    zdouble2 x;
    EXPECT_THROW(x.init(65536, 65537), std::length_error);
    EXPECT_FALSE(x.isAllocated());
    EXPECT_THROW(initRows(x, kInnerSize, 0, INT_MAX / kInnerSize + 1), std::length_error);
    EXPECT_FALSE(x.isAllocated());
}

TEST(Init0Test, TheObjectDoesNotGrow) {
    EXPECT_EQ(sizeof(zdouble1), sizeof(OldLayout<1>));
    EXPECT_EQ(sizeof(zdouble2), sizeof(OldLayout<2>));
    EXPECT_EQ(sizeof(zdouble3), sizeof(OldLayout<3>));
}

#ifdef JIARRAY_DEBUG
// Bounds are checked only with JIARRAY_DEBUG: a row outside the range traps.
TEST(Init0Test, ARowOutsideTheRangeTrapsInDebug) {
    zdouble2 x;
    initRows(x, 3, 10, 4);
    EXPECT_THROW(element(x, 0, 9), std::out_of_range);
    EXPECT_THROW(element(x, 0, 14), std::out_of_range);
    EXPECT_NO_THROW(element(x, 0, 13));
}
#endif
