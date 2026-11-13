// init0's range contract and large-offset addressing on the managed-memory JICudaArray, and the
// 32-bit unsigned device position of both JICudaArray and JIArray. Built for both index bases,
// with and without JIARRAY_DEBUG. JICudaArray is column-major only, so there is no row-major build.
#include <jiarray/JICudaArray.h>

#include <gtest/gtest.h>

#include <climits>
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

__global__ void cuReadCudaArray(cdouble2 array, int i, int row, double* out) {
    *out = array(i, row);
}

__global__ void cuReadHostTypeView(zdouble2 array, int i, int row, double* out) {
    *out = array(i, row);
}
} // namespace

TEST(CudaArrayInit0Test, Rank1RangeFollowsTheIndexBase) {
    cint1 a;
    a.init0(5, 9);
#if JIARRAY_OFFSET == 0
    EXPECT_EQ(a.getSize(), 4); // [5, 9)
#else
    EXPECT_EQ(a.getSize(), 5); // [5, 9]
#endif
    EXPECT_EQ(a.data(), &a(5));
    EXPECT_EQ(&a(8) - &a(5), 3);
}

TEST(CudaArrayInit0Test, Rank2RowStrideIsTheInnerSize) {
    cdouble2 x;
    x.init0(0, endOf(2), 10, endOf(13));
    EXPECT_EQ(x.getSize(), 12);
    EXPECT_EQ(x.data(), &x(0, 10));
    EXPECT_EQ(&x(0, 11) - &x(0, 10), 3);
    EXPECT_EQ(&x(2, 13) - x.data(), 11);
}

TEST(CudaArrayInit0Test, AnEqualFirstAndEndHoldsNoRowOnlyWhenZeroBased) {
    cdouble2 x;
    x.init0(0, 3, 7, 7);
#if JIARRAY_OFFSET == 0
    EXPECT_FALSE(x.isAllocated());
    EXPECT_EQ(x.getSize(), 0);
#else
    EXPECT_TRUE(x.isAllocated());
    EXPECT_EQ(x.getSize(), 4);
#endif
}

TEST(CudaArrayInit0Test, ALargeFirstRowAddressesTheSameElementOnHostAndDevice) {
    cdouble2 x;
    x.init0(0, endOf(kInnerSize - 1), kLargeFirst, endOf(kLargeFirst + 2));
    ASSERT_EQ(x.getSize(), 3 * kInnerSize);
    x                                  = 0.0;
    x(kInnerSize - 1, kLargeFirst + 2) = 7.0;
    EXPECT_EQ(&x(kInnerSize - 1, kLargeFirst + 2) - x.data(), 3 * kInnerSize - 1);

    cdouble1 out(1);
    out = 0.0;
    cuReadCudaArray<<<1, 1>>>(x, kInnerSize - 1, kLargeFirst + 2, out.data());
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
    EXPECT_EQ(out(kBase), 7.0);

    // The host type over the same managed buffer takes its own device path.
    zdouble2 view;
    view.init(kInnerSize, 3, x.data());
    view.setOffsets(0, kLargeFirst);
    out = 0.0;
    cuReadHostTypeView<<<1, 1>>>(view, kInnerSize - 1, kLargeFirst + 2, out.data());
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
    EXPECT_EQ(out(kBase), 7.0);
}

// 65536 x 65537 wraps an int to 65536: unchecked, it would allocate.
TEST(CudaArrayInit0Test, AnArrayPastAnIntPositionIsNotAllocated) {
    cdouble2 x;
    EXPECT_THROW(x.init(65536, 65537), std::length_error);
    EXPECT_FALSE(x.isAllocated());
}

#ifdef JIARRAY_DEBUG
TEST(CudaArrayInit0Test, ARowOutsideTheRangeTrapsInDebug) {
    cdouble2 x;
    x.init0(0, endOf(2), 10, endOf(13));
    EXPECT_THROW(x(0, 9), std::out_of_range);
    EXPECT_THROW(x(0, 14), std::out_of_range);
    EXPECT_NO_THROW(x(0, 13));
}
#endif
