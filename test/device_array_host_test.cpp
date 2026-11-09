// JIDeviceArray from a host-only compiler (no nvcc): the owner, transfers and
// view metadata must work with just cuda_runtime.h and CUDA::cudart.

#include <jiarray/JIDeviceArray.h>

#include <gtest/gtest.h>

#include <utility>

#ifdef __CUDACC__
#error "this test must be compiled by the host compiler"
#endif

using namespace dnegri::jiarray;

TEST(DeviceArrayHost, RoundTripAndViewsWithoutNvcc) {
    zfloat2 h;
    h.init0(0, 4, 1, 3);
    for (int n = 0; n < h.size(); ++n) {
        h[n] = 0.25f * n;
    }
    dfloat2 d;
    d = h;
    EXPECT_EQ(d.getOffset(1), 0);
    EXPECT_EQ(d.getSize(1), 5);

    dfloat2 moved = std::move(d);
    EXPECT_FALSE(d.isAllocated());

    zfloat2 back;
    moved.copyTo(back);
    for (int n = 0; n < h.size(); ++n) {
        ASSERT_EQ(back[n], h[n]);
    }

    dfloat1::View column = moved.slice(2);
    zfloat1 host_column = h.slice(2);
    EXPECT_EQ(column.size(), host_column.size());
    EXPECT_EQ(column.data() - moved.data(), host_column.data() - h.data());
}
