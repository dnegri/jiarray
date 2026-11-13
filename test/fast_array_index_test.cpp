// FastArray-indexed element access must honour the global JIARRAY_OFFSET.
// Built twice: JIARRAY_OFFSET=0 and the default 1, each in both storage orders.
#include <jiarray/JIArray.h>
#include <gtest/gtest.h>

using namespace dnegri::jiarray;

namespace {
constexpr int kBase = JIARRAY_OFFSET;
}

TEST(FastArrayIndexTest, Rank2MatchesScalarIndex) {
    zint2 arr(3, 4);
    for (int i = kBase; i < kBase + 3; ++i)
        for (int j = kBase; j < kBase + 4; ++j)
            arr(i, j) = 10 * i + j;

    FastArray<int, 2> idx({kBase + 1, kBase + 2});
    EXPECT_EQ(arr(idx), arr(kBase + 1, kBase + 2));

    const zint2& cref = arr;
    FastArray<int, 2> first({kBase, kBase});
    EXPECT_EQ(cref(first), arr(kBase, kBase));
}

TEST(FastArrayIndexTest, Rank3MatchesScalarIndex) {
    zint3 arr(2, 3, 4);
    int value = 0;
    for (int i = kBase; i < kBase + 2; ++i)
        for (int j = kBase; j < kBase + 3; ++j)
            for (int k = kBase; k < kBase + 4; ++k)
                arr(i, j, k) = value++;

    FastArray<int, 3> idx({kBase + 1, kBase + 2, kBase + 3});
    EXPECT_EQ(arr(idx), arr(kBase + 1, kBase + 2, kBase + 3));
}

TEST(FastArrayIndexTest, CustomOffsetsMatchScalarIndex) {
    zint2 arr;
    arr.init0(-1, 2 - kBase, 5, 8 - kBase); // [-1, 1] x [5, 7]: init0's end is exclusive when 0-based
    for (int i = -1; i <= 1; ++i)
        for (int j = 5; j <= 7; ++j)
            arr(i, j) = 100 * i + j;

    FastArray<int, 2> idx({0, 6});
    EXPECT_EQ(arr(idx), arr(0, 6));
}
