#include <gtest/gtest.h>
#include "../src/JIArray.h"  // Ensure this path points to your JIArray class

using namespace dnegri::jiarray;

TEST(JIArrayTestsWithOffset, SliceBasic) {
    JIArray<double, 2> array(3, 3); // Create a 3x3 array
    array = {1.0, 2.0, 3.0,
             4.0, 5.0, 6.0,
             7.0, 8.0, 9.0}; // Column-major order

    auto slice = array.slice(1); // Slice to get the first column
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 1.0);
    EXPECT_EQ(slice(2), 2.0);
    EXPECT_EQ(slice(3), 3.0);
}

TEST(JIArrayTestsWithOffset, SliceMiddleColumn) {
    JIArray<double, 2> array(3, 3); // Create a 3x3 array
    array = {1.0, 2.0, 3.0,
             4.0, 5.0, 6.0,
             7.0, 8.0, 9.0}; // Column-major order

    auto slice = array.slice(2); // Slice to get the second column
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 4.0);
    EXPECT_EQ(slice(2), 5.0);
    EXPECT_EQ(slice(3), 6.0);
}

TEST(JIArrayTestsWithOffset, SliceLastColumn) {
    JIArray<double, 2> array(3, 3); // Create a 3x3 array
    array = {1.0, 2.0, 3.0,
             4.0, 5.0, 6.0,
             7.0, 8.0, 9.0}; // Column-major order

    auto slice = array.slice(3); // Slice to get the third column
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 7.0);
    EXPECT_EQ(slice(2), 8.0);
    EXPECT_EQ(slice(3), 9.0);
}

TEST(JIArrayTestsWithOffset, SliceOn3DArray) {
    JIArray<double, 3> array(2, 2, 2); // Create a 2x2x2 array
    array = {1.0, 2.0,
             3.0, 4.0,
             5.0, 6.0,
             7.0, 8.0}; // Column-major order

    auto slice = array.slice(1); // Slice along the first dimension
    EXPECT_EQ(slice.size(), 4); // 2x2 slice

    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(2, 1), 2.0);
    EXPECT_EQ(slice(1, 2), 3.0);
    EXPECT_EQ(slice(2, 2), 4.0);
}

TEST(JIArrayTestsWithOffset, SliceOn3DArraySecondDimension) {
    JIArray<double, 3> array(2, 2, 2); // Create a 2x2x2 array
    array = {1.0, 2.0,
             3.0, 4.0,
             5.0, 6.0,
             7.0, 8.0}; // Column-major order

    auto slice = array.slice(2); // Slice along the second dimension
    EXPECT_EQ(slice.size(), 4); // 2x2 slice

    EXPECT_EQ(slice(1, 1), 5.0);
    EXPECT_EQ(slice(2, 1), 6.0);
    EXPECT_EQ(slice(1, 2), 7.0);
    EXPECT_EQ(slice(2, 2), 8.0);
}

TEST(JIArrayTestsWithOffset, MultiDimensionalSlicing) {
    JIArray<double, 3> array(3, 3, 3); // Create a 3x3x3 array
    array = {
        1.0,  2.0,  3.0,
        4.0,  5.0,  6.0,
        7.0,  8.0,  9.0,

        10.0, 11.0, 12.0,
        13.0, 14.0, 15.0,
        16.0, 17.0, 18.0,

        19.0, 20.0, 21.0,
        22.0, 23.0, 24.0,
        25.0, 26.0, 27.0
    }; // Column-major order

    // Slice the first 2D "plane"
    auto slice = array.slice(1);
    EXPECT_EQ(slice.size(), 9);
    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(2, 1), 2.0);
    EXPECT_EQ(slice(3, 1), 3.0);
    EXPECT_EQ(slice(1, 2), 4.0);
    EXPECT_EQ(slice(3, 3), 9.0);
}

TEST(JIArrayTestsWithOffset, SliceInvalidIndex) {
    JIArray<double, 2> array(3, 3); // Create a 3x3 array
    array = {1.0, 2.0, 3.0,
             4.0, 5.0, 6.0,
             7.0, 8.0, 9.0}; // Column-major order

    // Attempt to slice an invalid dimension index
    EXPECT_THROW(array.slice(4), std::out_of_range);
}

TEST(JIArrayTestsWithOffset, SliceNonSquareArray) {
    JIArray<double, 2> array(3, 2); // Create a 3x2 array
    array = {1.0, 2.0,
             3.0, 4.0,
             5.0, 6.0}; // Column-major order

    auto slice = array.slice(1); // Slice to get the first column
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 1.0);
    EXPECT_EQ(slice(2), 2.0);
    EXPECT_EQ(slice(3), 3.0);

    slice = array.slice(2); // Slice to get the second column
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 4.0);
    EXPECT_EQ(slice(2), 5.0);
    EXPECT_EQ(slice(3), 6.0);
}

TEST(JIArrayTestsWithOffset, SliceRectangular3DArray) {
    JIArray<double, 3> array(3, 2, 2); // Create a 3x2x2 array
    array = {
        1.0, 2.0,
        3.0, 4.0,
        5.0, 6.0,
        7.0, 8.0,
        9.0, 10.0,
        11.0, 12.0
    }; // Column-major order

    auto slice = array.slice(1); // Slice along the first "layer"
    EXPECT_EQ(slice.size(), 6); // Should return a 3x2 slice

    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(2, 1), 2.0);
    EXPECT_EQ(slice(3, 1), 3.0);
    EXPECT_EQ(slice(1, 2), 4.0);
    EXPECT_EQ(slice(2, 2), 5.0);
    EXPECT_EQ(slice(3, 2), 6.0);
}

TEST(JIArrayTestsWithOffset, SliceNonUniform3DArray) {
    JIArray<double, 3> array(4, 3, 2); // Create a 4x3x2 array
    array = {
        1.0,  2.0,
        3.0,  4.0,
        5.0,  6.0,
        7.0,  8.0,
        9.0, 10.0,
        11.0, 12.0,

        13.0, 14.0,
        15.0, 16.0,
        17.0, 18.0,
        19.0, 20.0,
        21.0, 22.0,
        23.0, 24.0
    }; // Column-major order

    // Slice the first "plane" along the last dimension
    auto slice = array.slice(1); 
    EXPECT_EQ(slice.size(), 12); // Should return a 4x3 slice

    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(2, 1), 2.0);
    EXPECT_EQ(slice(3, 1), 3.0);
    EXPECT_EQ(slice(4, 1), 4.0);
    EXPECT_EQ(slice(1, 2), 5.0);
    EXPECT_EQ(slice(2, 2), 6.0);
    EXPECT_EQ(slice(3, 2), 7.0);
    EXPECT_EQ(slice(4, 2), 8.0);
    EXPECT_EQ(slice(1, 3), 9.0);
    EXPECT_EQ(slice(2, 3), 10.0);
    EXPECT_EQ(slice(3, 3), 11.0);
    EXPECT_EQ(slice(4, 3), 12.0);
}

TEST(JIArrayTestsWithOffset, Slice4DArray) {
    JIArray<double, 4> array(2, 2, 2, 2); // Create a 2x2x2x2 array
    array = {
        1.0, 2.0,
        3.0, 4.0,
        5.0, 6.0,
        7.0, 8.0,
        
        9.0, 10.0,
        11.0, 12.0,
        13.0, 14.0,
        15.0, 16.0
    }; // Column-major order

    auto slice = array.slice(1); // Slice the first "3D layer"
    EXPECT_EQ(slice.size(), 8); // Should return a 2x2x2 slice

    EXPECT_EQ(slice(1, 1, 1), 1.0);
    EXPECT_EQ(slice(2, 1, 1), 2.0);
    EXPECT_EQ(slice(1, 2, 1), 3.0);
    EXPECT_EQ(slice(2, 2, 1), 4.0);
    EXPECT_EQ(slice(1, 1, 2), 5.0);
    EXPECT_EQ(slice(2, 1, 2), 6.0);
    EXPECT_EQ(slice(2, 2, 2), 8.0);
}

TEST(JIArrayTestsWithOffset, SliceIrregularSizedArray) {
    JIArray<double, 2> array(5, 4); // Create a 5x4 array
    array = {
        1.0,  2.0,  3.0,  4.0,
        5.0,  6.0,  7.0,  8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0
    }; // Column-major order

    auto slice = array.slice(3); // Get the third column
    EXPECT_EQ(slice.size(), 5);
    EXPECT_EQ(slice(1), 11.0);
    EXPECT_EQ(slice(2), 12.0);
    EXPECT_EQ(slice(3), 13.0);
    EXPECT_EQ(slice(4), 14.0);
    EXPECT_EQ(slice(5), 15.0);
}