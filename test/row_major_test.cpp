#include "../src/JIArray.h"
#include <algorithm>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

using namespace dnegri::jiarray;

// ============================================================================
// CONSTRUCTION AND INITIALIZATION TESTS (ROW-MAJOR)
// ============================================================================

TEST(JIArrayRowMajorTests, DefaultConstruction) {

    JIArray<double, 2> array;
    EXPECT_FALSE(array.isAllocated());
    EXPECT_EQ(array.size(), 0);
}

TEST(JIArrayRowMajorTests, SizeConstruction) {
    JIArray<double, 2> array(3, 4);
    EXPECT_TRUE(array.isAllocated());
    EXPECT_EQ(array.size(), 12);
    EXPECT_EQ(array.getSize(1), 3);
    EXPECT_EQ(array.getSize(2), 4);
}

TEST(JIArrayRowMajorTests, InitializerListConstruction) {
    JIArray<int, 1> array = {1, 2, 3, 4, 5};
    EXPECT_TRUE(array.isAllocated());
    EXPECT_EQ(array.size(), 5);
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(array(i + JIARRAY_OFFSET), i + 1);
    }
}

// ============================================================================
// ROW-MAJOR SPECIFIC ELEMENT ACCESS TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, ElementAccess2DRowMajor) {
    JIArray<double, 2> array(2, 3);
    // In row-major order, elements are stored row by row
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    // First row
    EXPECT_EQ(array(1, 1), 1.0);
    EXPECT_EQ(array(1, 2), 2.0);
    EXPECT_EQ(array(1, 3), 3.0);

    // Second row
    EXPECT_EQ(array(2, 1), 4.0);
    EXPECT_EQ(array(2, 2), 5.0);
    EXPECT_EQ(array(2, 3), 6.0);
}

TEST(JIArrayRowMajorTests, ElementAccess3DRowMajor) {
    JIArray<int, 3> array(2, 2, 2);
    // Row-major storage: elements vary fastest in last dimension
    array = {1, 2, 3, 4, 5, 6, 7, 8};

    EXPECT_EQ(array(1, 1, 1), 1);
    EXPECT_EQ(array(1, 1, 2), 2);
    EXPECT_EQ(array(1, 2, 1), 3);
    EXPECT_EQ(array(1, 2, 2), 4);
    EXPECT_EQ(array(2, 1, 1), 5);
    EXPECT_EQ(array(2, 1, 2), 6);
    EXPECT_EQ(array(2, 2, 1), 7);
    EXPECT_EQ(array(2, 2, 2), 8);
}

// ============================================================================
// ROW-MAJOR SLICING TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, SliceRowMajor2D) {
    JIArray<double, 2> array(3, 3);
    // Row-major order
    array = {1.0, 2.0, 3.0,  // First row
             4.0, 5.0, 6.0,  // Second row
             7.0, 8.0, 9.0}; // Third row

    auto slice = array.slice(1); // Slice first row
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 1.0);
    EXPECT_EQ(slice(2), 2.0);
    EXPECT_EQ(slice(3), 3.0);

    auto slice2 = array.slice(2); // Slice second row
    EXPECT_EQ(slice2.size(), 3);
    EXPECT_EQ(slice2(1), 4.0);
    EXPECT_EQ(slice2(2), 5.0);
    EXPECT_EQ(slice2(3), 6.0);
}

TEST(JIArrayRowMajorTests, SliceRowMajor3D) {
    JIArray<double, 3> array(2, 2, 2);
    // Row-major: first dimension varies slowest
    array = {1.0, 2.0,  // [1,1,:]
             3.0, 4.0,  // [1,2,:]
             5.0, 6.0,  // [2,1,:]
             7.0, 8.0}; // [2,2,:]

    auto slice = array.slice(1); // First "layer"
    EXPECT_EQ(slice.size(), 4);
    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(1, 2), 2.0);
    EXPECT_EQ(slice(2, 1), 3.0);
    EXPECT_EQ(slice(2, 2), 4.0);
}

// ============================================================================
// ROW-MAJOR ARITHMETIC OPERATIONS TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, ScalarAdditionRowMajor) {
    JIArray<double, 2> array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0}; // Row-major

    array += 5.0;
    EXPECT_EQ(array(1, 1), 6.0);
    EXPECT_EQ(array(1, 2), 7.0);
    EXPECT_EQ(array(2, 1), 8.0);
    EXPECT_EQ(array(2, 2), 9.0);
}

TEST(JIArrayRowMajorTests, ArrayMultiplicationRowMajor) {
    JIArray<double, 2> array1(2, 2);
    JIArray<double, 2> array2(2, 2);

    array1 = {2.0, 3.0, 4.0, 5.0}; // Row-major
    array2 = {1.0, 2.0, 3.0, 4.0}; // Row-major

    auto result = array1 * array2;
    EXPECT_EQ(result(1, 1), 2.0);  // 2*1
    EXPECT_EQ(result(1, 2), 6.0);  // 3*2
    EXPECT_EQ(result(2, 1), 12.0); // 4*3
    EXPECT_EQ(result(2, 2), 20.0); // 5*4
}

// ============================================================================
// ROW-MAJOR STATISTICAL OPERATIONS TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, AverageRowMajor) {
    JIArray<double, 2> array(2, 2);
    array = {1.0, 3.0, 5.0, 7.0}; // Row-major storage

    EXPECT_DOUBLE_EQ(array.average(), 4.0);
}

TEST(JIArrayRowMajorTests, SumRowMajor) {
    JIArray<int, 2> array(2, 3);
    array = {1, 2, 3, 4, 5, 6}; // Row-major

    EXPECT_EQ(array.sum(), 21);
}

// ============================================================================
// ROW-MAJOR MEMORY LAYOUT VERIFICATION TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, MemoryLayoutVerification) {
    JIArray<int, 2> array(2, 3);
    array = {1, 2, 3, 4, 5, 6}; // Should be stored as [1,2,3,4,5,6] in memory

    int *data = array.data();
    EXPECT_EQ(data[0], 1); // First element in memory
    EXPECT_EQ(data[1], 2); // Second element in memory
    EXPECT_EQ(data[2], 3); // Third element in memory
    EXPECT_EQ(data[3], 4); // Fourth element in memory
    EXPECT_EQ(data[4], 5); // Fifth element in memory
    EXPECT_EQ(data[5], 6); // Sixth element in memory
}

TEST(JIArrayRowMajorTests, MemoryLayout3DVerification) {
    JIArray<int, 3> array(2, 2, 2);
    // Fill in logical order for row-major
    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
            for (int k = 1; k <= 2; ++k) {
                array(i, j, k) = (i - 1) * 4 + (j - 1) * 2 + k;
            }
        }
    }

    int *data = array.data();
    // In row-major, memory should be arranged as:
    // [1,2,3,4,5,6,7,8] corresponding to [1,1,1], [1,1,2], [1,2,1], [1,2,2], [2,1,1], [2,1,2],
    // [2,2,1], [2,2,2]
    for (int i = 0; i < 8; ++i) {
        EXPECT_EQ(data[i], i + 1);
    }
}

// ============================================================================
// ROW-MAJOR ITERATOR TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, IteratorRowMajorOrder) {
    JIArray<int, 2> array(2, 3);
    // In row-major, initialization should follow row-wise pattern
    array = {1, 2, 3, 4, 5, 6};

    std::vector<int> expected = {1, 2, 3, 4, 5, 6}; // Row-major memory order
    std::vector<int> actual;

    for (const auto &value : array) {
        actual.push_back(value);
    }

    EXPECT_EQ(actual, expected);
}

// ============================================================================
// ROW-MAJOR RESHAPING TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, ReshapeRowMajor) {
    JIArray<int, 1> array1d(6);
    array1d = {1, 2, 3, 4, 5, 6};

    auto array2d = array1d.reshape(2, 3);
    EXPECT_EQ(array2d.size(), 6);

    // In row-major, reshape should preserve memory order
    EXPECT_EQ(array2d(1, 1), 1);
    EXPECT_EQ(array2d(1, 2), 2);
    EXPECT_EQ(array2d(1, 3), 3);
    EXPECT_EQ(array2d(2, 1), 4);
    EXPECT_EQ(array2d(2, 2), 5);
    EXPECT_EQ(array2d(2, 3), 6);
}

// ============================================================================
// ROW-MAJOR COMPARISON WITH COLUMN-MAJOR TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, RowMajorVsColumnMajorDifference) {
    JIArray<int, 2> array(2, 3);

    // Fill array using logical indices
    int counter = 1;
    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 3; ++j) {
            array(i, j) = counter++;
        }
    }

    // In row-major, memory should be: [1,2,3,4,5,6]
    // In column-major, it would be: [1,4,2,5,3,6]
    int *data = array.data();
    std::vector<int> memory_layout(data, data + 6);
    std::vector<int> expected_row_major = {1, 2, 3, 4, 5, 6};

    EXPECT_EQ(memory_layout, expected_row_major);
}

// ============================================================================
// ROW-MAJOR ASSIGNMENT TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, InitializerListAssignmentRowMajor) {
    JIArray<double, 2> array(2, 3);
    array = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};

    // Verify row-major storage pattern
    EXPECT_EQ(array(1, 1), 1.1);
    EXPECT_EQ(array(1, 2), 2.2);
    EXPECT_EQ(array(1, 3), 3.3);
    EXPECT_EQ(array(2, 1), 4.4);
    EXPECT_EQ(array(2, 2), 5.5);
    EXPECT_EQ(array(2, 3), 6.6);
}

// ============================================================================
// ROW-MAJOR COMPLEX OPERATIONS TESTS
// ============================================================================

TEST(JIArrayRowMajorTests, ComplexSlicingRowMajor) {
    JIArray<double, 4> array(2, 2, 2, 2);
    // Initialize with sequential values
    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
            for (int k = 1; k <= 2; ++k) {
                for (int l = 1; l <= 2; ++l) {
                    int value = (i - 1) * 8 + (j - 1) * 4 + (k - 1) * 2 + l;
                    array(i, j, k, l) = static_cast<double>(value);
                }
            }
        }
    }

    auto slice3d = array.slice(1); // Remove first dimension
    EXPECT_EQ(slice3d.size(), 8);

    auto slice2d = slice3d.slice(1); // Remove second dimension
    EXPECT_EQ(slice2d.size(), 4);

    auto slice1d = slice2d.slice(1); // Remove third dimension
    EXPECT_EQ(slice1d.size(), 2);
}

TEST(JIArrayRowMajorTests, DotProductRowMajor) {
    JIArray<double, 1> array1(4);
    JIArray<double, 1> array2(4);

    array1 = {1.0, 2.0, 3.0, 4.0};
    array2 = {5.0, 6.0, 7.0, 8.0};

    double result = dot(array1, array2);
    EXPECT_DOUBLE_EQ(result, 70.0); // 1*5 + 2*6 + 3*7 + 4*8 = 70
}

// Test to verify compilation with row-major preprocessor setting
TEST(JIArrayRowMajorTests, CompilationTest) {
    // This test simply verifies that the code compiles with row-major setting
    JIArray<int, 3> test_array(2, 2, 2);
    EXPECT_TRUE(test_array.isAllocated());
    EXPECT_EQ(test_array.size(), 8);
}

// ============================================================================
// ROW-MAJOR SLICING TESTS (Converted from Column-Major)
// ============================================================================

TEST(JIArrayRowMajorTests, SliceBasic) {
    JIArray<double, 2> array(3, 3);
    // Row-major: elements stored row by row
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    auto slice = array.slice(1); // Slice first row
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 1.0);
    EXPECT_EQ(slice(2), 2.0);
    EXPECT_EQ(slice(3), 3.0);
}

TEST(JIArrayRowMajorTests, SliceMiddleRow) {
    JIArray<double, 2> array(3, 3);                        // Create a 3x3 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; // Row-major order

    auto slice = array.slice(2); // Slice to get the second row
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 4.0);
    EXPECT_EQ(slice(2), 5.0);
    EXPECT_EQ(slice(3), 6.0);
}

TEST(JIArrayRowMajorTests, SliceLastRow) {
    JIArray<double, 2> array(3, 3);                        // Create a 3x3 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; // Row-major order

    auto slice = array.slice(3); // Slice to get the third row
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 7.0);
    EXPECT_EQ(slice(2), 8.0);
    EXPECT_EQ(slice(3), 9.0);
}

TEST(JIArrayRowMajorTests, SliceOn3DArray) {
    JIArray<double, 3> array(2, 2, 2);                // Create a 2x2x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // Row-major order

    auto slice = array.slice(1); // Slice along the first dimension
    EXPECT_EQ(slice.size(), 4);  // 2x2 slice

    // In row-major: [1,1,:] = [1,2], [1,2,:] = [3,4]
    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(1, 2), 2.0);
    EXPECT_EQ(slice(2, 1), 3.0);
    EXPECT_EQ(slice(2, 2), 4.0);
}

TEST(JIArrayRowMajorTests, SliceOn3DArraySecondDimension) {
    JIArray<double, 3> array(2, 2, 2);                // Create a 2x2x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // Row-major order

    auto slice = array.slice(2); // Slice along the second dimension
    EXPECT_EQ(slice.size(), 4);  // 2x2 slice

    // In row-major: [2,1,:] = [5,6], [2,2,:] = [7,8]
    EXPECT_EQ(slice(1, 1), 5.0);
    EXPECT_EQ(slice(1, 2), 6.0);
    EXPECT_EQ(slice(2, 1), 7.0);
    EXPECT_EQ(slice(2, 2), 8.0);
}

TEST(JIArrayRowMajorTests, MultiDimensionalSlicing) {
    JIArray<double, 3> array(3, 3, 3); // Create a 3x3x3 array
    // Row-major order: elements are laid out as [i][j][k] where k varies fastest
    array = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,   // [1,:,:]
             10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,  // [2,:,:]
             19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0}; // [3,:,:]

    // Slice the first 2D "plane" - first slice along dimension 1
    auto slice = array.slice(1);
    EXPECT_EQ(slice.size(), 9);
    EXPECT_EQ(slice(1, 1), 1.0); // [1,1]
    EXPECT_EQ(slice(1, 2), 2.0); // [1,2]
    EXPECT_EQ(slice(1, 3), 3.0); // [1,3]
    EXPECT_EQ(slice(2, 1), 4.0); // [2,1]
    EXPECT_EQ(slice(3, 3), 9.0); // [3,3]
}

TEST(JIArrayRowMajorTests, SliceInvalidIndex) {
    JIArray<double, 2> array(3, 3);                        // Create a 3x3 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; // Row-major order

    // Attempt to slice an invalid dimension index
    EXPECT_THROW(array.slice(4), std::out_of_range);
}

TEST(JIArrayRowMajorTests, SliceNonSquareArray) {
    JIArray<double, 2> array(3, 2);         // Create a 3x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // Row-major order: [1,2], [3,4], [5,6]

    auto slice = array.slice(1); // Slice to get the first row
    EXPECT_EQ(slice.size(), 2);
    EXPECT_EQ(slice(1), 1.0);
    EXPECT_EQ(slice(2), 2.0);

    slice = array.slice(2); // Slice to get the second row
    EXPECT_EQ(slice.size(), 2);
    EXPECT_EQ(slice(1), 3.0);
    EXPECT_EQ(slice(2), 4.0);
}

TEST(JIArrayRowMajorTests, SliceRectangular3DArray) {
    JIArray<double, 3> array(3, 2, 2);                                       // Create a 3x2x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}; // Row-major order

    auto slice = array.slice(1); // Slice along the first "layer"
    EXPECT_EQ(slice.size(), 4);  // Should return a 2x2 slice

    // In row-major: [1,:,:] = [[1,2], [3,4]]
    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(1, 2), 2.0);
    EXPECT_EQ(slice(2, 1), 3.0);
    EXPECT_EQ(slice(2, 2), 4.0);
}

TEST(JIArrayRowMajorTests, SliceNonUniform3DArray) {
    JIArray<double, 3> array(4, 3, 2); // Create a 4x3x2 array
    // Row-major order: rightmost index varies fastest
    array = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,   // [1,:,:]
             7.0,  8.0,  9.0,  10.0, 11.0, 12.0,  // [2,:,:]
             13.0, 14.0, 15.0, 16.0, 17.0, 18.0,  // [3,:,:]
             19.0, 20.0, 21.0, 22.0, 23.0, 24.0}; // [4,:,:]

    // Slice the 4th "plane" along the first dimension
    auto slice = array.slice(4);
    EXPECT_EQ(slice.size(), 6); // Should return a 3x2 slice

    // In row-major: [4,:,:] contains first 6 elements arranged as 3x2
    EXPECT_EQ(slice(1, 1), 19.0);
    EXPECT_EQ(slice(1, 2), 20.0);
    EXPECT_EQ(slice(2, 1), 21.0);
    EXPECT_EQ(slice(2, 2), 22.0);
    EXPECT_EQ(slice(3, 1), 23.0);
    EXPECT_EQ(slice(3, 2), 24.0);
}

TEST(JIArrayRowMajorTests, Slice4DArray) {
    JIArray<double, 4> array(2, 2, 2, 2); // Create a 2x2x2x2 array
    // Row-major order: last index varies fastest
    array = {1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,   // [1,:,:,:]
             9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0}; // [2,:,:,:]

    auto slice = array.slice(1); // Slice the first "3D layer"
    EXPECT_EQ(slice.size(), 8);  // Should return a 2x2x2 slice

    // In row-major: [1,:,:,:] = first 8 elements
    EXPECT_EQ(slice(1, 1, 1), 1.0);
    EXPECT_EQ(slice(1, 1, 2), 2.0);
    EXPECT_EQ(slice(1, 2, 1), 3.0);
    EXPECT_EQ(slice(1, 2, 2), 4.0);
    EXPECT_EQ(slice(2, 1, 1), 5.0);
    EXPECT_EQ(slice(2, 1, 2), 6.0);
    EXPECT_EQ(slice(2, 2, 2), 8.0);
}

TEST(JIArrayRowMajorTests, SliceIrregularSizedArray) {
    JIArray<double, 2> array(5, 4); // Create a 5x4 array
    // Row-major order: elements stored row by row
    array = {1.0,  2.0,  3.0,  4.0,   // Row 1
             5.0,  6.0,  7.0,  8.0,   // Row 2
             9.0,  10.0, 11.0, 12.0,  // Row 3
             13.0, 14.0, 15.0, 16.0,  // Row 4
             17.0, 18.0, 19.0, 20.0}; // Row 5

    auto slice = array.slice(3); // Get the third row
    EXPECT_EQ(slice.size(), 4);
    EXPECT_EQ(slice(1), 9.0);
    EXPECT_EQ(slice(2), 10.0);
    EXPECT_EQ(slice(3), 11.0);
    EXPECT_EQ(slice(4), 12.0);
}
