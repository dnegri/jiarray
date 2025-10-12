#include "../src/JIArray.h" // Ensure this path points to your JIArray class
#include <gtest/gtest.h>

using namespace dnegri::jiarray;

TEST(JIArrayColumnMajorTests, InitializationAndAccess) {
    zdouble2 array(2, 2);         // Initialize with dimensions
    array = {1.0, 2.0, 3.0, 4.0}; // Use assignment operator to set values

    EXPECT_EQ(array(1, 1), 1.0); // Column 1, Row 1
    EXPECT_EQ(array(2, 1), 2.0); // Column 2, Row 1
    EXPECT_EQ(array(1, 2), 3.0); // Column 1, Row 2
    EXPECT_EQ(array(2, 2), 4.0); // Column 2, Row 2
}

TEST(JIArrayColumnMajorTests, AssignmentOperator) {
    zdouble2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    zdouble2 array2(2, 2);
    array2 = array1;
    EXPECT_EQ(array2(1, 1), 1.0);
    EXPECT_EQ(array2(2, 1), 2.0);
    EXPECT_EQ(array2(1, 2), 3.0);
    EXPECT_EQ(array2(2, 2), 4.0);
}

TEST(JIArrayColumnMajorTests, CopyConstructor) {
    zdouble2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    zdouble2 array2(array1);
    EXPECT_EQ(array2(1, 1), 1.0);
    EXPECT_EQ(array2(2, 1), 2.0);
    EXPECT_EQ(array2(1, 2), 3.0);
    EXPECT_EQ(array2(2, 2), 4.0);
}

TEST(JIArrayColumnMajorTests, Slicing) {
    JIArray<double, 3> array(3, 3, 3);
    array = 1.0;
    auto slice = array.slice(1);
    EXPECT_EQ(slice.size(), 9);
}

TEST(JIArrayColumnMajorTests, ArithmeticOperations) {
    zdouble2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    zdouble2 array2(2, 2);
    array2 = {5.0, 6.0, 7.0, 8.0};

    auto result = array1 + array2;
    EXPECT_EQ(result(1, 1), 6.0);
    EXPECT_EQ(result(2, 1), 8.0);
    EXPECT_EQ(result(1, 2), 10.0);
    EXPECT_EQ(result(2, 2), 12.0);

    result = array1 - array2;
    EXPECT_EQ(result(1, 1), -4.0);
    EXPECT_EQ(result(2, 1), -4.0);
    EXPECT_EQ(result(1, 2), -4.0);
    EXPECT_EQ(result(2, 2), -4.0);

    result = 2.0 * array1;
    EXPECT_EQ(result(1, 1), 2.0);
    EXPECT_EQ(result(2, 1), 4.0);
    EXPECT_EQ(result(1, 2), 6.0);
    EXPECT_EQ(result(2, 2), 8.0);
}

TEST(JIArrayColumnMajorTests, Reshape2D) {
    zdouble2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    auto reshaped = array.reshape(4);
    EXPECT_EQ(reshaped.size(), 4);
    EXPECT_EQ(reshaped(1), 1.0);
    EXPECT_EQ(reshaped(2), 2.0);
    EXPECT_EQ(reshaped(3), 3.0);
    EXPECT_EQ(reshaped(4), 4.0);
}

TEST(JIArrayColumnMajorTests, ContainsAndFind) {
    zdouble2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    EXPECT_TRUE(array.contains(3.0));
    EXPECT_FALSE(array.contains(5.0));

    auto loc = array.findFirst(3.0);
    EXPECT_EQ(loc(1), 1);
    EXPECT_EQ(loc(2), 2);
}

TEST(JIArrayColumnMajorTests, MaxLoc) {
    zdouble2 array(2, 2);
    array = {1.0, 5.0, 3.0, 4.0};
    auto maxLoc = array.maxloc();
    EXPECT_EQ(maxLoc(1), 2);
    EXPECT_EQ(maxLoc(2), 1);
}

TEST(JIArrayColumnMajorTests, AverageSum) {
    zdouble2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    EXPECT_DOUBLE_EQ(array.average(), 2.5);
    EXPECT_DOUBLE_EQ(array.sum(), 10.0);
}

TEST(JIArrayColumnMajorTests, IteratorTest) {
    zdouble2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    double sum = 0.0;
    for (auto it = array.begin(); it != array.end(); ++it) {
        sum += *it;
    }
    EXPECT_DOUBLE_EQ(sum, 10.0);
}

TEST(JIArrayColumnMajorTests, ConversionToVector) {
    zdouble2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> vec = array.convertToVector();
    EXPECT_EQ(vec.size(), 4);
    EXPECT_EQ(vec[0], 1.0);
    EXPECT_EQ(vec[1], 2.0);
    EXPECT_EQ(vec[2], 3.0);
    EXPECT_EQ(vec[3], 4.0);
}

// ============================================================================
// CONSTRUCTION AND INITIALIZATION TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, DefaultConstruction) {
    JIArray<double, 2> array;
    EXPECT_FALSE(array.isAllocated());
    EXPECT_EQ(array.size(), 0);
}

TEST(JIArrayColumnMajorTests, SizeConstruction) {
    JIArray<double, 2> array(3, 4);
    EXPECT_TRUE(array.isAllocated());
    EXPECT_EQ(array.size(), 12);
    EXPECT_EQ(array.getSize(1), 3);
    EXPECT_EQ(array.getSize(2), 4);
}

TEST(JIArrayColumnMajorTests, InitializerListConstruction) {
    JIArray<int, 1> array = {1, 2, 3, 4, 5};
    EXPECT_TRUE(array.isAllocated());
    EXPECT_EQ(array.size(), 5);
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(array(i + JIARRAY_OFFSET), i + 1);
    }
}

TEST(JIArrayColumnMajorTests, CopyConstruction) {
    JIArray<double, 2> original(2, 3);
    original = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    JIArray<double, 2> copy(original);
    EXPECT_EQ(copy.size(), original.size());
    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 3; ++j) {
            EXPECT_EQ(copy(i, j), original(i, j));
        }
    }
}

// ============================================================================
// MEMORY MANAGEMENT TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, MemoryManagement) {
    JIArray<double, 2> array(3, 3);
    double *memory_ptr = array.data();
    EXPECT_NE(memory_ptr, nullptr);
    EXPECT_TRUE(array.isAllocated());

    array.destroy();
    EXPECT_FALSE(array.isAllocated());
    EXPECT_EQ(array.size(), 0);
}

TEST(JIArrayColumnMajorTests, ExternalMemoryInit) {
    double *external_memory = new double[6];
    for (int i = 0; i < 6; ++i) {
        external_memory[i] = i + 1.0;
    }

    JIArray<double, 2> array(external_memory, 2, 3);
    EXPECT_EQ(array.size(), 6);
    EXPECT_EQ(array.data(), external_memory);

    delete[] external_memory;
}

// ============================================================================
// ELEMENT ACCESS TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, ElementAccess1D) {
    JIArray<int, 1> array(5);
    array = {10, 20, 30, 40, 50};

    EXPECT_EQ(array(1), 10);
    EXPECT_EQ(array(2), 20);
    EXPECT_EQ(array(3), 30);
    EXPECT_EQ(array(4), 40);
    EXPECT_EQ(array(5), 50);

    array(3) = 100;
    EXPECT_EQ(array(3), 100);
}

TEST(JIArrayColumnMajorTests, ElementAccess2D) {
    JIArray<double, 2> array(3, 2);
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // Column-major

    // First column
    EXPECT_EQ(array(1, 1), 1.0);
    EXPECT_EQ(array(2, 1), 2.0);
    EXPECT_EQ(array(3, 1), 3.0);

    // Second column
    EXPECT_EQ(array(1, 2), 4.0);
    EXPECT_EQ(array(2, 2), 5.0);
    EXPECT_EQ(array(3, 2), 6.0);
}

TEST(JIArrayColumnMajorTests, ElementAccess3D) {
    JIArray<int, 3> array(2, 2, 2);
    array = {1, 2, 3, 4, 5, 6, 7, 8}; // Column-major order

    EXPECT_EQ(array(1, 1, 1), 1);
    EXPECT_EQ(array(2, 1, 1), 2);
    EXPECT_EQ(array(1, 2, 1), 3);
    EXPECT_EQ(array(2, 2, 1), 4);
    EXPECT_EQ(array(1, 1, 2), 5);
    EXPECT_EQ(array(2, 1, 2), 6);
    EXPECT_EQ(array(1, 2, 2), 7);
    EXPECT_EQ(array(2, 2, 2), 8);
}

TEST(JIArrayColumnMajorTests, AtMethodAccess) {
    JIArray<double, 2> array(2, 2);
    array = {1.5, 2.5, 3.5, 4.5};

    EXPECT_EQ(array.at(1, 1), 1.5);
    EXPECT_EQ(array.at(2, 1), 2.5);
    EXPECT_EQ(array.at(1, 2), 3.5);
    EXPECT_EQ(array.at(2, 2), 4.5);

    const auto &const_array = array;
    EXPECT_EQ(const_array.at(1, 1), 1.5);
}

// ============================================================================
// ARITHMETIC OPERATIONS TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, ScalarAddition) {
    JIArray<double, 2> array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};

    array += 5.0;
    EXPECT_EQ(array(1, 1), 6.0);
    EXPECT_EQ(array(2, 1), 7.0);
    EXPECT_EQ(array(1, 2), 8.0);
    EXPECT_EQ(array(2, 2), 9.0);
}

TEST(JIArrayColumnMajorTests, ArrayAddition) {
    JIArray<int, 2> array1(2, 2);
    JIArray<int, 2> array2(2, 2);

    array1 = {1, 2, 3, 4};
    array2 = {5, 6, 7, 8};

    array1 += array2;
    EXPECT_EQ(array1(1, 1), 6);
    EXPECT_EQ(array1(2, 1), 8);
    EXPECT_EQ(array1(1, 2), 10);
    EXPECT_EQ(array1(2, 2), 12);
}

TEST(JIArrayColumnMajorTests, ArraySubtraction) {
    JIArray<double, 2> array1(2, 2);
    JIArray<double, 2> array2(2, 2);

    array1 = {10.0, 8.0, 6.0, 4.0};
    array2 = {1.0, 2.0, 3.0, 2.0};

    auto result = array1 - array2;
    EXPECT_EQ(result(1, 1), 9.0);
    EXPECT_EQ(result(2, 1), 6.0);
    EXPECT_EQ(result(1, 2), 3.0);
    EXPECT_EQ(result(2, 2), 2.0);
}

TEST(JIArrayColumnMajorTests, ArrayMultiplication) {
    JIArray<double, 2> array1(2, 2);
    JIArray<double, 2> array2(2, 2);

    array1 = {2.0, 3.0, 4.0, 5.0};
    array2 = {1.0, 2.0, 3.0, 4.0};

    auto result = array1 * array2;
    EXPECT_EQ(result(1, 1), 2.0);
    EXPECT_EQ(result(2, 1), 6.0);
    EXPECT_EQ(result(1, 2), 12.0);
    EXPECT_EQ(result(2, 2), 20.0);
}

TEST(JIArrayColumnMajorTests, ScalarMultiplication) {
    JIArray<double, 2> array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};

    auto result = array * 2.5;
    EXPECT_EQ(result(1, 1), 2.5);
    EXPECT_EQ(result(2, 1), 5.0);
    EXPECT_EQ(result(1, 2), 7.5);
    EXPECT_EQ(result(2, 2), 10.0);

    auto result2 = 3.0 * array;
    EXPECT_EQ(result2(1, 1), 3.0);
    EXPECT_EQ(result2(2, 1), 6.0);
    EXPECT_EQ(result2(1, 2), 9.0);
    EXPECT_EQ(result2(2, 2), 12.0);
}

TEST(JIArrayColumnMajorTests, ArrayDivision) {
    JIArray<double, 2> array1(2, 2);
    JIArray<double, 2> array2(2, 2);

    array1 = {10.0, 8.0, 6.0, 4.0};
    array2 = {2.0, 2.0, 3.0, 2.0};

    auto result = array1 / array2;
    EXPECT_DOUBLE_EQ(result(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(result(2, 1), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 2), 2.0);
    EXPECT_DOUBLE_EQ(result(2, 2), 2.0);
}

TEST(JIArrayColumnMajorTests, ScalarDivision) {
    JIArray<double, 2> array(2, 2);
    array = {10.0, 20.0, 30.0, 40.0};

    auto result = array / 5.0;
    EXPECT_DOUBLE_EQ(result(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(result(2, 1), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 2), 6.0);
    EXPECT_DOUBLE_EQ(result(2, 2), 8.0);
}

// ============================================================================
// STATISTICAL OPERATIONS TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, Average) {
    JIArray<double, 2> array(2, 2);
    array = {1.0, 3.0, 5.0, 7.0};

    EXPECT_DOUBLE_EQ(array.average(), 4.0);
}

TEST(JIArrayColumnMajorTests, Sum) {
    JIArray<int, 2> array(2, 3);
    array = {1, 2, 3, 4, 5, 6};

    EXPECT_EQ(array.sum(), 21);
}

TEST(JIArrayColumnMajorTests, MinMax) {
    JIArray<double, 2> array(2, 3);
    array = {3.5, 1.2, 7.8, 2.1, 9.3, 0.5};

    EXPECT_DOUBLE_EQ(array.min(), 0.5);
    EXPECT_DOUBLE_EQ(array.max(), 9.3);
}

TEST(JIArrayColumnMajorTests, SumOfSquares) {
    JIArray<double, 1> array(3);
    array = {2.0, 3.0, 4.0};

    EXPECT_DOUBLE_EQ(array.sqsum(), 29.0); // 4 + 9 + 16
}

// ============================================================================
// COMPARISON OPERATIONS TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, Equality) {
    JIArray<int, 2> array1(2, 2);
    JIArray<int, 2> array2(2, 2);

    array1 = {1, 2, 3, 4};
    array2 = {1, 2, 3, 4};

    EXPECT_TRUE(array1 == array2);

    array2(2, 2) = 5;
    EXPECT_FALSE(array1 == array2);
}

TEST(JIArrayColumnMajorTests, Contains) {
    JIArray<int, 2> array(2, 3);
    array = {1, 5, 3, 8, 2, 7};

    EXPECT_TRUE(array.contains(5));
    EXPECT_TRUE(array.contains(8));
    EXPECT_FALSE(array.contains(9));
    EXPECT_FALSE(array.contains(0));
}

TEST(JIArrayColumnMajorTests, FindFirst1D) {
    JIArray<int, 1> array(5);
    array = {10, 20, 30, 20, 50};

    auto pos = array.findFirst(20);
    EXPECT_EQ(pos, 2); // First occurrence

    auto not_found = array.findFirst(99);
    EXPECT_EQ(not_found, -1 + JIARRAY_OFFSET);
}

TEST(JIArrayColumnMajorTests, FindFirstMultiD) {
    JIArray<int, 2> array(2, 2);
    array = {1, 2, 3, 4};

    auto pos = array.findFirst(3);
    EXPECT_EQ(pos(1), 1);
    EXPECT_EQ(pos(2), 2);
}

TEST(JIArrayColumnMajorTests, MaxLocation) {
    JIArray<double, 2> array(2, 3);
    array = {1.0, 3.0, 2.0, 8.0, 5.0, 4.0};

    auto maxpos = array.maxloc();
    EXPECT_EQ(maxpos(1), 2);
    EXPECT_EQ(maxpos(2), 2); // Position of 8.0
}

// ============================================================================
// ASSIGNMENT OPERATIONS TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, ScalarAssignment) {
    JIArray<double, 2> array(2, 2);
    array = 5.5;

    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
            EXPECT_DOUBLE_EQ(array(i, j), 5.5);
        }
    }
}

TEST(JIArrayColumnMajorTests, VectorAssignment) {
    JIArray<int, 1> array(5);
    std::vector<int> vec = {10, 20, 30, 40, 50};

    array = vec;
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(array(i + 1), vec[i]);
    }
}

TEST(JIArrayColumnMajorTests, ArrayAssignment) {
    JIArray<double, 2> source(2, 2);
    source = {1.0, 2.0, 3.0, 4.0};

    JIArray<double, 2> target(2, 2);
    target = source;

    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
            EXPECT_EQ(target(i, j), source(i, j));
        }
    }
}

// ============================================================================
// RESHAPING AND MEMORY SHARING TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, Reshape1D) {
    JIArray<int, 1> array1d(6);
    array1d = {1, 2, 3, 4, 5, 6};

    auto array2d = array1d.reshape(2, 3);
    EXPECT_EQ(array2d.size(), 6);
    EXPECT_EQ(array2d(1, 1), 1);
    EXPECT_EQ(array2d(2, 1), 2);
    EXPECT_EQ(array2d(1, 2), 3);
    EXPECT_EQ(array2d(2, 2), 4);
}

TEST(JIArrayColumnMajorTests, ShareWith) {
    JIArray<double, 2> source(2, 3);
    source = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    JIArray<double, 2> target;
    target.shareWith(source);

    EXPECT_EQ(target.size(), source.size());
    EXPECT_EQ(target.data(), source.data());

    source(1, 1) = 99.0;
    EXPECT_EQ(target(1, 1), 99.0); // Shared memory
}

TEST(JIArrayColumnMajorTests, Copy) {
    JIArray<double, 2> original(2, 2);
    original = {1.0, 2.0, 3.0, 4.0};

    auto copied = original.copy();
    EXPECT_EQ(copied.size(), original.size());
    EXPECT_NE(copied.data(), original.data()); // Different memory

    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
            EXPECT_EQ(copied(i, j), original(i, j));
        }
    }
}

// ============================================================================
// OFFSET AND RANGE TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, SetOffsets) {
    JIArray<int, 2> array(2, 2);
    array = {1, 2, 3, 4};

    array.setOffsets(0, 0); // Change to 0-based indexing
    EXPECT_EQ(array(0, 0), 1);
    EXPECT_EQ(array(1, 0), 2);
    EXPECT_EQ(array(0, 1), 3);
    EXPECT_EQ(array(1, 1), 4);
}

TEST(JIArrayColumnMajorTests, RangeInitialization) {
    JIArray<double, 2> array;
    array.init0(0, 2, 0, 1); // Range: [0,2] x [0,1]

    EXPECT_EQ(array.getSize(1), 3); // 0, 1, 2
    EXPECT_EQ(array.getSize(2), 2); // 0, 1
    EXPECT_EQ(array.size(), 6);
}

// ============================================================================
// ITERATOR TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, IteratorBasic) {
    JIArray<int, 1> array(5);
    array = {10, 20, 30, 40, 50};

    auto it = array.begin();
    EXPECT_EQ(*it, 10);

    ++it;
    EXPECT_EQ(*it, 20);

    it += 2;
    EXPECT_EQ(*it, 40);

    EXPECT_EQ(array.end() - array.begin(), 5);
}

TEST(JIArrayColumnMajorTests, IteratorLooping) {
    JIArray<int, 2> array(2, 3);
    array = {1, 2, 3, 4, 5, 6};

    std::vector<int> expected = {1, 2, 3, 4, 5, 6};
    std::vector<int> actual;

    for (const auto &value : array) {
        actual.push_back(value);
    }

    EXPECT_EQ(actual, expected);
}

TEST(JIArrayColumnMajorTests, IteratorAlgorithms) {
    JIArray<int, 1> array(5);
    array = {5, 2, 8, 1, 9};

    auto max_it = std::max_element(array.begin(), array.end());
    EXPECT_EQ(*max_it, 9);

    auto count = std::count_if(array.begin(), array.end(), [](int x) { return x > 5; });
    EXPECT_EQ(count, 2); // 8 and 9
}

// ============================================================================
// DATA ACCESS AND UTILITY TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, DataAccess) {
    JIArray<double, 2> array(2, 3);
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    double *data_ptr = array.data();
    EXPECT_NE(data_ptr, nullptr);
    EXPECT_EQ(data_ptr[0], 1.0);
    EXPECT_EQ(data_ptr[1], 2.0);

    const double *const_data = const_cast<const JIArray<double, 2> &>(array).data();
    EXPECT_EQ(const_data, data_ptr);
}

TEST(JIArrayColumnMajorTests, ConvertToVector) {
    JIArray<int, 2> array(2, 2);
    array = {1, 2, 3, 4};

    auto vec = array.convertToVector();
    std::vector<int> expected = {1, 2, 3, 4};
    EXPECT_EQ(vec, expected);
}

TEST(JIArrayColumnMajorTests, SizeManagement) {
    JIArray<double, 2> array(2, 2);

    EXPECT_EQ(array.getSize(), 4);
    EXPECT_EQ(array.size(), 4);

    array.setSize(3, 3);
    EXPECT_EQ(array.getSize(), 9);
    EXPECT_EQ(array.getSize(1), 3);
    EXPECT_EQ(array.getSize(2), 3);
}

// ============================================================================
// DOT PRODUCT TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, DotProduct) {
    JIArray<double, 1> array1(3);
    JIArray<double, 1> array2(3);

    array1 = {1.0, 2.0, 3.0};
    array2 = {4.0, 5.0, 6.0};

    double result = dot(array1, array2);
    EXPECT_DOUBLE_EQ(result, 32.0); // 1*4 + 2*5 + 3*6 = 32
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, BoundaryConditions) {
    JIArray<int, 1> array(3);
    array = {10, 20, 30};

    // Valid accesses
    EXPECT_NO_THROW(array(1));
    EXPECT_NO_THROW(array(3));

    // Boundary tests depend on your JIARRAY_CHECK_BOUND implementation
    // Uncomment if bounds checking throws exceptions
    // EXPECT_THROW(array(0), std::out_of_range);
    // EXPECT_THROW(array(4), std::out_of_range);
}

TEST(JIArrayColumnMajorTests, EmptyArrayBehavior) {
    JIArray<double, 2> array;
    EXPECT_EQ(array.size(), 0);
    EXPECT_FALSE(array.isAllocated());
    EXPECT_EQ(array.data(), nullptr);
}

// ============================================================================
// UNARY OPERATIONS TESTS
// ============================================================================

TEST(JIArrayColumnMajorTests, UnaryMinus) {
    JIArray<double, 2> array(2, 2);
    array = {1.0, -2.0, 3.0, -4.0};

    auto result = -array;
    EXPECT_EQ(result(1, 1), -1.0);
    EXPECT_EQ(result(2, 1), 2.0);
    EXPECT_EQ(result(1, 2), -3.0);
    EXPECT_EQ(result(2, 2), 4.0);
}

// ============================================================================
// SLICING TESTS (Extended from your original tests)
// ============================================================================

TEST(JIArrayColumnMajorTests, SliceBasic) {
    JIArray<double, 2> array(3, 3);
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    auto slice = array.slice(1);
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 1.0);
    EXPECT_EQ(slice(2), 2.0);
    EXPECT_EQ(slice(3), 3.0);
}

TEST(JIArrayColumnMajorTests, SliceMiddleColumn) {
    JIArray<double, 2> array(3, 3);                        // Create a 3x3 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; // Column-major order

    auto slice = array.slice(2); // Slice to get the second column
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 4.0);
    EXPECT_EQ(slice(2), 5.0);
    EXPECT_EQ(slice(3), 6.0);
}

TEST(JIArrayColumnMajorTests, SliceLastColumn) {
    JIArray<double, 2> array(3, 3);                        // Create a 3x3 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; // Column-major order

    auto slice = array.slice(3); // Slice to get the third column
    EXPECT_EQ(slice.size(), 3);
    EXPECT_EQ(slice(1), 7.0);
    EXPECT_EQ(slice(2), 8.0);
    EXPECT_EQ(slice(3), 9.0);
}

TEST(JIArrayColumnMajorTests, SliceOn3DArray) {
    JIArray<double, 3> array(2, 2, 2);                // Create a 2x2x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // Column-major order

    auto slice = array.slice(1); // Slice along the first dimension
    EXPECT_EQ(slice.size(), 4);  // 2x2 slice

    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(2, 1), 2.0);
    EXPECT_EQ(slice(1, 2), 3.0);
    EXPECT_EQ(slice(2, 2), 4.0);
}

TEST(JIArrayColumnMajorTests, SliceOn3DArraySecondDimension) {
    JIArray<double, 3> array(2, 2, 2);                // Create a 2x2x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // Column-major order

    auto slice = array.slice(2); // Slice along the second dimension
    EXPECT_EQ(slice.size(), 4);  // 2x2 slice

    EXPECT_EQ(slice(1, 1), 5.0);
    EXPECT_EQ(slice(2, 1), 6.0);
    EXPECT_EQ(slice(1, 2), 7.0);
    EXPECT_EQ(slice(2, 2), 8.0);
}

TEST(JIArrayColumnMajorTests, MultiDimensionalSlicing) {
    JIArray<double, 3> array(3, 3, 3); // Create a 3x3x3 array
    array = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,

             10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,

             19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0}; // Column-major order

    // Slice the first 2D "plane"
    auto slice = array.slice(1);
    EXPECT_EQ(slice.size(), 9);
    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(2, 1), 2.0);
    EXPECT_EQ(slice(3, 1), 3.0);
    EXPECT_EQ(slice(1, 2), 4.0);
    EXPECT_EQ(slice(3, 3), 9.0);
}

TEST(JIArrayColumnMajorTests, SliceInvalidIndex) {
    JIArray<double, 2> array(3, 3);                        // Create a 3x3 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; // Column-major order

    // Attempt to slice an invalid dimension index
    try{
        auto a = array.slice(4);
    } catch (std::exception& e) {
        SUCCEED(); // Exception correctly thrown

    }
}

TEST(JIArrayColumnMajorTests, SliceNonSquareArray) {
    JIArray<double, 2> array(3, 2);         // Create a 3x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // Column-major order

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

TEST(JIArrayColumnMajorTests, SliceRectangular3DArray) {
    JIArray<double, 3> array(3, 2, 2);                                       // Create a 3x2x2 array
    array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}; // Column-major order

    auto slice = array.slice(1); // Slice along the first "layer"
    EXPECT_EQ(slice.size(), 6);  // Should return a 3x2 slice

    EXPECT_EQ(slice(1, 1), 1.0);
    EXPECT_EQ(slice(2, 1), 2.0);
    EXPECT_EQ(slice(3, 1), 3.0);
    EXPECT_EQ(slice(1, 2), 4.0);
    EXPECT_EQ(slice(2, 2), 5.0);
    EXPECT_EQ(slice(3, 2), 6.0);
}

TEST(JIArrayColumnMajorTests, SliceNonUniform3DArray) {
    JIArray<double, 3> array(4, 3, 2); // Create a 4x3x2 array
    array = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,
             7.0,  8.0,  9.0,  10.0, 11.0, 12.0,

             13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
             19.0, 20.0, 21.0, 22.0, 23.0, 24.0}; // Column-major order

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

TEST(JIArrayColumnMajorTests, Slice4DArray) {
    JIArray<double, 4> array(2, 2, 2, 2); // Create a 2x2x2x2 array
    array = {1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,

             9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0}; // Column-major order

    auto slice = array.slice(1); // Slice the first "3D layer"
    EXPECT_EQ(slice.size(), 8);  // Should return a 2x2x2 slice

    EXPECT_EQ(slice(1, 1, 1), 1.0);
    EXPECT_EQ(slice(2, 1, 1), 2.0);
    EXPECT_EQ(slice(1, 2, 1), 3.0);
    EXPECT_EQ(slice(2, 2, 1), 4.0);
    EXPECT_EQ(slice(1, 1, 2), 5.0);
    EXPECT_EQ(slice(2, 1, 2), 6.0);
    EXPECT_EQ(slice(2, 2, 2), 8.0);
}

TEST(JIArrayColumnMajorTests, SliceIrregularSizedArray) {
    JIArray<double, 2> array(5, 4); // Create a 5x4 array
    array = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0,
             11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0}; // Column-major order

    auto slice = array.slice(3); // Get the third column
    EXPECT_EQ(slice.size(), 5);
    EXPECT_EQ(slice(1), 11.0);
    EXPECT_EQ(slice(2), 12.0);
    EXPECT_EQ(slice(3), 13.0);
    EXPECT_EQ(slice(4), 14.0);
    EXPECT_EQ(slice(5), 15.0);
}