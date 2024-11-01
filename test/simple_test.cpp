#include <gtest/gtest.h>
#include "../src/JIArray.h"  // Ensure this path points to your JIArray class

using namespace dnegri::jiarray;

TEST(JIArrayTestsWithOffset, InitializationAndAccess) {
    double2 array(2, 2);  // Initialize with dimensions
    array = {1.0, 2.0, 3.0, 4.0};  // Use assignment operator to set values

    EXPECT_EQ(array(1, 1), 1.0);  // Column 1, Row 1
    EXPECT_EQ(array(2, 1), 2.0);  // Column 2, Row 1
    EXPECT_EQ(array(1, 2), 3.0);  // Column 1, Row 2
    EXPECT_EQ(array(2, 2), 4.0);  // Column 2, Row 2
}

TEST(JIArrayTestsWithOffset, AssignmentOperator) {
    double2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    double2 array2(2, 2);
    array2 = array1;
    EXPECT_EQ(array2(1, 1), 1.0);
    EXPECT_EQ(array2(2, 1), 2.0);
    EXPECT_EQ(array2(1, 2), 3.0);
    EXPECT_EQ(array2(2, 2), 4.0);
}

TEST(JIArrayTestsWithOffset, CopyConstructor) {
    double2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    double2 array2(array1);
    EXPECT_EQ(array2(1, 1), 1.0);
    EXPECT_EQ(array2(2, 1), 2.0);
    EXPECT_EQ(array2(1, 2), 3.0);
    EXPECT_EQ(array2(2, 2), 4.0);
}

TEST(JIArrayTestsWithOffset, Slicing) {
    JIArray<double, 3> array(3, 3, 3);
    array = 1.0;
    auto slice = array.slice(1);
    EXPECT_EQ(slice.size(), 9);
}

TEST(JIArrayTestsWithOffset, ArithmeticOperations) {
    double2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    double2 array2(2, 2);
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

TEST(JIArrayTestsWithOffset, Reshape) {
    double2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    auto reshaped = array.reshape(4);
    EXPECT_EQ(reshaped.size(), 4);
    EXPECT_EQ(reshaped(1), 1.0);
    EXPECT_EQ(reshaped(2), 2.0);
    EXPECT_EQ(reshaped(3), 3.0);
    EXPECT_EQ(reshaped(4), 4.0);
}

TEST(JIArrayTestsWithOffset, ContainsAndFind) {
    double2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    EXPECT_TRUE(array.contains(3.0));
    EXPECT_FALSE(array.contains(5.0));

    auto loc = array.findFirst(3.0);
    EXPECT_EQ(loc(1), 1);
    EXPECT_EQ(loc(2), 2);
}

TEST(JIArrayTestsWithOffset, MaxLoc) {
    double2 array(2, 2);
    array = {1.0, 5.0, 3.0, 4.0};
    auto maxLoc = array.maxloc();
    EXPECT_EQ(maxLoc(1), 2);
    EXPECT_EQ(maxLoc(2), 1);
}

TEST(JIArrayTestsWithOffset, AverageSum) {
    double2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    EXPECT_DOUBLE_EQ(array.average(), 2.5);
    EXPECT_DOUBLE_EQ(array.sum(), 10.0);
}

TEST(JIArrayTestsWithOffset, IteratorTest) {
    double2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    double sum = 0.0;
    for (auto it = array.begin(); it != array.end(); ++it) {
        sum += *it;
    }
    EXPECT_DOUBLE_EQ(sum, 10.0);
}

TEST(JIArrayTestsWithOffset, ShareWith) {
    double2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    double2 array2;
    array2.shareWith(array1);

    EXPECT_EQ(array2(1, 1), 1.0);
    EXPECT_EQ(array2(2, 1), 2.0);
    EXPECT_EQ(array2(1, 2), 3.0);
    EXPECT_EQ(array2(2, 2), 4.0);
}

TEST(JIArrayTestsWithOffset, ConversionToVector) {
    double2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> vec = array.convertToVector();
    EXPECT_EQ(vec.size(), 4);
    EXPECT_EQ(vec[0], 1.0);
    EXPECT_EQ(vec[1], 2.0);
    EXPECT_EQ(vec[2], 3.0);
    EXPECT_EQ(vec[3], 4.0);
}

TEST(JIArrayTestsWithOffset, ScalarMultiplication) {
    double2 array(2, 2);
    array = {1.0, 2.0, 3.0, 4.0};
    auto result = 2.0 * array;
    EXPECT_EQ(result(1, 1), 2.0);
    EXPECT_EQ(result(2, 1), 4.0);
    EXPECT_EQ(result(1, 2), 6.0);
    EXPECT_EQ(result(2, 2), 8.0);
}

TEST(JIArrayTestsWithOffset, DotProduct) {
    double2 array1(2, 2);
    array1 = {1.0, 2.0, 3.0, 4.0};
    double2 array2(2, 2);
    array2 = {5.0, 6.0, 7.0, 8.0};
    double dot_product = dot(array1, array2);
    EXPECT_DOUBLE_EQ(dot_product, 70.0);
}