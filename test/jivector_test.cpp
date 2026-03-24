#include <jiarray/JIVector.h>
#include <gtest/gtest.h>
#include <string>

using namespace dnegri::jiarray;

// ============================================================================
// JIVector Construction
// ============================================================================

TEST(JIVectorTest, DefaultConstructorEmpty) {
    JIVector<int> v;
    EXPECT_EQ(v.size(), 0u);
}

TEST(JIVectorTest, InitializerListConstructor) {
    JIVector<int> v{10, 20, 30};
    EXPECT_EQ(v.size(), 3u);
    EXPECT_EQ(v[1], 10);
    EXPECT_EQ(v[2], 20);
    EXPECT_EQ(v[3], 30);
}

TEST(JIVectorTest, FromStdVector) {
    std::vector<double> sv{1.5, 2.5, 3.5};
    JIVector<double> v(sv);
    EXPECT_EQ(v.size(), 3u);
    EXPECT_DOUBLE_EQ(v[1], 1.5);
    EXPECT_DOUBLE_EQ(v[3], 3.5);
}

// ============================================================================
// JIVector 1-based Access
// ============================================================================

TEST(JIVectorTest, BracketOperator1Based) {
    JIVector<int> v{100, 200, 300};
    // operator[] is 1-based
    EXPECT_EQ(v[1], 100);
    EXPECT_EQ(v[2], 200);
    EXPECT_EQ(v[3], 300);
    v[2] = 999;
    EXPECT_EQ(v[2], 999);
}

TEST(JIVectorTest, ConstBracketOperator) {
    const JIVector<int> v{10, 20, 30};
    EXPECT_EQ(v[1], 10);
    EXPECT_EQ(v[3], 30);
}

TEST(JIVectorTest, AtMethod1Based) {
    JIVector<int> v{5, 10, 15};
    EXPECT_EQ(v.at(1), 5);
    EXPECT_EQ(v.at(3), 15);
    v.at(2) = 42;
    EXPECT_EQ(v.at(2), 42);
}

TEST(JIVectorTest, ConstAtMethod) {
    const JIVector<int> v{5, 10, 15};
    EXPECT_EQ(v.at(1), 5);
    EXPECT_EQ(v.at(3), 15);
}

TEST(JIVectorTest, AtMethodOutOfRange) {
    JIVector<int> v{1, 2, 3};
    EXPECT_THROW(v.at(0), std::out_of_range);
    EXPECT_THROW(v.at(4), std::out_of_range);
}

TEST(JIVectorTest, ParenthesisOperator) {
    JIVector<int> v{7, 8, 9};
    EXPECT_EQ(v(1), 7);
    EXPECT_EQ(v(3), 9);
    v(2) = 77;
    EXPECT_EQ(v(2), 77);
}

TEST(JIVectorTest, ConstParenthesisOperator) {
    const JIVector<int> v{7, 8, 9};
    EXPECT_EQ(v(1), 7);
    EXPECT_EQ(v(3), 9);
}

// ============================================================================
// JIVector Memory Access
// ============================================================================

TEST(JIVectorTest, GetMemory) {
    JIVector<int> v{10, 20, 30};
    int* p = v.getMemory(0);
    EXPECT_EQ(*p, 10);
    p = v.getMemory(2);
    EXPECT_EQ(*p, 30);
}

TEST(JIVectorTest, GetPointer) {
    JIVector<int> v{1, 2, 3};
    int* p = v.get_pointer();
    EXPECT_EQ(p[0], 1);
    EXPECT_EQ(p[1], 2);
    EXPECT_EQ(p[2], 3);
}

// ============================================================================
// JIVector Assignment
// ============================================================================

TEST(JIVectorTest, AssignFromStdVector) {
    JIVector<int> v;
    std::vector<int> sv{4, 5, 6};
    v = sv;
    EXPECT_EQ(v.size(), 3u);
    EXPECT_EQ(v[1], 4);
    EXPECT_EQ(v[3], 6);
}

TEST(JIVectorTest, AssignFromJIArray1D) {
    zint1 arr{10, 20, 30};
    JIVector<int> v;
    v = arr;
    EXPECT_EQ(v.size(), 3u);
    EXPECT_EQ(v[1], 10);
    EXPECT_EQ(v[2], 20);
    EXPECT_EQ(v[3], 30);
}

// ============================================================================
// JIVector Search
// ============================================================================

TEST(JIVectorTest, ContainsFound) {
    JIVector<int> v{10, 20, 30, 40};
    EXPECT_TRUE(v.contains(20));
    EXPECT_TRUE(v.contains(40));
}

TEST(JIVectorTest, ContainsNotFound) {
    JIVector<int> v{10, 20, 30};
    EXPECT_FALSE(v.contains(99));
}

TEST(JIVectorTest, FindFound) {
    JIVector<int> v{10, 20, 30, 20};
    // find returns 1-based index
    EXPECT_EQ(v.find(10), 1);
    EXPECT_EQ(v.find(20), 2); // first occurrence
    EXPECT_EQ(v.find(30), 3);
}

TEST(JIVectorTest, FindNotFound) {
    JIVector<int> v{1, 2, 3};
    // not found returns JIARRAY_OFFSET - 1 = 0
    EXPECT_EQ(v.find(99), 0);
}

// ============================================================================
// JIVector Insert
// ============================================================================

TEST(JIVectorTest, InsertInitializerList) {
    JIVector<int> v{1, 2, 3};
    v.insert(v.end(), {4, 5});
    EXPECT_EQ(v.size(), 5u);
    EXPECT_EQ(v[4], 4);
    EXPECT_EQ(v[5], 5);
}

TEST(JIVectorTest, InsertFromJIVector) {
    JIVector<int> v{1, 2};
    JIVector<int> other{3, 4};
    v.insert(v.end(), other);
    EXPECT_EQ(v.size(), 4u);
    EXPECT_EQ(v[1], 1);
    EXPECT_EQ(v[3], 3);
    EXPECT_EQ(v[4], 4);
}

TEST(JIVectorTest, InsertFromStdVector) {
    JIVector<int> v{10};
    std::vector<int> sv{20, 30};
    v.insert(v.end(), sv);
    EXPECT_EQ(v.size(), 3u);
    EXPECT_EQ(v[1], 10);
    EXPECT_EQ(v[2], 20);
    EXPECT_EQ(v[3], 30);
}

TEST(JIVectorTest, InsertFromJIArray) {
    JIVector<int> v{1};
    zint1 arr{2, 3};
    v.insert(v.end(), arr);
    EXPECT_EQ(v.size(), 3u);
    EXPECT_EQ(v[2], 2);
    EXPECT_EQ(v[3], 3);
}

TEST(JIVectorTest, InsertIteratorRange) {
    JIVector<int> v{1, 2};
    JIVector<int> src{10, 20, 30};
    v.insert(v.end(), src.begin(), src.begin() + 2);
    EXPECT_EQ(v.size(), 4u);
    EXPECT_EQ(v[3], 10);
    EXPECT_EQ(v[4], 20);
}

TEST(JIVectorTest, InsertAtBeginning) {
    JIVector<int> v{3, 4};
    v.insert(v.begin(), {1, 2});
    EXPECT_EQ(v.size(), 4u);
    EXPECT_EQ(v[1], 1);
    EXPECT_EQ(v[2], 2);
    EXPECT_EQ(v[3], 3);
    EXPECT_EQ(v[4], 4);
}

// ============================================================================
// JIVector Iteration
// ============================================================================

TEST(JIVectorTest, BeginEnd) {
    JIVector<int> v{10, 20, 30};
    std::vector<int> collected(v.begin(), v.end());
    EXPECT_EQ(collected, (std::vector<int>{10, 20, 30}));
}

TEST(JIVectorTest, ConstBeginEnd) {
    const JIVector<int> v{10, 20, 30};
    int sum = 0;
    for (auto it = v.begin(); it != v.end(); ++it) {
        sum += *it;
    }
    EXPECT_EQ(sum, 60);
}

TEST(JIVectorTest, RangeBasedFor) {
    JIVector<int> v{1, 2, 3, 4};
    int sum = 0;
    for (auto x : v) {
        sum += x;
    }
    EXPECT_EQ(sum, 10);
}

// ============================================================================
// JIVector String specialization
// ============================================================================

TEST(JIVectorTest, StringVector) {
    JIVector<std::string> v{"hello", "world"};
    EXPECT_EQ(v[1], "hello");
    EXPECT_EQ(v[2], "world");
    EXPECT_TRUE(v.contains("hello"));
    EXPECT_FALSE(v.contains("missing"));
    EXPECT_EQ(v.find("world"), 2);
}

// ============================================================================
// Type alias tests
// ============================================================================

TEST(JIVectorAliasTest, ZVectorAlias) {
    zvector<int> v{1, 2, 3};
    EXPECT_EQ(v[1], 1);
    EXPECT_EQ(v.size(), 3u);
}

TEST(JIVectorAliasTest, ZIntsAlias) {
    zints v{10, 20, 30};
    EXPECT_EQ(v[2], 20);
}

TEST(JIVectorAliasTest, ZDoublesAlias) {
    zdoubles v{1.5, 2.5};
    EXPECT_DOUBLE_EQ(v[1], 1.5);
}

TEST(JIVectorAliasTest, ZStringsAlias) {
    zstrings v{"a", "b"};
    EXPECT_EQ(v[1], "a");
}
