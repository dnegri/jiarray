// Tests for the ffor / ffor_back / zfor macros — in particular the
// optional `step` argument, which used to silently expand to 1 because
// of a position-counting bug in the GET_STEP macro.

#include <jiarray/JIArray.h>
#include <gtest/gtest.h>

TEST(Ffor, DefaultStep_OneBased) {
    int sum = 0;
    ffor(i, 1, 5) sum += i;
    EXPECT_EQ(sum, 15);  // 1+2+3+4+5
}

TEST(Ffor, ExplicitStepOne_OneBased) {
    int sum = 0;
    ffor(i, 1, 5, 1) sum += i;
    EXPECT_EQ(sum, 15);
}

TEST(Ffor, ExplicitStepTwo_OneBased) {
    // Pre-fix bug: the user-supplied "2" was dropped and the loop
    // walked at step 1, producing 2+3+4+5+6+7+8+9+10 = 54.
    int sum = 0;
    ffor(i, 2, 10, 2) sum += i;
    EXPECT_EQ(sum, 30);  // 2+4+6+8+10
}

TEST(Ffor, NegativeBoundsStepTwo) {
    int sum = 0;
    ffor(i, -4, 4, 2) sum += i;
    EXPECT_EQ(sum, 0);  // -4-2+0+2+4
}

TEST(Ffor, ExpressionInBounds) {
    int n = 5;
    int sum = 0;
    ffor(i, 1, n - 1) sum += i;
    EXPECT_EQ(sum, 10);  // 1+2+3+4
}

TEST(FforBack, DefaultStep) {
    int sum = 0;
    ffor_back(i, 5, 1) sum += i;
    EXPECT_EQ(sum, 15);
}

TEST(FforBack, ExplicitStepTwo) {
    int sum = 0;
    ffor_back(i, 10, 2, 2) sum += i;
    EXPECT_EQ(sum, 30);  // 10+8+6+4+2
}

TEST(Zfor, DefaultStep_OneBased) {
    int sum = 0;
    zfor(i, 5) sum += i;
    EXPECT_EQ(sum, 15);  // 1..5 inclusive in JIARRAY_OFFSET=1 mode
}

TEST(Ffor, IterationCountMatchesStep) {
    // Sanity: the number of iterations actually equals
    // floor((end - begin) / step) + 1 for inclusive ranges.
    int count = 0;
    ffor(i, 0, 100, 7) ++count;
    EXPECT_EQ(count, 15);  // 0,7,14,...,98 → 15 values
}
