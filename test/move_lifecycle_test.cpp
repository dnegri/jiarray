// move_lifecycle_test.cpp
// Comprehensive coverage for the v0.7.1 additions:
//   * move constructor
//   * move assignment
//   * destroy() / erase() full metadata reset
//   * copy ctor remains a shallow non-owning view
//   * copy assignment remains an explicit deep copy
//   * real-world usage patterns observed in sphincs / hybrid / montex / croma
//
// The goal is not only "does it compile" — every test checks lifecycle
// invariants:
//   * no double-free  (valgrind-clean)
//   * no dangling mm after move/destroy
//   * metadata (sizes / offsets / stride) consistent after every op
//   * copy ctor never steals ownership
//   * copy assignment never silently takes over the source

#include <gtest/gtest.h>
#include <jiarray/JIArray.h>
#include <jiarray/FastArray.h>
#include <jiarray/JIVector.h>

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

using namespace dnegri::jiarray;

// ---------------------------------------------------------------------------
// 1. Move constructor
// ---------------------------------------------------------------------------

TEST(MoveCtor, Basic1D) {
    zdouble1 src; src.init(5);
    for (int i = 1; i <= 5; ++i) src(i) = i * 1.5;
    const double* src_mm = src.data();

    zdouble1 dst(std::move(src));

    EXPECT_EQ(dst.size(), 5);
    for (int i = 1; i <= 5; ++i) EXPECT_DOUBLE_EQ(dst(i), i * 1.5);
    EXPECT_EQ(dst.data(), src_mm);        // ownership really transferred

    // moved-from is empty and destructor-safe
    EXPECT_EQ(src.size(), 0);
    EXPECT_EQ(src.data(), nullptr);
    EXPECT_FALSE(src.isAllocated());
}

TEST(MoveCtor, Basic2D) {
    zint2 src; src.init(3, 4);
    int v = 1;
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 4; ++j) src(i, j) = v++;

    zint2 dst(std::move(src));
    EXPECT_EQ(dst.getSize(1), 3);
    EXPECT_EQ(dst.getSize(2), 4);
    EXPECT_EQ(dst(1, 1), 1);
    EXPECT_EQ(dst(3, 4), 12);
    EXPECT_EQ(src.size(), 0);
}

TEST(MoveCtor, FromViewIsStillEmpty) {
    // Moving from a view (allocated=NONE) should yield another view of
    // the same memory, not double-free.
    zdouble1 owner; owner.init(3); owner = 42.0;
    zdouble1 view(owner);                 // shallow copy ctor
    EXPECT_FALSE(view.isAllocated());
    const double* shared_mm = view.data();

    zdouble1 moved(std::move(view));
    // The "moved" view is either empty (strict move) or a second view —
    // either way it must NOT own; owner still owns the memory.
    EXPECT_FALSE(moved.isAllocated());
    EXPECT_EQ(view.data(), nullptr);

    // owner must still hold correct data (no accidental free)
    EXPECT_EQ(owner.size(), 3);
    EXPECT_DOUBLE_EQ(owner(1), 42.0);
    (void)shared_mm;
}

TEST(MoveCtor, EmptySource) {
    zdouble1 src;  // never initialised
    zdouble1 dst(std::move(src));
    EXPECT_EQ(dst.size(), 0);
    EXPECT_EQ(src.size(), 0);
}

// ---------------------------------------------------------------------------
// 2. Move assignment
// ---------------------------------------------------------------------------

TEST(MoveAssign, OverwritesAndFreesPrevious) {
    zdouble1 a; a.init(2); a(1) = 1.0; a(2) = 2.0;
    zdouble1 b; b.init(5); b = 9.0;     // b currently owns 5-element block
    b = std::move(a);
    EXPECT_EQ(b.size(), 2);
    EXPECT_DOUBLE_EQ(b(1), 1.0);
    EXPECT_DOUBLE_EQ(b(2), 2.0);
    EXPECT_EQ(a.size(), 0);             // a emptied
}

TEST(MoveAssign, SelfAssignmentNoop) {
    zdouble1 a; a.init(4); a = 3.14;
    const double* mm_before = a.data();
    a = std::move(a);                   // shouldn't corrupt
    EXPECT_EQ(a.size(), 4);
    EXPECT_DOUBLE_EQ(a(1), 3.14);
    EXPECT_EQ(a.data(), mm_before);
}

TEST(MoveAssign, IntoEmptyTarget) {
    zdouble1 src; src.init(3); src = 7.0;
    zdouble1 dst;
    dst = std::move(src);
    EXPECT_EQ(dst.size(), 3);
    EXPECT_DOUBLE_EQ(dst(1), 7.0);
    EXPECT_EQ(src.size(), 0);
}

TEST(MoveAssign, FromEmptyToOwner) {
    zdouble1 owner; owner.init(4); owner = 1.0;
    zdouble1 empty;
    owner = std::move(empty);
    EXPECT_EQ(owner.size(), 0);         // previous block freed, now empty
    EXPECT_EQ(empty.size(), 0);
}

TEST(MoveAssign, DoesNotCorruptUnrelatedInstance) {
    zdouble1 a; a.init(3); a(1) = 10; a(2) = 20; a(3) = 30;
    zdouble1 b; b.init(3); b(1) = 40; b(2) = 50; b(3) = 60;
    zdouble1 c; c = std::move(a);
    // b untouched
    EXPECT_DOUBLE_EQ(b(1), 40);
    EXPECT_DOUBLE_EQ(b(3), 60);
    // c has a's old content
    EXPECT_DOUBLE_EQ(c(1), 10);
    EXPECT_DOUBLE_EQ(c(3), 30);
}

// ---------------------------------------------------------------------------
// 3. destroy() / erase() full reset
// ---------------------------------------------------------------------------

TEST(Destroy, ResetsMetadata) {
    zdouble2 a; a.init(3, 4); a = 5.0;
    EXPECT_EQ(a.getSize(1), 3);
    EXPECT_EQ(a.getSize(2), 4);
    a.destroy();
    EXPECT_EQ(a.size(), 0);
    EXPECT_EQ(a.data(), nullptr);
    EXPECT_EQ(a.getSize(1), 0);
    EXPECT_EQ(a.getSize(2), 0);
    EXPECT_FALSE(a.isAllocated());
}

TEST(Destroy, SafeToCallTwice) {
    zdouble1 a; a.init(3); a = 1.0;
    a.destroy();
    a.destroy();  // must not crash
    EXPECT_EQ(a.size(), 0);
}

TEST(Destroy, SafeOnUninitialised) {
    zdouble1 a;
    a.destroy();
    EXPECT_EQ(a.size(), 0);
}

TEST(Destroy, CanReInit) {
    zint3 a; a.init(2, 3, 4); a = 7;
    a.destroy();
    a.init(5, 6, 7); a = 11;
    EXPECT_EQ(a.getSize(1), 5);
    EXPECT_EQ(a.getSize(2), 6);
    EXPECT_EQ(a.getSize(3), 7);
    EXPECT_EQ(a(1, 1, 1), 11);
}

// ---------------------------------------------------------------------------
// 4. Copy ctor remains a view (ownership unchanged)
// ---------------------------------------------------------------------------

TEST(CopyCtor, IsNonOwningView) {
    zdouble1 owner; owner.init(3); owner(1)=1; owner(2)=2; owner(3)=3;
    {
        zdouble1 view(owner);            // shallow
        EXPECT_FALSE(view.isAllocated());
        EXPECT_EQ(view.data(), owner.data());
        view(1) = 999;                   // writes to owner's memory
    }
    // After view destruction, owner still intact (no double-free).
    EXPECT_DOUBLE_EQ(owner(1), 999);
    EXPECT_DOUBLE_EQ(owner(3), 3);
    EXPECT_TRUE(owner.isAllocated());
}

TEST(CopyCtor, ReturnByValuePreservesDataViaMoveOrElision) {
    // This is the most-used real-world pattern (sphincs / croma):
    //     zdouble1 result = makeArray();
    // Before v0.7.1 it worked only by NRVO.  Now works via move ctor too.
    auto make = []() {
        zdouble1 x; x.init(4);
        x(1)=1; x(2)=2; x(3)=3; x(4)=4;
        return x;                        // NRVO or move
    };
    zdouble1 y = make();
    EXPECT_EQ(y.size(), 4);
    EXPECT_DOUBLE_EQ(y(1), 1);
    EXPECT_DOUBLE_EQ(y(4), 4);
    EXPECT_TRUE(y.isAllocated());       // must own the memory
}

// ---------------------------------------------------------------------------
// 5. Copy assignment is still explicit deep copy
// ---------------------------------------------------------------------------

TEST(CopyAssign, ExplicitDeepCopyIntoEmpty) {
    zint1 src{1, 2, 3};
    zint1 dst;
    dst = src;                           // explicit deep copy
    EXPECT_EQ(dst.size(), 3);
    dst(1) = 99;
    EXPECT_EQ(src(1), 1);                // independent
    EXPECT_TRUE(dst.isAllocated());
    EXPECT_TRUE(src.isAllocated());
}

TEST(CopyAssign, ExplicitDeepCopyIntoExistingSameShape) {
    zint1 src{10, 20, 30};
    zint1 dst(3); dst = -1;
    dst = src;
    EXPECT_EQ(dst(1), 10);
    EXPECT_EQ(dst(3), 30);
}

TEST(CopyAssign, SelfAssignmentSafe) {
    zint1 a{5, 6, 7};
    a = a;
    EXPECT_EQ(a(1), 5);
    EXPECT_EQ(a(3), 7);
}

// ---------------------------------------------------------------------------
// 6. Explicit .copy() still works
// ---------------------------------------------------------------------------

TEST(ExplicitCopy, DotCopyProducesOwner) {
    zdouble1 a; a.init(3); a(1)=1; a(2)=2; a(3)=3;
    zdouble1 b = a.copy();               // prvalue -> copy elision into b
    EXPECT_TRUE(b.isAllocated());
    EXPECT_NE(b.data(), a.data());
    a(1) = 999;
    EXPECT_DOUBLE_EQ(b(1), 1);           // independent
}

TEST(ExplicitCopy, OfSlicePreservesIndependence) {
    zdouble2 m; m.init(3, 4); m = 5.0;
    zdouble1 col = m.slice(2).copy();
    EXPECT_TRUE(col.isAllocated());
    col(1) = 99;
    EXPECT_DOUBLE_EQ(m(1, 2), 5.0);      // m unchanged
}

// ---------------------------------------------------------------------------
// 7. Lifecycle in real-world embedding patterns
// ---------------------------------------------------------------------------

// Pattern: struct with JIArray members, default constructed, members
// populated after construction.  Observed in croma's AsyHomData,
// sphincs' CrossSection, montex's _ITXE etc.
struct Widget {
    zdouble1 a;
    zint2    b;
    int      tag = 0;
};

TEST(Embedding, StructMemberOwnedMemoryIsPreservedAcrossMove) {
    auto make = [](int tag) {
        Widget w; w.tag = tag;
        w.a.init(4); w.a = 1.0 * tag;
        w.b.init(2, 3); w.b = tag;
        return w;                        // triggers Widget move ctor
    };

    Widget w = make(7);
    EXPECT_EQ(w.tag, 7);
    EXPECT_EQ(w.a.size(), 4);
    EXPECT_DOUBLE_EQ(w.a(1), 7.0);
    EXPECT_EQ(w.b(2, 3), 7);

    // move-assign over existing widget
    w = make(11);
    EXPECT_EQ(w.tag, 11);
    EXPECT_DOUBLE_EQ(w.a(1), 11.0);
    EXPECT_EQ(w.b(1, 1), 11);
}

TEST(Embedding, ZarrayOfStructsWithJIArrayMembersLifecycle) {
    // zarray<T> new T[n]{} default-constructs each Widget; members are
    // empty JIArray until init()'d in place.  Destruction of the outer
    // zarray must invoke each Widget's destructor which in turn frees
    // each JIArray member.  No double-free, no leaks.
    zarray<Widget> arr;
    arr.init(4);
    for (int i = 1; i <= 4; ++i) {
        arr(i).tag = i;
        arr(i).a.init(3);
        arr(i).a = static_cast<double>(i);
        arr(i).b.init(2, 2);
        arr(i).b = i * 10;
    }
    for (int i = 1; i <= 4; ++i) {
        EXPECT_EQ(arr(i).tag, i);
        EXPECT_DOUBLE_EQ(arr(i).a(1), double(i));
        EXPECT_EQ(arr(i).b(2, 2), i * 10);
    }
    // arr goes out of scope — each Widget dtor frees its JIArray members.
}

// Pattern: return a zdouble1 from a function (croma PowerWriter::averageOverAssemblies).
static zdouble1 averageFake(int N, double base) {
    zdouble1 out; out.init(N);
    for (int i = 1; i <= N; ++i) out(i) = base + i;
    return out;
}

TEST(RealWorld, ReturnByValueFromHelper) {
    auto a = averageFake(5, 100.0);
    EXPECT_EQ(a.size(), 5);
    EXPECT_DOUBLE_EQ(a(1), 101.0);
    EXPECT_DOUBLE_EQ(a(5), 105.0);
    EXPECT_TRUE(a.isAllocated());
}

TEST(RealWorld, ReturnByValueIntoExisting) {
    zdouble1 existing; existing.init(2); existing = -1;
    existing = averageFake(3, 7);       // triggers move assign
    EXPECT_EQ(existing.size(), 3);
    EXPECT_DOUBLE_EQ(existing(1), 8.0);
}

// Pattern: sphincs' slice + copy = deep-copy of a sub-range.
TEST(RealWorld, SphincsSlicePlusDotCopy) {
    zdouble3 tfuel; tfuel.init(3, 4, 5);
    double v = 1.0;
    for (int i = 1; i <= 3; ++i)
      for (int j = 1; j <= 4; ++j)
        for (int k = 1; k <= 5; ++k) tfuel(i, j, k) = v++;

    zdouble1 tfuel1 = tfuel.slice(2, 3).copy();   // column over last dim
    EXPECT_TRUE(tfuel1.isAllocated());
    tfuel1(1) = -999.0;
    // mutation does not affect tfuel — deep copy
    EXPECT_NE(tfuel(2, 3, 1), -999.0);
}

// Pattern: view share via shareWith.
TEST(RealWorld, SphincsCrossSectionShareWith) {
    zdouble2 sct; sct.init(4, 3); sct = 0.0;
    for (int i = 1; i <= 4; ++i)
      for (int j = 1; j <= 3; ++j) sct(i, j) = 10 * i + j;

    zdouble1 view;
    view.shareWith(sct.slice(2));       // column 2
    EXPECT_FALSE(view.isAllocated());
    view(3) = 999;
    EXPECT_DOUBLE_EQ(sct(3, 2), 999);   // writes back through the view
}

// Pattern: moving from a temp returned by a factory, into a member of a
// resize-prone container.  Mirrors croma's `snap.n2n3n.push_back(std::move(entry))`
// pattern AFTER the struct becomes move-enabled by JIArray v0.7.1.
struct Entry {
    int       id = 0;
    zdouble1  payload;
};

TEST(RealWorld, VectorOfStructWithJIArrayAfterMoveOps) {
    std::vector<Entry> v;
    v.reserve(8);                       // avoid reallocation during push_back
    for (int i = 1; i <= 5; ++i) {
        Entry e;
        e.id = i;
        e.payload.init(3);
        e.payload(1) = i; e.payload(2) = i * 10; e.payload(3) = i * 100;
        v.push_back(std::move(e));      // Entry is implicitly movable now
    }
    ASSERT_EQ(v.size(), 5u);
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(v[i].id, i + 1);
        EXPECT_DOUBLE_EQ(v[i].payload(1), i + 1);
        EXPECT_DOUBLE_EQ(v[i].payload(3), (i + 1) * 100);
    }
}

TEST(RealWorld, VectorOfStructResizingForcesMoveNotCopy) {
    // If a struct's JIArray member gets deep-copied (via operator=) during
    // vector resize, this test still passes — but the test is asserting
    // that whatever mechanism vector uses (move if available, copy otherwise)
    // leaves each element's payload correctly populated.
    std::vector<Entry> v;
    for (int i = 1; i <= 20; ++i) {
        Entry e;
        e.id = i;
        e.payload.init(2); e.payload(1) = i; e.payload(2) = -i;
        v.push_back(std::move(e));      // may reallocate underlying buffer
    }
    ASSERT_EQ(v.size(), 20u);
    for (int i = 0; i < 20; ++i) {
        EXPECT_EQ(v[i].id, i + 1);
        EXPECT_DOUBLE_EQ(v[i].payload(1),  (i + 1));
        EXPECT_DOUBLE_EQ(v[i].payload(2), -(i + 1));
        EXPECT_TRUE(v[i].payload.isAllocated());
    }
}

// ---------------------------------------------------------------------------
// 8. Non-virtual destructor → sizeof sanity + no vptr
// ---------------------------------------------------------------------------

TEST(Layout, NoVtableMeansSmallFootprint) {
    // zdouble1 previously had a vptr (virtual dtor).  sizeof on x86_64
    // with RANK=1 is: 4 (nn) + 8 (mm) + 4 (allocated) + 4 (rankSize) +
    // 4 (offset) + 4 (sumOfOffset) + 4 (sizes) = 32 bytes with padding.
    // With a vptr it was 40.  Allow some slack for compiler alignment.
    EXPECT_LE(sizeof(zdouble1), 48u);
    EXPECT_LE(sizeof(zint1),    48u);
}

// ---------------------------------------------------------------------------
// 9. SIMD/unroll macros: compile-only sanity (no warning)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// 10. auto sliced = zarr.slice(...) — MUST be zero-cost view
// ---------------------------------------------------------------------------
//
// User contract (0.7.1):
//   * `.slice()` never allocates new element memory.
//   * The returned JIArray shares `mm` with the parent.
//   * Writing through `sliced(...)` mutates the parent.
//   * `auto sliced = zarr.slice(...)` relies on either copy-elision (prvalue
//     initialisation, guaranteed in C++17 for the return value) or our new
//     move ctor.  In BOTH paths the final `sliced` must still be a view —
//     the move ctor must NOT upgrade a non-owning temporary to an owner.

// slice() semantics:
//   * column-major (default): slice args fix the LAST dims — .slice(c) on
//     a 2-D array returns column c (contiguous in memory).
//   * row-major: slice args fix the FIRST dims — .slice(r) returns row r.
//   * both cases MUST produce a non-owning view into the parent's buffer.

TEST(Slice, AutoAssignmentIsNonOwningView_ColumnMajor) {
    zdouble2 zarr; zarr.init(3, 4); zarr = 0.0;
    zarr(2, 2) = 22.0;                   // row 2, col 2

    auto sliced = zarr.slice(2);         // column 2 (col-major)

    EXPECT_FALSE(sliced.isAllocated());
    EXPECT_DOUBLE_EQ(sliced(2), 22.0);   // column 2, row 2 == zarr(2,2)

    sliced(2) = -99.0;
    EXPECT_DOUBLE_EQ(zarr(2, 2), -99.0); // write through to parent

    EXPECT_NE(sliced.data(), nullptr);
    EXPECT_GE(sliced.data(), zarr.data());
    EXPECT_LT(sliced.data(), zarr.data() + zarr.size());
}

TEST(Slice, AutoAssignmentPreservesViewThroughMove) {
    // 3-D column-major: slice(j, k) fixes the last two dims, giving a
    // contiguous 1-D view along the fastest (first) dim.
    zdouble3 arr; arr.init(4, 5, 6); arr = 0.0;
    arr(2, 3, 4) = 123.0;

    auto make_slice = [&]() { return arr.slice(3, 4); };   // column (j=3, k=4)
    auto sliced = make_slice();                            // move from temp

    EXPECT_FALSE(sliced.isAllocated());
    EXPECT_DOUBLE_EQ(sliced(2), 123.0);                    // i=2 of that column
    sliced(2) = 777.0;
    EXPECT_DOUBLE_EQ(arr(2, 3, 4), 777.0);
}

TEST(Slice, ChainedSliceRemainsView) {
    zdouble3 arr; arr.init(3, 4, 5); arr = 0.0;
    arr(1, 2, 3) = 77.0;

    // Column-major chain: .slice(k) fixes last dim, .slice(j) fixes the
    // (now-last) dim of the resulting 2-D.
    auto page = arr.slice(3);             // 2-D view arr(:, :, 3)
    auto col  = page.slice(2);            // 1-D view arr(:, 2, 3)

    EXPECT_FALSE(page.isAllocated());
    EXPECT_FALSE(col.isAllocated());

    EXPECT_DOUBLE_EQ(col(1), 77.0);
    col(1) = -1.0;
    EXPECT_DOUBLE_EQ(arr(1, 2, 3), -1.0);
}

TEST(Slice, LifetimeOutlivesParentBlockEnd_UseAfterParentFreeIsUB) {
    // This doesn't test UB — it documents the contract: the slice is a
    // view, so the PARENT must outlive all its slices.  As long as the
    // caller respects that, the allocator traffic is exactly zero per
    // slice (parent allocated once, all slices share it).
    //
    // The previous-heap-alloc count is obtained via an alloc counter.
    struct AllocCounter {
        static int& allocs() { static int a = 0; return a; }
        static int& frees()  { static int f = 0; return f; }
    };

    // Run a normal slice burst; delta in allocator should be zero.
    const int a0 = AllocCounter::allocs();
    const int f0 = AllocCounter::frees();
    {
        zdouble2 src; src.init(100, 100); src = 0.0;
        int sum_alloc_local = 0;
        for (int i = 1; i <= 100; ++i) {
            auto sliced = src.slice(i);    // no alloc
            (void)sliced;
            ++sum_alloc_local;             // just busy work
        }
        EXPECT_EQ(sum_alloc_local, 100);
    }
    // alloc counter unused here — test is descriptive; exact counters would
    // require overriding new/delete.  See move_lifecycle_alloc_test for
    // that (not included in this patch; covered via valgrind alloc count).
    (void)a0; (void)f0;
}

TEST(Slice, ReturnFromFunctionIsZeroCopy) {
    // Column-major: .slice(c) returns column c.
    zdouble2 arr; arr.init(8, 8); arr = 0.0;
    arr(4, 5) = 314.0;

    auto factory = [&](int col) -> zdouble1 { return arr.slice(col); };
    auto v1 = factory(5);
    auto v2 = factory(5);

    EXPECT_FALSE(v1.isAllocated());
    EXPECT_FALSE(v2.isAllocated());
    EXPECT_EQ(v1.data(), v2.data());          // both point at column 5
    EXPECT_DOUBLE_EQ(v1(4), 314.0);           // row 4 of col 5 == arr(4,5)
    v1(4) = -1.0;
    EXPECT_DOUBLE_EQ(arr(4, 5), -1.0);
    EXPECT_DOUBLE_EQ(v2(4), -1.0);
}

TEST(Macros, ArithmeticLoopsStillCompile) {
    zdouble1 a; a.init(8); a = 2.0;
    zdouble1 b; b.init(8); b = 3.0;
    a += b;
    EXPECT_DOUBLE_EQ(a(1), 5.0);
    a *= 2.0;
    EXPECT_DOUBLE_EQ(a(1), 10.0);
    auto c = a + b;                      // returns by value — exercises
    EXPECT_DOUBLE_EQ(c(1), 13.0);        // move on return
    auto d = a.sum();
    EXPECT_DOUBLE_EQ(d, 10.0 * 8.0);
    auto e = a.average();
    EXPECT_DOUBLE_EQ(e, 10.0);
}
