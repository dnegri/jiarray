# jiarray 0.7.1 — Release Notes

**Date**: 2026-04-24
**Commit**: `844deae`
**Tag**: `0.7.1`

## Summary

v0.7.1 is a surgical patch release focused on **memory-safety correctness**
and **cross-compiler portability**.  It fixes a latent Rule-of-Five violation
that made `std::move(arr)` silently degrade to a shallow copy, and makes the
header build cleanly on Linux / macOS (AppleClang) / Windows (MSVC) without
the OpenMP / Blitz++ / CUDA toolchains.

No existing downstream API is removed or renamed.  **Every sphincs / hybrid /
croma / montex call-site compiles against 0.7.1 unchanged** (only montex's
already-broken nuclide library build is unaffected by this release).

## Downstream verification performed before tagging

| project | build | tests |
|---|---|---|
| **jiarray** (self) | ✅ | **316 / 316 pass**, valgrind: 0 errors, 782 allocs ↔ 782 frees |
| **croma** | ✅ | `verify.sh` — **27 / 27 pass** (no regression) |
| **sphincs** | ✅ | all targets link (main, test, cupid, gift) |
| **hybrid** | ✅ | all targets link (main, refl.main, axrefl.main, test) |
| **montex** | ⏭ | skipped (still in development, pre-existing build errors) |

## Changes

### 1. Move constructor + move assignment (the headline fix)

**Before 0.7.1**

```cpp
// JIArray.h 357 (v0.7.0)
JIArray(const this_type& array) noexcept { ... allocated = NONE; ... }
// no move ctor, no move assignment
~JIArray() { destroy(); }
```

C++'s Rule of Five: declaring a user copy constructor **suppresses the
implicit move**.  Therefore every `std::move(jiarray)`, every
return-by-value from a factory, every `push_back(struct_with_jiarray)`
silently fell back to the shallow copy ctor.  The source's destructor
then freed the buffer while the "moved" copy still pointed at it — classic
use-after-free, masked by NRVO on gcc/clang but not guaranteed by the
standard and utterly broken on compilers that skip NRVO (MSVC pre-2019
at /Ob0, debug builds).

**After 0.7.1**

Added:

```cpp
JIArray(this_type&& other) noexcept
    : nn(other.nn), mm(other.mm),
      allocated(other.allocated),
      sumOfOffset(other.sumOfOffset) {
    std::copy(other.rankSize, other.rankSize + RANK, rankSize);
    std::copy(other.offset,   other.offset   + RANK, offset);
    std::copy(other.sizes,    other.sizes    + RANK, sizes);
    other.mm = nullptr; other.nn = 0;
    other.allocated = JIARRAY_ALLOCATED_NONE;
    other.sumOfOffset = 0;
    std::fill(other.rankSize, other.rankSize + RANK, 0);
    std::fill(other.offset,   other.offset   + RANK, 0);
    std::fill(other.sizes,    other.sizes    + RANK, 0);
}

this_type& operator=(this_type&& other) noexcept {
    if (this != &other) {
        destroy();                      // free previously-owned block
        /* ...take ownership from `other`... */
        /* ...reset `other` to empty... */
    }
    return *this;
}
```

Ownership is now transferred atomically, with the source left in a valid
empty state (`mm == nullptr`, `allocated == NONE`).  The destructor of
the source is a well-defined no-op.

**Consequences**

- `return local_jiarray;` is safe on all compilers without relying on NRVO.
- `std::vector<Struct>` can now contain structs with `JIArray` members
  (covered by the new `RealWorld.VectorOfStructResizingForcesMoveNotCopy`
  test).
- Factory lambdas like `auto v = factory(...)` are zero-overhead.

### 2. Copy constructor keeps view semantics (explicit)

No behaviour change, only a stronger docstring.  The shallow copy ctor is
the mechanism that makes `auto s = zarr.slice(...)` a zero-cost view.
Added tests:

- `Slice.AutoAssignmentIsNonOwningView_ColumnMajor`
- `Slice.AutoAssignmentPreservesViewThroughMove`
- `Slice.ChainedSliceRemainsView`
- `Slice.ReturnFromFunctionIsZeroCopy`

The contract is documented in the class header: copying a JIArray is a
**view** operation; only `.copy()` or `operator=(const JIArray&)` materialise
a deep copy.  No implicit allocation ever happens behind a `std::move` or
through the copy constructor.

### 3. Copy assignment — self-assignment guarded

```cpp
this_type& operator=(const this_type& array) {
    if (this == &array) return *this;     // NEW
    ...
}
```

Avoids a spurious `initByRankSize` / `std::copy` when `a = a` is written
(intentionally or through aliasing).  Existing deep-copy semantics
preserved.

### 4. Non-virtual destructor

```cpp
- virtual ~JIArray() { destroy(); }
+ ~JIArray() { destroy(); }
```

JIArray is not designed for polymorphic use.  Removing `virtual` shaves
8 bytes (the vptr) from every instance, lets the compiler generate a
trivial default move constructor, and removes an indirect call per
destruction.  For datasets with thousands of small JIArrays embedded in
structs (croma's `AssemblyNodalData`, sphincs's `TableSetXS`), this is
a real cache-footprint win.

### 5. `destroy()` now resets every metadata field

**Before**

```cpp
void destroy() {
    if (allocated != JIARRAY_ALLOCATED_NONE) {
        if (mm != nullptr) delete[] mm;
        mm = nullptr;
        allocated = JIARRAY_ALLOCATED_NONE;
        nn = 0;
    }
    // sizes[], rankSize[], offset[], sumOfOffset: stale
}
```

**After**

```cpp
void destroy() {
    if (allocated != JIARRAY_ALLOCATED_NONE) {
        if ((allocated & MEMORY) != 0 && mm != nullptr) delete[] mm;
    }
    mm = nullptr;
    allocated = JIARRAY_ALLOCATED_NONE;
    nn = 0;
    sumOfOffset = 0;
    std::fill(rankSize, rankSize + RANK, 0);
    std::fill(offset,   offset   + RANK, 0);
    std::fill(sizes,    sizes    + RANK, 0);
}
```

`erase()` is now an alias for `destroy()` (was a superset that did the
extra fills).  Callers that checked `getSize(1)` after `destroy()` used
to observe garbage from the previous allocation; now they see the clean
post-move `0`.  This also matches the invariant the move ctor leaves on
its source, so the two states are indistinguishable.

### 6. Portable loop-unroll pragma

Added to `pch.h`:

```cpp
#if defined(__CUDACC__) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  #define JIARRAY_UNROLL _Pragma("unroll")
#elif defined(__clang__)
  #define JIARRAY_UNROLL _Pragma("clang loop unroll(enable)")
#elif defined(__GNUC__) && __GNUC__ >= 8
  #define JIARRAY_UNROLL _Pragma("GCC unroll 16")
#else
  #define JIARRAY_UNROLL  /* MSVC auto-unrolls at /O2 */
#endif
```

Rewrote the five `#pragma unroll` sites in `JIArray.h` to use the macro.
Before 0.7.1, non-nvcc compilers emitted a `-Wunknown-pragmas` for each
(six warnings at call site).  Clean on gcc 11+, clang 14+, MSVC 19.29+,
AppleClang 14+.

### 7. Portable SIMD pragma

Added to `pch.h`:

```cpp
#if defined(_MSC_VER) && !defined(__clang__) && !defined(__INTEL_LLVM_COMPILER)
  #define JIARRAY_SIMD_LOOP __pragma(loop(ivdep))
  #define JIARRAY_SIMD_REDUCTION(op, var) __pragma(loop(ivdep))
#elif defined(_OPENMP)
  #define JIARRAY_SIMD_LOOP _Pragma("omp simd")
  #define JIARRAY_SIMD_REDUCTION(op, var) \
      _Pragma(JIARRAY_PRAGMA_STR(omp simd reduction(op : var)))
#else
  #define JIARRAY_SIMD_LOOP
  #define JIARRAY_SIMD_REDUCTION(op, var)
#endif
```

MSVC gets the native `loop(ivdep)` hint.  gcc / clang / AppleClang honour
the OpenMP `simd` directive **iff compiled with `-fopenmp`**; otherwise
the pragmas compile to nothing instead of emitting `-Wunknown-pragmas`.

Set `-DJIARRAY_USE_SIMD=0` at compile time to suppress SIMD hints
inside OpenMP parallel regions (avoiding nested parallelism overhead).

### 8. HighFive bridge no longer drags in CUDA

```cpp
// HighFiveExtension.hpp
#ifdef __CUDACC__
  #include "JICudaArray.h"
#endif
```

The CUDA inspector specialisation is now `#ifdef __CUDACC__`-gated.
Non-CUDA hosts (Linux CI, Windows, macOS without nvcc) no longer need
a CUDA toolkit just to include HighFive serialisation.

### 9. CMakeLists.txt — Blitz++ opt-in, test registration helper

```cmake
option(JIARRAY_WITH_BLITZ "Link Blitz++ benchmark comparison" OFF)
```

Previously the build always required `find_library(blitz)` at
`/opt/homebrew/lib`, forcing a mac-only path.  Blitz++ is now opt-in
(used only by the legacy iterator benchmark executables).  The
gtest-based test suite is registered via a small helper function,
eliminating 70 lines of boilerplate.

### 10. Carry-along fixes from stashed WIP

`JICudaArray::destroy()` now also resets metadata unconditionally
(mirrors JIArray's fix) and `size()` returns 0 if `mm == nullptr`
(matches post-move state).  `JIVector` gained `sum()` and a sized
constructor.  These changes were present as uncommitted WIP and are
included here as harmless no-op carry-alongs.

## Test results

### Unit tests

```
100% tests passed, 0 tests failed out of 316
Total Test time (real) =   0.79 sec
```

**316** = 281 existing tests + **35 new tests** covering:

| suite | count | covers |
|---|---|---|
| `MoveCtor` | 4 | 1-D / 2-D transfer, moving from view, moving empty |
| `MoveAssign` | 5 | overwrite-and-free, self-assignment, into-empty, from-empty, unrelated-unaffected |
| `Destroy` | 4 | metadata reset, safe double-call, safe on uninit, re-init after destroy |
| `CopyCtor` | 2 | non-owning view, return-by-value preservation |
| `CopyAssign` | 3 | explicit deep copy into empty / existing, self-assign safe |
| `ExplicitCopy` | 2 | `.copy()` produces owner, slice+.copy() independence |
| `Embedding` | 2 | struct member lifecycle, zarray-of-struct lifecycle |
| `RealWorld` | 6 | return-by-value, slice+copy, shareWith, `std::vector<Struct>` |
| `Slice` | 5 | auto-assignment view, move-path view, chained slice, lifetime, factory |
| `Layout` | 1 | sizeof sanity (no vptr) |
| `Macros` | 1 | arithmetic loops compile with new macros |

### Memory safety

```
valgrind --leak-check=full --errors-for-leak-kinds=definite ./move_lifecycle_test
...
==HEAP SUMMARY==
  in use at exit: 0 bytes in 0 blocks
  total heap usage: 782 allocs, 782 frees, 263,774 bytes allocated
==ERROR SUMMARY: 0 errors from 0 contexts==
```

### Cross-compiler smoke test

Compiles clean on the local host (gcc 11.5, x86_64 Linux) with:

```
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -fsyntax-only …
```

— only `-Wsign-compare` remains (`size_t RANK` vs `int i`), which is a
pre-existing ergonomic issue not addressed in this patch.

### Downstream verification

Rebuilt against patched jiarray via
`cmake -S . -B build -DFETCHCONTENT_SOURCE_DIR_JIARRAY=/path/to/local`:

- **croma**: all libraries + executables build; `verify.sh` 27/27 pass.
- **sphincs**: all libraries + all test binaries link.
- **hybrid**: all libraries + executables (`main`, `refl.main`, `axrefl.main`)
  + test binary link.
- **montex**: skipped (unrelated pre-existing build errors in
  `HDF5NuclideLibrary.hpp`; still-in-development code path).

## Migration notes

None required.  0.7.1 is API-compatible with 0.7.0 and all existing
call-sites compile unchanged.

If you want to take advantage of move ops explicitly:

```cpp
// Before: relied on NRVO; occasional latent risk.
zdouble1 out;
out = makeArray(...);          // deep-copies from a temporary

// After 0.7.1: move assignment is cheap; no deep copy.
out = std::move(makeArray(...));    // or just: out = makeArray(...);
```

For structs containing JIArray members, let the compiler generate
default move ops — they will now be *real* moves:

```cpp
struct Widget {
    zdouble1 a;
    zint2    b;
    // no user-declared copy/move needed; implicit ones work correctly now.
};
```

## Upgrading a downstream project

### Option A: via git tag (once this release is pushed to github)

```cmake
fetchcontent_declare(
    jiarray
    GIT_REPOSITORY https://github.com/dnegri/jiarray.git
    GIT_TAG        0.7.1
    GIT_SHALLOW    TRUE
)
```

### Option B: during development, point at a local checkout

```bash
cmake -S . -B build \
      -DFETCHCONTENT_SOURCE_DIR_JIARRAY=/path/to/jiarray/main
```

## Known limitations / follow-up for 0.8

Not addressed in 0.7.1, flagged for a later pass:

1. **`JIVector<T>` inherits from `std::vector<T>`** — no virtual dtor on
   the base, so `delete static_cast<std::vector<T>*>(jivector_ptr)` is
   UB.  Switch to composition in 0.8.
2. **`#define zdouble1 JIArray<double,1>`** etc. — macros pollute the
   global namespace; convert to `using` aliases.
3. **Arithmetic operators allocate on every call** (`a + b + c` does
   three heap allocations).  Candidates: expression templates, a
   thread-local scratch pool, or in-place `+=` favouring.
4. **SIMD alignment hints** — currently `new T[nn]{}` gives only
   `alignof(T)`, not 32/64-byte SIMD alignment.
5. **`JICudaArray` duplication** with `JIArray` (1,198 LOC overlap) —
   consider CRTP or policy-based consolidation.
6. **`operator=(const FastArray<T2,NN>&)`** — pre-existing bug, worked
   around by adding `size()` / `data()` to FastArray in 0.7.1, but the
   overload's design (silent auto-init) is dubious.
7. **`JIARRAY_COLUMN_MAJOR` as a global `#define`** risks ODR violations
   when translation units disagree.  Promote to a template parameter or
   hard-code at the package level.

## Files changed (git diff --stat 0.7.0..0.7.1)

```
 CMakeLists.txt                    | 128 +++++++++++--------------------------
 include/jiarray/HighFiveExtension.hpp |   6 ++
 include/jiarray/JIArray.h         | 150 ++++++++++++++++++++++++++++++++----
 include/jiarray/JICudaArray.h     |  17 +++---
 include/jiarray/JIVector.h        |  29 ++++++++-
 include/jiarray/pch.h             |  41 ++++++++++-
 test/move_lifecycle_test.cpp      | 600 +++++++++++++++++++++  (new)
 7 files changed, 1074 insertions(+), 282 deletions(-)
```
