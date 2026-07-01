# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

JIArray is a **header-only C++17 N-dimensional array library** for scientific computing, with Fortran-style column-major layout and 1-based indexing by default. Everything lives under `include/jiarray/` in namespace `dnegri::jiarray`. There is no compiled library — the CMake target `jiarray` (alias of `dnegri.jiarray`) is `INTERFACE`-only.

The full public API is documented in `docs/manual.md`. Read it before adding or changing API surface; keep it in sync with code changes.

## Build & test

```bash
mkdir -p build && cd build
cmake ..              # Debug build defines JIARRAY_DEBUG (bounds checks on)
cmake --build .
ctest                 # or: ctest --output-on-failure
```

Tests require GoogleTest (`find_package(GTest)`); they are silently skipped if it is absent. Tests only build when JIArray is the top-level project (`PROJECT_IS_TOP_LEVEL`); override with `-DJIARRAY_BUILD_TESTS=ON`.

Run a single test executable or filter:

```bash
./build/dnegri.jiarray.column_major_test                       # one suite
./build/dnegri.jiarray.column_major_test --gtest_filter='*Slice*'
ctest -R column_major                                          # by ctest name
```

The Blitz++ benchmark suites (`column_major_iter_test`, `row_major_iter_test`) are opt-in via `-DJIARRAY_WITH_BLITZ=ON` (Mac Homebrew path) and are off by default — don't expect them in CI.

## Compile-time configuration (the central design idea)

Behavior is selected by macros defined **before including any header** (or via `-D`). Defaults live in `pch.h`:

- `JIARRAY_COLUMN_MAJOR` — `1` (default, Fortran/first-index-fastest) or `0` (row-major/C). **Changes slicing direction**: column-major `slice()` peels the last dimension; row-major peels the first.
- `JIARRAY_OFFSET` — `1` (default, 1-based) or `0` (0-based). Affects indexing *and* the `ffor`/`zfor` loop macros' bounds.
- `JIARRAY_DEBUG` — enables bounds/rank/size checks (throw `std::out_of_range` / `std::invalid_argument`); no-ops otherwise. Set automatically in Debug builds.
- `JIARRAY_USE_SIMD` — `1` (default, emits `#pragma omp simd`) or `0`. **Set to `0` when calling JIArray arithmetic inside `#pragma omp parallel` regions** to avoid nested parallelism.

Because these are compile-time, the test suite compiles the *same* source twice — e.g. `column_major_test` vs `row_major_test`, `fast_array_test` vs `fast_array_rowmajor_test` — with the row-major variant adding `-DJIARRAY_COLUMN_MAJOR=0`. **Any change to layout-dependent logic must pass both variants.**

## Headers and their roles

- `pch.h` — shared macros: config defaults, `JIARRAY_HD` (CUDA `__host__ __device__`), `JIARRAY_SIMD_LOOP`/`JIARRAY_SIMD_REDUCTION`/`JIARRAY_UNROLL` portable pragmas, `JIARRAY_CHECK_*` debug guards, and the `JIARRAY_ALLOCATED_*` bitmask constants.
- `JIArray.h` — the main `JIArray<T, RANK>` template plus all `z*` type aliases (`zint1`, `zdouble2`, …) and the `ffor`/`ffor_back`/`zfor` loop macros.
- `JIArrayExpr.h` — lazy expression templates backing element-wise `+ - * /` and unary `-` (see "Expression templates" below). Included by `JIArray.h`.
- `FastArray.h` — `FastArray<T, Dims...>`, a variadic, stack-allocated, compile-time-sized array (1D/2D/ND, also backs string arrays). CUDA-safe (`JIARRAY_HD`).
- `JIVector.h` — `JIVector<T>`, a thin `std::vector` subclass giving 1-based `operator()`/`at()` and `find`/`contains`.
- `JICudaArray.h` — `JICudaArray<T, RANK>`, a CUDA device-memory analogue of JIArray (only relevant under `nvcc`/`__CUDACC__`).
- `HighFiveExtension.hpp` — `HighFive::inspector` specializations enabling direct HDF5 read/write of JIArray/JIVector (needs HighFive available).

## Memory ownership model — the key gotcha

`JIArray` mixes owning and non-owning instances, tracked by the `allocated` bitmask (`JIARRAY_ALLOCATED_NONE/MEMORY/RANKSIZE/OFFSET/ALL`). Get this wrong and you get use-after-free or double-free:

- **Copy constructor (`JIArray(const JIArray&)`) is a SHALLOW, non-owning view** — it shares `mm` and sets `allocated = NONE`, so its destructor is a no-op. This exists so slices and return-by-value views are cheap.
- **Copy assignment (`operator=(const JIArray&)`) is a DEEP copy** — allocates if the target is empty, else asserts matching shape and copies elements. A source that is empty resets the target to empty (guards against UAF).
- **Move ctor/assignment transfer ownership** and null out the source.
- `slice()` and `reshape()` return views (share memory). `.copy()` is the explicit deep copy; `.shareWith()` makes an explicit view.

When writing tests or examples, prefer `.copy()` for an independent array. `auto x = arr;` invokes the copy constructor, so `x` is a non-owning view aliasing `arr`'s storage — the original stays the owner and must outlive `x`.

## Expression templates — the second `auto` gotcha

Element-wise `+ - * /` and unary `-` are **lazy expression templates** (in `JIArrayExpr.h`): each operator returns an expression node, and a chain like `d = a + b + c` evaluates in one fused, allocation-free pass on assignment/construction (no per-operator temporaries — this is the F4 fix). Consequences:

- **Never capture an operator result with `auto`** — `auto r = a + b;` binds a lazy node holding references to `a`/`b` (it can't be indexed/reassigned and dangles if the operands die). Write the concrete type: `zdouble1 r = a + b;` materialises it. Same rule as Eigen/xtensor. Tests must follow this.
- Nodes evaluate **linearly over flat storage** (`eval(i)`), so they are layout-agnostic — no column/row-major branching in the node code. Keep it that way for any new op.
- Nodes/operators are `JIARRAY_HD`, but the materializing `JIArray(expr)` ctor and `operator=(expr)` allocate, so they are **host-only** (not `JIARRAY_HD`). Don't add `JIARRAY_HD` to anything that allocates, or NVCC warns (`20011-D`).
- In-place `+= -= *= /=` stay eager single-pass (the other fast path); don't route them through expressions.

## Other implementation notes

- `ffor(i, b, e, step)` uses C++20 `__VA_OPT__`. On MSVC this needs the conforming preprocessor — the CMake `INTERFACE` target propagates `/Zc:preprocessor` to consumers (no-op on GCC/Clang). Bounds are inclusive when `JIARRAY_OFFSET=1`.
- "Not found" sentinels return `JIARRAY_OFFSET - 1` (i.e. `0` in 1-based mode) from `findFirst`/`find`.
- Arithmetic reductions/loops route through the `JIARRAY_SIMD_*` macros; keep new element-wise ops consistent with that pattern.

## Conventions

- Match the existing column-major/1-based assumptions and the heavy Doxygen comment style on public members.
- This is a header-only template library: no `.cpp` implementation files for the library itself — all code goes in the headers.
