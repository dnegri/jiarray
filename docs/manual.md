# JIArray Reference Manual

## Table of Contents

- [Configuration](#configuration)
  - [Storage Order](#storage-order)
  - [Index Offset](#index-offset)
  - [Debug Mode](#debug-mode)
  - [SIMD](#simd)
- [JIArray](#jiarray)
  - [Type Aliases](#type-aliases)
  - [Construction](#construction)
  - [Initialization](#initialization)
  - [Memory Management](#memory-management)
  - [Element Access](#element-access)
  - [Size and Offset Queries](#size-and-offset-queries)
  - [Slicing](#slicing)
  - [Reshape](#reshape)
  - [Assignment](#assignment)
  - [Arithmetic](#arithmetic)
  - [Statistics](#statistics)
  - [Search](#search)
  - [Utility](#utility)
  - [Iterators](#iterators)
  - [Loop Macros](#loop-macros)
- [Multi-Dimensional Arrays](#multi-dimensional-arrays)
- [FastArray](#fastarray)
  - [FastArray2D](#fastarray2d)
  - [StringFastArray](#stringfastarray)
  - [FastArray Type Aliases](#fastarray-type-aliases)
- [JIVector](#jivector)
  - [JIVector Type Aliases](#jivector-type-aliases)
- [JIDeviceArray](#jidevicearray)
  - [Device Type Aliases](#device-type-aliases)
  - [Allocation](#allocation)
  - [Host-Device Transfer](#host-device-transfer)
  - [Kernel Views and Indexing](#kernel-views-and-indexing)
  - [Indexing Cost](#indexing-cost)
  - [Debug Checks](#debug-checks)
  - [CUDA Tests](#cuda-tests)
- [Integration](#integration)
  - [CMake](#cmake)
  - [HighFive HDF5 Extension](#highfive-hdf5-extension)

---

## Configuration

JIArray behavior is controlled by compile-time macros. Define them **before** including any JIArray header, or pass via compiler flags (`-D`).

### Storage Order

```cpp
#define JIARRAY_COLUMN_MAJOR 1   // default: column-major (Fortran-style)
#define JIARRAY_COLUMN_MAJOR 0   // row-major (C-style)
```

**Column-major** (default): the first index varies fastest in memory.

```
2D array (3 rows x 4 cols), column-major memory layout:

  (1,1) (1,2) (1,3) (1,4)
  (2,1) (2,2) (2,3) (2,4)
  (3,1) (3,2) (3,3) (3,4)

Memory: [(1,1),(2,1),(3,1), (1,2),(2,2),(3,2), (1,3),(2,3),(3,3), (1,4),(2,4),(3,4)]
         ---- col 1 -----  ---- col 2 -----  ---- col 3 -----  ---- col 4 -----
```

**Row-major**: the last index varies fastest in memory.

```
Same 3x4 array, row-major memory layout:

Memory: [(1,1),(1,2),(1,3),(1,4), (2,1),(2,2),(2,3),(2,4), (3,1),(3,2),(3,3),(3,4)]
         ------ row 1 --------  ------ row 2 --------  ------ row 3 --------
```

**3D example** -- array(2, 3, 4):

| Order | Stride (rankSize) | Fastest-varying index |
|---|---|---|
| Column-major | `[1, 2, 6]` | 1st dimension |
| Row-major | `[12, 4, 1]` | 3rd (last) dimension |

**Performance tip**: iterate in memory order for cache efficiency. For column-major, the innermost loop should be the first index. For row-major, the innermost loop should be the last index.

```cpp
// Column-major (default) -- fast iteration
zdouble2 mat(1000, 1000);
for (int j = 1; j <= 1000; ++j)      // outer: 2nd index
    for (int i = 1; i <= 1000; ++i)   // inner: 1st index (contiguous)
        mat(i, j) = i + j;

// Row-major -- fast iteration
// compile with -DJIARRAY_COLUMN_MAJOR=0
for (int i = 1; i <= 1000; ++i)      // outer: 1st index
    for (int j = 1; j <= 1000; ++j)   // inner: 2nd index (contiguous)
        mat(i, j) = i + j;
```

### Index Offset

```cpp
#define JIARRAY_OFFSET 1   // default: 1-based indexing (Fortran-style)
#define JIARRAY_OFFSET 0   // 0-based indexing (C-style)
```

With `JIARRAY_OFFSET=1` (default):
```cpp
zint1 arr(5);          // indices: 1, 2, 3, 4, 5
arr(1) = 10;           // first element
arr(5) = 50;           // last element
```

With `JIARRAY_OFFSET=0`:
```cpp
zint1 arr(5);          // indices: 0, 1, 2, 3, 4
arr(0) = 10;           // first element
arr(4) = 50;           // last element
```

**Custom offsets per dimension** with `init0(first, end, ...)`. The range follows the
index base, as `ffor` does: `[first, end]` with `JIARRAY_OFFSET=1`, `[first, end)` with
`JIARRAY_OFFSET=0` (so `first == end` allocates nothing when 0-based):
```cpp
zdouble2 arr;
arr.init0(5, 10, 20, 25);   // 1-based: dim1 [5..10], dim2 [20..25]  (6 x 6)
                            // 0-based: dim1 [5..9],  dim2 [20..24]  (5 x 5)
arr(5, 20) = 1.0;           // valid in both
arr(9, 24) = 2.0;           // valid in both; arr(10, 25) only when 1-based
```

**Custom offsets after init** with `setOffsets()`:
```cpp
zint2 arr(3, 4);
arr.setOffsets(0, 0);        // switch to 0-based
arr(0, 0) = 42;
```

### Debug Mode

```cpp
#define JIARRAY_DEBUG   // enable bounds checking (define before including headers)
```

When active, out-of-bounds access throws `std::out_of_range`, size mismatches throw `std::invalid_argument`, and rank mismatches throw `std::invalid_argument`.

When not defined, all checks are no-ops for maximum performance.

### SIMD

```cpp
#define JIARRAY_USE_SIMD 1   // default: enable OpenMP SIMD directives
#define JIARRAY_USE_SIMD 0   // disable (use inside OpenMP parallel regions)
```

When enabled, arithmetic and reduction operations use `#pragma omp simd` for auto-vectorization. Disable when calling JIArray methods from within `#pragma omp parallel` regions to avoid nested parallelism overhead.

```cpp
// Safe: SIMD enabled (default)
zdouble1 arr(10000);
arr += 1.0;   // auto-vectorized

// Inside OpenMP parallel region: disable SIMD
#define JIARRAY_USE_SIMD 0
#include <jiarray/JIArray.h>

#pragma omp parallel for
for (int i = 0; i < n; ++i) {
    arrays[i] += 1.0;   // no nested SIMD
}
```

---

## JIArray

The main class template:

```cpp
template <class T, size_t RANK = 1, class SEQ = std::make_index_sequence<RANK>>
class JIArray;
```

- `T` -- element type (must be default constructible)
- `RANK` -- number of dimensions (must be > 0)

### Type Aliases

| Alias | Type |
|---|---|
| `zint1` .. `zint5` | `JIArray<int, 1>` .. `JIArray<int, 5>` |
| `zdouble1` .. `zdouble6` | `JIArray<double, 1>` .. `JIArray<double, 6>` |
| `zfloat1` .. `zfloat5` | `JIArray<float, 1>` .. `JIArray<float, 5>` |
| `zbool1` .. `zbool5` | `JIArray<bool, 1>` .. `JIArray<bool, 5>` |
| `zstring1` .. `zstring5` | `JIArray<std::string, 1>` .. `JIArray<std::string, 5>` |
| `zarray<T, N>` | `JIArray<T, N>` |

These are namespaced `using` aliases in `dnegri::jiarray` — reach them via
`using namespace dnegri::jiarray;` or qualification. (They were global preprocessor
macros in prior releases; the change is source-compatible for code that already imports
the namespace, and no longer collides with unrelated identifiers named e.g. `zint1`.)

### Construction

```cpp
// Default (empty, no allocation)
zint2 arr;

// With dimensions (allocates zero-initialized memory)
zint2 mat(3, 4);           // 3x4 matrix
zint3 cube(2, 3, 4);       // 2x3x4 cube

// 1D from initializer list
zint1 v{10, 20, 30};

// Copy constructor (shallow -- shares memory)
zint1 view(v);             // view shares memory with v

// External memory (array does not own the buffer)
double buf[12];
zdouble2 wrapper(buf, 3, 4);   // wraps buf as 3x4 matrix
```

### Initialization

```cpp
// init() -- allocate with dimensions
zint2 arr;
arr.init(3, 4);

// init0() -- allocate with custom index ranges
zdouble2 arr;
arr.init0(5, 10, 20, 25);     // dim1=[5..10], dim2=[20..25] (1-based); [5..9], [20..24] (0-based)

// init() with external memory
double buf[12];
zdouble2 arr;
arr.init(3, 4, buf);          // wraps buf, does NOT own it

// initByRankSize() -- advanced: specify strides, offsets, memory
zint1 arr;
int rankSizes[] = {1};
int offsets[] = {1};
arr.initByRankSize(5, rankSizes, offsets);   // allocates 5 elements

// initByRankSize() with external memory
int buf[6];
zint2 arr;
int rankSizes[] = {1, 2};     // column-major 2x3
int offsets[] = {1, 1};
arr.initByRankSize(6, rankSizes, offsets, buf);
```

### Memory Management

```cpp
arr.destroy();    // deallocate memory (if owned), reset size to 0
arr.erase();      // destroy + reset all metadata (dimensions, offsets, strides)
arr.isAllocated();  // true if array owns or references memory
```

After `destroy()` or `erase()`, the array can be re-initialized with `init()`.

### Memory Ownership & View Lifetime

A `JIArray` is either an **owner** (it allocated its buffer and frees it on destruction)
or a **non-owning view** (it points into memory owned by someone else). Getting this wrong
is the main source of dangling reads, so the rules are explicit:

| Operation | Result | Owns memory? |
|---|---|---|
| `zint1 a(n);` / `a.init(n)` | fresh allocation | **owner** |
| `zint1 b(a);` (copy constructor) | **shallow view** aliasing `a` | no — view |
| `b = a;` (copy assignment) | **deep copy** into `b` (allocates if empty) | owner |
| `auto c = std::move(a);` | ownership transferred; `a` left empty | owner |
| `a.slice(i)` / `a.reshape(...)` | **view** sharing `a`'s buffer | no — view |
| `a.copy()` | explicit **deep copy** | owner |
| `b.shareWith(a)` | explicit **view** of `a` | no — view |
| `zdouble2 w(buf, r, c)` | wraps external `buf` | no — view |

**A view must not outlive the owner of its memory.** Using a view after the owner is
destroyed (or reallocated via `init`/`setSize`) is a use-after-free.

```cpp
zint1 make() { zint1 t{1,2,3}; return t.slice(0); }  // ✗ view into a buffer freed on return
zint1 owner{1,2,3};
zint1 v = owner.slice(1);   // ok while `owner` is alive
owner.init(10);             // ✗ reallocates — `v` now dangles
```

**Two `auto` traps to know:**

```cpp
auto x = arr;        // ✗ copy constructor → x is a VIEW aliasing arr (not a copy!)
zint1 x = arr.copy();// ✅ independent owner

auto r = a + b;      // ✗ r is a lazy expression node holding refs to a,b (see Arithmetic)
zint1 r = a + b;     // ✅ materialised owner
```

Rule of thumb: when you want an **independent array**, take a concrete type and use `.copy()`
(or assign an expression into a concrete type). Reserve `auto` for genuine views whose owner
you know outlives them (e.g. a local `auto col = mat.slice(2);` used within the same scope).

**Catching lifetime bugs.** A dangling-view read is a heap-use-after-free, so build your tests
with **AddressSanitizer** (`-fsanitize=address`) to detect it precisely — ASan reports the
read site, the free site (the owner's `destroy()`), and the original allocation, with no
special library flag. Destroying the *owner* itself and then indexing it is additionally
caught by `JIARRAY_DEBUG` bounds checks (its sizes reset to zero).

### Element Access

```cpp
zint2 mat(3, 4);

// operator() -- primary access method (1-based by default)
mat(2, 3) = 42;
int val = mat(2, 3);

// at() -- same as operator()
mat.at(2, 3) = 42;

// Const access
const zint2& cref = mat;
int val = cref(2, 3);

// Access with FastArray index
FastArray<int, 2> idx({2, 3});
int val = mat(idx);

// data() -- raw pointer to first element
int* p = mat.data();

// data(INDEX...) -- pointer to subarray start
// Column-major: data(j) gives pointer to column j
const int* col2 = mat.data(2);   // points to start of column 2
// col2[0] = mat(1,2), col2[1] = mat(2,2), col2[2] = mat(3,2)

// getMemory(offset) -- pointer offset from start
int* p = mat.getMemory(0);   // same as data()
int* p3 = mat.getMemory(3);  // 4th element in linear storage
```

### Size and Offset Queries

```cpp
zint3 arr(2, 3, 4);

arr.size();           // 24 (total elements)
arr.getSize();        // 24 (same)
arr.getSize(1);       // 2 (dimension 1 size, 1-based rank index)
arr.getSize(2);       // 3
arr.getSize(3);       // 4

arr.getSizeOfRank();  // int[3] = {2, 3, 4}
arr.getRankSize();    // int[3] = strides, e.g. {1, 2, 6} for column-major

arr.getOffset();      // int[3] = {1, 1, 1} for default JIARRAY_OFFSET=1
arr.getOffset(1);     // 1 (offset for dimension 1)
```

### Slicing

Slicing creates a **zero-copy view** into the original array, reducing the rank by the number of indices provided.

**Column-major** slicing removes dimensions from the **right** (last dimension first):

```cpp
// Column-major (default)
zint3 cube(2, 3, 4);    // 2x3x4

zint2 plane = cube.slice(2);      // cube(:, :, 2) -> 2x3 view
zint1 line  = cube.slice(2, 3);   // cube(:, 2, 3) -> size-2 view
```

In column-major, `slice(k)` selects the k-th "page" along the last dimension.

**Row-major** slicing removes dimensions from the **left** (first dimension first):

```cpp
// Row-major (compile with -DJIARRAY_COLUMN_MAJOR=0)
zint3 cube(2, 3, 4);    // 2x3x4

zint2 plane = cube.slice(1);      // cube(1, :, :) -> 3x4 view
zint1 line  = cube.slice(1, 2);   // cube(1, 2, :) -> size-4 view
```

**Slicing is a view** -- modifications through the slice affect the original:

```cpp
zdouble2 mat(3, 4);
mat = 0.0;
auto col = mat.slice(2);   // column 2
col(1) = 99.0;
// mat(1, 2) is now 99.0
```

**Const slicing** works on const arrays:

```cpp
const zdouble2& cmat = mat;
auto col = cmat.slice(2);   // const view, read-only
```

### Reshape

Creates a view with different dimensionality (same total size):

```cpp
zdouble1 flat(12);
for (int i = 1; i <= 12; ++i) flat(i) = i;

auto mat = flat.reshape(3, 4);   // 3x4 view
// Column-major: mat(1,1)=1, mat(2,1)=2, mat(3,1)=3, mat(1,2)=4, ...
```

### Assignment

```cpp
zint1 arr(5);

// Scalar -- fill all elements
arr = 42;

// Initializer list
arr = {1, 2, 3, 4, 5};

// std::vector
std::vector<int> v{10, 20, 30, 40, 50};
arr = v;

// C-style array pointer
int raw[] = {1, 2, 3, 4, 5};
arr = raw;

// Another JIArray (deep copy, allocates if empty)
zint1 src{1, 2, 3};
zint1 dst;
dst = src;         // dst is allocated and filled

// FastArray
FastArray<int, 5> fa({10, 20, 30, 40, 50});
zint1 arr2;
arr2 = fa;         // auto-allocates
```

### Arithmetic

Element-wise `+ - * /` (array/array and array/scalar) and unary `-` are **expression
templates**: an operator returns a lazy expression node, and a whole chain such as
`d = a + b + c` is evaluated in a **single, allocation-free pass** when it is assigned to (or
used to construct) a `JIArray`. This eliminates the per-operator temporaries that the earlier
eager operators created — on the streaming benchmark the operator chain went from ~7.9× slower
than a hand-fused loop to running **at the memory-bandwidth ceiling** (see `docs/EVALUATION.md`).

> ⚠️ **Do not capture an operator result with `auto`.** Because operators return an expression
> node (not a `JIArray`), write the concrete array type so the result materialises:
> ```cpp
> zdouble1 c = a + b;   // ✅ materialises into a JIArray
> auto     c = a + b;   // ✗ c is a lazy node holding references to a and b;
>                       //    it cannot be indexed/reassigned and dangles if a or b die
> ```
> This is the same rule Eigen/xtensor use. SIMD vectorization of the fused pass is still
> controlled by `JIARRAY_USE_SIMD`.

**Array-array operations** (element-wise, same size required):

```cpp
zdouble1 a{1.0, 2.0, 3.0};
zdouble1 b{4.0, 5.0, 6.0};

zdouble1 c = a + b;     // {5.0, 7.0, 9.0}
zdouble1 d = a - b;     // {-3.0, -3.0, -3.0}
zdouble1 e = a * b;     // {4.0, 10.0, 18.0}
zdouble1 f = a / b;     // {0.25, 0.4, 0.5}
zdouble1 g = a + b + c; // fused: one pass, no temporaries

a += b;             // a = {5.0, 7.0, 9.0}
a -= b;             // a = {1.0, 2.0, 3.0}
a *= b;             // a = {4.0, 10.0, 18.0}
a /= b;             // a = {1.0, 2.0, 3.0}
```

**Scalar operations** (both orderings):

```cpp
zdouble1 arr{1.0, 2.0, 3.0};

zdouble1 r1 = arr + 10.0;    // {11.0, 12.0, 13.0}
zdouble1 r2 = 10.0 + arr;    // {11.0, 12.0, 13.0}
zdouble1 r3 = arr * 3.0;     // {3.0, 6.0, 9.0}
zdouble1 r4 = 3.0 * arr;     // {3.0, 6.0, 9.0}
zdouble1 r5 = arr / 2.0;     // {0.5, 1.0, 1.5}
zdouble1 r6 = 6.0 / arr;     // {6.0, 3.0, 2.0}

arr += 10.0;              // {11.0, 12.0, 13.0}
arr *= 2.0;               // {22.0, 24.0, 26.0}
arr /= 2.0;               // {11.0, 12.0, 13.0}
```

**Unary negation**:

```cpp
zdouble1 arr{1.0, -2.0, 3.0};
zdouble1 neg = -arr;          // {-1.0, 2.0, -3.0}
```

**Notes**:
- In-place `+= -= *= /=` remain eager single-pass operations (no expression templates).
- For floating-point `/=` by scalar, the library uses reciprocal multiplication (`*= 1/val`)
  for performance. (The expression-template `array / scalar` divides per element.)

### Statistics

```cpp
zdouble1 arr{1.0, 2.0, 3.0, 4.0, 5.0};

arr.sum();        // 15.0
arr.average();    // 3.0
arr.min();        // 1.0
arr.max();        // 5.0
arr.sqsum();      // 55.0  (sum of squares, returns double)

zdouble1 a{1.0, 2.0, 3.0};
zdouble1 b{4.0, 5.0, 6.0};
dot(a, b);        // 32.0  (1*4 + 2*5 + 3*6)
```

### Search

```cpp
zint1 arr{10, 20, 30, 20, 50};

// contains
arr.contains(30);     // true
arr.contains(99);     // false

// findFirst (1D) -- returns 1-based index
arr.findFirst(20);    // 2 (first occurrence)
arr.findFirst(99);    // 0 (JIARRAY_OFFSET - 1, not found)

// findFirst (multi-dimensional) -- returns FastArray<int, RANK>
zint2 mat(3, 3);
mat = 0;
mat(2, 3) = 42;
auto loc = mat.findFirst(42);
// loc(1) = 2, loc(2) = 3

// maxloc -- location of maximum element
zint1 v{5, 30, 10, 25};
auto ml = v.maxloc();
// ml(1) = 2  (30 is at index 2)

// maxloc with range (1-based linear range)
auto ml2 = v.maxloc(2, 4);   // search indices 2..4 -> max is 30 at index 2

// Multi-dimensional maxloc
zint2 mat(3, 4);
mat = 0;
mat(2, 3) = 99;
auto loc = mat.maxloc();
// loc(1) = 2, loc(2) = 3
```

### Utility

```cpp
// Deep copy (independent of original)
zint1 original{1, 2, 3};
auto copied = original.copy();
copied(1) = 999;        // original(1) is still 1

// Share memory (creates a view)
zint1 view;
view.shareWith(original);
view(2) = 42;           // original(2) is now 42

// Convert to std::vector
auto vec = original.convertToVector();   // std::vector<int>{1, 42, 3}

// Comparison
zint1 a{1, 2, 3};
zint1 b{1, 2, 3};
zint1 c{1, 2, 4};
a == b;   // true
a != c;   // true

// setSize -- resize only if dimensions differ
zint1 arr(5);
arr(1) = 42;
arr.setSize(5);   // no-op, data preserved
arr.setSize(3);   // reallocates, data lost
```

### Iterators

JIArray provides `Iterator` (mutable) and `ConstIterator` (read-only), both random-access.

```cpp
zint1 arr{10, 20, 30, 40};

// Range-based for
for (auto& val : arr) {
    val *= 2;
}

// STL algorithms
auto it = std::max_element(arr.begin(), arr.end());
int count = std::count_if(arr.begin(), arr.end(), [](int v){ return v > 30; });

// Const iteration
const zint1& cref = arr;
for (auto it = cref.cbegin(); it != cref.cend(); ++it) {
    // *it is const
}

// Iterator arithmetic
auto it = arr.begin();
it += 2;           // advance by 2
auto val = it[1];  // subscript access
auto diff = arr.end() - arr.begin();  // 4
```

### Loop Macros

```cpp
// ffor -- inclusive range loop (adapts to JIARRAY_OFFSET)
// JIARRAY_OFFSET=1: for (int i = begin; i <= end; ...)
// JIARRAY_OFFSET=0: for (int i = begin; i < end; ...)
ffor(i, 1, 10) {
    // i = 1, 2, ..., 10  (1-based)
}

// ffor_back -- reverse loop
ffor_back(i, 10, 1) {
    // i = 10, 9, ..., 1
}

// zfor -- shorthand for ffor(i, JIARRAY_OFFSET, end)
zfor(i, n) {
    // i = 1..n (1-based) or 0..(n-1) (0-based)
}
```

`end` is evaluated as one complete expression (`ffor(i, 1, c ? a : b)` loops
to `c ? a : b`), and must not be `bool`: the former break-condition idiom
`zfor(i, n && !found)` is a compile error since 0.8.3 -- test the condition
and `break` inside the loop instead:

```cpp
zfor(i, n) {
    if (found) {
        break;
    }
    // ...
}
```

---

## Multi-Dimensional Arrays

### 1D

```cpp
zint1 v(10);
ffor(i, 1, 10) v(i) = i * i;
```

### 2D

```cpp
zdouble2 mat(3, 4);   // 3 rows, 4 columns
mat(2, 3) = 3.14;

// Slice column (column-major)
zdouble1 col = mat.slice(2);   // all rows of column 2
```

### 3D

```cpp
zdouble3 cube(2, 3, 4);
cube(1, 2, 3) = 1.5;

// Slice to 2D
zdouble2 plane = cube.slice(3);      // cube(:, :, 3)

// Slice to 1D
zdouble1 line = cube.slice(2, 3);    // cube(:, 2, 3)
```

### 4D and beyond

```cpp
zint4 arr(2, 3, 4, 5);       // 120 elements
arr(1, 2, 3, 4) = 42;

zint3 sub = arr.slice(4);    // 2x3x4 view
zint2 sub2 = arr.slice(3, 4); // 2x3 view
```

### Memory Layout

**Column-major 3D** array(2, 3, 4) -- strides `[1, 2, 6]`:

```
Linear index = (i-1)*1 + (j-1)*2 + (k-1)*6

Memory order: (1,1,1),(2,1,1), (1,2,1),(2,2,1), (1,3,1),(2,3,1),
              (1,1,2),(2,1,2), (1,2,2),(2,2,2), ...
```

**Row-major 3D** array(2, 3, 4) -- strides `[12, 4, 1]`:

```
Linear index = (i-1)*12 + (j-1)*4 + (k-1)*1

Memory order: (1,1,1),(1,1,2),(1,1,3),(1,1,4), (1,2,1),(1,2,2),(1,2,3),(1,2,4),
              (1,3,1),(1,3,2),(1,3,3),(1,3,4), (2,1,1),(2,1,2), ...
```

### Iteration Order Performance

For column-major, iterate with the first index innermost:

```cpp
zdouble3 arr(NI, NJ, NK);

// Fast (column-major)
for (int k = 1; k <= NK; ++k)
    for (int j = 1; j <= NJ; ++j)
        for (int i = 1; i <= NI; ++i)   // innermost = first dim
            arr(i, j, k) = compute(i, j, k);

// Slow (strided access)
for (int i = 1; i <= NI; ++i)
    for (int j = 1; j <= NJ; ++j)
        for (int k = 1; k <= NK; ++k)
            arr(i, j, k) = compute(i, j, k);
```

For row-major, the opposite: iterate with the last index innermost.

---

## FastArray

Compile-time fixed-size N-dimensional array, stack-allocated, with offset-based indexing.

```cpp
template <typename T, int Offset, std::size_t... Dims>
class FastArrayBase;                       // the implementation

template <typename T, std::size_t... Dims>
using FastArray = FastArrayBase<T, JIARRAY_OFFSET, Dims...>;   // the usual form
```

### Index base

Every dimension is indexed from `Offset`. `FastArray` and the `f<type>` aliases
follow the global `JIARRAY_OFFSET`; `farray0` and the `f<type>0` aliases are
0-based regardless of it; any other base is `FastArrayBase<T, Offset, Dims...>`.
`operator[]` is always 0-based flat access, and `findFirst` returns
`Offset - 1` when not found.

```cpp
fdouble0<8> coef;               // coef(0) .. coef(7): series coefficient by power
FastArrayBase<int, 5, 3> a;     // a(5) .. a(7)
fdouble1d<8, 0> h;              // legacy spelling, also 0-based
```

The trailing `OFFSET` argument of the legacy aliases (`fdouble1d<N, OFFSET>`,
`FastArray2D<T, I, J, OFFSET>`, `StringFastArray<N, OFFSET>`, ...) was silently
ignored in 0.8.0-0.8.1 and is honoured again from 0.8.2.

### Construction

```cpp
FastArray<int, 5> a;                    // zero-initialized
FastArray<double, 3> b(3.14);           // fill with scalar
FastArray<int, 4> c({10, 20, 30, 40});  // initializer list
FastArray<int, 3> d(c);                 // copy (independent)
int raw[] = {7, 8, 9};
FastArray<int, 3> e(raw);               // from C-array
```

### Access

```cpp
FastArray<int, 3> arr({10, 20, 30});

arr(1);       // 10  (1-based via operator())
arr(3);       // 30
arr[0];       // 10  (0-based via operator[])
arr[2];       // 30
arr.front();  // 10
arr.back();   // 30
arr.data();   // int* to underlying array
```

### Operations

```cpp
FastArray<int, 3> arr({10, 20, 30});

arr.findFirst(20);   // 2 (1-based index, or OFFSET-1 if not found)
arr.min();           // 10
arr.max();           // 30

arr += 5;            // {15, 25, 35}
arr *= 2;            // {30, 50, 70}
arr /= 10;          // {3, 5, 7}
auto r = arr * 3;   // {9, 15, 21} (new array)

arr = 42;            // fill: {42, 42, 42}
```

### Iterators

FastArray provides full STL iterator support:

```cpp
FastArray<int, 5> arr({5, 3, 1, 4, 2});

// Forward
for (auto& v : arr) { /* ... */ }

// STL algorithms
std::sort(arr.begin(), arr.end());

// Reverse
for (auto it = arr.rbegin(); it != arr.rend(); ++it) { /* ... */ }

// Const
for (auto it = arr.cbegin(); it != arr.cend(); ++it) { /* ... */ }

arr.size();      // 5
arr.empty();     // false
arr.max_size();  // 5
```

### FastArray2D

2D fixed-size array with column-major layout:

```cpp
template <typename T, std::size_t I, std::size_t J, std::size_t OFFSET = JIARRAY_OFFSET>
using FastArray2D = FastArrayBase<T, OFFSET, I, J>;   // compatibility alias
```

```cpp
FastArray2D<int, 3, 4> mat;          // 3x4, zero-initialized
FastArray2D<double, 2, 2> m(1.0);   // fill with 1.0
FastArray2D<int, 2, 2> m2({1, 2, 3, 4});  // column-major init

m2(1, 1);   // 1
m2(2, 1);   // 2
m2(1, 2);   // 3
m2(2, 2);   // 4

m2 = 0;     // fill with 0
```

### StringFastArray

Fixed-size string array:

```cpp
StringFastArray<3> names({"Alice", "Bob", "Charlie"});
names(1);       // "Alice"
names(2) = "Dave";
names = std::string("default");   // fill all with "default"
```

### FastArray Type Aliases

| Alias | Type | Index base |
|---|---|---|
| `fint<Dims...>`, `fdouble<...>`, `ffloat<...>`, `fbool<...>`, `fchar<...>`, `fshort<...>`, `flong<...>`, `fstring<...>` | `FastArray<type, Dims...>` | `JIARRAY_OFFSET` |
| `farray<T, Dims...>` | `FastArray<T, Dims...>` | `JIARRAY_OFFSET` |
| `fint0<Dims...>`, `fdouble0<...>`, `ffloat0<...>`, `fbool0<...>`, `fchar0<...>`, `fshort0<...>`, `flong0<...>`, `fstring0<...>` | `FastArrayBase<type, 0, Dims...>` | 0 |
| `farray0<T, Dims...>` | `FastArrayBase<T, 0, Dims...>` | 0 |
| `fint1d<N, O>`, `fdouble1d<N, O>`, `ffloat1d<N, O>`, `fbool1d<N, O>`, `fstring1d<N, O>` (deprecated) | `FastArrayBase<type, O, N>` | `O` (default `JIARRAY_OFFSET`) |
| `fint2d<I, J, O>`, `fdouble2d<I, J, O>`, `ffloat2d<I, J, O>`, `fbool2d<I, J, O>` (deprecated) | `FastArrayBase<type, O, I, J>` | `O` (default `JIARRAY_OFFSET`) |
| `FastArray2D<T, I, J, O>`, `StringFastArray<N, O>` (compatibility) | `FastArrayBase<T, O, ...>` | `O` (default `JIARRAY_OFFSET`) |

---

## JIVector

Dynamic vector with 1-based indexing, inheriting from `std::vector`.

```cpp
template <typename T, typename Alloc = std::allocator<T>>
class JIVector : public std::vector<T, Alloc>;
```

### Construction and Access

```cpp
JIVector<int> v;                        // empty
JIVector<int> v2{10, 20, 30};           // initializer list
JIVector<double> v3(std::vector<double>{1.0, 2.0});  // from std::vector

v2[1];     // 10  (1-based)
v2[3];     // 30
v2.at(2);  // 20  (1-based, throws std::out_of_range if invalid)
v2(2);     // 20  (same as at())

// Memory access
v2.getMemory(0);    // pointer to first element
v2.get_pointer();   // same as std::vector::data()
```

### Search

```cpp
JIVector<int> v{10, 20, 30, 20};

v.contains(20);   // true
v.contains(99);   // false

v.find(20);       // 2 (1-based index of first occurrence)
v.find(99);       // 0 (JIARRAY_OFFSET - 1, not found)
```

### Insert

```cpp
JIVector<int> v{1, 2, 3};

v.insert(v.end(), {4, 5});           // from initializer list
v.insert(v.end(), other_jivector);   // from JIVector
v.insert(v.end(), std_vector);       // from std::vector
v.insert(v.end(), jiarray_1d);       // from JIArray<T, 1>
v.insert(v.begin(), first, last);    // from iterator range
```

### Assignment

```cpp
JIVector<int> v;

std::vector<int> sv{1, 2, 3};
v = sv;                 // from std::vector

zint1 arr{10, 20, 30};
v = arr;                // from JIArray (appends elements via push_back)
```

### JIVector Type Aliases

| Alias | Type |
|---|---|
| `zvector<T>` | `JIVector<T>` |
| `zints` | `JIVector<int>` |
| `zdoubles` | `JIVector<double>` |
| `zfloats` | `JIVector<float>` |
| `zstrings` | `JIVector<std::string>` |

---

## JIDeviceArray

`JIDeviceArray<T, RANK, Index>` (`#include <jiarray/JIDeviceArray.h>`) is a CUDA
device-memory array with JIArray shape, offsets and indexing. It has two halves:

- **`JIDeviceArray`** -- the host-side owner. Allocates device-only memory
  (`cudaMalloc`, or stream-ordered `cudaMallocAsync`), RAII, **move-only** (copy is
  deleted), no virtual. It moves data to and from host `JIArray`s. It never
  dereferences device memory on the host.
- **`JIDeviceView`** (`JIDeviceArray::View`, from `view()`) -- a non-owning,
  **trivially copyable** view passed to kernels **by value**. In kernels it is
  indexed with JIArray syntax: `a(i, j, k)`, `slice`, `reshape`, `zfor`/`ffor`.

Unlike `JICudaArray` (managed memory, used by montex), there is no host element
access and no page migration. Element arithmetic, statistics, search and
expression assignment are not provided -- write kernels.

The header works under nvcc and also under a host-only compiler that sees
`cuda_runtime.h` (link `CUDA::cudart`). `fillAsync` launches a kernel, so it
needs nvcc. If the CUDA runtime header is not found, the header is empty.

### Device Type Aliases

| Alias | Type |
|---|---|
| `dint1` .. `dint5` | `JIDeviceArray<int, 1>` .. `JIDeviceArray<int, 5>` |
| `dfloat1` .. `dfloat5` | `JIDeviceArray<float, 1>` .. `JIDeviceArray<float, 5>` |
| `ddouble1` .. `ddouble6` | `JIDeviceArray<double, 1>` .. `JIDeviceArray<double, 6>` |
| `darray<T, N, Index>` | `JIDeviceArray<T, N, Index>` |
| `X::View` / `X::ConstView` | `JIDeviceView<T, N, Index>` / `JIDeviceView<const T, N, Index>` |

`d` stands for device. The `c*` names (`cdouble2`, ...) are already global macros
of `JICudaArray.h`, and `f*` are FastArray, so a `c`/`f` prefix would collide.
These are namespaced `using` aliases, like the `z*` ones.

### Allocation

```cpp
ddouble3 flux(ng, nxy, nz);       // cudaMalloc, zero-filled, indices start at JIARRAY_OFFSET
dint2 map;
map.init0(0, 2, -1, 3);           // map(0..2, -1..3), as JIArray::init0
map.setOffsets(1, 1);             // re-base; views taken afterwards see the new offsets
dfloat2 work;
work.initAsync(stream, 216, n);   // cudaMallocAsync + cudaMemsetAsync on `stream`;
                                  // released with cudaFreeAsync on the same stream
ddouble3 owner = std::move(flux); // ownership moves; flux is empty
owner.destroy();                  // explicit release (the destructor does the same)
```

New memory is zero-filled, matching `JIArray::init`. Sizes and the offset-folded
index range must fit in `Index`. If they do not, `init` throws `std::length_error`;
use a 64-bit `Index` instead.

### Host-Device Transfer

```cpp
zdouble3 host(ng, nxy, nz);
ddouble3 dev;
dev = host;                   // allocate to host's sizes AND offsets, upload (blocking)
dev.copyTo(host);             // download (blocking); an empty host array is allocated

dfloat3 dev32;
dev32.uploadAs(host);         // FP64 -> FP32 conversion on upload (host staging buffer)
zdouble3 back;
dev32.copyTo(back);           // FP32 -> FP64 on download (copyTo converts when types differ)

dev.copyFromAsync(host, stream);  // cudaMemcpyAsync; shapes must already match
dev.copyToAsync(host, stream);    // read `host` only after the stream is synchronized
dev.fillAsync(1.0, stream);       // grid-stride fill kernel (nvcc only)
```

- If the target is already allocated, its sizes and offsets must equal those of the
  source. Otherwise `std::invalid_argument` is thrown. This check always runs; it
  does not depend on `JIARRAY_DEBUG`.
- An empty host source empties the device array, as JIArray copy assignment does.
- The async copies overlap with host work only when the host array lives in pinned
  memory.

### Kernel Views and Indexing

```cpp
__global__ void scale(ddouble3::View flux, ddouble1::ConstView factor) {
    const int ixy = blockIdx.x * blockDim.x + threadIdx.x + 1;  // 1-based
    if (ixy > flux.getSize(2)) return;
    zfor (iz, flux.getSize(3)) {
        zfor (ig, flux.getSize(1)) {
            flux(ig, ixy, iz) *= factor(ig);
        }
    }
}
scale<<<grid, block>>>(flux.view(), factor.view());
```

Views use the same names and meanings as JIArray:

| Member | Meaning |
|---|---|
| `a(i, j, ...)`, `a.at(...)` | element at offset-relative indices |
| `a(FastArray<int, RANK>)` | element at a FastArray of indices |
| `a[n]` | 0-based flat element |
| `a.data()`, `a.data(i, ...)` | device pointer to the start / to the partial-index position (as `JIArray::data(INDEX...)`) |
| `size()`, `getSize()`, `getSize(r)` | element count / extent of dimension `r` (numbered from `JIARRAY_OFFSET`) |
| `getSizeOfRank()`, `getRankSize()`, `getOffset()`, `getOffset(r)` | extents, strides, first indices |
| `slice(i, ...)` | lower-rank view. Column-major fixes the **trailing** dimensions (`a.slice(k)` = `a(:, :, k)`), row-major the **leading** ones, exactly as JIArray |
| `reshape(n1, n2, ...)` | same elements, new shape, default offsets |
| `setOffsets(...)` | re-base this view |

`slice`, `reshape` and `setOffsets` only change metadata. They also work inside
kernels (`ddouble1::View line = a.slice(j, k);`). The owner has the same queries,
plus `slice`/`reshape`/`data(i, ...)`, which return views or device pointers. To
wrap an external device buffer, use `ddouble2::View v(device_ptr, n1, n2);`.

**Constness follows `std::span`, not JIArray.** A `const` view object still writes
its elements. For read-only access, use `ConstView` (`JIDeviceView<const T>`).
That is what `view()` on a `const` owner returns, and `View` converts to it
implicitly.

**Lifetime.** A view never frees memory. It must not outlive its owner, and it must
not be used after `init`/`destroy` on the owner reallocates or frees the memory
(the same rule as JIArray views).

### Indexing Cost

`a(i, j, k)` compiles to the same address arithmetic as a hand-written 0-based raw
index `p[i + nx * (j + ny * k)]`:

- When the view is built, the offsets are folded into a shifted base pointer
  (`origin = data - Σ offset·stride`). No `- offset` term remains in the access.
- Every view is contiguous: allocation, `slice` and `reshape` all preserve this.
  So the fastest dimension has unit stride, and its multiply is dropped at compile
  time.
- The linear index is computed in `Index` (`int` by default, 32-bit
  multiply-adds). For arrays beyond 2^31 elements, use
  `JIDeviceArray<T, RANK, long long>`.

Measurements on an RTX 4090 (CUDA 12.6, sm_89). Kernel times are the medians of
51 interleaved launches, taken over 3 runs of `dnegri.jiarray.device_array_bench`.
Instruction counts come from `cuobjdump -sass`:

| Kernel (FP64) | raw `__restrict__` | raw | view | view + `__ldg` |
|---|---|---|---|---|
| 7-point stencil 256³ | 0.287 ms / 50 SASS | 0.290 ms / 50 | 0.292 ms / 53 | 0.290 ms / 57 |
| 216×216 GEMV × 20000 nodes | 1.459 ms / 191 | 1.459 ms / 191 | 1.457 ms / 200 | -- |

A single access costs the same as raw (two `IMAD` + one `IMAD.WIDE` for rank 3).
In the stencil, the remaining few instructions come from `out` carrying its own
strides. The raw kernel reuses `in`'s index `c` for `out[c]`; the view kernel
cannot know that the two shapes match. Kernel arguments cannot carry
`__restrict__` on struct members, so loads do not take the read-only path by
default. Where that matters, write `__ldg(&a(i, j, k))`.

**Compile-time extents (not implemented).** A fixed extent such as the 216 face
moments could be given as a template argument
(`JIDeviceView<T, 2, int, Extent<216, kDynamic>>`), with strides that are
compile-time constants. The measured gain was small. `gemvRawFixed` (stride 216
as a literal) runs 0.3 % faster than the same kernel with the stride as a
runtime argument. On the GPU a runtime stride is a constant-bank operand of the
same `IMAD`, so the template is deferred until a kernel shows a real need.

### Debug Checks

With `JIARRAY_DEBUG`, indices, `getSize(r)`/`getOffset(r)`, `slice` indices and
`reshape` sizes are checked:

- On the host, these checks throw the same exceptions as JIArray
  (`std::out_of_range`, `std::invalid_argument`).
- In a kernel, a failed check prints the index and its range, then calls
  `__trap()`. The kernel aborts, and the next synchronization returns
  `cudaErrorLaunchFailure`.

Without `JIARRAY_DEBUG`, these checks compile to nothing.

### CUDA Tests

```bash
cmake -S . -B build-cuda -DJIARRAY_WITH_CUDA_TESTS=ON   # CMAKE_CUDA_ARCHITECTURES defaults to native
cmake --build build-cuda -j 8
ctest --test-dir build-cuda -R DeviceArray --output-on-failure
./build-cuda/dnegri.jiarray.device_array_bench          # raw vs view timing
```

This builds three test executables:

- `device_array_test`: the column-major build.
- `device_array_rowmajor_test`: the same source built row-major.
- `device_array_host_test`: the header built by the host compiler only.

The default (non-CUDA) build is unchanged.

---

## Integration

### CMake

**FetchContent** (recommended):

```cmake
include(FetchContent)
FetchContent_Declare(
  jiarray
  GIT_REPOSITORY https://github.com/dnegri/jiarray.git
  GIT_TAG main
)
FetchContent_MakeAvailable(jiarray)

target_link_libraries(your_target PRIVATE jiarray)
```

**add_subdirectory**:

```cmake
add_subdirectory(external/jiarray)
target_link_libraries(your_target PRIVATE jiarray)
```

Tests are automatically disabled when jiarray is not the top-level project. Override with:

```cmake
set(JIARRAY_BUILD_TESTS ON CACHE BOOL "" FORCE)
```

**Manual include path**:

```cmake
target_include_directories(your_target PRIVATE path/to/jiarray/include)
```

Or in your compiler command:

```bash
g++ -std=c++17 -Ipath/to/jiarray/include your_code.cpp
```

### HighFive HDF5 Extension

`HighFiveExtension.hpp` provides `HighFive::inspector` specializations for JIArray and JIVector, enabling direct HDF5 file I/O:

```cpp
#include <jiarray/HighFiveExtension.hpp>

HighFive::File file("data.h5", HighFive::File::ReadWrite | HighFive::File::Create);

// Write
zdouble2 mat(100, 200);
// ... fill mat ...
file.createDataSet("matrix", mat);

// Read
zdouble2 loaded(100, 200);
file.getDataSet("matrix").read(loaded);
```

Requires [HighFive](https://github.com/BlueBrain/HighFive) to be available.
