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

**Custom offsets per dimension** with `init0()`:
```cpp
zdouble2 arr;
arr.init0(5, 10, 20, 25);   // dim1: [5..10], dim2: [20..25]
arr(5, 20) = 1.0;           // valid
arr(10, 25) = 2.0;          // valid
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
arr.init0(5, 10, 20, 25);     // dim1=[5..10], dim2=[20..25]

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

All arithmetic operations support optional SIMD vectorization via `JIARRAY_USE_SIMD`.

**Array-array operations** (element-wise, same size required):

```cpp
zdouble1 a{1.0, 2.0, 3.0};
zdouble1 b{4.0, 5.0, 6.0};

auto c = a + b;     // {5.0, 7.0, 9.0}
auto d = a - b;     // {-3.0, -3.0, -3.0}
auto e = a * b;     // {4.0, 10.0, 18.0}
auto f = a / b;     // {0.25, 0.4, 0.5}

a += b;             // a = {5.0, 7.0, 9.0}
a -= b;             // a = {1.0, 2.0, 3.0}
a *= b;             // a = {4.0, 10.0, 18.0}
a /= b;             // a = {1.0, 2.0, 3.0}
```

**Scalar operations** (both orderings):

```cpp
zdouble1 arr{1.0, 2.0, 3.0};

auto r1 = arr + 10.0;    // {11.0, 12.0, 13.0}
auto r2 = 10.0 + arr;    // {11.0, 12.0, 13.0}
auto r3 = arr * 3.0;     // {3.0, 6.0, 9.0}
auto r4 = 3.0 * arr;     // {3.0, 6.0, 9.0}
auto r5 = arr / 2.0;     // {0.5, 1.0, 1.5}
auto r6 = 6.0 / arr;     // {6.0, 3.0, 2.0}

arr += 10.0;              // {11.0, 12.0, 13.0}
arr *= 2.0;               // {22.0, 24.0, 26.0}
arr /= 2.0;               // {11.0, 12.0, 13.0}
```

**Unary negation**:

```cpp
zdouble1 arr{1.0, -2.0, 3.0};
auto neg = -arr;          // {-1.0, 2.0, -3.0}
```

**Note**: For floating-point `/=` by scalar, the library uses reciprocal multiplication (`*= 1/val`) for performance.

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

Compile-time fixed-size array, stack-allocated, with offset-based indexing.

```cpp
template <class T, std::size_t SIZE, std::size_t OFFSET = JIARRAY_OFFSET>
class FastArray;
```

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
template <class T, std::size_t SIZE1, std::size_t SIZE2, std::size_t OFFSET = JIARRAY_OFFSET>
class FastArray2D;
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

| Alias | Type |
|---|---|
| `fint1d<N>` | `FastArray<int, N>` |
| `fdouble1d<N>` | `FastArray<double, N>` |
| `ffloat1d<N>` | `FastArray<float, N>` |
| `fbool1d<N>` | `FastArray<bool, N>` |
| `fstring1d<N>` | `StringFastArray<N>` |
| `fint2d<I, J>` | `FastArray2D<int, I, J>` |
| `fdouble2d<I, J>` | `FastArray2D<double, I, J>` |
| `ffloat2d<I, J>` | `FastArray2D<float, I, J>` |
| `fbool2d<I, J>` | `FastArray2D<bool, I, J>` |
| `farray<T, N>` | `FastArray<T, N>` |

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
