# JIArray

High-performance N-dimensional array library for C++ scientific computing.

Header-only, C++17, with Fortran-style column-major layout and 1-based indexing by default.

## Features

- **Header-only** -- just add `include/` to your include path
- **N-dimensional** -- 1D through 6D+ with compile-time rank safety
- **Flexible storage order** -- column-major (Fortran) or row-major (C), compile-time selectable
- **Flexible indexing** -- 1-based (default) or 0-based, custom offsets per dimension
- **Slicing** -- zero-copy view semantics, slice ND to (N-1)D
- **Arithmetic** -- element-wise `+`, `-`, `*`, `/` with arrays and scalars, SIMD-accelerated
- **Statistics** -- `sum()`, `average()`, `min()`, `max()`, `sqsum()`, `dot()`
- **Search** -- `contains()`, `findFirst()`, `maxloc()`
- **STL-compatible** -- random-access iterators, range-based for loops, `std::sort` etc.
- **Fixed-size arrays** -- `FastArray<T, N>` for stack-allocated compile-time-sized arrays
- **Dynamic vector** -- `JIVector<T>` with 1-based indexing
- **HDF5 integration** -- HighFive extension for file I/O

## Quick Start

### CMake FetchContent

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

Tests are automatically disabled when fetched as a dependency. To force enable/disable:

```cmake
set(JIARRAY_BUILD_TESTS OFF)  # or ON
FetchContent_MakeAvailable(jiarray)
```

### add_subdirectory

```cmake
add_subdirectory(external/jiarray)
target_link_libraries(your_target PRIVATE jiarray)
```

### Manual

Copy the `include/jiarray/` directory into your project and add the parent to your include path.

## Usage

```cpp
#include <jiarray/JIArray.h>
using namespace dnegri::jiarray;

// 1D array with initializer list
zint1 v{10, 20, 30, 40, 50};
v(3) = 99;                    // 1-based indexing
int s = v.sum();              // 199

// 2D array (3x4 matrix)
zdouble2 mat(3, 4);
mat = 0.0;
mat(2, 3) = 3.14;

// Slicing -- zero-copy view
zdouble1 col = mat.slice(2);  // column 2 (column-major)

// Arithmetic
zdouble1 a{1.0, 2.0, 3.0};
zdouble1 b{4.0, 5.0, 6.0};
auto c = a + b;               // {5.0, 7.0, 9.0}
double d = dot(a, b);         // 32.0

// Iteration
for (auto& val : mat) {
    val *= 2.0;
}

// Fortran-style loop macro (inclusive bounds)
ffor(i, 1, 5) {
    v(i) += 1;
}
```

See [docs/manual.md](docs/manual.md) for the complete API reference.

## Build & Test

```bash
mkdir build && cd build
cmake ..
cmake --build .
ctest
```

## Requirements

- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.21+ (for `PROJECT_IS_TOP_LEVEL`)
- Google Test (optional, for tests)

## License

MIT
