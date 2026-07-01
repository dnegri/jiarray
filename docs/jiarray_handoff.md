# jiarray — Evaluation & Upgrade Handoff

Owner: dnegri (Jooil). Header-only C++17 N-dimensional array library for scientific
computing (Fortran ergonomics: 1-based indexing, column-major, `ffor`, zero-copy slicing,
SIMD element-wise ops, HDF5 I/O, optional CUDA path). This doc is the working brief for a
Claude Code session that will **officially evaluate and then upgrade** the library.

> Line refs below are approximate (captured from a prior review). **Re-verify against the
> current tree before editing.** Every performance claim here was measured, not assumed —
> hold yourself to the same standard (see Methodology).

---

## 1. Repository facts (as reviewed)

- Header-only, C++17. `include/jiarray/`: `pch.h` (~89 L), `JIArray.h` (~2068 L),
  `FastArray.h` (~301 L, stack fixed-size), `JIVector.h` (~160 L, 1-based dynamic vector),
  `JICudaArray.h` (~1199 L, CUDA path via `JIARRAY_HD`), `HighFiveExtension.hpp` (HDF5).
- Compile-time config (in `pch.h`): `JIARRAY_OFFSET=1` (1-based default),
  `JIARRAY_COLUMN_MAJOR=1` (Fortran default), `JIARRAY_HD` = `__host__ __device__` under NVCC.
- Type names are **macros**: `#define zdouble1 JIArray<double,1>` … `zint1..5`, `zbool1..5`
  (JIArray.h ~L1963+).
- Tests: Google Test + `ctest`; CMake (`FetchContent`, `add_subdirectory`). Optional cereal
  serialization (`JIARRAY_CEREAL`).
- **Version inconsistency to reconcile:** release tag `v0.4.0` (Dec 2024) but
  `RELEASE_NOTES_0.7.1.md` present in tree.

---

## 2. Findings already established this session (source-verified + benchmarked)

### F1 — Element-wise arithmetic is EAGER with materialized temporaries (no expression templates)
Binary `operator+ - * /` (JIArray.h ~L1246+) each allocate a fresh array via
`initByRankSize(...)`, run one `#pragma omp simd` pass, and return by value:
```cpp
this_type operator+(const this_type& a) const {
    this_type result; result.initByRankSize(...);   // heap alloc (temporary)
    JIARRAY_SIMD_LOOP for (int i=0;i<nn;++i) result.mm[i] = mm[i] + a.mm[i];
    return result;
}
```
There is **no cross-operator fusion**. Compound-assignment `+= -= *= /=` is in-place
(single pass, no alloc) and is the current fast path.

### F2 — Bounds checking is OPT-IN via `JIARRAY_DEBUG`, decoupled from `NDEBUG`
`pch.h` ~L63:
```cpp
#if defined(JIARRAY_DEBUG) && !defined(__CUDA_ARCH__)
  #define JIARRAY_CHECK_BOUND(i,beg,end) if (i<beg||i>end) throw std::out_of_range(...)
#else
  #define JIARRAY_CHECK_BOUND(i,beg,end)      // no-op
#endif
```
Opposite polarity from `assert`: **off by default** (zero overhead in release regardless of
NDEBUG), on only with `-DJIARRAY_DEBUG`. Uses `throw` (not `assert`). Device code is never
checked. Consequence: correctness safety only exists if a debug config explicitly defines it.

### F3 — SIMD requires `-fopenmp`
`JIARRAY_SIMD_LOOP` = `#pragma omp simd` only when `_OPENMP` is defined; MSVC falls back to
`loop(ivdep)`; otherwise the pragma is empty (auto-vec only). `JIARRAY_USE_SIMD` defaults to
1; set `-DJIARRAY_USE_SIMD=0` to suppress (e.g. when calling ops inside an existing OpenMP
parallel region). To get advertised SIMD: build `-fopenmp -O3 -march=native`.

### F4 — Benchmark: the operator sugar is a 4.6× hot-path trap; the container itself is zero-overhead
`d = a + b + c`, n=16M doubles (128 MB/array), `-O3 -march=native -fopenmp`, bounds off,
single thread (~13 GB/s VM ceiling), best-of-7:

| variant | time (ms) | eff GB/s | vs fused |
|---|---|---|---|
| (A) operator chain `d=a+b+c` | 3698.8 | 2.8 | — |
| (B) `ffor` fused | 808.5 | 12.7 | **4.57× faster** |
| (C) in-place `d=a; d+=b; d+=c` | 1288.4 | 7.9 | 2.87× faster |
| (D) raw ptr fused (ceiling) | 791.9 | 12.9 | 4.67× faster |

Key reads: (B) ≈ (D) within 2% → **1-based `operator()`→`at()` compiles to raw-pointer code
when bounds are off; the container/indexing abstraction is genuinely zero-overhead.** The
entire 4.6× gap is the operator temporaries — dominated by **repeated 128 MB heap alloc +
page faults per expression**, not merely extra memory passes (naive bandwidth model predicts
~1.5×; measured 4.6×). The benchmark harness is `bench.cpp` (ship it into `benchmark/`).

---

## 3. Evaluation plan (audit before upgrading)

Produce `docs/EVALUATION.md` scoring each area with evidence.

1. **Correctness & safety** — test coverage %; rank/type safety (already good via
   `static_assert`); const-correctness; **view/slice lifetime** (slices are non-owning raw
   pointers → dangling-view risk: audit and document the ownership model); behavior on
   empty/degenerate/zero-size arrays; offset and row/col-major variants all tested?
2. **Performance** — reproduce F4; audit every element-wise op for temporaries; SIMD
   effectiveness (inspect asm on Compiler Explorer / `-fopt-info-vec`); **alignment**
   (default allocator gives 16 B; AVX-512 wants 64 B → aligned load opportunity); slice
   **aliasing** inhibiting vectorization; allocation strategy & first-touch/NUMA behavior.
3. **API design** — macro type names (`#define zdouble1`) leak into the global namespace and
   can't be scoped/forward-declared → hazard; `ffor` macro hygiene; naming/consistency;
   documentation of the 1-based + column-major contract.
4. **Modernization gap** — C++17 today; opportunities: `std::mdspan` interop (ties directly
   to slicing/`submdspan`/layout), `std::span`, concepts vs SFINAE, `std::assume_aligned`,
   ranges, three-way comparison.
5. **Portability & build** — compiler matrix (GCC/Clang/MSVC/NVCC) actually green? CMake
   modernity; presence of CI, sanitizers (ASan/UBSan), warnings-as-errors.
6. **GPU story** — `JICudaArray.h` is CUDA-only and partly parallel to `JIArray.h`. Audit
   duplication/divergence; decide long-term direction (bespoke CUDA vs `mdspan`/Kokkos
   interop vs unified host/device type).
7. **Release hygiene** — reconcile v0.4.0 vs 0.7.1; adopt semver + `CHANGELOG.md`; API docs
   (the headers already have Doxygen comments — wire up doc generation).

---

## 4. Upgrade roadmap (prioritized; evidence first)

**P0 — measured, highest impact**
- **Fix the arithmetic temporary problem (F4).** Evaluate two routes and pick with data:
  (a) **expression templates** so `d=a+b+c` fuses to one allocation-free pass — proper fix,
  but adds template complexity and must stay `JIARRAY_HD`-clean for the CUDA path and
  debuggable; (b) keep eager but add a fused-expression/`eval()` layer + preallocated
  workspace ops + clear docs steering hot loops to `ffor`/`+=`. Acceptance test: operator
  chain must come within ~1.2× of `ffor` fused on `bench.cpp`. Do **not** regress the
  in-place path.
- **Alignment**: optional 64 B aligned allocator + `std::assume_aligned` in kernels; measure
  AVX-512 delta.

**P1**
- Replace `#define zdouble1 …` with namespaced `using` aliases (keep macros as a
  deprecated back-compat shim for one release).
- **`std::mdspan` interop**: `.to_mdspan()` and construct-from-mdspan; evaluate backing
  slices on `layout_stride` / exposing `submdspan`.
- **Slice lifetime safety**: document ownership; optional debug-mode guard against dangling
  views; consider owning vs non-owning slice distinction.
- Ensure element-wise kernels carry `__restrict`/no-alias guarantees so slice operands
  vectorize.

**P2**
- C++20/23 modernization (concepts, `std::span`, ranges, `[[likely]]`) behind a standard
  detection so C++17 users keep working.
- CI: GitHub Actions matrix (GCC/Clang/MSVC + NVCC build-only), ASan/UBSan test runs,
  `-Wall -Wextra -Werror`, and a benchmark-regression job seeded from `bench.cpp`.
- GPU direction decision + de-duplication of `JICudaArray.h` vs `JIArray.h`.
- Release hygiene: version reconcile, semver, `CHANGELOG.md`, Doxygen site.
- Coverage: raise and add edge cases (offset variants, row vs column major, slice
  correctness, degenerate shapes).

---

## 5. Methodology & guardrails (non-negotiable)

- **Measure, don't assert.** Every perf change is validated with `bench.cpp` (and asm
  inspection where relevant); record before/after numbers in the PR. No perf regressions on
  the existing fast paths (`+=`, fused, indexing).
- **Baseline first.** Before editing: build all configs, run `ctest`, run ASan+UBSan, run
  `bench.cpp`, and capture the numbers as the reference. 
- **Incremental branches.** One concern per branch; each lands with tests green on the
  compiler matrix + a benchmark delta note.
- **Preserve the public API.** 1-based + column-major defaults, existing type names, and
  operator semantics stay working. Deprecate before removing; bump semver accordingly.
- **CUDA-clean.** Anything touched in shared code paths must stay `JIARRAY_HD`-compatible
  and compile under NVCC.
- Build/bench command of record:
  `g++ -std=c++17 -O3 -march=native -fopenmp -DNDEBUG -I include bench.cpp -o bench`
  (add `-DJIARRAY_DEBUG` only for the safety/bounds test config).

---

## 6. Kickoff prompts for Claude Code

**Phase 1 — evaluation**
```
Read this handoff and the jiarray source. Reproduce the F4 benchmark with bench.cpp on this
machine and record the numbers. Then produce docs/EVALUATION.md scoring the seven areas in
section 3 with concrete evidence (file:line, asm, or measurements) — findings and severity
only, no code changes yet. Confirm or correct F1–F4 against the current tree.
```

**Phase 2 — upgrade (after we review the evaluation)**
```
Implement the P0 items from the upgrade roadmap on a branch. Establish the baseline first
(ctest + ASan/UBSan + bench.cpp). For the arithmetic fix, prototype both routes (expression
templates vs fused/workspace), benchmark each against bench.cpp, and recommend one with
data. Every change: tests green on GCC+Clang, CUDA still compiles, benchmark delta recorded,
public API preserved. Do not proceed to P1 until P0 is reviewed.
```

