# jiarray — Evaluation (Phase 1, findings only)

Scope: audit the seven areas in `docs/jiarray_handoff.md` §3 and confirm/correct F1–F4.
**No code changes** — findings + severity + evidence only. Companion to the handoff brief.

Severity scale: **Critical** (UAF/data loss possible) · **High** (major perf/correctness trap
hit on normal use) · **Medium** (real hazard, workaround exists) · **Low** (polish/enhancement)
· **Info** (no action, recorded for completeness).

---

## Environment of record

| | |
|---|---|
| CPU | AMD Ryzen Threadripper PRO 7985WX (64c), AVX-512 (F/DQ/BW/VL/VBMI/IFMA/CD) |
| Compiler | GCC 11.5.0 (Red Hat 11.5.0-11) |
| Build | `g++ -std=c++17 -O3 -march=native -fopenmp -DNDEBUG -I include` |
| Bench run | `taskset -c 0 ./bench` (pinned single core), best-of-7, n=16M doubles |
| Tree | HEAD `0d07b8a`, `git describe` = **0.8.0** |

This is **bare-metal**, not the "~13 GB/s VM" of the handoff. Absolute GB/s below are **not
comparable** to the handoff's table; only the **ratios** carry over (and they are larger here).

---

## F1–F4 confirmation

### F1 — Eager arithmetic with materialized temporaries — **CONFIRMED (exact)**
Every binary/scalar `+ - * /` allocates a fresh result via `initByRankSize(...)`, runs one
SIMD pass, returns by value; no cross-operator fusion.
- `operator+` alloc+loop: `JIArray.h:1246-1255` (alloc `:1249`)
- `operator-` `:1278`, `operator*` `:1340`, `operator/` (friend) `:1396`
- scalar `+ * /`: `:1419`, `:1447`, `:1473`, `:1494` — all `initByRankSize` + one pass
- fast path `+= -= *= /=` is in-place, single pass, no alloc: `:1217`, `:1231`, `:1263`, `:1325`, `:1380`
Line numbers in the handoff (~L1246) are accurate against the current tree.

### F2 — Bounds checks opt-in via `JIARRAY_DEBUG`, decoupled from `NDEBUG` — **CONFIRMED**
`pch.h:63-81`: `JIARRAY_CHECK_BOUND/RANK/SIZE/NOT_ALLOCATED` expand to `throw` only under
`JIARRAY_DEBUG && !__CUDA_ARCH__`; no-ops otherwise. Off by default regardless of `NDEBUG`;
device code never checked. Tests compile **with** `-DJIARRAY_DEBUG` (`CMakeLists.txt:63`), so
the checked path is exercised. Accurate as written.

### F3 — SIMD requires `-fopenmp` — **CONFIRMED**
`pch.h:44-61`: `JIARRAY_SIMD_LOOP` = `#pragma omp simd` only under `_OPENMP`; MSVC →
`loop(ivdep)`; otherwise empty (auto-vec only). `JIARRAY_USE_SIMD` defaults to 1
(`JIArray.h:47-56`); `-DJIARRAY_USE_SIMD=0` neutralizes the pragmas. Accurate.

### F4 — Operator sugar is a hot-path trap; container is zero-overhead — **REPRODUCED, ratios larger**
`d=a+b+c`, n=16M, single core, best-of-7 (two runs shown; stable):

| variant | run1 ms | run2 ms | eff GB/s | vs fused |
|---|---|---|---|---|
| (A) operator chain `d=a+b+c` | 1633.7 | 1639.6 | 6.3 | — |
| (B) `ffor` fused | 207.6 | 210.1 | 49.3 | **7.87×** |
| (C) in-place `d=a; d+=b; d+=c` | 372.4 | 373.7 | 27.5 | 4.39× |
| (D) raw ptr fused (ceiling) | 207.2 | 207.2 | 49.4 | 7.88× |

Qualitative claims **hold and are stronger** here:
- **(B) ≈ (D) within <1.5%** → 1-based `operator()`→`at()` compiles to raw-pointer code with
  bounds off; the container/indexing abstraction is genuinely **zero-overhead**. Confirmed.
- The operator-chain gap is **7.9× vs the handoff's 4.6×**. On faster hardware the fused kernel
  gets cheaper while the eager path stays dominated by **per-expression 128 MB alloc + first-touch
  page faults** (fixed cost that doesn't shrink with bandwidth) — so the relative penalty grows.
- **(C) in-place is 4.4× faster than (A)** and remains the correct guidance for hot loops that
  must use array ops.

**Correction to the handoff:** the "single thread (~13 GB/s VM ceiling)" note is environment-
specific; here single-thread streaming tops ~49 GB/s. Report ratios, not absolute GB/s, when
comparing machines.

---

## Area 1 — Correctness & safety — **overall good; one Medium lifetime hazard**

- **Test suite: 413 tests, 100% pass** in ~1.0 s (`ctest`), built with `JIARRAY_DEBUG` on.
  Both layouts covered (`column_major_test` + `row_major_test`, `fast_array_test` +
  `fast_array_rowmajor_test`; `CMakeLists.txt:70-84`). Suites: coverage, slice, move-lifecycle,
  ffor, jivector, iter. **[Info — healthy]**
- **ASan clean**: `move_lifecycle_test` (40 tests) rebuilt `-fsanitize=address` — no
  leaks/UAF/double-free. (UBSan not runnable in this env: `libubsan.so.1.0.0` missing from the
  toolchain — a *host* gap, not a code finding.) **[Info]**
- **Rank/type safety**: strong via `static_assert` on rank match and `all_integral_v`
  (`JIArray.h:219-220`, `126`). **[Info — good]**
- **View/slice lifetime (Medium).** Copy ctor is a **non-owning shallow view**
  (`JIArray.h:349-356`, `allocated=NONE`), and `slice()`/`reshape()` return views
  (`:676`, `:702`, `:957`). A view (incl. `auto x = arr;`) outliving its owner is a silent
  **use-after-free** with **no debug-mode guard**. Documented in `CLAUDE.md` but unenforced at
  runtime. **Severity: Medium.**
- **Empty/degenerate handling (Low, good).** `init()` early-returns on any zero dimension
  (`:470-473`); copy-assign resets to empty when source is null/zero (`:1059-1062`);
  `operator=(const T&)`, `max()`, `min()` guard `assert(nn>0)` (`:992`, `:1170`, `:1179`) — but
  these are `assert`, i.e. compiled out under `NDEBUG` → calling `max()` on an empty array in
  release is UB. **Severity: Low.**
- **const-correctness (Low).** `const` `slice()` returns a **`const` value**
  (`JIArray.h:702`) which inhibits move and is an unusual signature; it also relies on an internal
  `const_cast` of `mm` (`:713`). Works, but the const view does not actually prevent writing
  through the returned object's non-const members. **Severity: Low.**

## Area 2 — Performance — **High (the operator-temporary trap); alignment payoff over-stated**

- **Headline (High):** F4 above — eager operator temporaries cost **7.9×** vs fused. This is the
  P0 target and the finding is confirmed with stronger numbers.
- **SIMD works.** `-fopt-info-vec` shows the op loops vectorize: `JIArray.h:1235` (`+=`),
  `:1252` (`+`), and the fused bench loop (`bench.cpp:42`). **[Info — good]**
- **Vector width / alignment — payoff is ~0 for streaming (corrects handoff P0).** With
  `-march=native` on AVX-512 hardware GCC 11 still emits **32-byte (AVX2/256-bit)** vectors
  (znver default `-mprefer-vector-width=256`). Forcing `-mprefer-vector-width=512` changed (B)
  from 210.1 → 210.9 ms — **no measurable gain**. `new double[]` returns **16-byte** alignment
  (measured; AVX-512 aligned load wants 64 B). **But** these ops are **memory-bandwidth-bound**,
  so neither 512-bit vectors nor 64-byte alignment help the streaming case. The handoff's
  "measure AVX-512 delta" for alignment should expect **≈0** here; alignment only pays off for
  **compute-bound / cache-resident** kernels, which this library's element-wise ops are not.
  **Severity: Low** (do it for correctness/other archs, not for streaming speed).
- **Aliasing inhibits nothing critical but is real (Low).** The fused loop through `operator()`
  triggers *"loop versioned for vectorization because of possible aliasing"* (`bench.cpp:42`):
  GCC emits a runtime overlap check + two code paths because `mm` carries no `__restrict`. Hits
  the ceiling here (long, distinct buffers) but will cost on short loops and slice operands.
  **Severity: Low.**

## Area 3 — API design — **Medium (macro identifiers leak globally)**

- **Type-name macros (Medium).** `#define zdouble1 JIArray<double,1>` … (`JIArray.h:1963-2005`,
  `zbool/zint/zdouble/zfloat/zstring` × ranks) are **preprocessor macros**: they leak into the
  global namespace, cannot be scoped or forward-declared, and will silently clobber any
  identifier named `zint1` etc. A clean alternative — `zarray<T,N>` template alias — already
  exists (`:2065`). **Severity: Medium.**
- **Loop-macro hygiene (Medium).** `ffor`, `ffor_back`, `zfor`, and especially the helpers
  `GET_STEP` / `GET_STEP_IMPL` (`JIArray.h:2035-2058`) are unqualified global macros with
  collision-prone names (`GET_STEP_IMPL` in particular). No `JIARRAY_` prefix. **Severity: Medium.**
- **1-based + column-major contract is well documented** (README, `docs/manual.md`, `CLAUDE.md`).
  Per-rank `operator()` overloads for debugger watch windows (`:864-875`) are a sensible
  pragmatic touch. **[Info — good]**

## Area 4 — Modernization gap — **Low (enhancement)**

- C++17 throughout; heavy SFINAE (`std::enable_if_t` on nearly every template). No `std::mdspan`,
  `std::span`, concepts, `std::assume_aligned`, ranges, or `<=>`. Clear opportunities:
  `.to_mdspan()`/construct-from-mdspan (ties to slicing + `layout_stride`), concepts to replace
  the `enable_if` walls, `assume_aligned` once an aligned allocator lands. All additive behind a
  standard-detection guard so C++17 users keep working. **Severity: Low.**

## Area 5 — Portability & build — **Medium (no CI; matrix "green" is unverified)**

- **No CI.** No `.github/` (confirmed absent). The GCC/Clang/MSVC/NVCC matrix asserted in the
  handoff is **not actually verified anywhere**; only GCC 11.5 is exercised here. **Severity: Medium.**
- **No sanitizer targets, no `-Wall -Wextra -Werror`** in `CMakeLists.txt` (grep: none). ASan
  had to be wired by hand for this eval. **Severity: Medium.**
- CMake is clean and modern-ish (INTERFACE target, `PROJECT_IS_TOP_LEVEL` test gating, Blitz
  opt-in, `/Zc:preprocessor` propagated for MSVC `__VA_OPT__`). GTest via `find_package` →
  suites silently skipped if absent. **[Info — reasonable]**

## Area 6 — GPU story — **Medium (parallel duplication drifting out of sync)**

- `JICudaArray.h` (1199 L) is a **near-parallel reimplementation** of `JIArray` with the same
  operator surface and the **same eager-temporary pattern** (`initByRankSize` allocations;
  `operator+` `:636/:646`, `*` `:693/:719`, `/` `:756`, scalar ops `:795-841`).
- **Active divergence** already visible: CUDA `operator+=(array)` **returns by value**
  (`JICudaArray.h:627`) whereas the host `operator+=` returns `this_type&`
  (`JIArray.h:1231`) — same-named operators with different semantics across the two types.
  Every fix to shared logic must be applied twice. **Severity: Medium** (maintenance/consistency,
  not a live bug on the host path).

## Area 7 — Release hygiene — **Low-Medium (handoff's version claim is stale)**

- **Correction:** the handoff's "release tag v0.4.0 (Dec 2024) but RELEASE_NOTES_0.7.1.md
  present" is **out of date**. `git describe` = **0.8.0**; tags run `v0.1…v0.5`, then
  `0.6, 0.7.0-0.7.3, 0.8.0` (the `v` prefix was dropped after `v0.5` — itself an inconsistency).
- Remaining real issues: **no `CHANGELOG.md`**; `RELEASE_NOTES_0.7.1.md` is stale (no 0.8.0
  notes); `JIArray.h:21` still says `@version 0.6`; **no single source of version truth** (CMake
  `project()` sets no VERSION). **Severity: Low-Medium.**

---

## Severity roll-up

| Area | Top finding | Severity |
|---|---|---|
| 2 Performance | Eager operator temporaries — 7.9× vs fused (F4) | **High** |
| 1 Correctness | Non-owning view/slice can dangle; no runtime guard | Medium |
| 3 API | Global macro type-names + loop macros leak/clobber | Medium |
| 5 Build | No CI, no sanitizers, no `-Werror`; matrix unverified | Medium |
| 6 GPU | `JICudaArray` duplicates `JIArray` and is diverging | Medium |
| 7 Release | Stale version notes/`@version`; no CHANGELOG/VERSION | Low-Medium |
| 2 Performance | 64 B alignment / AVX-512 width — ~0 gain for streaming | Low |
| 4 Modernization | No mdspan/span/concepts (all additive) | Low |

**Corrections to the handoff to carry into Phase 2:** (a) F4 gap is ~7.9× on this box, not 4.6×,
and only ratios transfer across machines; (b) the P0 "alignment → AVX-512 delta" item should be
expected to yield ~0 for streaming element-wise ops (measured), so it is not a perf lever here;
(c) the version state is 0.8.0, not v0.4.0/0.7.1.
