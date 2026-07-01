# P0 — Arithmetic fusion: prototype results, decision & outcome

Branch: `p0-arithmetic-fusion`.
**Decision (owner): Option 1 — full expression templates** (type named `JIArrayExpr`), accepting
the "receive operator results into a concrete type, not `auto`" contract.
**Status: implemented and verified — see [Outcome](#outcome-implemented) below.**

## Baseline (this machine, `taskset -c 0`, best-of-7, n=16M doubles)

| | time (ms) | GB/s | vs `ffor` |
|---|---|---|---|
| eager operator chain `d=a+b+c` (status quo) | 1634 | 6.3 | **7.9×** |
| `ffor` fused (target) | 208 | 49.3 | 1.0× |
| in-place `d=a; d+=b; d+=c` | 372 | 27.5 | 1.8× |
| raw ptr ceiling | 207 | 49.4 | 1.0× |

`ctest`: 413/413 green. ASan (lifecycle suite): clean. (UBSan runtime unavailable in this env.)

## Prototype: both fusion routes reach the ceiling

Standalone prototype (`scratchpad/proto.cpp`), same harness:

| route | time (ms) | GB/s | vs `ffor` |
|---|---|---|---|
| Route A — expression templates (lazy, fused on assign) | 211.7 | 48.4 | **1.00×** |
| Route B — named fused op (`add_into(d,a,b,c)`) | 211.5 | 48.5 | **1.00×** |
| `ffor` fused reference | 212.5 | 48.2 | 1.00× |

**Performance is a tie: both crush the ≤1.2× acceptance bar (both at 1.0×).**
The decision is therefore purely **API preservation + complexity + CUDA cost**, not speed.

## The real constraint: what fusing the literal `+` costs

To make the *literal* `d = a + b + c` fast, `operator+` must become **lazy** (return an
expression node, not a `JIArray`). Measured blast radius on the current suite:

- **1 hard breakage** — `auto result = array1 + array2; … result = array1 - array2;`
  (`column_major_test.cpp:50` then `:56`): reassigning an `auto`-deduced expression variable with
  a *different* expression type cannot compile (its type is fixed at declaration). No ET scheme
  fixes this without type erasure (which reintroduces the overhead).
- **~16 sites** index an operator result as an array — `auto result = a+b; result(1,1)`. These
  compile only if expression nodes **reimplement offset-aware `operator()` / `size()`** for *both*
  column-major and row-major, and that machinery must be **mirrored into `JICudaArray.h`**
  (`JIARRAY_HD`-clean) — a large, error-prone surface across two files and two layouts.

## Three implementable options (all reach ceiling)

| | Fixes bare `d=a+b+c`? | API break | Complexity / CUDA mirror | Ergonomics |
|---|---|---|---|---|
| **1. Full lazy ETs** (`operator+` lazy) | ✅ yes | 1 test pattern + must reimpl expr indexing | High (2 files, 2 layouts) | best (`d=a+b+c`) |
| **2. Opt-in ETs via `fuse()`** (`d = fuse(a)+b+c`) | ⚠️ opt-in only | **none** | Low (small, additive) | `d = fuse(a)+b+c` |
| **3. Named fused op** (`add_into(d,a,b,c)`) | ❌ no | **none** | Lowest | least natural |

Option 2 keeps eager `Arr+Arr` untouched (API 100% preserved, all 413 tests green); the lazy
overloads only engage once an operand is wrapped by `fuse(...)`, so `d = fuse(a)+b+c` fuses to a
single allocation-free pass (measured 1.0×) while `auto r = a+b` stays an ordinary materialized
array. It reads almost like the natural syntax and its CUDA mirror is a handful of overloads.

## Recommendation

**Option 2 (opt-in expression templates seeded by `fuse()`).** It hits the memory ceiling,
preserves the entire public API and all tests, stays cheap to mirror for CUDA, and reads close
to natural syntax. Full lazy ETs (Option 1) are the textbook fix but impose an API regression
plus a large two-file/two-layout complexity burden for a marginal ergonomic gain over `fuse()`;
that cost is not justified given `fuse()` already reaches the ceiling.

Regardless of route, keep the eager operators (back-compat) and **document `ffor` / `+=` / the
new fused path as the hot-path guidance** — that closes the F4 footgun for the common case.

## Outcome (implemented)

Owner chose **Option 1**. Implementation:

- **`include/jiarray/JIArrayExpr.h`** (new): lazy nodes `JIArrayLeaf` / `JIArrayScalar` /
  `JIArrayBinaryExpr` / `JIArrayUnaryExpr` under a CRTP base, `to_expr()`, and the lazy
  `+ - * /` + unary `-` operators (SFINAE-gated so only JIArray/expr operands engage — FastArray
  etc. are untouched). All nodes are `JIARRAY_HD`; evaluation is linear over flat storage, so it
  is layout-agnostic.
- **`JIArray.h`**: removed the eager array/array + scalar operators and unary `-`; added a
  materializing constructor and `operator=` from an expression (single fused pass). These two
  allocate, so they are **host-only** (not `JIARRAY_HD`) — verified no NVCC `20011-D` warning.
  In-place `+= -= *= /=` unchanged.
- **Tests**: 19 `auto result = <op>` sites converted to concrete types (the new contract).

### Measured delta (`bench.cpp`, `taskset -c 0`, best-of-7, quiet box)

| variant | before | after | note |
|---|---|---|---|
| (A) operator chain `d=a+b+c` | 1634 ms / 6.3 GB/s | **205 ms / 49.9 GB/s** | **7.9× → 1.00× of ceiling** |
| (B) `ffor` fused | 208 ms | 220 ms | reference |
| (C) in-place `d=a;d+=b;d+=c` | 372 ms | 371 ms | fast path preserved |
| (D) raw ptr ceiling | 207 ms | 205 ms | (A) now == ceiling |

Acceptance bar (operator chain ≤ 1.2× of `ffor`): **met — 0.93× (faster than `ffor`).**

### Verification

- `ctest`: **413/413** green (GCC, `JIARRAY_DEBUG`, column + row major).
- Clang 20: coverage/colmajor/rowmajor/move/fastarray suites all green.
- NVCC: `JIArray.h` (host ET) and `JICudaArray.h` compile, no warnings.
- ASan: `move_lifecycle` (40) + ET correctness snippet — clean.
- Public API preserved except the documented `auto`→concrete-type contract (manual.md updated).

### Residual (documented, not fixed here)

`auto r = a + b;` still compiles (binds a node) and can dangle — the contract is convention, not
compiler-enforced. A future opt-in debug guard could catch it. Not in P0 scope.
