// bench.cpp — jiarray element-wise expression: operator chain vs fused
#include <jiarray/JIArray.h>
#include <chrono>
#include <cstdio>
#include <vector>
#include <algorithm>
using namespace dnegri::jiarray;
using clk = std::chrono::steady_clock;

static double sink = 0.0;  // defeat dead-code elimination

template <class F>
double best_ms(F&& f, int reps) {
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        auto t0 = clk::now();
        f();
        auto t1 = clk::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        best = std::min(best, ms);
    }
    return best;
}

int main() {
    const int    n     = 16'000'000;          // 128 MB per double array -> memory bound
    const int    iters = 20;                    // passes fused inside each timed call
    const int    reps  = 7;                     // take best of reps
    const double GB    = 1e9;
    const double bytes_per_elem_move = 8.0;

    zdouble1 a(n), b(n), c(n), d(n);
    ffor(i, 1, n) { a(i) = 1.0 + (i & 7); b(i) = 2.0; c(i) = 0.5; }

    // (A) operator chain: d = a + b + c   (allocates 2 temporaries per expression)
    double msA = best_ms([&]{
        for (int k = 0; k < iters; ++k) { d = a + b + c; sink += d(1); }
    }, reps);

    // (B) fused ffor loop: single pass, no temporary
    double msB = best_ms([&]{
        for (int k = 0; k < iters; ++k) { ffor(i, 1, n) d(i) = a(i) + b(i) + c(i); sink += d(1); }
    }, reps);

    // (C) in-place accumulation: d = a; d += b; d += c;
    double msC = best_ms([&]{
        for (int k = 0; k < iters; ++k) { d = a; d += b; d += c; sink += d(1); }
    }, reps);

    // (D) raw pointer fused baseline (hardware ceiling)
    std::vector<double> ra(n), rb(n), rc(n), rd(n);
    for (int i = 0; i < n; ++i) { ra[i]=1.0+(i&7); rb[i]=2.0; rc[i]=0.5; }
    double msD = best_ms([&]{
        for (int k = 0; k < iters; ++k) {
            const double* __restrict pa=ra.data(); const double* __restrict pb=rb.data();
            const double* __restrict pc=rc.data(); double* __restrict pd=rd.data();
            #pragma omp simd
            for (int i = 0; i < n; ++i) pd[i] = pa[i] + pb[i] + pc[i];
            sink += rd[0];
        }
    }, reps);

    // effective bandwidth: fused expression touches 4 arrays (read a,b,c + write d) per pass
    auto bw = [&](double ms){ return (4.0 * n * bytes_per_elem_move * iters) / (ms/1e3) / GB; };

    printf("n=%d  iters=%d  best-of-%d   (sink=%.3g)\n\n", n, iters, reps, sink);
    printf("%-38s %10s %12s %8s\n", "variant", "time(ms)", "eff GB/s", "vs fused");
    printf("%-38s %10.1f %12.1f %8s\n", "(A) operator chain  d=a+b+c", msA, bw(msA), "");
    printf("%-38s %10.1f %12.1f %7.2fx\n", "(B) ffor fused      d(i)=..", msB, bw(msB), msA/msB);
    printf("%-38s %10.1f %12.1f %7.2fx\n", "(C) in-place        d=a;d+=b;d+=c", msC, bw(msC), msA/msC);
    printf("%-38s %10.1f %12.1f %7.2fx\n", "(D) raw ptr fused (ceiling)", msD, bw(msD), msA/msD);
    return 0;
}

