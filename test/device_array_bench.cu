// Raw-pointer vs JIDeviceView kernels: timing (median of interleaved runs) and
// identical results.  Compare instruction counts with
//   cuobjdump -sass dnegri.jiarray.device_array_bench
// Column-major layout is assumed by the raw kernels (default build).

#include <jiarray/JIDeviceArray.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

using namespace dnegri::jiarray;

namespace {

static_assert(JIARRAY_COLUMN_MAJOR == 1 && JIARRAY_OFFSET == 1, "raw kernels assume the default layout");

constexpr double kCenter = -6.0;    ///< 7-point Laplacian centre weight.
constexpr double kNeighbour = 1.0;  ///< 7-point Laplacian neighbour weight.
constexpr int kFace = 216;          ///< SPX P3 face moments per node (venus GEMV size).
constexpr int kThreads = 224;       ///< 7 warps: one thread per row.

// ---------------------------------------------------------------------------
// 7-point stencil on an nx × ny × nz grid (interior points only)
// ---------------------------------------------------------------------------

__global__ void stencilRawRestrict(const double* __restrict__ in, double* __restrict__ out, int nx, int ny, int nz) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z;
    if (i < 1 || i > nx - 2 || j < 1 || j > ny - 2 || k < 1 || k > nz - 2) {
        return;
    }
    const int c = i + nx * (j + ny * k);
    const int plane = nx * ny;
    out[c] = kCenter * in[c] +
             kNeighbour * (in[c - 1] + in[c + 1] + in[c - nx] + in[c + nx] + in[c - plane] + in[c + plane]);
}

__global__ void stencilRaw(const double* in, double* out, int nx, int ny, int nz) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z;
    if (i < 1 || i > nx - 2 || j < 1 || j > ny - 2 || k < 1 || k > nz - 2) {
        return;
    }
    const int c = i + nx * (j + ny * k);
    const int plane = nx * ny;
    out[c] = kCenter * in[c] +
             kNeighbour * (in[c - 1] + in[c + 1] + in[c - nx] + in[c + nx] + in[c - plane] + in[c + plane]);
}

__global__ void stencilView(ddouble3::ConstView in, ddouble3::View out) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const int k = blockIdx.z + 1;
    if (i < 2 || i > in.getSize(1) - 1 || j < 2 || j > in.getSize(2) - 1 || k < 2 || k > in.getSize(3) - 1) {
        return;
    }
    out(i, j, k) = kCenter * in(i, j, k) + kNeighbour * (in(i - 1, j, k) + in(i + 1, j, k) + in(i, j - 1, k) +
                                                         in(i, j + 1, k) + in(i, j, k - 1) + in(i, j, k + 1));
}

/// @brief stencilView with the read-only data path requested explicitly (`__ldg` on the element address).
__global__ void stencilViewLdg(ddouble3::ConstView in, ddouble3::View out) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const int k = blockIdx.z + 1;
    if (i < 2 || i > in.getSize(1) - 1 || j < 2 || j > in.getSize(2) - 1 || k < 2 || k > in.getSize(3) - 1) {
        return;
    }
    out(i, j, k) = kCenter * __ldg(&in(i, j, k)) +
                   kNeighbour * (__ldg(&in(i - 1, j, k)) + __ldg(&in(i + 1, j, k)) + __ldg(&in(i, j - 1, k)) +
                                 __ldg(&in(i, j + 1, k)) + __ldg(&in(i, j, k - 1)) + __ldg(&in(i, j, k + 1)));
}

// ---------------------------------------------------------------------------
// Node-batched GEMV: y(:, n) = A x(:, n), A is kFace × kFace (column-major)
// ---------------------------------------------------------------------------

__global__ void gemvRawFixed(const double* __restrict__ a, const double* __restrict__ x, double* __restrict__ y) {
    const int i = threadIdx.x;
    const int n = blockIdx.x;
    if (i >= kFace) {
        return;
    }
    double sum = 0.0;
    for (int j = 0; j < kFace; ++j) {
        sum += a[i + kFace * j] * x[j + kFace * n];
    }
    y[i + kFace * n] = sum;
}

__global__ void gemvRaw(const double* __restrict__ a, const double* __restrict__ x, double* __restrict__ y,
                        int num_face) {
    const int i = threadIdx.x;
    const int n = blockIdx.x;
    if (i >= num_face) {
        return;
    }
    double sum = 0.0;
    for (int j = 0; j < num_face; ++j) {
        sum += a[i + num_face * j] * x[j + num_face * n];
    }
    y[i + num_face * n] = sum;
}

__global__ void gemvRawNoRestrict(const double* a, const double* x, double* y, int num_face) {
    const int i = threadIdx.x;
    const int n = blockIdx.x;
    if (i >= num_face) {
        return;
    }
    double sum = 0.0;
    for (int j = 0; j < num_face; ++j) {
        sum += a[i + num_face * j] * x[j + num_face * n];
    }
    y[i + num_face * n] = sum;
}

__global__ void gemvView(ddouble2::ConstView a, ddouble2::ConstView x, ddouble2::View y) {
    const int i = threadIdx.x + 1;
    const int n = blockIdx.x + 1;
    if (i > a.getSize(1)) {
        return;
    }
    double sum = 0.0;
    zfor(j, a.getSize(2)) {
        sum += a(i, j) * x(j, n);
    }
    y(i, n) = sum;
}

void check(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        std::fprintf(stderr, "%s: %s\n", what, cudaGetErrorString(status));
        std::exit(1);
    }
}

/// @brief Median kernel time [ms] of each variant over interleaved rounds.
std::vector<double> timeInterleaved(const std::vector<std::function<void()>>& variants, int num_rounds) {
    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    check(cudaEventCreate(&start), "event");
    check(cudaEventCreate(&stop), "event");
    std::vector<std::vector<float>> samples(variants.size());
    for (const auto& launch : variants) {
        launch();  // warm-up
    }
    for (int round = 0; round < num_rounds; ++round) {
        for (size_t v = 0; v < variants.size(); ++v) {
            check(cudaEventRecord(start), "record");
            variants[v]();
            check(cudaEventRecord(stop), "record");
            check(cudaEventSynchronize(stop), "sync");
            float elapsed = 0.0f;
            check(cudaEventElapsedTime(&elapsed, start, stop), "elapsed");
            samples[v].push_back(elapsed);
        }
    }
    std::vector<double> medians;
    for (auto& s : samples) {
        std::sort(s.begin(), s.end());
        medians.push_back(s[s.size() / 2]);
    }
    check(cudaEventDestroy(start), "event");
    check(cudaEventDestroy(stop), "event");
    return medians;
}

double computeMaxDifference(const zdouble1& a, const zdouble1& b) {
    double worst = 0.0;
    zfor(n, a.size()) {
        worst = std::max(worst, std::abs(a(n) - b(n)));
    }
    return worst;
}

void runStencil(int edge, int num_rounds) {
    zdouble3 host_in(edge, edge, edge);
    for (int n = 0; n < host_in.size(); ++n) {
        host_in[n] = 1.0 / (1.0 + n % 97);
    }
    ddouble3 in;
    in = host_in;
    ddouble3 out_raw(edge, edge, edge);
    ddouble3 out_view(edge, edge, edge);

    const dim3 block(32, 8);
    const dim3 grid((edge + 31) / 32, (edge + 7) / 8, edge);
    const std::vector<std::function<void()>> variants = {
        [&] { stencilRawRestrict<<<grid, block>>>(in.data(), out_raw.data(), edge, edge, edge); },
        [&] { stencilRaw<<<grid, block>>>(in.data(), out_raw.data(), edge, edge, edge); },
        [&] { stencilView<<<grid, block>>>(in.view(), out_view.view()); },
        [&] { stencilViewLdg<<<grid, block>>>(in.view(), out_view.view()); },
    };
    const std::vector<double> ms = timeInterleaved(variants, num_rounds);
    check(cudaDeviceSynchronize(), "stencil");

    zdouble1 raw_flat(out_raw.size());
    zdouble1 view_flat(out_view.size());
    check(cudaMemcpy(raw_flat.data(), out_raw.data(), raw_flat.size() * sizeof(double), cudaMemcpyDeviceToHost), "d2h");
    check(cudaMemcpy(view_flat.data(), out_view.data(), view_flat.size() * sizeof(double), cudaMemcpyDeviceToHost),
          "d2h");
    std::printf("stencil %d^3 FP64: raw+restrict %.4f ms | raw %.4f ms | view %.4f ms | view+__ldg %.4f ms | "
                "max|diff| %.1e\n",
                edge, ms[0], ms[1], ms[2], ms[3], computeMaxDifference(raw_flat, view_flat));
}

void runGemv(int num_nodes, int num_rounds) {
    zdouble2 host_a(kFace, kFace);
    zdouble2 host_x(kFace, num_nodes);
    for (int n = 0; n < host_a.size(); ++n) {
        host_a[n] = 1.0 / (1.0 + n % 31);
    }
    for (int n = 0; n < host_x.size(); ++n) {
        host_x[n] = 1.0 / (2.0 + n % 17);
    }
    ddouble2 a;
    a = host_a;
    ddouble2 x;
    x = host_x;
    ddouble2 y_raw(kFace, num_nodes);
    ddouble2 y_view(kFace, num_nodes);

    const std::vector<std::function<void()>> variants = {
        [&] { gemvRawFixed<<<num_nodes, kThreads>>>(a.data(), x.data(), y_raw.data()); },
        [&] { gemvRaw<<<num_nodes, kThreads>>>(a.data(), x.data(), y_raw.data(), kFace); },
        [&] { gemvRawNoRestrict<<<num_nodes, kThreads>>>(a.data(), x.data(), y_raw.data(), kFace); },
        [&] { gemvView<<<num_nodes, kThreads>>>(a.view(), x.view(), y_view.view()); },
    };
    const std::vector<double> ms = timeInterleaved(variants, num_rounds);
    check(cudaDeviceSynchronize(), "gemv");

    zdouble1 raw_flat(y_raw.size());
    zdouble1 view_flat(y_view.size());
    check(cudaMemcpy(raw_flat.data(), y_raw.data(), raw_flat.size() * sizeof(double), cudaMemcpyDeviceToHost), "d2h");
    check(cudaMemcpy(view_flat.data(), y_view.data(), view_flat.size() * sizeof(double), cudaMemcpyDeviceToHost),
          "d2h");
    std::printf("gemv %dx%d x %d nodes FP64: raw fixed-216 %.4f ms | raw %.4f ms | raw no-restrict %.4f ms | "
                "view %.4f ms | max|diff| %.1e\n",
                kFace, kFace, num_nodes, ms[0], ms[1], ms[2], ms[3], computeMaxDifference(raw_flat, view_flat));
}

}  // namespace

int main(int argc, char** argv) {
    const int num_rounds = argc > 1 ? std::stoi(argv[1]) : 31;
    runStencil(256, num_rounds);
    runGemv(20000, num_rounds);
    return 0;
}
