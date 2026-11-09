// JIDeviceArray / JIDeviceView tests.  Built twice: column-major (default) and
// row-major (-DJIARRAY_COLUMN_MAJOR=0); every layout check compares against the
// host JIArray compiled with the same macros, so both orders are covered by the
// same assertions.

#include <jiarray/JIDeviceArray.h>

#include <gtest/gtest.h>

#include <cstdlib>
#include <type_traits>
#include <utility>

using namespace dnegri::jiarray;

namespace {

static_assert(std::is_trivially_copyable_v<ddouble3::View>, "view must be a by-value kernel argument");
static_assert(std::is_trivially_copyable_v<ddouble3::ConstView>, "const view must be a by-value kernel argument");
static_assert(!std::is_copy_constructible_v<ddouble3>, "owner must be move-only");
static_assert(std::is_nothrow_move_constructible_v<ddouble3>, "owner must be movable");

/// @brief Value unique to an index triple (offsets may be negative).
__host__ __device__ double computeCode(int i, int j, int k) {
    return 10000.0 * i + 100.0 * j + k;
}

void synchronize() {
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
}

/// @brief One thread writes computeCode over the whole index box with JIArray loop macros.
__global__ void writeCodeSerial(ddouble3::View a) {
    ffor(k, a.getOffset(3), a.getOffset(3) + a.getSize(3) - 1) {
        ffor(j, a.getOffset(2), a.getOffset(2) + a.getSize(2) - 1) {
            ffor(i, a.getOffset(1), a.getOffset(1) + a.getSize(1) - 1) {
                a(i, j, k) = computeCode(i, j, k);
            }
        }
    }
}

/// @brief One thread per element on a 3-D grid; thread (x, y, z) owns index (offset + x, ...).
__global__ void writeCodeParallel(ddouble3::View a) {
    const int i = a.getOffset(1) + static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int j = a.getOffset(2) + static_cast<int>(blockIdx.y);
    const int k = a.getOffset(3) + static_cast<int>(blockIdx.z);
    if (i > a.getOffset(1) + a.getSize(1) - 1) {
        return;
    }
    a(i, j, k) = computeCode(i, j, k);
}

/// @brief One thread fills a default-offset array with zfor, the 1-based kernel idiom.
__global__ void writeCodeZfor(ddouble3::View a) {
    zfor(k, a.getSize(3)) {
        zfor(j, a.getSize(2)) {
            zfor(i, a.getSize(1)) {
                a(i, j, k) = computeCode(i, j, k);
            }
        }
    }
}

/// @brief Index through a FastArray, as JIArray::operator()(const FastArray<int, RANK>&).
__global__ void writeCodeFastArray(ddouble3::View a) {
    FastArray<int, 3> idx;
    zfor(k, a.getSize(3)) {
        zfor(j, a.getSize(2)) {
            zfor(i, a.getSize(1)) {
                idx[0] = i;
                idx[1] = j;
                idx[2] = k;
                a(idx) = computeCode(i, j, k);
            }
        }
    }
}

/// @brief Negate every element of a 2-D view over its own index box.
__global__ void negate2(ddouble2::View s) {
    ffor(j, s.getOffset(2), s.getOffset(2) + s.getSize(2) - 1) {
        ffor(i, s.getOffset(1), s.getOffset(1) + s.getSize(1) - 1) {
            s(i, j) = -s(i, j);
        }
    }
}

/// @brief Slice inside the kernel and add 0.5 to the 1-D line.
__global__ void addHalfToLine(ddouble3::View a, int first, int second) {
    ddouble1::View line = a.slice(first, second);
    ffor(i, line.getOffset(1), line.getOffset(1) + line.getSize(1) - 1) {
        line(i) += 0.5;
    }
}

/// @brief Reduce a read-only view into out[0].
__global__ void sumConst(ddouble3::ConstView a, double* out) {
    double total = 0.0;
    for (int n = 0; n < a.size(); ++n) {
        total += a[n];
    }
    out[0] = total;
}

/// @brief Write 10*i + j through a reshaped 2-D view built on the device.
__global__ void writeReshaped(ddouble1::View flat, int rows, int cols) {
    ddouble2::View m = flat.reshape(rows, cols);
    zfor(j, cols) {
        zfor(i, rows) {
            m(i, j) = 10.0 * i + j;
        }
    }
}

__global__ void writeCode64(JIDeviceView<float, 2, long long> a) {
    ffor(j, a.getOffset(2), a.getOffset(2) + a.getSize(2) - 1) {
        ffor(i, a.getOffset(1), a.getOffset(1) + a.getSize(1) - 1) {
            a(i, j) = static_cast<float>(100 * i + j);
        }
    }
}

__global__ void readOutOfBounds(ddouble1::View a, double* out) {
    out[0] = a(a.getSize(1) + 1);
}

/// @brief Launch an out-of-bounds read and exit with 3 if the kernel failed (death-test child body).
[[noreturn]] void runOutOfBoundsAndExit() {
    ddouble1 a(8);
    ddouble1 out(1);
    readOutOfBounds<<<1, 1>>>(a.view(), out.data());
    std::exit(cudaDeviceSynchronize() == cudaSuccess ? 0 : 3);
}

/// @brief Host JIArray filled with computeCode over its index box.
void fillCode(zdouble3& h) {
    ffor(k, h.getOffset(3), h.getOffset(3) + h.getSize(3) - 1) {
        ffor(j, h.getOffset(2), h.getOffset(2) + h.getSize(2) - 1) {
            ffor(i, h.getOffset(1), h.getOffset(1) + h.getSize(1) - 1) {
                h(i, j, k) = computeCode(i, j, k);
            }
        }
    }
}

/// @brief Every element of `h` equals computeCode at its JIArray index.
void expectCode(const zdouble3& h) {
    ffor(k, h.getOffset(3), h.getOffset(3) + h.getSize(3) - 1) {
        ffor(j, h.getOffset(2), h.getOffset(2) + h.getSize(2) - 1) {
            ffor(i, h.getOffset(1), h.getOffset(1) + h.getSize(1) - 1) {
                ASSERT_EQ(h(i, j, k), computeCode(i, j, k)) << i << "," << j << "," << k;
            }
        }
    }
}

template <class View, class Host>
void expectSameLayout(const View& v, const Host& h) {
    ASSERT_EQ(v.size(), h.size());
    for (size_t d = 0; d < View::rank; ++d) {
        EXPECT_EQ(v.getSizeOfRank()[d], h.getSizeOfRank()[d]) << "dim " << d;
        EXPECT_EQ(v.getOffset()[d], h.getOffset()[d]) << "dim " << d;
        EXPECT_EQ(v.getRankSize()[d], h.getRankSize()[d]) << "dim " << d;
    }
}

void expectEqualArrays(const zdouble3& a, const zdouble3& b) {
    ASSERT_EQ(a.size(), b.size());
    for (int n = 0; n < a.size(); ++n) {
        ASSERT_EQ(a[n], b[n]) << "flat " << n;
    }
}

}  // namespace

TEST(DeviceArray, UploadDownloadRoundTrip) {
    zdouble3 h(3, 4, 5);
    fillCode(h);
    ddouble3 d;
    d = h;
    expectSameLayout(d, h);

    zdouble3 fresh;
    d.copyTo(fresh);
    expectSameLayout(fresh, h);
    expectEqualArrays(fresh, h);

    zdouble3 preallocated(3, 4, 5);
    d.copyTo(preallocated);
    expectEqualArrays(preallocated, h);
}

TEST(DeviceArray, InitZeroFills) {
    {
        ddouble2 dirty(64, 64);
        dirty.fillAsync(7.0, nullptr);
        synchronize();
    }
    ddouble2 d(64, 64);
    zdouble2 h;
    d.copyTo(h);
    for (int n = 0; n < h.size(); ++n) {
        ASSERT_EQ(h[n], 0.0);
    }
}

TEST(DeviceArray, KernelIndexingMatchesHost) {
    ddouble3 d(5, 6, 7);
    writeCodeZfor<<<1, 1>>>(d.view());
    synchronize();
    zdouble3 h;
    d.copyTo(h);
    expectCode(h);

    ddouble3 p(37, 6, 7);
    writeCodeParallel<<<dim3(2, 6, 7), 32>>>(p.view());
    synchronize();
    zdouble3 hp;
    p.copyTo(hp);
    expectCode(hp);
}

TEST(DeviceArray, FastArrayIndex) {
    ddouble3 d(4, 3, 2);
    writeCodeFastArray<<<1, 1>>>(d.view());
    synchronize();
    zdouble3 h;
    d.copyTo(h);
    expectCode(h);
}

TEST(DeviceArray, Init0AndSetOffsets) {
    ddouble3 d;
    d.init0(0, 2, -1, 3, 2, 4);
    EXPECT_EQ(d.getOffset(1), 0);
    EXPECT_EQ(d.getOffset(2), -1);
    EXPECT_EQ(d.getOffset(3), 2);
    EXPECT_EQ(d.getSize(1), 3);
    EXPECT_EQ(d.getSize(2), 5);
    EXPECT_EQ(d.getSize(3), 3);
    writeCodeSerial<<<1, 1>>>(d.view());
    synchronize();
    zdouble3 h;
    d.copyTo(h);
    EXPECT_EQ(h.getOffset(2), -1);
    expectCode(h);

    d.setOffsets(5, 6, 7);
    EXPECT_EQ(d.getOffset(1), 5);
    writeCodeSerial<<<1, 1>>>(d.view());
    synchronize();
    zdouble3 h2;
    d.copyTo(h2);
    expectCode(h2);

    zdouble3 host_offset;
    host_offset.init0(-2, 0, 1, 2, 3, 6);
    fillCode(host_offset);
    ddouble3 adopted;
    adopted = host_offset;
    expectSameLayout(adopted, host_offset);
    writeCodeParallel<<<dim3(1, 2, 4), 32>>>(adopted.view());
    synchronize();
    zdouble3 back;
    adopted.copyTo(back);
    expectCode(back);
}

TEST(DeviceArray, SliceMatchesJIArray) {
    zdouble3 h;
    h.init0(0, 3, -1, 2, 1, 5);
    fillCode(h);
    ddouble3 d;
    d = h;

    // Host-built slice: same layout and base address as JIArray::slice.
    const int fix = 2;
    ddouble2::View s = d.slice(fix);
    zdouble2 hs = h.slice(fix);
    expectSameLayout(s, hs);
    EXPECT_EQ(s.data() - d.data(), hs.data() - h.data());

    negate2<<<1, 1>>>(s);
    synchronize();
    ffor(j, hs.getOffset(2), hs.getOffset(2) + hs.getSize(2) - 1) {
        ffor(i, hs.getOffset(1), hs.getOffset(1) + hs.getSize(1) - 1) {
            hs(i, j) = -hs(i, j);
        }
    }

    // Two-index slice built inside the kernel.
    const int first = JIARRAY_COLUMN_MAJOR ? 1 : 2;
    const int second = JIARRAY_COLUMN_MAJOR ? 3 : 0;
    ddouble1::View line = d.slice(first, second);
    zdouble1 hline = h.slice(first, second);
    expectSameLayout(line, hline);
    EXPECT_EQ(line.data() - d.data(), hline.data() - h.data());
    addHalfToLine<<<1, 1>>>(d.view(), first, second);
    synchronize();
    ffor(i, hline.getOffset(1), hline.getOffset(1) + hline.getSize(1) - 1) {
        hline(i) += 0.5;
    }

    zdouble3 back;
    d.copyTo(back);
    expectEqualArrays(back, h);
}

TEST(DeviceArray, ReshapeMatchesJIArray) {
    zdouble1 h(12);
    ddouble1 d(12);
    ddouble2::View m = d.reshape(3, 4);
    zdouble2 hm = h.reshape(3, 4);
    expectSameLayout(m, hm);

    writeReshaped<<<1, 1>>>(d.view(), 3, 4);
    synchronize();
    zfor(j, 4) {
        zfor(i, 3) {
            hm(i, j) = 10.0 * i + j;
        }
    }
    zdouble1 back;
    d.copyTo(back);
    zfor(n, 12) {
        ASSERT_EQ(back(n), h(n)) << n;
    }
}

TEST(DeviceArray, DataPartialIndexMatchesJIArray) {
    zdouble3 h(3, 4, 5);
    ddouble3 d;
    d = h;
    EXPECT_EQ(d.data(2) - d.data(), h.data(2) - h.data());
    EXPECT_EQ(d.data(2, 3) - d.data(), h.data(2, 3) - h.data());
    EXPECT_EQ(d.data(2, 3, 4) - d.data(), &h(2, 3, 4) - h.data());
    EXPECT_EQ(&d.view()(3, 2, 1) - d.data(), &h(3, 2, 1) - h.data());
}

TEST(DeviceArray, ConstViewReads) {
    zdouble3 h(4, 4, 4);
    fillCode(h);
    ddouble3 d;
    d = h;
    const ddouble3& cd = d;
    ddouble1 total(1);
    sumConst<<<1, 1>>>(cd.view(), total.data());
    synchronize();
    zdouble1 htotal;
    total.copyTo(htotal);
    EXPECT_DOUBLE_EQ(htotal(1), h.sum());

    ddouble3::ConstView converted = d.view();
    EXPECT_EQ(converted.data(), d.data());
}

TEST(DeviceArray, MoveTransfersOwnershipWithoutLeak) {
    size_t free_before = 0;
    size_t total = 0;
    ASSERT_EQ(cudaMemGetInfo(&free_before, &total), cudaSuccess);
    constexpr int kCount = 32 << 20;  // 256 MiB of doubles
    for (int repeat = 0; repeat < 8; ++repeat) {
        ddouble1 a(kCount);
        const double* raw = a.data();
        ddouble1 b(std::move(a));
        EXPECT_FALSE(a.isAllocated());
        EXPECT_EQ(a.data(), nullptr);
        EXPECT_EQ(a.size(), 0);
        EXPECT_EQ(b.data(), raw);

        ddouble1 c(kCount);
        c = std::move(b);  // releases c's own buffer
        EXPECT_FALSE(b.isAllocated());
        EXPECT_EQ(c.data(), raw);
        c.destroy();
        EXPECT_FALSE(c.isAllocated());
    }
    size_t free_after = 0;
    ASSERT_EQ(cudaMemGetInfo(&free_after, &total), cudaSuccess);
    EXPECT_GE(free_after + (64u << 20), free_before);
}

TEST(DeviceArray, StreamOrderedAllocation) {
    cudaStream_t stream = nullptr;
    ASSERT_EQ(cudaStreamCreate(&stream), cudaSuccess);
    {
        ddouble2 d;
        d.initAsync(stream, 100, 50);
        EXPECT_EQ(d.getAllocation(), DeviceAllocation::StreamOrdered);
        zdouble2 h(100, 50);
        d.copyToAsync(h, stream);
        ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
        for (int n = 0; n < h.size(); ++n) {
            ASSERT_EQ(h[n], 0.0);
        }
        ddouble2 moved = std::move(d);
        EXPECT_EQ(moved.getAllocation(), DeviceAllocation::StreamOrdered);
    }
    ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
    ASSERT_EQ(cudaStreamDestroy(stream), cudaSuccess);
}

TEST(DeviceArray, UploadAsConvertsPrecision) {
    zdouble2 h(5, 7);
    zfor(n, h.size()) {
        h[n - 1] = 1.0 / (3.0 + n);
    }
    dfloat2 f;
    f.uploadAs(h);
    expectSameLayout(f, h);

    zfloat2 hf;
    f.copyTo(hf);
    zdouble2 hd;
    f.copyTo(hd);
    for (int n = 0; n < h.size(); ++n) {
        ASSERT_EQ(hf[n], static_cast<float>(h[n]));
        ASSERT_EQ(hd[n], static_cast<double>(static_cast<float>(h[n])));
    }
}

TEST(DeviceArray, AsyncCopyAndFill) {
    cudaStream_t stream = nullptr;
    ASSERT_EQ(cudaStreamCreate(&stream), cudaSuccess);
    constexpr int kCount = 100000;
    ddouble1 d(kCount);
    d.fillAsync(2.5, stream);
    zdouble1 h(kCount);
    d.copyToAsync(h, stream);
    ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
    zfor(n, kCount) {
        ASSERT_EQ(h(n), 2.5);
    }

    zdouble1 source(kCount);
    zfor(n, kCount) {
        source(n) = n;
    }
    d.copyFromAsync(source, stream);
    zdouble1 back(kCount);
    d.copyToAsync(back, stream);
    ASSERT_EQ(cudaStreamSynchronize(stream), cudaSuccess);
    zfor(n, kCount) {
        ASSERT_EQ(back(n), n);
    }
    ASSERT_EQ(cudaStreamDestroy(stream), cudaSuccess);
}

TEST(DeviceArray, Index64) {
    JIDeviceArray<float, 2, long long> d;
    d.init0(-3, 4, 2, 6);
    writeCode64<<<1, 1>>>(d.view());
    synchronize();
    zfloat2 h;
    d.copyTo(h);
    EXPECT_EQ(h.getOffset(1), -3);
    ffor(j, 2, 6) {
        ffor(i, -3, 4) {
            ASSERT_EQ(h(i, j), static_cast<float>(100 * i + j));
        }
    }
}

TEST(DeviceArray, ShapeMismatchThrows) {
    ddouble2 d(3, 4);
    zdouble2 wrong_size(4, 3);
    EXPECT_THROW(d = wrong_size, std::invalid_argument);
    zdouble2 wrong_offset(3, 4);
    wrong_offset.setOffsets(0, 0);
    EXPECT_THROW(d = wrong_offset, std::invalid_argument);
    EXPECT_THROW(d.copyToAsync(wrong_size, nullptr), std::invalid_argument);

    ddouble2 empty;
    zdouble2 h(3, 4);
    EXPECT_THROW(empty.copyFromAsync(h, nullptr), std::invalid_argument);

    dint1 huge;
    EXPECT_THROW(huge.init(3000000000LL), std::length_error);
}

#ifdef JIARRAY_DEBUG
TEST(DeviceArray, HostBoundsChecks) {
    ddouble3 d(3, 4, 5);
    EXPECT_THROW(d.getSize(4), std::out_of_range);
    EXPECT_THROW(d.getOffset(0), std::out_of_range);
    const int too_far = JIARRAY_COLUMN_MAJOR ? 6 : 4;
    EXPECT_THROW(d.slice(too_far), std::out_of_range);
    EXPECT_THROW(d.view()(4, 1, 1), std::out_of_range);
    EXPECT_THROW(d.view()(1, 1, 0), std::out_of_range);
    ddouble1 flat(12);
    EXPECT_THROW(flat.reshape(5, 5), std::invalid_argument);
}

TEST(DeviceArrayDeathTest, DeviceBoundsTrap) {
    ::testing::GTEST_FLAG(death_test_style) = "threadsafe";
    EXPECT_EXIT(runOutOfBoundsAndExit(), ::testing::ExitedWithCode(3), "");
}
#endif
