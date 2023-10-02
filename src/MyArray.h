#pragma once

#include <assert.h>

#include <algorithm>
#include <cstddef>  // For std::ptrdiff_t
#include <iterator>
#include <iterator>  // For std::forward_iterator_tag
#include <tuple>
#include <utility>

#ifndef MYARRAY_OFFSET
#define MYARRAY_OFFSET 1
#endif

#ifdef MYARRAY_DEBUG
#define MYARRAY_CHECK_NOT_ALLOCATED() assert(allocated == false)
#define MYARRAY_CHECK_BOUND(i, beg, end) assert(i >= beg && i <= end)
#define MYARRAY_CHECK_RANK(r1, r2) assert(r1 == r2)
#define MYARRAY_CHECK_SIZE(r1, r2) assert(r1 == r2)
#else
#define MYARRAY_CHECK_NOT_ALLOCATED()
#define MYARRAY_CHECK_BOUND(i, beg, end)
#define MYARRAY_CHECK_RANK(r1, r2)
#define MYARRAY_CHECK_SIZE(r1, r2)
#endif

#define MYARRAY_ALLOCATED_NONE 0
#define MYARRAY_ALLOCATED_MEMORY 1
#define MYARRAY_ALLOCATED_RANKSIZE 2
#define MYARRAY_ALLOCATED_OFFSET 4
#define MYARRAY_ALLOCATED_ALL 7
#define MYARRAY_ALLOCATED_RANKSIZE_OFFSET 6

template <class T, std::size_t SIZE, std::size_t OFFSET = MYARRAY_OFFSET>
class FastArray {
   public:
    T memory[SIZE]{};

   public:
    inline T &operator()(int i) {
        MYARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return memory[i - OFFSET];
    }
    inline const T &operator()(int i) const {
        MYARRAY_CHECK_BOUND(i, OFFSET, OFFSET + SIZE - 1);
        return memory[i - OFFSET];
    }
};

#ifndef SOURCE1_H_
#define SOURCE1_H_
extern int64_t sizeOfMyArray;
#endif

template <class T, std::size_t RANK, class = std::make_index_sequence<RANK>>
class MyArray;

template <class T, std::size_t RANK, std::size_t... INTS>
class MyArray<T, RANK, std::index_sequence<INTS...>> {
   public:
    T *memory = nullptr;
    int allocated = MYARRAY_ALLOCATED_NONE;
    int size = 0;
    int rankSize[RANK]{};
    int offset[RANK]{};
    int sumOfOffset = 0;

    static const int RANK_1 = RANK - 1;
    static const int RANK_2 = RANK - 2;
    static const int RANK_3 = RANK - 3;
    static const int RANK_4 = RANK - 4;
    static const int RANK_5 = RANK - 5;

   private:
#ifdef MYARRAY_DEBUG
    int sizeOfRank[RANK]{};
#endif
   private:
    void init(int i) {
        size = i;
        sizeOfMyArray += sizeof(T) * size;
        memory = new T[size]{};
        for (int ioffset = 0; ioffset < 1; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        allocated = MYARRAY_ALLOCATED_ALL;
        sumOfOffset = offset[0];
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
#endif
    };

    void init(int i, int j) {
        size = i * j;
        sizeOfMyArray += sizeof(T) * size;
        memory = new T[size]{};
        for (int ioffset = 0; ioffset < 2; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        sumOfOffset = offset[0] + rankSize[1] * offset[1];
        allocated = MYARRAY_ALLOCATED_ALL;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
#endif
    };

    void init(int i, int j, int k) {
        size = i * j * k;
        sizeOfMyArray += sizeof(T) * size;
        memory = new T[size]{};
        for (int ioffset = 0; ioffset < 3; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2];
        allocated = MYARRAY_ALLOCATED_ALL;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
#endif
    };

    void init(int i, int j, int k, int l) {
        size = i * j * k * l;
        sizeOfMyArray += sizeof(T) * size;
        memory = new T[size]{};
        for (int ioffset = 0; ioffset < 4; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3];
        allocated = MYARRAY_ALLOCATED_ALL;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
#endif
    };

    void init(int i, int j, int k, int l, int m) {
        size = i * j * k * l * m;
        sizeOfMyArray += sizeof(T) * size;
        memory = new T[size]{};
        for (int ioffset = 0; ioffset < 5; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        rankSize[4] = rankSize[3] * l;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] + rankSize[4] * offset[4];
        allocated = MYARRAY_ALLOCATED_ALL;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
        sizeOfRank[4] = m;
#endif
    };

    void init(int i, T *memory_) {
        size = i;
        sizeOfMyArray += sizeof(T) * size;
        memory = memory_;
        for (int ioffset = 0; ioffset < 1; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        sumOfOffset = offset[0];
        allocated = MYARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
#endif
    };

    void init(int i, int j, T *memory_) {
        size = i * j;
        sizeOfMyArray += sizeof(T) * size;
        memory = memory_;
        for (int ioffset = 0; ioffset < 2; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        sumOfOffset = offset[0] + rankSize[1] * offset[1];
        allocated = MYARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
#endif
    };

    void init(int i, int j, int k, T *memory_) {
        size = i * j * k;
        sizeOfMyArray += sizeof(T) * size;
        memory = memory_;
        for (int ioffset = 0; ioffset < 3; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2];
        allocated = MYARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
#endif
    };

    void init(int i, int j, int k, int l, T *memory_) {
        size = i * j * k * l;
        sizeOfMyArray += sizeof(T) * size;
        memory = memory_;
        for (int ioffset = 0; ioffset < 4; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3];
        allocated = MYARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
#endif
    };

    void init(int i, int j, int k, int l, int m, T *memory_) {
        size = i * j * k * l * m;
        sizeOfMyArray += sizeof(T) * size;
        memory = memory_;
        for (int ioffset = 0; ioffset < 5; ioffset++)
            offset[ioffset] = MYARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        rankSize[4] = rankSize[3] * l;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] + rankSize[4] * offset[4];
        allocated = MYARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef MYARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
        sizeOfRank[4] = m;
#endif
    };

    void setOffset(int i) {
        this->offset[0] = i;
        sumOfOffset = offset[0];
    }

    void setOffset(int i, int j) {
        this->offset[0] = i;
        this->offset[1] = j;
        sumOfOffset = offset[0] + rankSize[1] * offset[1];
    }

    void setOffset(int i, int j, int k) {
        this->offset[0] = i;
        this->offset[1] = j;
        this->offset[2] = k;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2];
    }

    void setOffset(int i, int j, int k, int l) {
        this->offset[0] = i;
        this->offset[1] = j;
        this->offset[2] = k;
        this->offset[3] = l;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3];
    }

    void setOffset(int i, int j, int k, int l, int m) {
        this->offset[0] = i;
        this->offset[1] = j;
        this->offset[2] = k;
        this->offset[3] = l;
        this->offset[4] = m;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] + rankSize[4] * offset[4];
    }

   public:
    void destroy() {
        if (allocated != MYARRAY_ALLOCATED_NONE) {
            if ((allocated & MYARRAY_ALLOCATED_MEMORY) != 0 && memory != nullptr)
                delete[] memory;
            memory = nullptr;
            allocated = MYARRAY_ALLOCATED_NONE;
            size = 0;
        }
    }

    MyArray() {}

    MyArray(const MyArray<T, RANK> &array) {
        size = array.size;
        memory = array.memory;
        std::copy((int *)array.rankSize, ((int *)array.rankSize) + RANK, rankSize);
        std::copy((int *)array.offset, ((int *)array.offset) + RANK, offset);
        allocated = MYARRAY_ALLOCATED_NONE;
#ifdef MYARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sizeOfRank[i] = array.sizeOfRank[i];
#endif
    }

    virtual ~MyArray() {
        destroy();
    }

    MyArray(decltype(INTS)... args) {
        init(args...);
    };

    MyArray(T *memory_, decltype(INTS)... args) {
        init(args..., memory_);
    };

    MyArray(T *memory_) : memory(memory_) {
    }

    MyArray(const int &size_, int *rankSize_, int *offset_, const int &sumOfOffset_, T *memory_
#ifdef MYARRAY_DEBUG
            ,
            int sizeOfRank_[]) {
#else
    ) {
#endif

        size = size_;
        std::copy(rankSize_, rankSize_ + RANK, rankSize);
        std::copy(offset_, offset_ + RANK, offset);
        memory = memory_;
        sumOfOffset = sumOfOffset_;
#ifdef MYARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sizeOfRank[i] = sizeOfRank_[i];
#endif
    };

    //    MyArray(int l, int m, MyArray<T,RANK>* a) {
    //        size = this->rankSize[RANK_2];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l;
    //    };
    //
    //    MyArray(int k, int l, int m, MyArray* a) {
    //        size = this->rankSize[RANK_3];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l + this->rankSize[RANK_3] * k;
    //    };
    //
    //    MyArray(int j, int k, int l, int m, MyArray* a) {
    //        size = this->rankSize[RANK_4];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l + this->rankSize[RANK_3] * k + this->rankSize[RANK_4] * j;
    //    };
    //
    //    MyArray(int i, int j, int k, int l, int m, MyArray* a) {
    //        size = this->rankSize[RANK_5];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l + this->rankSize[RANK_3] * k + this->rankSize[RANK_4] * j +
    //                 this->rankSize[RANK_5] * i;
    //    }

    template <std::size_t = RANK, class = std::index_sequence<INTS...>>
    inline void setOffsets(decltype(INTS)... args) {
        setOffset(args...);
    }

    void setSize(decltype(INTS)... args) {
        destroy();
        init(args...);
    };

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int m) {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 1);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);

        int sumOfOffset = offset[0];

        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this->rankSize[RANK_1],
                                       this->rankSize,
                                       this->offset,
                                       sumOfOffset,
                                       p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int l, int m) {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 2);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this->rankSize[RANK_2], this->rankSize, this->offset, sumOfOffset,
                                       this->memory + this->rankSize[RANK_1] * (m - offset[RANK_1]) + this->rankSize[RANK_2] * (l - offset[RANK_2])
#ifdef MYARRAY_DEBUG
                                           ,
                                       this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int k, int l, int m) {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 3);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);
        MyArray<T, RANK2> array = MyArray<T, RANK2>(this->rankSize[RANK_3], this->rankSize, this->offset, sumOfOffset, p_memory
#ifdef MYARRAY_DEBUG
                                                    ,
                                                    this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int j, int k, int l, int m) {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 4);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this->rankSize[RANK_4], this->rankSize, this->offset, sumOfOffset, p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int i, int j, int k, int l, int m) {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 5);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        MYARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_5] * (i - offset[RANK_5]);
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this->rankSize[RANK_5],
                                       this->rankSize,
                                       this->offset,
                                       sumOfOffset,
                                       p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int m) const {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 1);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        auto this_noconst = const_cast<MyArray<T, RANK> *>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this_noconst->rankSize[RANK_1], this_noconst->rankSize, this_noconst->offset, sumOfOffset, p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int l, int m) const {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 2);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        auto this_noconst = const_cast<MyArray<T, RANK> *>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this_noconst->rankSize[RANK_2], this_noconst->rankSize, this_noconst->offset, sumOfOffset, p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int k, int l, int m) const {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 3);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        auto this_noconst = const_cast<MyArray<T, RANK> *>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this_noconst->rankSize[RANK_3], this_noconst->rankSize, this_noconst->offset, sumOfOffset, p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int j, int k, int l, int m) const {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 4);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        auto this_noconst = const_cast<MyArray<T, RANK> *>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this_noconst->rankSize[RANK_4], this_noconst->rankSize, this_noconst->offset, sumOfOffset, p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline MyArray<T, RANK2> slice(int i, int j, int k, int l, int m) const {
        MYARRAY_CHECK_RANK(RANK, RANK2 + 5);
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        MYARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);
        auto this_noconst = const_cast<MyArray<T, RANK> *>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto *p_memory = memory;
        p_memory += this->rankSize[RANK_5] * (i - offset[RANK_5]);
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = MyArray<T, RANK2>(this_noconst->rankSize[RANK_5], this_noconst->rankSize, this_noconst->offset, sumOfOffset, p_memory
#ifdef MYARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    inline T *getMemory(int idx) {
        return memory + idx;
    }

    inline T *data() {
        return memory;
    }

    inline const T *data() const {
        return memory;
    }

    inline T *data(int m) {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[RANK_1]);
    }

    inline T *data(int l, int m) {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]);
    }

    inline T *data(int k, int l, int m) {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) + this->rankSize[RANK_3] * (k - offset[3]);
    }

    inline T *data(int j, int k, int l, int m) {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) +
               this->rankSize[RANK_4] * (j - offset[4]);
    }

    inline T *data(int i, int j, int k, int l, int m) {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        MYARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) +
               this->rankSize[RANK_4] * (j - offset[4]) + this->rankSize[RANK_5] * (i - offset[5]);
    }

    inline T *data(int m) const {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[RANK_1]);
    }

    inline T *data(int l, int m) const {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]);
    }

    inline T *data(int k, int l, int m) const {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) + this->rankSize[RANK_3] * (k - offset[3]);
    }

    inline T *data(int j, int k, int l, int m) const {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) +
               this->rankSize[RANK_4] * (j - offset[4]);
    }

    inline T *data(int i, int j, int k, int l, int m) const {
        MYARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        MYARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        MYARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        MYARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        MYARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) +
               this->rankSize[RANK_4] * (j - offset[4]) + this->rankSize[RANK_5] * (i - offset[5]);
    }

    inline T sum() const {
        T sum = 0.0;
        for (int i = 0; i < size; ++i) {
            sum += memory[i];
        }
        return sum;
    }

    inline T max() const {
        T mx = 0.0;
        for (int i = 0; i < size; ++i) {
            if (memory[i] > mx)
                mx = memory[i];
        }
        return mx;
    }

    inline T min() const {
        T mx = INT_MAX;
        for (int i = 0; i < size; ++i) {
            if (memory[i] < mx)
                mx = memory[i];
        }
        return mx;
    }

    inline int getSize() const {
        return size;
    }

    inline const int *getRankSize() const {
        return rankSize;
    }

#ifdef MYARRAY_DEBUG
    inline const int *getSizeOfRank() const {
        return sizeOfRank;
    }
#endif

    inline const int *getOffset() const {
        return offset;
    }

    inline T &operator()(int i) {
        MYARRAY_CHECK_RANK(RANK, 1);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline T &operator()(int i, int j) {
        MYARRAY_CHECK_RANK(RANK, 2);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline T &operator()(int i, int j, int k) {
        MYARRAY_CHECK_RANK(RANK, 3);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline T &operator()(int i, int j, int k, int l) {
        MYARRAY_CHECK_RANK(RANK, 4);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline T &operator()(int i, int j, int k, int l, int m) {
        MYARRAY_CHECK_RANK(RANK, 5);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        MYARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }
    inline T &operator()(const FastArray<int, RANK> &idx) {
        int pos = idx(1) - sumOfOffset;
        MYARRAY_CHECK_BOUND(idx(1), offset[0], offset[0] + sizeOfRank[0] - 1);
        for (int i = 1; i < RANK; i++) {
            MYARRAY_CHECK_BOUND(idx(i + 1), offset[i], offset[i] + sizeOfRank[i] - 1);
            pos += rankSize[i] * idx(i + 1);
        }

        return memory[pos];
    }

    inline const T &operator()(const FastArray<int, RANK> &idx) const {
        int pos = idx(1) - sumOfOffset;
        MYARRAY_CHECK_BOUND(idx(1), offset[0], offset[0] + sizeOfRank[0] - 1);
        for (int i = 1; i < RANK; i++) {
            MYARRAY_CHECK_BOUND(idx(i + 1), offset[i], offset[i] + sizeOfRank[i] - 1);
            pos += rankSize[i] * idx(i + 1);
        }

        return memory[pos];
    }

    inline const T &operator()(int i) const {
        MYARRAY_CHECK_RANK(RANK, 1);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline const T &operator()(int i, int j) const {
        MYARRAY_CHECK_RANK(RANK, 2);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline const T &operator()(int i, int j, int k) const {
        MYARRAY_CHECK_RANK(RANK, 3);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline const T &operator()(int i, int j, int k, int l) const {
        MYARRAY_CHECK_RANK(RANK, 4);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline const T &operator()(int i, int j, int k, int l, int m) const {
        MYARRAY_CHECK_RANK(RANK, 5);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        MYARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }

    inline const T &at(int i) const {
        MYARRAY_CHECK_RANK(RANK, 1);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline const T &at(int i, int j) const {
        MYARRAY_CHECK_RANK(RANK, 2);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline const T &at(int i, int j, int k) const {
        MYARRAY_CHECK_RANK(RANK, 3);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline const T &at(int i, int j, int k, int l) const {
        MYARRAY_CHECK_RANK(RANK, 4);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline const T &at(int i, int j, int k, int l, int m) const {
        MYARRAY_CHECK_RANK(RANK, 5);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        MYARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }

    inline T &at(int i) {
        MYARRAY_CHECK_RANK(RANK, 1);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline T &at(int i, int j) {
        MYARRAY_CHECK_RANK(RANK, 2);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline T &at(int i, int j, int k) {
        MYARRAY_CHECK_RANK(RANK, 3);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline T &at(int i, int j, int k, int l) {
        MYARRAY_CHECK_RANK(RANK, 4);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline T &at(int i, int j, int k, int l, int m) {
        MYARRAY_CHECK_RANK(RANK, 5);
        MYARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        MYARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        MYARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        MYARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        MYARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }

    inline MyArray<T, RANK> &operator=(const std::initializer_list<T> &list) {
        MYARRAY_CHECK_SIZE(this->size, list.size());
        int idx = 0;
        for (const T &val : list) {
            memory[idx++] = val;
        }
        return *this;
    }

    inline MyArray<T, RANK> &operator=(const T &val) {
        for (int i = 0; i < size; ++i) {
            memory[i] = val;
        }

        return *this;
    }

    inline MyArray<T, RANK> &operator=(const T *val) {
        for (int i = 0; i < size; ++i) {
            memory[i] = val[i];
        }

        return *this;
    }

    inline MyArray<T, RANK> &operator=(const MyArray<T, RANK> &array) {
        if (allocated == MYARRAY_ALLOCATED_NONE && memory == nullptr) {
            size = array.getSize();

            auto arrayRankSize = array.getRankSize();
            auto arrayOffset = array.getOffset();
            for (int i = 0; i < RANK; ++i) {
                rankSize[i] = arrayRankSize[i];
                offset[i] = arrayOffset[i];
            }

            memory = new T[size];
            auto arrayMemory = array.data();
            for (int i = 0; i < size; ++i) {
                memory[i] = arrayMemory[i];
            }

            sumOfOffset = array.sumOfOffset;

            allocated = MYARRAY_ALLOCATED_ALL;

#ifdef MYARRAY_DEBUG
            auto arraySizeOfRank = array.getSizeOfRank();
            for (int i = 0; i < RANK; i++)
                sizeOfRank[i] = arraySizeOfRank[i];
#endif

            return *this;
        } else {
#ifdef MYARRAY_DEBUG
            MYARRAY_CHECK_SIZE(size, array.getSize());

            auto arraySizeOfRank = array.getSizeOfRank();
            auto arrayOffset = array.getOffset();
            for (int rank = 0; rank < RANK; ++rank) {
                assert(sizeOfRank[rank] == arraySizeOfRank[rank]);
                assert(offset[rank] == arrayOffset[rank]);
            }
#endif
            auto arrayMemory = array.data();
            for (int i = 0; i < size; ++i) {
                memory[i] = arrayMemory[i];
            }
            return *this;
        }
    }

    inline MyArray<T, RANK> operator+(const MyArray<T, RANK> &array) {
        MYARRAY_CHECK_SIZE(size, array.size);
        MyArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < size; ++i) {
            out.memory[i] = this->memory[i] + array.memory[i];
        }
        return out;
    }

    inline MyArray<T, RANK> operator+(const MyArray<T, RANK> &array) const {
        MYARRAY_CHECK_SIZE(size, array.size);
        MyArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < size; ++i) {
            out.memory[i] = this->memory[i] + array.memory[i];
        }
        return out;
    }

    inline MyArray<T, RANK> operator-(const MyArray<T, RANK> &array) {
        MYARRAY_CHECK_SIZE(size, array.size);
        MyArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < size; ++i) {
            out.memory[i] = this->memory[i] - array.memory[i];
        }
        return out;
    }

    inline MyArray<T, RANK> friend operator-(const MyArray<T, RANK> &array) {
        MyArray<T, RANK> out;
        out = array;
        for (int i = 0; i < out.size; ++i) {
            out.memory[i] = -out.memory[i];
        }
        return out;
    }

    inline MyArray<T, RANK> operator-(const MyArray<T, RANK> &array) const {
        MYARRAY_CHECK_SIZE(size, array.size);
        MyArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < size; ++i) {
            out.memory[i] = this->memory[i] - array.memory[i];
        }
        return out;
    }

    inline MyArray<T, RANK> operator*(const MyArray<T, RANK> &array) {
        MYARRAY_CHECK_SIZE(size, array.size);
        MyArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < size; ++i) {
            out.memory[i] = this->memory[i] * array.memory[i];
        }
        return out;
    }

    inline MyArray<T, RANK> operator*(const T &val) {
        MyArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < size; ++i) {
            out.memory[i] = this->memory[i] * val;
        }
        return out;
    }

    inline MyArray<T, RANK> &operator*=(const T &val) {
        for (int i = 0; i < size; ++i) {
            this->memory[i] *= val;
        }
        return *this;
    }

    inline MyArray<T, RANK> operator*(const MyArray<T, RANK> &array) const {
        MyArray<T, RANK> sum;

        sum.size = array.size;
        sum.memory = new T[sum.size];
        sum.allocated = MYARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < size; ++i) {
            sum.memory[i] = this->memory[i] * array.memory[i];
        }
#ifdef MYARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline void operator/=(const T &val) {
        auto rval = 1.0 / val;
        for (int i = 0; i < size; ++i) {
            memory[i] *= rval;
        }
    }

    inline MyArray<T, RANK> operator/(const T &val) {
        MyArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < size; ++i) {
            out.memory[i] = this->memory[i] / val;
        }
        return out;
    }

    inline double sqsum() const {
        double sum = 0;
        for (int i = 0; i < size; ++i) {
            sum += memory[i] * memory[i];
        }

        return sum;
    }

    inline MyArray<T, RANK> copy() const {
        MyArray<T, RANK> array;
        array = *this;
        return array;
    }

    inline MyArray<T, RANK> friend operator+(const T &val, const MyArray<T, RANK> &array) {
        MyArray<T, RANK> sum;

        sum.size = array.size;
        sum.memory = new T[sum.size];
        sum.allocated = MYARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < sum.size; ++i) {
            sum.memory[i] = val + array.memory[i];
        }
#ifdef MYARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline MyArray<T, RANK> friend operator*(const T &val, const MyArray<T, RANK> &array) {
        MyArray<T, RANK> sum;

        sum.size = array.size;
        sum.memory = new T[sum.size];
        sum.allocated = MYARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < sum.size; ++i) {
            sum.memory[i] = val * array.memory[i];
        }
#ifdef MYARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline MyArray<T, RANK> friend operator/(const T &val, const MyArray<T, RANK> &array) {
        MyArray<T, RANK> sum;

        sum.size = array.size;
        sum.memory = new T[sum.size];
        sum.allocated = MYARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < sum.size; ++i) {
            sum.memory[i] = val / array.memory[i];
        }
#ifdef MYARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline T friend dot(const MyArray<T, RANK> &array1, const MyArray<T, RANK> &array2) {
#ifdef MYARRAY_DEBUG
        for (int rank = 0; rank < RANK; ++rank) {
            assert(array1.sizeOfRank[rank] == array2.sizeOfRank[rank]);
            assert(array1.offset[rank] == array2.offset[rank]);
            MYARRAY_CHECK_SIZE(array1.size, array2.size);
        }
#endif

        T sum = 0.0;
        for (int i = 0; i < array1.size; ++i) {
            sum += array1.memory[i] * array2.memory[i];
        }
        return sum;
    }

    inline FastArray<int, RANK> maxloc() {
        int maxloc = 0;
        T maxval = memory[0];
        for (int i = 1; i < this->size; ++i) {
            if (memory[i] > maxval) {
                maxloc = i;
                maxval = memory[i];
            }
        }

        FastArray<int, RANK> maxLocation;

        for (int rank = RANK - 1; rank > 0; --rank) {
            maxLocation.memory[rank] = maxloc / rankSize[rank];
            maxloc -= (maxLocation.memory[rank]) * rankSize[rank];
        }
        maxLocation.memory[0] = maxloc;

        for (int rank = 0; rank < RANK; ++rank) {
            maxLocation.memory[rank] += offset[rank];
        }

        return maxLocation;
    }

    inline FastArray<int, RANK> maxloc(int from, int to) {
        int maxloc = 0;
        T maxval = memory[from - MYARRAY_OFFSET];
        for (int i = from - MYARRAY_OFFSET + 1; i < to - MYARRAY_OFFSET; ++i) {
            if (memory[i] > maxval) {
                maxloc = i;
                maxval = memory[i];
            }
        }

        FastArray<int, RANK> maxLocation;

        for (int rank = RANK - 1; rank > 0; ++rank) {
            maxLocation.memory[rank] = maxloc / rankSize[rank] + MYARRAY_OFFSET;
            maxloc -= maxLocation.memory[rank] * rankSize[rank];
        }
        maxLocation.memory[0] = maxloc + MYARRAY_OFFSET;
        return maxLocation;
    }

    struct Iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;
        using pointer = T *;    // or also value_type*
        using reference = T &;  // or also value_type&

        Iterator(pointer ptr) : m_ptr(ptr){};
        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Iterator &operator++() {
            m_ptr++;
            return *this;
        }

        // Postfix increment
        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const Iterator &a, const Iterator &b) { return a.m_ptr == b.m_ptr; };
        friend bool operator!=(const Iterator &a, const Iterator &b) { return a.m_ptr != b.m_ptr; };

       private:
        pointer m_ptr;
    };

   public:
    Iterator begin() { return Iterator(memory); };
    Iterator end() { return Iterator(memory + size); };  // 200 is out of bounds
};

#define alloc1(array, i) array.setSize(i);
#define alloc2(array, i, j) array.setSize(i, j);
#define alloc3(array, i, j, k) array.setSize(i, j, k);
#define alloc4(array, m, i, j, k) array.setSize(m, i, j, k);
#define alloc5(array, m, i, j, k, l) array.setSize(m, i, j, k, l);

#define allocx(array, m, i, j, k, l, FUNC, ...) FUNC

#define alloc(...) allocx(__VA_ARGS__, \
                          alloc5,      \
                          alloc4,      \
                          alloc3,      \
                          alloc2,      \
                          alloc1)(__VA_ARGS__)

#define alloc02(array, i1, i2)  \
    array.setSize(i2 - i1 + 1); \
    array.setOffsets(i1)
#define alloc03(array, i1, i2, j1) ERROR
#define alloc04(array, i1, i2, j1, j2)       \
    array.setSize(i2 - i1 + 1, j2 - j1 + 1); \
    array.setOffsets(i1, j1)
#define alloc05(array, i1, i2, j1, j2, k1) ERROR
#define alloc06(array, i1, i2, j1, j2, k1, k2)            \
    array.setSize(i2 - i1 + 1, j2 - j1 + 1, k2 - k1 + 1); \
    array.setOffsets(i1, j1, k1)
#define alloc07(array, m1, i1, i2, j1, j2, k1, k2) ERROR
#define alloc08(array, m1, m2, i1, i2, j1, j2, k1, k2)                 \
    array.setSize(m2 - m1 + 1, i2 - i1 + 1, j2 - j1 + 1, k2 - k1 + 1); \
    array.setOffsets(m1, i1, j1, k1)
#define alloc09(array, m1, m2, i1, i2, j1, j2, k1, k2, l1) ERROR
#define alloc10(array, m1, m2, i1, i2, j1, j2, k1, k2, l1, l2)                      \
    array.setSize(m2 - m1 + 1, i2 - i1 + 1, j2 - j1 + 1, k2 - k1 + 1, l2 - l1 + 1); \
    array.setOffsets(m1, i1, j1, k1, l1)

#define alloc0x(array, m1, m2, i1, i2, j1, j2, k1, k2, l1, l2, FUNC, ...) FUNC

#define alloc0(...) alloc0x(__VA_ARGS__, \
                            alloc10,     \
                            alloc09,     \
                            alloc08,     \
                            alloc07,     \
                            alloc06,     \
                            alloc05,     \
                            alloc04,     \
                            alloc03,     \
                            alloc02)(__VA_ARGS__)
