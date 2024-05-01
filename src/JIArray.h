#pragma once

#include "pch.h"
#include "FastArray.h"
#include <vector>

namespace dnegri::jiarray {
    
template <class T, std::size_t RANK, class = std::make_index_sequence<RANK>>
class JIArray;

template <class T, std::size_t RANK, std::size_t... INTS>
class JIArray<T, RANK, std::index_sequence<INTS...>> {
private:
    int n_data      = 0;

public:
    T*  memory    = nullptr;
    int allocated = JIARRAY_ALLOCATED_NONE;
    int rankSize[RANK]{};
    int offset[RANK]{};
    int sumOfOffset = 0;

    static const int RANK_1 = RANK - 1;
    static const int RANK_2 = RANK - 2;
    static const int RANK_3 = RANK - 3;
    static const int RANK_4 = RANK - 4;
    static const int RANK_5 = RANK - 5;

private:
#ifdef JIARRAY_DEBUG
    int sizeOfRank[RANK]{};
#endif
private:
    void init(int i) {
        n_data = i;
        memory = new T[n_data]{};
        for (int ioffset = 0; ioffset < 1; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        allocated   = JIARRAY_ALLOCATED_ALL;
        sumOfOffset = offset[0];
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
#endif
    };

    void init(int i, int j) {
        n_data = i * j;
        memory = new T[n_data]{};
        for (int ioffset = 0; ioffset < 2; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        sumOfOffset = offset[0] + rankSize[1] * offset[1];
        allocated   = JIARRAY_ALLOCATED_ALL;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
#endif
    };

    void init(int i, int j, int k) {
        n_data = i * j * k;
        memory = new T[n_data]{};
        for (int ioffset = 0; ioffset < 3; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2];
        allocated   = JIARRAY_ALLOCATED_ALL;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
#endif
    };

    void init(int i, int j, int k, int l) {
        n_data = i * j * k * l;
        memory = new T[n_data]{};
        for (int ioffset = 0; ioffset < 4; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3];
        allocated   = JIARRAY_ALLOCATED_ALL;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
#endif
    };

    void init(int i, int j, int k, int l, int m) {
        n_data = i * j * k * l * m;
        memory = new T[n_data]{};
        for (int ioffset = 0; ioffset < 5; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        rankSize[4] = rankSize[3] * l;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] +
                      rankSize[4] * offset[4];
        allocated = JIARRAY_ALLOCATED_ALL;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
        sizeOfRank[4] = m;
#endif
    };

    void init(int i, int j, int k, int l, int m, int n) {
        n_data = i * j * k * l * m * n;
        memory = new T[n_data]{};
        for (int ioffset = 0; ioffset < 6; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        rankSize[4] = rankSize[3] * l;
        rankSize[5] = rankSize[4] * m;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] +
                      rankSize[4] * offset[4] + rankSize[5] * offset[5];
        allocated = JIARRAY_ALLOCATED_ALL;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
        sizeOfRank[4] = m;
        sizeOfRank[5] = n;
#endif
    };
    void init(int i, T* memory_) {
        n_data   = i;
        memory = memory_;
        for (int ioffset = 0; ioffset < 1; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        sumOfOffset = offset[0];
        allocated   = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
#endif
    };

    void init(int i, int j, T* memory_) {
        n_data   = i * j;
        memory = memory_;
        for (int ioffset = 0; ioffset < 2; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        sumOfOffset = offset[0] + rankSize[1] * offset[1];
        allocated   = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
#endif
    };

    void init(int i, int j, int k, T* memory_) {
        n_data   = i * j * k;
        memory = memory_;
        for (int ioffset = 0; ioffset < 3; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2];
        allocated   = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
#endif
    };

    void init(int i, int j, int k, int l, T* memory_) {
        n_data   = i * j * k * l;
        memory = memory_;
        for (int ioffset = 0; ioffset < 4; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3];
        allocated   = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
#endif
    };

    void init(int i, int j, int k, int l, int m, T* memory_) {
        n_data   = i * j * k * l * m;
        memory = memory_;
        for (int ioffset = 0; ioffset < 5; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        rankSize[4] = rankSize[3] * l;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] +
                      rankSize[4] * offset[4];
        allocated = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
        sizeOfRank[4] = m;
#endif
    };

    void init(int i, int j, int k, int l, int m, int n, T* memory_) {
        n_data   = i * j * k * l * m * n;
        memory = memory_;
        for (int ioffset = 0; ioffset < 6; ioffset++)
            offset[ioffset] = JIARRAY_OFFSET;
        rankSize[0] = 1;
        rankSize[1] = rankSize[0] * i;
        rankSize[2] = rankSize[1] * j;
        rankSize[3] = rankSize[2] * k;
        rankSize[4] = rankSize[3] * l;
        rankSize[5] = rankSize[4] * m;
        sumOfOffset = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] +
                      rankSize[4] * offset[4] + rankSize[5] * offset[5];
        allocated = JIARRAY_ALLOCATED_RANKSIZE_OFFSET;
#ifdef JIARRAY_DEBUG
        sizeOfRank[0] = i;
        sizeOfRank[1] = j;
        sizeOfRank[2] = k;
        sizeOfRank[3] = l;
        sizeOfRank[4] = m;
        sizeOfRank[5] = n;
#endif
    };

    void setOffset(int i) {
        this->offset[0] = i;
        sumOfOffset     = offset[0];
    }

    void setOffset(int i, int j) {
        this->offset[0] = i;
        this->offset[1] = j;
        sumOfOffset     = offset[0] + rankSize[1] * offset[1];
    }

    void setOffset(int i, int j, int k) {
        this->offset[0] = i;
        this->offset[1] = j;
        this->offset[2] = k;
        sumOfOffset     = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2];
    }

    void setOffset(int i, int j, int k, int l) {
        this->offset[0] = i;
        this->offset[1] = j;
        this->offset[2] = k;
        this->offset[3] = l;
        sumOfOffset     = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3];
    }

    void setOffset(int i, int j, int k, int l, int m) {
        this->offset[0] = i;
        this->offset[1] = j;
        this->offset[2] = k;
        this->offset[3] = l;
        this->offset[4] = m;
        sumOfOffset     = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] +
                      rankSize[4] * offset[4];
    }

    void setOffset(int i, int j, int k, int l, int m, int n) {
        this->offset[0] = i;
        this->offset[1] = j;
        this->offset[2] = k;
        this->offset[3] = l;
        this->offset[4] = m;
        this->offset[5] = n;
        sumOfOffset     = offset[0] + rankSize[1] * offset[1] + rankSize[2] * offset[2] + rankSize[3] * offset[3] +
                      rankSize[4] * offset[4] + rankSize[5] * offset[5];
    }

public:
    void destroy() {
        if (allocated != JIARRAY_ALLOCATED_NONE) {
            if ((allocated & JIARRAY_ALLOCATED_MEMORY) != 0 && memory != nullptr) {
                delete[] memory;
            }
            memory    = nullptr;
            allocated = JIARRAY_ALLOCATED_NONE;
            n_data      = 0;
        }
    }

    void erase() {
        destroy();
        std::fill(rankSize, rankSize + RANK, 0);
        std::fill(offset, offset + RANK, 0);
        sumOfOffset = 0;
#ifdef JIARRAY_DEBUG
        std::fill(sizeOfRank, sizeOfRank + RANK, 0);
#endif
    }

    JIArray() {
    }

    JIArray(std::initializer_list<T> il) : JIArray(il.size()) {
        std::copy(il.begin(), il.end(), memory);
    }


    JIArray(const JIArray<T, RANK>& array) {
        n_data   = array.n_data;
        memory = array.memory;
        std::copy((int*)array.rankSize, ((int*)array.rankSize) + RANK, rankSize);
        std::copy((int*)array.offset, ((int*)array.offset) + RANK, offset);
        allocated = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
#ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sizeOfRank[i] = array.sizeOfRank[i];
#endif
    }

    virtual ~JIArray() {
        destroy();
    }

    JIArray(decltype(INTS)... args) {
        init(args...);
    };

    JIArray(T* memory_, decltype(INTS)... args) {
        init(args..., memory_);
    };

    JIArray(T* memory_)
        : memory(memory_) {
    }

    JIArray(const int& size_, int* rankSize_, int* offset_, const int& sumOfOffset_, T* memory_
#ifdef JIARRAY_DEBUG
            ,
            int sizeOfRank_[]) {
#else
    ) {
#endif

        n_data = size_;
        std::copy(rankSize_, rankSize_ + RANK, rankSize);
        std::copy(offset_, offset_ + RANK, offset);
        memory      = memory_;
        sumOfOffset = sumOfOffset_;
#ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sizeOfRank[i] = sizeOfRank_[i];
#endif
    };

    //    JIArray(int l, int m, JIArray<T,RANK>* a) {
    //        size = this->rankSize[RANK_2];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l;
    //    };
    //
    //    JIArray(int k, int l, int m, JIArray* a) {
    //        size = this->rankSize[RANK_3];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l + this->rankSize[RANK_3] *
    //        k;
    //    };
    //
    //    JIArray(int j, int k, int l, int m, JIArray* a) {
    //        size = this->rankSize[RANK_4];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l + this->rankSize[RANK_3] *
    //        k + this->rankSize[RANK_4] * j;
    //    };
    //
    //    JIArray(int i, int j, int k, int l, int m, JIArray* a) {
    //        size = this->rankSize[RANK_5];
    //        memory = this->memory + this->rankSize[RANK_1] * m + this->rankSize[RANK_2] * l + this->rankSize[RANK_3] *
    //        k + this->rankSize[RANK_4] * j +
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
    inline JIArray<T, RANK2> slice(int m) {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 1);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);

        int sumOfOffset = offset[0];

        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this->rankSize[RANK_1], this->rankSize, this->offset, sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline JIArray<T, RANK2> slice(int l, int m) {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 2);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(
            this->rankSize[RANK_2], this->rankSize, this->offset, sumOfOffset,
            this->memory + this->rankSize[RANK_1] * (m - offset[RANK_1]) + this->rankSize[RANK_2] * (l - offset[RANK_2])
#ifdef JIARRAY_DEBUG
                ,
            this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline JIArray<T, RANK2> slice(int k, int l, int m) {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 3);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);
        JIArray<T, RANK2> array =
            JIArray<T, RANK2>(this->rankSize[RANK_3], this->rankSize, this->offset, sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                              ,
                              this->sizeOfRank
#endif
            );
        return array;
    }

    template <int RANK2>
    inline JIArray<T, RANK2> slice(int j, int k, int l, int m) {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 4);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this->rankSize[RANK_4], this->rankSize, this->offset, sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline JIArray<T, RANK2> slice(int i, int j, int k, int l, int m) {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 5);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        JIARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);
        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_5] * (i - offset[RANK_5]);
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this->rankSize[RANK_5], this->rankSize, this->offset, sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline const JIArray<T, RANK2> slice(int m) const {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 1);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        auto this_noconst = const_cast<JIArray<T, RANK>*>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this_noconst->rankSize[RANK_1], this_noconst->rankSize, this_noconst->offset,
                                       sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline const JIArray<T, RANK2> slice(int l, int m) const {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 2);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        auto this_noconst = const_cast<JIArray<T, RANK>*>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this_noconst->rankSize[RANK_2], this_noconst->rankSize, this_noconst->offset,
                                       sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline const JIArray<T, RANK2> slice(int k, int l, int m) const {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 3);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        auto this_noconst = const_cast<JIArray<T, RANK>*>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this_noconst->rankSize[RANK_3], this_noconst->rankSize, this_noconst->offset,
                                       sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline const JIArray<T, RANK2> slice(int j, int k, int l, int m) const {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 4);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        auto this_noconst = const_cast<JIArray<T, RANK>*>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this_noconst->rankSize[RANK_4], this_noconst->rankSize, this_noconst->offset,
                                       sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    template <int RANK2>
    inline const JIArray<T, RANK2> slice(int i, int j, int k, int l, int m) const {
        JIARRAY_CHECK_RANK(RANK, RANK2 + 5);
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        JIARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);
        auto this_noconst = const_cast<JIArray<T, RANK>*>(this);

        int sumOfOffset = offset[0];
        for (int irank = 1; irank < RANK2; irank++)
            sumOfOffset += this->rankSize[irank] * this->offset[irank];

        auto* p_memory = memory;
        p_memory += this->rankSize[RANK_5] * (i - offset[RANK_5]);
        p_memory += this->rankSize[RANK_4] * (j - offset[RANK_4]);
        p_memory += this->rankSize[RANK_3] * (k - offset[RANK_3]);
        p_memory += this->rankSize[RANK_2] * (l - offset[RANK_2]);
        p_memory += this->rankSize[RANK_1] * (m - offset[RANK_1]);

        auto array = JIArray<T, RANK2>(this_noconst->rankSize[RANK_5], this_noconst->rankSize, this_noconst->offset,
                                       sumOfOffset, p_memory
#ifdef JIARRAY_DEBUG
                                       ,
                                       this_noconst->sizeOfRank
#endif
        );
        return array;
    }

    inline T* getMemory(int idx = 0) {
        return memory + idx;
    }

    inline const T* getMemory(int idx = 0) const {
        return memory + idx;
    }

    inline T* data() {
        return memory;
    }

    inline const T* data() const {
        return memory;
    }

    inline T* data(int m) {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[RANK_1]);
    }

    inline T* data(int l, int m) {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]);
    }

    inline T* data(int k, int l, int m) {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]);
    }

    inline T* data(int j, int k, int l, int m) {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) + this->rankSize[RANK_4] * (j - offset[4]);
    }

    inline T* data(int i, int j, int k, int l, int m) {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        JIARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) + this->rankSize[RANK_4] * (j - offset[4]) +
               this->rankSize[RANK_5] * (i - offset[5]);
    }

    inline T* data(int m) const {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[RANK_1]);
    }

    inline T* data(int l, int m) const {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]);
    }

    inline T* data(int k, int l, int m) const {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]);
    }

    inline T* data(int j, int k, int l, int m) const {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) + this->rankSize[RANK_4] * (j - offset[4]);
    }

    inline T* data(int i, int j, int k, int l, int m) const {
        JIARRAY_CHECK_BOUND(m, offset[RANK_1], offset[RANK_1] + sizeOfRank[RANK_1] - 1);
        JIARRAY_CHECK_BOUND(l, offset[RANK_2], offset[RANK_2] + sizeOfRank[RANK_2] - 1);
        JIARRAY_CHECK_BOUND(k, offset[RANK_3], offset[RANK_3] + sizeOfRank[RANK_3] - 1);
        JIARRAY_CHECK_BOUND(j, offset[RANK_4], offset[RANK_4] + sizeOfRank[RANK_4] - 1);
        JIARRAY_CHECK_BOUND(i, offset[RANK_5], offset[RANK_5] + sizeOfRank[RANK_5] - 1);

        return this->memory + this->rankSize[RANK_1] * (m - offset[1]) + this->rankSize[RANK_2] * (l - offset[2]) +
               this->rankSize[RANK_3] * (k - offset[3]) + this->rankSize[RANK_4] * (j - offset[4]) +
               this->rankSize[RANK_5] * (i - offset[5]);
    }

    const int& size() const {
        return n_data;
    }

    inline T average() const {
        T sum = 0.0;
        for (int i = 0; i < n_data; ++i) {
            sum += memory[i];
        }
        return sum / n_data;
    }

    inline T sum() const {
        T sum = 0.0;
        for (int i = 0; i < n_data; ++i) {
            sum += memory[i];
        }
        return sum;
    }

    inline T max() const {
        T mx = memory[0];
        for (int i = 1; i < n_data; ++i) {
            if (memory[i] > mx)
                mx = memory[i];
        }
        return mx;
    }

    inline T min() const {
        T mx = memory[0];
        for (int i = 1; i < n_data; ++i) {
            if (memory[i] < mx)
                mx = memory[i];
        }
        return mx;
    }

    inline int getSize() const {
        return n_data;
    }

    inline const int* getRankSize() const {
        return rankSize;
    }

#ifdef JIARRAY_DEBUG
    inline const int* getSizeOfRank() const {
        return sizeOfRank;
    }
#endif

    inline const int* getOffset() const {
        return offset;
    }

    inline T& operator()(int i) {
        JIARRAY_CHECK_RANK(RANK, 1);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline T& operator()(int i, int j) {
        JIARRAY_CHECK_RANK(RANK, 2);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline T& operator()(int i, int j, int k) {
        JIARRAY_CHECK_RANK(RANK, 3);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline T& operator()(int i, int j, int k, int l) {
        JIARRAY_CHECK_RANK(RANK, 4);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline T& operator()(int i, int j, int k, int l, int m) {
        JIARRAY_CHECK_RANK(RANK, 5);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        JIARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }

    inline T& operator()(int i, int j, int k, int l, int m, int n) {
        JIARRAY_CHECK_RANK(RANK, 6);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        JIARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        JIARRAY_CHECK_BOUND(n, offset[5], offset[5] + sizeOfRank[5] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m + rankSize[5] * n - sumOfOffset];
    }

    inline T& operator()(const FastArray<int, RANK>& idx) {
        int pos = idx(1) - sumOfOffset;
        JIARRAY_CHECK_BOUND(idx(1), offset[0], offset[0] + sizeOfRank[0] - 1);
        for (int i = 1; i < RANK; i++) {
            JIARRAY_CHECK_BOUND(idx(i + 1), offset[i], offset[i] + sizeOfRank[i] - 1);
            pos += rankSize[i] * idx(i + 1);
        }

        return memory[pos];
    }

    inline const T& operator()(const FastArray<int, RANK>& idx) const {
        int pos = idx(1) - sumOfOffset;
        JIARRAY_CHECK_BOUND(idx(1), offset[0], offset[0] + sizeOfRank[0] - 1);
        for (int i = 1; i < RANK; i++) {
            JIARRAY_CHECK_BOUND(idx(i + 1), offset[i], offset[i] + sizeOfRank[i] - 1);
            pos += rankSize[i] * idx(i + 1);
        }

        return memory[pos];
    }

    inline const T& operator()(int i) const {
        JIARRAY_CHECK_RANK(RANK, 1);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline const T& operator()(int i, int j) const {
        JIARRAY_CHECK_RANK(RANK, 2);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline const T& operator()(int i, int j, int k) const {
        JIARRAY_CHECK_RANK(RANK, 3);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline const T& operator()(int i, int j, int k, int l) const {
        JIARRAY_CHECK_RANK(RANK, 4);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline const T& operator()(int i, int j, int k, int l, int m) const {
        JIARRAY_CHECK_RANK(RANK, 5);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        JIARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }

    inline const T& at(int i) const {
        JIARRAY_CHECK_RANK(RANK, 1);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline const T& at(int i, int j) const {
        JIARRAY_CHECK_RANK(RANK, 2);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline const T& at(int i, int j, int k) const {
        JIARRAY_CHECK_RANK(RANK, 3);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline const T& at(int i, int j, int k, int l) const {
        JIARRAY_CHECK_RANK(RANK, 4);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline const T& at(int i, int j, int k, int l, int m) const {
        JIARRAY_CHECK_RANK(RANK, 5);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        JIARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }

    inline T& at(int i) {
        JIARRAY_CHECK_RANK(RANK, 1);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        return memory[i - sumOfOffset];
    }

    inline T& at(int i, int j) {
        JIARRAY_CHECK_RANK(RANK, 2);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        return memory[i + rankSize[1] * j - sumOfOffset];
    }

    inline T& at(int i, int j, int k) {
        JIARRAY_CHECK_RANK(RANK, 3);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k - sumOfOffset];
    }

    inline T& at(int i, int j, int k, int l) {
        JIARRAY_CHECK_RANK(RANK, 4);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l - sumOfOffset];
    }

    inline T& at(int i, int j, int k, int l, int m) {
        JIARRAY_CHECK_RANK(RANK, 5);
        JIARRAY_CHECK_BOUND(i, offset[0], offset[0] + sizeOfRank[0] - 1);
        JIARRAY_CHECK_BOUND(j, offset[1], offset[1] + sizeOfRank[1] - 1);
        JIARRAY_CHECK_BOUND(k, offset[2], offset[2] + sizeOfRank[2] - 1);
        JIARRAY_CHECK_BOUND(l, offset[3], offset[3] + sizeOfRank[3] - 1);
        JIARRAY_CHECK_BOUND(m, offset[4], offset[4] + sizeOfRank[4] - 1);
        return memory[i + rankSize[1] * j + rankSize[2] * k + rankSize[3] * l + rankSize[4] * m - sumOfOffset];
    }

    inline JIArray<T, RANK>& operator=(const std::initializer_list<T>& list) {
        JIARRAY_CHECK_SIZE(this->n_data, list.size());
        int idx = 0;
        for (const T& val : list) {
            memory[idx++] = val;
        }
        return *this;
    }

    inline JIArray<T, RANK>& operator=(const T& val) {
        assert(memory != nullptr);

        for (int i = 0; i < n_data; ++i) {
            memory[i] = val;
        }

        return *this;
    }

    inline JIArray<T, RANK>& operator=(const std::vector<T>& val) {
        JIARRAY_CHECK_SIZE(n_data, val.size());
        
        for (int i = 0; i < n_data; ++i) {
            memory[i] = val[i];
        }

        return *this;
    }

    inline JIArray<T, RANK>& operator=(const T* array) {
        assert(memory != nullptr);

        for (int i = 0; i < n_data; ++i) {
            memory[i] = array[i];
        }

        return *this;
    }

    inline JIArray<T, RANK>& operator=(const JIArray<T, RANK>& array) {
        if (allocated == JIARRAY_ALLOCATED_NONE && memory == nullptr) {
            n_data = array.getSize();

            auto arrayRankSize = array.getRankSize();
            auto arrayOffset   = array.getOffset();
            for (int i = 0; i < RANK; ++i) {
                rankSize[i] = arrayRankSize[i];
                offset[i]   = arrayOffset[i];
            }

            memory           = new T[n_data];
            auto arrayMemory = array.data();
            for (int i = 0; i < n_data; ++i) {
                memory[i] = arrayMemory[i];
            }

            sumOfOffset = array.sumOfOffset;

            allocated = JIARRAY_ALLOCATED_ALL;

#ifdef JIARRAY_DEBUG
            auto arraySizeOfRank = array.getSizeOfRank();
            for (int i = 0; i < RANK; i++)
                sizeOfRank[i] = arraySizeOfRank[i];
#endif

            return *this;
        } else {
#ifdef JIARRAY_DEBUG
            JIARRAY_CHECK_SIZE(n_data, array.getSize());

            auto arraySizeOfRank = array.getSizeOfRank();
            auto arrayOffset     = array.getOffset();
            for (int rank = 0; rank < RANK; ++rank) {
                assert(sizeOfRank[rank] == arraySizeOfRank[rank]);
                assert(offset[rank] == arrayOffset[rank]);
            }
#endif
            auto arrayMemory = array.data();
            for (int i = 0; i < n_data; ++i) {
                memory[i] = arrayMemory[i];
            }
            return *this;
        }
    }

    inline bool operator==(const JIArray<T, RANK>& array) const {
        if (n_data != array.n_data)
            return false;

        for (int i = 0; i < n_data; ++i) {
            if (array.memory[i] != this->memory[i])
                return false;
        }

        return true;
    }

    inline JIArray<T, RANK>& operator+=(const T& val) {
        for (int i = 0; i < n_data; ++i) {
            this->memory[i] += val;
        }
        return *this;
    }

    inline JIArray<T, RANK> operator+=(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);

        for (int i = 0; i < n_data; ++i) {
            this->memory[i] += array.memory[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> operator+(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] + array.memory[i];
        }
        return out;
    }

    inline JIArray<T, RANK> operator+(const JIArray<T, RANK>& array) const {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] + array.memory[i];
        }
        return out;
    }

    inline JIArray<T, RANK>& operator-=(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);
        for (int i = 0; i < n_data; ++i) {
            memory[i] -= array.memory[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> operator-(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] - array.memory[i];
        }
        return out;
    }

    inline JIArray<T, RANK> friend operator-(const JIArray<T, RANK>& array) {
        JIArray<T, RANK> out;
        out = array;
        for (int i = 0; i < out.n_data; ++i) {
            out.memory[i] = -out.memory[i];
        }
        return out;
    }

    inline JIArray<T, RANK> operator-(const JIArray<T, RANK>& array) const {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] - array.memory[i];
        }
        return out;
    }

    inline JIArray<T, RANK> operator*(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] * array.memory[i];
        }
        return out;
    }

    inline JIArray<T, RANK> operator*(const T& val) {
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] * val;
        }
        return out;
    }

    inline JIArray<T, RANK>& operator*=(const T& val) {
        for (int i = 0; i < n_data; ++i) {
            this->memory[i] *= val;
        }
        return *this;
    }

    inline JIArray<T, RANK>& operator*=(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);

        for (int i = 0; i < n_data; ++i) {
            this->memory[i] *= array.memory[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> operator*(const JIArray<T, RANK>& array) const {
        JIArray<T, RANK> sum;

        sum.n_data      = array.n_data;
        sum.memory    = new T[sum.n_data];
        sum.allocated = JIARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < n_data; ++i) {
            sum.memory[i] = this->memory[i] * array.memory[i];
        }
#ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline void operator/=(const T& val) {
        auto rval = 1.0 / val;
        for (int i = 0; i < n_data; ++i) {
            memory[i] *= rval;
        }
    }

    inline JIArray<T, RANK> operator/=(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);

        for (int i = 0; i < n_data; ++i) {
            this->memory[i] = this->memory[i] / array.memory[i];
        }
        return *this;
    }

    inline JIArray<T, RANK> operator/(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_SIZE(n_data, array.n_data);
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] / array.memory[i];
        }
        return out;
    }


    inline JIArray<T, RANK> operator/(const T& val) {
        JIArray<T, RANK> out;
        out = *this;
        for (int i = 0; i < n_data; ++i) {
            out.memory[i] = this->memory[i] / val;
        }
        return out;
    }

    inline double sqsum() const {
        double sum = 0;
        for (int i = 0; i < n_data; ++i) {
            sum += memory[i] * memory[i];
        }

        return sum;
    }

    inline JIArray<T, RANK> copy() const {
        JIArray<T, RANK> array;
        array = *this;
        return array;
    }

    inline JIArray<T, RANK> friend operator+(const T& val, const JIArray<T, RANK>& array) {
        JIArray<T, RANK> sum;

        sum.n_data      = array.n_data;
        sum.memory    = new T[sum.n_data];
        sum.allocated = JIARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < sum.n_data; ++i) {
            sum.memory[i] = val + array.memory[i];
        }
#ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline JIArray<T, RANK> friend operator*(const T& val, const JIArray<T, RANK>& array) {
        JIArray<T, RANK> sum;

        sum.n_data      = array.n_data;
        sum.memory    = new T[sum.n_data];
        sum.allocated = JIARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < sum.n_data; ++i) {
            sum.memory[i] = val * array.memory[i];
        }
#ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline JIArray<T, RANK> friend operator/(const T& val, const JIArray<T, RANK>& array) {
        JIArray<T, RANK> sum;

        sum.n_data      = array.n_data;
        sum.memory    = new T[sum.n_data];
        sum.allocated = JIARRAY_ALLOCATED_ALL;

        std::copy(array.rankSize, array.rankSize + RANK, sum.rankSize);
        std::copy(array.offset, array.offset + RANK, sum.offset);
        for (int i = 0; i < sum.n_data; ++i) {
            sum.memory[i] = val / array.memory[i];
        }
#ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sum.sizeOfRank[i] = array.sizeOfRank[i];
#endif

        return sum;
    }

    inline T friend dot(const JIArray<T, RANK>& array1, const JIArray<T, RANK>& array2) {
#ifdef JIARRAY_DEBUG
        for (int rank = 0; rank < RANK; ++rank) {
            assert(array1.sizeOfRank[rank] == array2.sizeOfRank[rank]);
            assert(array1.offset[rank] == array2.offset[rank]);
            JIARRAY_CHECK_SIZE(array1.n_data, array2.n_data);
        }
#endif

        T sum = 0.0;
        for (int i = 0; i < array1.n_data; ++i) {
            sum += array1.memory[i] * array2.memory[i];
        }
        return sum;
    }

    inline bool contains(const T& item) const {
        for (int i = 0; i < n_data; ++i) {
            if (memory[i] == item)
                return true;
        }
        return false;
    }

    inline FastArray<int, RANK> findFirst(const T& item) const {
        int loc = -1;
        T   value;
        for (int i = 0; i < this->n_data; ++i) {
            if (memory[i] == item) {
                loc   = i;
                value = memory[i];
                break;
            }
        }

        FastArray<int, RANK> location;
        location = (-1 + JIARRAY_OFFSET);

        if (loc == -1)
            return location;

        for (int rank = RANK - 1; rank > 0; --rank) {
            location.memory[rank] = loc / rankSize[rank];
            loc -= (location.memory[rank]) * rankSize[rank];
        }
        location.memory[0] = loc;

        for (int rank = 0; rank < RANK; ++rank) {
            location.memory[rank] += offset[rank];
        }

        return location;
    }

    inline FastArray<int, RANK> maxloc() {
        int maxloc = 0;
        T   maxval = memory[0];
        for (int i = 1; i < this->n_data; ++i) {
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
        T   maxval = memory[from - JIARRAY_OFFSET];
        for (int i = from - JIARRAY_OFFSET + 1; i < to - JIARRAY_OFFSET; ++i) {
            if (memory[i] > maxval) {
                maxloc = i;
                maxval = memory[i];
            }
        }

        FastArray<int, RANK> maxLocation;

        for (int rank = RANK - 1; rank > 0; ++rank) {
            maxLocation.memory[rank] = maxloc / rankSize[rank] + JIARRAY_OFFSET;
            maxloc -= maxLocation.memory[rank] * rankSize[rank];
        }
        maxLocation.memory[0] = maxloc + JIARRAY_OFFSET;
        return maxLocation;
    }

    inline void shareWith(const JIArray<T, RANK>& array) {
        JIARRAY_CHECK_NOT_ALLOCATED();

        n_data   = array.n_data;
        memory = array.memory;
        std::copy((int*)array.rankSize, ((int*)array.rankSize) + RANK, rankSize);
        std::copy((int*)array.offset, ((int*)array.offset) + RANK, offset);
        allocated = JIARRAY_ALLOCATED_NONE;
        sumOfOffset = array.sumOfOffset;
#ifdef JIARRAY_DEBUG
        for (int i = 0; i < RANK; i++)
            sizeOfRank[i] = array.sizeOfRank[i];
#endif
    }

    struct Iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*; // or also value_type*
        using reference         = T&; // or also value_type&

        Iterator(pointer ptr)
            : m_ptr(ptr){};

        reference operator*() {
            return *m_ptr;
        }

        pointer operator->() {
            return m_ptr;
        }

        const pointer operator->() const {
            return m_ptr;
        }

        // Prefix increment
        Iterator& operator++() {
            m_ptr++;
            return *this;
        }

        const Iterator& operator++() const {
            m_ptr++;
            return *this;
        }

        // Postfix increment
        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        const Iterator operator++(int) const {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const Iterator& a, const Iterator& b) {
            return a.m_ptr == b.m_ptr;
        };
        friend bool operator!=(const Iterator& a, const Iterator& b) {
            return a.m_ptr != b.m_ptr;
        };

    private:
        pointer m_ptr;
    };

public:
    Iterator begin() {
        return Iterator(memory);
    };
    Iterator end() {
        return Iterator(memory + n_data);
    }; // 200 is out of bounds

    const Iterator begin() const {
        return Iterator(memory);
    };
    const Iterator end() const {
        return Iterator(memory + n_data);
    }; // 200 is out of bounds


public:
    std::vector<T> convertToVector() {
        std::vector<T> vec(memory, memory + n_data);
        return vec;
    }
};

#define alloc1(array, i)                array.setSize(i);
#define alloc2(array, i, j)             array.setSize(i, j);
#define alloc3(array, i, j, k)          array.setSize(i, j, k);
#define alloc4(array, m, i, j, k)       array.setSize(m, i, j, k);
#define alloc5(array, m, i, j, k, l)    array.setSize(m, i, j, k, l);
#define alloc6(array, m, i, j, k, l, n) array.setSize(m, i, j, k, l, n);

#define allocx(array, m, i, j, k, l, n, FUNC, ...) FUNC

#define alloc(...) allocx(__VA_ARGS__, alloc6, alloc5, alloc4, alloc3, alloc2, alloc1)(__VA_ARGS__)

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

#define alloc0(...) \
    alloc0x(__VA_ARGS__, alloc10, alloc09, alloc08, alloc07, alloc06, alloc05, alloc04, alloc03, alloc02)(__VA_ARGS__)

#define bool1   JIArray<bool, 1>
#define bool2   JIArray<bool, 2>
#define bool3   JIArray<bool, 3>
#define bool4   JIArray<bool, 4>
#define bool5   JIArray<bool, 5>
#define int1    JIArray<int, 1>
#define int2    JIArray<int, 2>
#define int3    JIArray<int, 3>
#define int4    JIArray<int, 4>
#define int5    JIArray<int, 5>
#define double1 JIArray<double, 1>
#define double2 JIArray<double, 2>
#define double3 JIArray<double, 3>
#define double4 JIArray<double, 4>
#define double5 JIArray<double, 5>
#define double6 JIArray<double, 6>

#define float1 JIArray<float, 1>
#define float2 JIArray<float, 2>
#define float3 JIArray<float, 3>
#define float4 JIArray<float, 4>
#define float5 JIArray<float, 5>

#define string1 JIArray<string, 1>
#define string2 JIArray<string, 2>
#define string3 JIArray<string, 3>
#define string4 JIArray<string, 4>
#define string5 JIArray<string, 5>

#define ffor(i, begin, end)      for (int i = begin; i <= end; ++i)
#define ffor_back(i, begin, end) for (int i = begin; i >= end; --i)
} // namespace dnegri::jiarray
