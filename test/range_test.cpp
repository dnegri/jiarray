#include <iostream>
#include <ctime>
#include "blitz/array.h"
#include <numeric>

using namespace blitz;

#define m_phif(ig,l,k) d_phif[(k*nxy+l)*ng+ig]
#define m_xsnf(ig,l,k) d_xsnf[(k*nxy+l)*ng+ig]

#define op +

int ng=8;
int nx=16*30;
int ny=16*30;
int nz=300;
int nxy=nx*nz;

inline double abs_rate(int ig, int l, int k, Array<float,3>& xs,  Array<float,3>& flux) {
    return flux(ig,l,k) op xs(ig,l,k);
}

double abs_rate(int l, int k, Array<float,3>& xs,  Array<float,3>& flux) {
    double sum=0.0;
    for (int ig = 0; ig < ng; ++ig) {
        sum += (flux(ig,l,k) op xs(ig,l,k));
    }
    return sum;
}
double abs_rate( int k, Array<float,3>& xs,  Array<float,3>& flux) {
    double sum=0.0;
    for (int l = 0; l < nxy; ++l) {
        for (int ig = 0; ig < ng; ++ig) {
            sum += (flux(ig,l,k) op xs(ig,l,k));
        }
    }
    return sum;
}

 double abs_rate(const Array<float,2>& xs,const  Array<float,2>& flux) {
    double sum=0.0;

    for (int l = 0; l < nxy; ++l) {
        for (int ig = 0; ig < ng; ++ig) {
            sum += (flux(ig,l) op xs(ig,l));
        }
    }

    return sum;
}
double abs_rate(const Array<float,1>& xs,const  Array<float,1>& flux) {
    double sum=0.0;

        for (int ig = 0; ig < ng; ++ig) {
            sum += (flux(ig) op xs(ig));
        }

    return sum;
}

double abs_rate(float* xs,float* flux) {
    double sum=0.0;

    for (int ig = 0; ig < ng; ++ig) {
        sum += *xs++ op *flux++;
    }

    return sum;
}

void range_test() {


    Array<float, 3> phif(ng,nxy,nz,ColumnMajorArray<3>()), xsnf(ng,nxy,nz,ColumnMajorArray<3>());
    Array<float, 2> psi(nxy,nz,ColumnMajorArray<2>());

    Array<float, 3> phif3(ng,nxy,nz), xsnf3(ng,nxy,nz);
    Array<float, 2> psi3(nxy,nz);

    Array<float, 3> phif2(nz,nxy,ng), xsnf2(nz,nxy,ng);
    Array<float, 2> psi2(nz,nxy);

    float * d_phif = new float[ng*nxy*nz];
    float * d_xsnf = new float[ng*nxy*nz];
    float * d_psi = new float[nxy*nz];

    float * p_phif = d_phif, *p_xsnf=d_xsnf;

    auto begin = clock();

    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                *p_phif++ = r;
                *p_xsnf++ = r*0.5;
            }
        }
    }
    auto end = clock();
    auto elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by legacy is %12.5f\n", elapsed_secs);


    begin = clock();
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                phif(ig,l,k) = r;
                xsnf(ig,l,k) = r*0.5;
            }
        }
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by blitz column major is %12.5f\n", elapsed_secs);

    begin = clock();
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                phif2(k,l,ig) = r;
                xsnf2(k,l,ig) = r*0.5;
            }
        }
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by blitz row major is %12.5f\n", elapsed_secs);

    begin = clock();
    for (int ig = 0; ig < ng; ++ig) {
        for (int l = 0; l < nxy; ++l) {
            for (int k = 0; k < nz; ++k) {
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                phif3(ig,l,k) = r;
                xsnf3(ig,l,k) = r*0.5;
            }
        }
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by blitz group major is %12.5f\n", elapsed_secs);

    {
        begin = clock();
        Range all = Range::all();
        auto sum=0.0;
        for (int k = 0; k < nz; ++k) {
            for (int l = 0; l < nxy; ++l) {
//                sum += abs_rate(l,k,xsnf, phif);
                for (int ig = 0; ig < ng; ++ig) {
                    sum+= abs_rate(ig,l,k,xsnf, phif);
                }
            }
        }
        end = clock();

        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi.data() , psi.data()+nxy*nz , sum);

        printf("Elapsed time with blitz range 1 is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
    }

    {
        begin = clock();
        Range all = Range::all();
        auto sum=0.0;
        for (int k = 0; k < nz; ++k) {
            sum+= abs_rate(k,xsnf, phif);
        }
        end = clock();

        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi.data() , psi.data()+nxy*nz , sum);

        printf("Elapsed time with blitz range 1 is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
    }

    {
        begin = clock();
        Range all = Range::all();
        auto sum=0.0;
        for (int k = 0; k < nz; ++k) {
            for (int l = 0; l < nxy; ++l) {
                sum += abs_rate(l,k,xsnf, phif);
            }
        }
        end = clock();

        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi.data() , psi.data()+nxy*nz , sum);

        printf("Elapsed time with blitz range 1 is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
    }

    begin = clock();
    p_phif = d_phif, p_xsnf=d_xsnf;
    auto sum=0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                sum += *(p_phif) op *(p_xsnf);
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(d_psi , d_psi+nxy*nz , sum);
    printf("Elapsed time with pointer is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    int i = 0;
    sum=0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                sum += d_phif[i] op d_xsnf[i];
                i++;
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(d_psi , d_psi+nxy*nz , sum);

    printf("Elapsed time with [i] is %12.5f seconds and sum is %e\n", elapsed_secs, sum);


    begin = clock();
    sum=0.0;
    i = 0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            auto lk0 = k*nxy+l;
            for (int ig = 0; ig < ng; ++ig) {
                sum += d_phif[lk0*ng+ig] op d_xsnf[lk0*ng+ig];
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(d_psi , d_psi+nxy*nz , sum);

    printf("Elapsed time with (ig,lk) index is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum=0.0;
    for (int ig = 0; ig < ng; ++ig) {
        for (int l = 0; l < nxy; ++l) {
            for (int k = 0; k < nz; ++k) {
                sum += (phif3(ig,l,k) op xsnf3(ig,l,k));
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi2.data() , psi2.data()+nxy*nz , sum);

    printf("Elapsed time with blitz group major is %12.5f seconds and sum is %e\n", elapsed_secs, sum);


    begin = clock();
    sum=0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                sum += d_phif[(k*nxy+l)*ng+ig] op d_xsnf[(k*nxy+l)*ng+ig];
            }
        }
    }
    end = clock();
//    accumulate(d_psi , d_psi+nxy*nz , sum);

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time with (ig,l,k) index is %12.5f seconds and sum is %e\n", elapsed_secs, sum);


    begin = clock();
    sum=0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                sum += m_phif(ig,l,k) op m_xsnf(ig,l,k);
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(d_psi , d_psi+nxy*nz , sum);

    printf("Elapsed time with (ig,l,k) macro is %12.5f seconds and sum is %e\n", elapsed_secs, sum);


    begin = clock();
    sum=0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                sum += (phif(ig,l,k) op xsnf(ig,l,k));
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi.data() , psi.data()+nxy*nz , sum);

    printf("Elapsed time with blitz column major is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum=0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                sum += (phif2(k,l,ig) op xsnf2(k,l,ig));
//                psi2(k,l) += (phif(k,l,ig) op xsnf(k,l,ig));
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi2.data() , psi2.data()+nxy*nz , sum);

    printf("Elapsed time with blitz row major is %12.5f seconds and sum is %e\n", elapsed_secs, sum);


    {
        begin = clock();
        Range all = Range::all();
        auto sum=0.0;
        for (int k = 0; k < nz; ++k) {
            sum += abs_rate(xsnf(all,all,k), phif(all,all,k));
        }
        end = clock();

        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi.data() , psi.data()+nxy*nz , sum);

        printf("Elapsed time with blitz range 2 is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
    }

    {
        begin = clock();
        Range all = Range::all();
        auto sum=0.0;
        float* d_xsnf = xsnf.data();
        float* d_phif = phif.data();
        for (int k = 0; k < nz; ++k) {
            for (int l = 0; l < nxy; ++l) {
//                sum += abs_rate(xsnf(all,l,k), phif(all,l,k));
//                const Array<float,1>& xs = xsnf(all,l,k);
//                const Array<float,1>& flux = phif(all,l,k);
//                for (int ig = 0; ig < ng; ++ig) {
//                    sum += m_phif(ig,l,k) op m_xsnf(ig,l,k);
//                }
                    sum += abs_rate(&m_xsnf(0,l,k), &m_phif(0,l,k));

            }
        }
        end = clock();

        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    accumulate(psi.data() , psi.data()+nxy*nz , sum);

        printf("Elapsed time with blitz range 2 is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
    }
}