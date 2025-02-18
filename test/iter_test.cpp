#include "../src/JIArray.h"
#include "blitz/array.h"
#include <ctime>
#include <iostream>
#include <numeric>

using namespace blitz;
using namespace dnegri::jiarray;

#define m_phif(ig, l, k) d_phif[(k * nxy + l) * ng + ig]
#define m_xsnf(ig, l, k) d_xsnf[(k * nxy + l) * ng + ig]
#define m_psi(l, k)      d_psi[k * nxy + l]

#define op /

int iter_test() {

    int ng  = 8;
    int nx  = 16 * 30;
    int ny  = 16 * 30;
    int nz  = 300;
    int nxy = nx * nz;

    Array<float, 3> phif(ng, nxy, nz, ColumnMajorArray<3>()), xsnf(ng, nxy, nz, ColumnMajorArray<3>());
    Array<float, 2> psi(nxy, nz, ColumnMajorArray<2>());

    Array<float, 3> phif3(ng, nxy, nz), xsnf3(ng, nxy, nz);
    Array<float, 2> psi3(nxy, nz);

    Array<float, 3> phif2(nz, nxy, ng), xsnf2(nz, nxy, ng);
    Array<float, 2> psi2(nz, nxy);

    float* d_phif = new float[ng * nxy * nz];
    float* d_xsnf = new float[ng * nxy * nz];
    float* d_psi  = new float[nxy * nz];

    float *p_phif = d_phif, *p_xsnf = d_xsnf;

    zfloat3 ji_phif3(ng, nxy, nz);
    zfloat3 ji_xsnf3(ng, nxy, nz);
    zfloat2 ji_psi3(nxy, nz);

    auto begin = clock();

    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                float r   = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                *p_phif++ = r;
                *p_xsnf++ = r * 0.5;
            }
        }
    }
    auto end          = clock();
    auto elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by legacy is %12.5f\n", elapsed_secs);

    begin = clock();
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                float r        = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                phif(ig, l, k) = r;
                xsnf(ig, l, k) = r * 0.5;
            }
        }
    }
    end          = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by blitz column major is %12.5f\n", elapsed_secs);

    begin = clock();
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            for (int ig = 0; ig < ng; ++ig) {
                float r         = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                phif2(k, l, ig) = r;
                xsnf2(k, l, ig) = r * 0.5;
            }
        }
    }
    end          = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by blitz row major is %12.5f\n", elapsed_secs);

    begin = clock();
    for (int ig = 0; ig < ng; ++ig) {
        for (int l = 0; l < nxy; ++l) {
            for (int k = 0; k < nz; ++k) {
                float r         = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                phif3(ig, l, k) = r;
                xsnf3(ig, l, k) = r * 0.5;
            }
        }
    }
    end          = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by blitz group major is %12.5f\n", elapsed_secs);

    begin = clock();
    zfor(k, nz) {
        zfor(l, nxy) {
            zfor(ig, ng) {
                float r            = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                ji_phif3(ig, l, k) = r;
                ji_xsnf3(ig, l, k) = r * 0.5;
            }
        }
    }
    end          = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time for initializing values by jiarray is %12.5f\n", elapsed_secs);

    begin  = clock();
    p_phif = d_phif, p_xsnf = d_xsnf;
    auto sum = 0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            m_psi(l, k) = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                sum += *(p_phif)op * (p_xsnf);
                //                m_psi(l,k) += *p_phif++ op *p_xsnf++;
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(d_psi , d_psi+nxy*nz , sum);
    printf("Elapsed time with pointer is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    int i = 0;
    sum   = 0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            m_psi(l, k) = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                sum += d_phif[i] op d_xsnf[i];
                //                m_psi(l,k) += d_phif[i] op d_xsnf[i];
                i++;
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(d_psi , d_psi+nxy*nz , sum);

    printf("Elapsed time with [i] is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum   = 0.0;
    i     = 0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            auto lk0    = k * nxy + l;
            m_psi(l, k) = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                sum += d_phif[lk0 * ng + ig] op d_xsnf[lk0 * ng + ig];
                //                m_psi(l,k) += d_phif[lk0*ng+ig] op d_xsnf[lk0*ng+ig];
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(d_psi , d_psi+nxy*nz , sum);

    printf("Elapsed time with (ig,lk) index is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum   = 0.0;
    for (int ig = 0; ig < ng; ++ig) {
        for (int l = 0; l < nxy; ++l) {
            for (int k = 0; k < nz; ++k) {
                //                psi3(l,k) = 0.0;
                sum += (phif3(ig, l, k) op xsnf3(ig, l, k));
                //                psi2(k,l) += (phif(k,l,ig) op xsnf(k,l,ig));
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(psi2.data() , psi2.data()+nxy*nz , sum);

    printf("Elapsed time with blitz group major is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum   = 0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            m_psi(l, k) = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                auto                                      lk0 = k * nxy + l;
                sum += d_phif[(k * nxy + l) * ng + ig] op d_xsnf[(k * nxy + l) * ng + ig];
                //                m_psi(l,k) += d_phif[(k*nxy+l)*ng+ig] op d_xsnf[(k*nxy+l)*ng+ig];
            }
        }
    }
    end = clock();
    //    accumulate(d_psi , d_psi+nxy*nz , sum);

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    printf("Elapsed time with (ig,l,k) index is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum   = 0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            m_psi(l, k) = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                auto                       lk0 = k * nxy + l;
                sum += m_phif(ig, l, k) op m_xsnf(ig, l, k);
                //                m_psi(l,k) += m_phif(ig,l,k) op m_xsnf(ig,l,k);
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(d_psi , d_psi+nxy*nz , sum);

    printf("Elapsed time with (ig,l,k) macro is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum   = 0.0;
    zfor(k, nz) {
        zfor(l, nxy) {
            ji_psi3(l, k) = 0.0;
            zfor(ig, ng) {
                sum += (ji_phif3(ig, l, k) op ji_xsnf3(ig, l, k));
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(psi2.data() , psi2.data()+nxy*nz , sum);

    printf("Elapsed time with jiarray is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum   = 0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            psi(l, k) = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                sum += (phif(ig, l, k) op xsnf(ig, l, k));
                //                psi(l,k) += (phif(ig,l,k) op xsnf(ig,l,k));
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(psi.data() , psi.data()+nxy*nz , sum);

    printf("Elapsed time with blitz column major is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    begin = clock();
    sum   = 0.0;
    for (int k = 0; k < nz; ++k) {
        for (int l = 0; l < nxy; ++l) {
            psi2(k, l) = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                sum += (phif2(k, l, ig) op xsnf2(k, l, ig));
                //                psi2(k,l) += (phif(k,l,ig) op xsnf(k,l,ig));
            }
        }
    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    accumulate(psi2.data() , psi2.data()+nxy*nz , sum);

    printf("Elapsed time with blitz row major is %12.5f seconds and sum is %e\n", elapsed_secs, sum);

    return 0;
}

// #include <iostream>
// #include <ctime>
// #include <numeric>
// #include "blitz/array.h"

// #define ffor(i,begin,end) for(int i = begin; i<=end; ++i)

// #define m_xs(ig,lp,la,k) xs[ig+lp*ng+la*ngp+k*ngpa]
// #define m_flux(ig,lp,la,k) flux[ig+lp*ng+la*ngp+k*ngpa]
// #define m_rate(ig,lp,la,k) rate[ig+lp*ng+la*ngp+k*ngpa]

// int ng = 4;
// int np = 256;
// int na = 241;
// int nz = 30;

// int ngp = ng*np;
// int ngpa = ng*np*na;
// int ngpaz = ng*np*na*nz;

// int nloop = 100;

// using namespace blitz;

// void multiply(int k, Array<double,4>& xs, Array<double,4>& flux, Array<double,4>& rate, double& sum) {

//     ffor(la,1,na) {
//         ffor(lp,1,np) {
//             ffor(ig,1,ng) {
//                 rate(ig,lp,la,k) = xs(ig,lp,la,k) * flux(ig,lp,la,k);
//                 sum+=rate(ig,lp,la,k);
//             }
//         }
//     }
// }

// void multiply(Array<double,3> xs, Array<double,3> flux, Array<double,3> rate, double& sum) {

//     ffor(la,1,na) {
//         ffor(lp,1,np) {
//             ffor(ig,1,ng) {
//                 rate(ig,lp,la) = xs(ig,lp,la) * flux(ig,lp,la);
//                 sum+=rate(ig,lp,la);
//             }
//         }
//     }
// }

// void multiply(double* xs, double* flux, double* rate, double& sum) {
//     int idx=-1;
//     ffor(la,1,na) {
//         ffor(lp,1,np) {
//             ffor(ig,1,ng) {
//                 idx++;
//                 rate[idx] = xs[idx] * flux[idx];
//                 sum+=rate[idx];
//             }
//         }
//     }
// }

// void m_multiply(int k, double* xs, double* flux, double* rate, double& sum) {
//     int idx=-1;
//     ffor(la,1,na) {
//         ffor(lp,1,np) {
//             ffor(ig,1,ng) {
//                 idx++;
//                 m_rate(ig,lp,la,k) = m_xs(ig,lp,la,k) * m_flux(ig,lp,la,k);
//                 sum+=m_rate(ig,lp,la,k);
//             }
//         }
//     }
// }

// int main() {

//     Range all = Range::all();

//     Array<double,4> xs(ng, np, na, nz);
//     Array<double,4> flux(ng, np, na, nz);
//     Array<double,4> rate(ng, np, na, nz);

//     auto size = ng*np*na;

//     ffor(i, 1, 10) {

//         {
//             auto xs_ = xs.data();
//             auto flux_ = flux.data();
//             auto rate_ = flux.data();

//             auto begin = clock();
//             auto sum=0.0;
//             ffor(k, 1, nz) {
//                 ffor(il, 1, nloop)
//                 m_multiply(k-1, xs_, flux_, rate_,sum);
//             }
//             auto end = clock();
//             auto elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//             printf("Elapsed time with macro is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
//         }

//         {
//             auto begin = clock();
//             auto sum=0.0;
//             ffor(k, 1, nz) {
//                 ffor(il, 1, nloop)
// //                multiply(xs(all, all, all, k),
// //                         flux(all, all, all, k),
// //                         rate(all, all, all, k), sum);
//                 multiply(k,xs,
//                          flux,
//                          rate,sum);
// //                ffor(la,0,na-1) {
// //                    ffor(lp,0,np-1) {
// //                        ffor(ig,0,ng-1) {
// //                            rate(ig,lp,la,k) = xs(ig,lp,la,k) * flux(ig,lp,la,k);
// //                            sum += rate(ig,lp,la,k);
// //                        }
// //                    }
// //                }
//             }
//             auto end = clock();
//             auto elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

//             printf("Elapsed time with blitz is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
//         }

//         {
//             auto xs_ = xs.data();
//             auto flux_ = flux.data();
//             auto rate_ = flux.data();

//             auto begin = clock();
//             auto sum = 0.0;
//             ffor(k, 1, nz) {
//                 ffor(il, 1, nloop)
//                 multiply(xs_, flux_, rate_, sum);
//                 xs_ += size;
//                 flux_ += size;
//                 rate_ += size;
//             }
//             auto end = clock();
//             auto elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//             printf("Elapsed time with legacy is %12.5f seconds and sum is %e\n", elapsed_secs, sum);
//         }

//     }

// //    subarray(xs, sub);
//     return 0;
// }
