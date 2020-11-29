// example.cpp


#include "mcmc.h"
#include "fwd.h"


// struct MultiModal {
//     f64 m;
//     f64 s;

//     template <typename T>
//     INLINE
//     T normalpdf( T x, f64 m, f64 s ) {
//         return 1.0 / (std::sqrt(2.0*M_PI) * s) * std::exp( - 0.5 * (x - m)*(x - m) / (s*s) );
//     }

//     template <typename T>
//     T f( T *x, size_t n )
//     {
//         f64 res = 0.0;

//         for ( u32 i = 0; i < n; ++i ) {
//             res += std::log( normalpdf( x[i], m, s ) + normalpdf( x[i], -m, s ) );
//         }

//         return res;
//     }

//     void grad( f64 *g, f64 *x, size_t n )
//     {
//         FGradient( g, x, n, f<Fwd<f64>> );
//     }
// };



f64 _M = 2.5;
f64 _S = 1.0;

template <typename T>
INLINE
T normalpdf( T x, f64 m, f64 s ) {
    return 1.0 / (std::sqrt(2.0*M_PI) * s) * std::exp( - 0.5 * (x - m)*(x - m) / (s*s) );
}


template <typename T>
T target( T *x, size_t n )
{
    T res = 0.0;

    for ( u32 i = 0; i < n; ++i ) {
        res += std::log( normalpdf( x[i], _M, _S ) + normalpdf( x[i], -_M, _S ) );
    }

    return res;
}


void grad( f64 *g, f64 *x, size_t n )
{
    FGradient( g, x, n, target<Fwd<f64>> );
}


void gradFD( f64 *g, f64 *x, size_t n )
{
    FiniteDifferences( g, x, n, 1e-8, target<f64> );
}





i32 main( i32 argn, const char **argv )
{
    u64 seed = 34875901;

    f64 eps     = 0.9;
    u32 intTime = 1;
    u32 nIter   = 1e3;

#define n 2

    f64 iValsD[n] = { 0 };
    Mat(iVals);
    iVals.data = iValsD; iVals.rows = n; iVals.cols = 1;

    f64 iMassD[n*n] = {
        1.0, 0.0,
        0.0, 1.0
    };
    Mat(iMass);
    iMass.data = iMassD; iMass.rows = n; iMass.cols = n;

    // MultiModal mm = { .m = 2.5, .s = 1.0 };

    // printf("%.4f\n", target( iValsD, n ) );


    Results res = RunHMC(
        _ALLOCATOR_DEFAULT,
        eps,
        intTime,
        iMass,
        nIter,
        iVals,
        seed,
        target<f64>,
        grad
    );

    // MPrint(res.sample, "sample");
    MPrintRLiteral(res.sample, "sample_ad");


    Results resfd = RunHMC(
        _ALLOCATOR_DEFAULT,
        eps,
        intTime,
        iMass,
        nIter,
        iVals,
        seed,
        target<f64>,
        gradFD
    );

    MPrintRLiteral(resfd.sample, "sample_fd");


    return 0;
}

