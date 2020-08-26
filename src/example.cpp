// example.cpp


#include "mcmc.h"
#include "fwd.h"



f64 normalpdf( f64 x, f64 mu, f64 s ) {
    return 1.0 / (sqrt(2.0*M_PI) * s) * exp( - 0.5 * (x - mu)*(x - mu) / (s*s) );
}


#define M  2.5
#define S  1.0

f64 MultipleNormals( f64 *x, size_t n )
{
    f64 mu = M;
    f64 s  = S;

    f64 res = 0.0;

    for ( u32 i = 0; i < n; ++i ) {
        res += log( normalpdf( x[i], mu, s ) + normalpdf( x[i], -mu, s ) );
    }

    return res;
}

void MultipleNormalsGrad( f64 *g, f64 *x, size_t n )
{
#define E 1e-10

    f64 f0 = MultipleNormals( x, n );

    for ( u32 i = 0; i < n; ++i ) {
        x[i] += E;
        g[i] = ( MultipleNormals( x, n ) - f0 ) / E;
        x[i] -= E;
    }

#undef E
}


// struct FGrad {
//     f64      eps;
//     mcmc_f   f;

//     INLINE
//     void operator() ( f64 *g, f64 *x, size_t n ) {
//         FiniteDifferences( g, x, n, eps, f );
//     }
// };



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


    Results res = RunHMC(
        _ALLOCATOR_DEFAULT,
        eps,
        intTime,
        iMass,
        nIter,
        iVals,
        seed,
        MultipleNormals,
        MultipleNormalsGrad
    );

    // // MPrint(res.sample, "sample");
    // MPrintRLiteral(res.sample, "sample");


    return 0;
}

