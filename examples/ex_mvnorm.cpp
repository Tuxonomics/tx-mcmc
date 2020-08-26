// ex_mvnorm.cpp

#include "../src/mcmc.h"


i32 main( i32 argn, const char **args )
{
    prng::Xoshiro256StarStar state;
    prng::Generator g = prng::InitXoshiro256StarStar( &state, 64883763 );

    // for ( u32 i = 0; i < 10; ++i ) {
    //     printf("%.10f\n", prng::Normal(g));
    // }


    Mat(x);
    MInit( _ALLOCATOR_DEFAULT, &x, 10000, 3 );

    Mat(cov);
    MInit( _ALLOCATOR_DEFAULT, &cov, x.cols, x.cols );

    Mat(chol);
    MInit( _ALLOCATOR_DEFAULT, &chol, cov.rows, cov.cols );

    Mat(cholT);
    MInit( _ALLOCATOR_DEFAULT, &cholT, cov.rows, cov.cols );


    // cov[0] = 504; cov[1] = 360; cov[2] = 180;
    // cov[3] = 360; cov[4] = 360; cov[5] = 0;
    // cov[6] = 180; cov[7] =   0; cov[8] = 720;

    cov[0] = 1.0; cov[1] = 0.2; cov[2] = 0.5;
    cov[3] = 0.2; cov[4] = 1.0; cov[5] = 0.0;
    cov[6] = 0.5; cov[7] = 0.0; cov[8] = 1.0;

    // MPrintRLiteral( cov, "covariance" );
    // MPrint( cov, "covariance" );


    MCopy( chol, cov );
    MCholeski( chol, 'U' );

    MTranspose( cholT, chol );

    // MPrint( chol, "choleski" );

    Mat(xSlice);
    xSlice.rows = x.cols; xSlice.cols = 1;

    // MPrint( xSlice, "xSlice" );

    for ( u32 i = 0; i < x.rows; ++i ) {
        xSlice.data = x.data + i * x.cols;
        MvNormal( g, cholT, xSlice );
    }


    MPrintRLiteral( x );


    return 0;
}

