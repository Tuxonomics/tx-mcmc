// mvnormal.h

// NOTE(jonas):
//  - expects transpose of upper triangular matrix from Choleski decomposition
//    so that cholCov is a lower triangular matrix
//  - overwrites results in x

void MvNormal( prng::Generator g, Matrix cholCovT, Matrix x )
{
	ASSERT_MSG( x.cols == 1,
				"MvNormal: x needs to be col vector." );

	for ( u32 i = 0; i < x.rows; ++i ) {
		x[i] = prng::Normal( g );
	}

	MMulTriV( cholCovT, x, 'L' );
}

