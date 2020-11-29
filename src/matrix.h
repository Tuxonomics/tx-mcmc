// matrix.h


// Row-major matrix interface option for using BLAS and LAPACK.


// Forward declare BLAS/LAPACK stuff
#if 1 // && defined(USE_BLAS) // 0 if BLAS/LAPACK header included
extern "C" {

    typedef enum {CblasRowMajor=101, CblasColMajor=102} CBLAS_LAYOUT;
    typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;
    typedef enum {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
    typedef enum {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
    typedef enum {CblasLeft=141, CblasRight=142} CBLAS_SIDE;

    void cblas_dgemm(
        const int Layout, const int transa, const int transb, const int m,
        const int n, const int k, const double alpha, const double *a,
        const int lda, const double *b, const int ldb, const double beta,
        double *c, const int ldc
    );

    void cblas_dtrmv(
        const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag,
        const int n, const double *a, const int lda, double *x, const int incx
    );


    i32 LAPACKE_dpotrf( i32 matrix_layout, char uplo, i32 n, f64 *a, i32 lda);

}
#endif



struct Matrix {
    f64  *data;
    u32   rows;
    u32   cols;

    f64 &operator[]( u32 idx ) {
        ASSERT_MSG_VA(
            0 <= idx && idx < rows*cols,
            "Index %lu >= %lu, access is out of bounds.\n", idx, rows*cols
        );
        return data[idx];
    }

    f64 const &operator[]( u32 idx ) const {
        ASSERT_MSG_VA(
            0 <= idx && idx < rows*cols,
            "Index %lu >= %lu, access is out of bounds.\n", idx, rows*cols
        );
        return data[idx];
    }
};


#define Mat(x) Matrix x = { 0 };


INLINE
void MInit( Allocator al, Matrix *m, u32 rows, u32 cols )
{
    m->rows  = rows;
    m->cols  = cols;
    m->data  = (f64 *) Alloc( al, (rows * cols * sizeof(f64)) );

    ASSERT_MSG( m->data, "MatMake: allocation failed!" );
}


INLINE
void MInitZero( Allocator al, Matrix *m, u32 rows, u32 cols )
{
    m->rows  = rows;
    m->cols  = cols;
    m->data  = (f64 *) Calloc( al, rows * cols, sizeof(f64) );

    ASSERT_MSG( m->data, "MatMake: allocation failed!" );
}


INLINE
void MInitIdentity( Allocator al, Matrix *m, u32 i )
{
    m->rows  = i;
    m->cols  = i;
    m->data  = (f64 *) Calloc( al, i*i, sizeof(f64) );

    ASSERT_MSG( m->data, "MatMake: allocation failed!" );

    for ( u32 j = 0; j < i; ++j ) {
        (*m)[j*i + j] = 1.0;
    }
}


INLINE
void MFree( Allocator al, Matrix *m )
{
    if ( m->data ) {
        Free(al, m->data);
    }
    m->rows = 0;
    m->cols = 0;
}


void MPrint( Matrix m, const char *name = NULL )
{
    printf("Matrix (%s) = {\n", name);
    printf("\t.rows = %u\n", m.rows);
    printf("\t.cols = %u\n", m.cols);
    printf("\t.data = {\n\n");
    if ( m.data ) {
        u32 dim;
        u32 prec = 1 + (u32) std::log10( (f64) (m.rows * m.cols) );
        for (u32 i=0; i<m.rows; ++i) {
            for (u32 j=0; j<m.cols; ++j) {
                dim = i*m.cols + j;
                printf("[%.*d] %15.9f  ", prec, dim, m[dim]);
            }
            printf("\n");
        }
        printf("\n\t}\n");
    }
    else {
        printf("\t- NULL -\n");
    }
    printf("}\n");
}

void MPrintRLiteral( Matrix m, const char *name = "data" )
{
    u32 j=0, i=0;

    printf("%s = cbind(\n", name);
    for ( j = 0; j < m.cols; ++j ) {
        printf("c( %.5f", m[i * m.cols + j]);
        for ( i = 1; i < m.rows; ++i ) {
            printf(", %.5f", m[i * m.cols + j]);
        }
        i = 0;
        printf(")");
        if ( j < m.cols-1 ) {
            printf(",\n");
        }
    }
    printf(")\n");
}


INLINE
void MCopy( Matrix z, Matrix x )
{
    ASSERT_MSG( z.rows == x.rows, "MCopy: rows of input don't match." );
    ASSERT_MSG( z.cols == x.cols, "MCopy: cols of input don't match." );
    memcpy( z.data, x.data, x.rows*x.cols*sizeof(f64) );
}


INLINE
void MAdd( Matrix c, Matrix a, Matrix b )
{
    ASSERT_MSG( c.rows == a.rows && c.rows == b.rows,
                "MAdd: rows of input don't match." );
    ASSERT_MSG( c.cols == a.cols && c.cols == b.cols,
                "MAdd: cols of input don't match." );
    for ( u32 i = 0; i < a.rows*a.cols; ++i ) {
        c[i] = a[i] + b[i];
    }
}


INLINE
void MScale( Matrix c, Matrix x, f64 a )
{
    ASSERT_MSG( c.rows == x.rows,
                "MAdd: rows of input don't match." );
    ASSERT_MSG( c.cols == x.cols,
                "MAdd: cols of input don't match." );
    for ( u32 i = 0; i < x.rows*x.cols; ++i ) {
        c[i] = x[i] * a;
    }
}


INLINE
void MSet( Matrix m, f64 val )
{
    for ( u32 i=0; i<(m.rows * m.cols); ++i ) {
        m[i] = val;
    }
}


INLINE
void MZero( Matrix m )
{
    memset( m.data, 0, m.rows*m.cols*sizeof(f64) );
}


INLINE
void MSetRow( Matrix m, u32 i, Matrix row )
{
    ASSERT_MSG( m.cols == row.rows*row.cols,
                "MSetRow: Size of row needs to be equal to columns in m." );
    u32 idx = i * m.cols;
    memcpy( m.data + idx, row.data, row.rows * row.cols * sizeof(f64) );
}


INLINE
void MMul( Matrix c, Matrix a, Matrix b )
{
    ASSERT_MSG( c.data != a.data && c.data != b.data,
                "MMul: can't be used in-place!" );
    ASSERT_MSG(
        a.rows == c.rows && a.cols == b.rows && b.cols == c.cols,
        "MatMul: dimensions of input don't match!" );

#if defined(USE_BLAS)
    cblas_dgemm (
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        (i32)a.rows, (i32)b.cols, (i32)a.cols,
        1.0,
        a.data, (i32)a.cols,
        b.data, (i32)b.cols,
        0,
        c.data, (i32)b.cols
    );
#else
    /* A is n x m, B is m x p, C is n x p */

    u32 n = a.rows;
    u32 m = a.cols;
    u32 p = b.cols;

    i32 i, j, k;

    f64 *ra  = a.data;
    f64 *rb  = b.data;
    f64 *rc  = c.data;

    memset(c.data, 0, n * p * sizeof(f64));

    for ( i = 0; i < n; ++i ) {
        for ( k = 0; k < m; ++k ) {
            for ( j = 0; j < p; ++j ) {
                rc[j] = rc[j] + ra[k] * rb[j];
            }
            rb += p;
        }
        ra += m;
        rc += p;
        rb -= m*p;
    }
#endif
}


INLINE
void MTranspose( Matrix z, Matrix x )
{
    ASSERT_MSG( z.rows == x.cols && z.cols == x.rows,
                "MTranspose: dimensions of input don't match." );
    ASSERT_MSG( z.data != x.data,
                "MTranspose: can't be used in-place." );

    for ( u32 i = 0; i < x.rows; ++i ) {
        for ( u32 j = 0; j < x.cols; ++j ) {
            z[ j*z.cols + i ] = x[ i*x.cols + j ];
        }
    }
}




#if !defined(USE_BLAS)

// TODO: use algorithm for upper too

b32 mCholeskiLower( f64 *a, i32 n )
// Numerical Recipes in C
// Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky T
// decomposition, A = L · L'. On input, only the upper triangle of a need be given; it is not
// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
// elements which are returned in p[1..n].
{
#define A(i,j) a[i*n + j]

    i32 i, j, k;
    f64 sum;

    for ( i = 0; i < n; ++i ) {
        for ( j = i; j < n; ++j ) {
            for ( sum = A(i,j), k = i-1; k >= 0; --k ){
                sum -= A(i,k) * A(j,k);
            }
            if (i == j) {
                if ( sum <= 0.0 ) {
                    return 0;
                }
                A(i,i) = std::sqrt(sum);
            }
            else {
                A(j,i) = sum / A(i,i);
            }
        }
    }
    return 1;

#undef A
}

#endif


// TODO: use Choleski decomposition to check if symmetric matrix is positive definite


/* Choleski decomposition in-place, uplo is 'U' for upper, else lower */
b32 MCholeski( Matrix x, char uplo )
{
    ASSERT_MSG( x.rows == x.cols, "MCholeski: dst and x need to be square matrices!" );

    u32 n = x.rows;

    Mat(row); row.data = x.data; row.rows = 1;

#if defined(USE_BLAS)

    i32 info = LAPACKE_dpotrf( CblasRowMajor, uplo, (i32) n, x.data, (i32) n );

    if ( uplo == 'U' ) {
        for ( u32 r=1; r<n; ++r ) {
            row.data += n;
            row.cols  = r;
            MZero( row );
        }
    }
    else {
        for ( u32 r=0; r<(n-1); ++r ) {
            row.data = x.data + r*n + r+1;
            row.cols = n - r - 1;
            MZero( row );
        }
    }

    return (info == 0);

#else

    b32 success = mCholeskiLower( x.data, (i32) n );

    if ( !success ) {
        return 0;
    }

    if ( uplo == 'U' ) {
        Mat(tmp);
        MInit( _ALLOCATOR_DEFAULT, &tmp, n, n );
        defer({ MFree(_ALLOCATOR_DEFAULT, &tmp ); });

        MTranspose( tmp, x );

        row.data = tmp.data;

        for ( u32 r=1; r<n; ++r ) {
            row.data += n;
            row.cols  = r;
            MZero( row );
        }

        MCopy(x, tmp);
    }
    else {
        for ( u32 r=0; r<(n-1); ++r ) {
            row.data = x.data + r*n + r+1;
            row.cols = n - r - 1;
            MZero( row );
        }
    }

    return 1;
#endif
}


#if !defined(USE_BLAS)

// Lower Triangular Multiply with Vector
//
// | x  0  0  0 |  | x |
// | x  x  0  0 |  | x |
// | x  x  x  0 |  | x |
// | x  x  x  x |  | x |

void mTrmvLower( f64 *a, f64 *x, ssize_t n )
{
#define A(i,j) a[i*n + j]

    for ( ssize_t i = n-1; i >= 0; --i ) {
        x[i] = A(i,i) * x[i];
        for ( ssize_t j = i-1; j >= 0; --j ) {
            x[i] += A(i,j) * x[j];
        }
    }

#undef A
}



// Upper Triangular Multiply with Vector
//
// | x  x  x  x |  | x |
// | 0  x  x  x |  | x |
// | 0  0  x  x |  | x |
// | 0  0  0  x |  | x |

void mTrmvUpper( f64 *a, f64 *x, ssize_t n )
{
#define A(i,j) a[i*n + j]

    for ( ssize_t i = 0; i < n; ++i ) {
        x[i] = A(i,i) * x[i];
        for ( ssize_t j = i+1; j < n; ++j ) {
            x[i] += A(i,j) * x[j];
        }
    }

#undef A
}

#endif

/* Multiply triangular matrix with vector, result is written to x */
void MMulTriV( Matrix chol, Matrix x, char uplo )
{
    ASSERT_MSG(
        chol.rows == chol.cols && chol.rows == x.rows * x.cols,
        "MTriV: dimensions of input don't match!");

    u32 n = chol.rows;

#if defined(USE_BLAS)

    if ( uplo == 'U' ) {
        cblas_dtrmv(
            CblasRowMajor, CblasUpper, CblasNoTrans,
            CblasNonUnit, (i32) n, chol.data, (i32) n, x.data, 1
        );
    }
    else {
        cblas_dtrmv(
            CblasRowMajor, CblasLower, CblasNoTrans,
            CblasNonUnit, (i32) n, chol.data, (i32) n, x.data, 1
        );
    }

#else

    if ( uplo == 'U' ) {
        mTrmvUpper( chol.data, x.data, (ssize_t) n );
    }
    else {
        mTrmvLower( chol.data, x.data, (ssize_t) n );
    }

#endif

}


b32 MEqual( Matrix a, Matrix b, f64 eps )
{
    if ( a.rows != b.rows || a.cols != b.cols ) {
        return 0;
    }

    for ( u32 i = 0; i < a.rows*a.cols; ++i ) {
        if ( ! Equal(a[i], b[i], eps) ) {
            return 0;
        }
    }

    return 1;
}


#if TEST

void test_init( void )
{
    Mat(x);

    MInit( _ALLOCATOR_DEFAULT, &x, 3, 3 );
    defer({ MFree( _ALLOCATOR_DEFAULT, &x ); });
}


void test_init_zero( void )
{
    Mat(x);

    MInitZero( _ALLOCATOR_DEFAULT, &x, 3, 3);
    defer({ MFree(_ALLOCATOR_DEFAULT, &x ); });

    for ( u32 i = 0; i < x.rows*x.cols; ++i ) {
        ASSERT( Equal( x[i], 0.0, 1e-12 ) );
    }
}


void test_init_identity( void )
{
    Mat(x);

    MInitIdentity( _ALLOCATOR_DEFAULT, &x, 3);
    defer({ MFree(_ALLOCATOR_DEFAULT, &x); });

    f64 test = 0.0;

    for ( u32 i = 0; i < x.rows; ++i ) {
        for ( u32 j = 0; j < x.cols; ++j ) {
            if ( i == j ) {
                test = 1.0;
            }
            else {
                test = 0.0;
            }
            ASSERT( Equal( x[i*x.cols + j], test, 1e-12 ) );
        }
    }
}


void test_copy( void )
{
    Mat(x);
    Mat(y);

    MInitIdentity(_ALLOCATOR_DEFAULT, &x, 3);
    MInit(_ALLOCATOR_DEFAULT, &y, 3, 3);
    defer({
        MFree(_ALLOCATOR_DEFAULT, &x);
        MFree(_ALLOCATOR_DEFAULT, &y);
    });

    MCopy(y, x);

    MEqual( x, y, 1e-12 );
}


void test_add_scale( void )
{
    Mat(x);
    Mat(y);
    Mat(z);

    MInitIdentity(_ALLOCATOR_DEFAULT, &z, 3);
    MInitIdentity(_ALLOCATOR_DEFAULT, &x, 3);
    MInit(_ALLOCATOR_DEFAULT, &y, 3, 3);
    defer({
        MFree(_ALLOCATOR_DEFAULT, &x);
        MFree(_ALLOCATOR_DEFAULT, &y);
        MFree(_ALLOCATOR_DEFAULT, &z);
    });

    MAdd( y, x, x );
    MScale( z, z, 2.0 );

    MEqual( z, y, 1e-12 );
}


void test_set_zero( void )
{
    Mat(x);
    Mat(y);

    MInit(_ALLOCATOR_DEFAULT, &x, 3, 3);
    MInit(_ALLOCATOR_DEFAULT, &y, 3, 3);
    defer({
        MFree(_ALLOCATOR_DEFAULT, &x);
        MFree(_ALLOCATOR_DEFAULT, &y);
    });

    MSet( x, 0.0 );
    MZero( y );

    MEqual( y, x, 1e-12 );
}


void test_set_row( void )
{
    Mat(x);
    Mat(y);
    Mat(z);

    MInitIdentity(_ALLOCATOR_DEFAULT, &x, 3);
    MInitZero(_ALLOCATOR_DEFAULT, &y, 1, 3);
    defer({
        MFree(_ALLOCATOR_DEFAULT, &x);
        MFree(_ALLOCATOR_DEFAULT, &y);
    });

    MSetRow( x, 1, y );

    z.data = x.data + x.cols; z.rows = 1; z.cols = x.cols;

    MEqual( z, y, 1e-12 );
}


void test_mmul( void )
{
    Mat(x);
    Mat(y);
    Mat(z);

    MInitIdentity(_ALLOCATOR_DEFAULT, &x, 3);
    MInit(_ALLOCATOR_DEFAULT, &y, 3, 3);
    MInit(_ALLOCATOR_DEFAULT, &z, 3, 3);
    defer({
        MFree(_ALLOCATOR_DEFAULT, &x);
        MFree(_ALLOCATOR_DEFAULT, &y);
        MFree(_ALLOCATOR_DEFAULT, &z);
    });

    y[0] = 1.0; y[1] = 2.0; y[2] = 3.0;
    y[3] = 4.0; y[4] = 5.0; y[5] = 6.0;
    y[6] = 7.0; y[7] = 8.0; y[8] = 9.0;

    MMul( z, y, x );

    MEqual( z, x, 1e-12 );
}


void test_choleski( void )
{
#define n 3

    f64 xData[n*n] = {
        4.0,  12.0, -16.0,
        12.0, 37.0, -43.0,
        -16.0, -43.0, 98.0
    };

    f64 uData[n*n] = { 0 };
    f64 lData[n*n] = { 0 };
    f64 tData[n*n] = { 0 };

    Matrix x = { .data = xData, .rows = n, .cols = n };
    Matrix u = { .data = uData, .rows = n, .cols = n };
    Matrix l = { .data = lData, .rows = n, .cols = n };
    Matrix t = { .data = tData, .rows = n, .cols = n };


    MCopy( u, x );

    MCholeski( u, 'U' );

    MTranspose( l, u );

    MMul( t, l, u );

    ASSERT( MEqual( t, x, 1e-10 ) );


    MCopy( l, x );

    MCholeski( l, 'L' );

    MTranspose( u, l );

    MMul( t, l, u );

    ASSERT( MEqual( t, x, 1e-10 ) );


#undef n
}


void test_mmultriv_u( void )
{
#define n 3
    f64 xData[] = {
        4.0,  12.0, -16.0,
        12.0, 37.0, -43.0,
        -16.0, -43.0, 98.0
    };
    Matrix x = { .data = xData, .rows = n, .cols = n };

    f64 lData[n*n] = { 0 };
    Matrix l = { .data = lData, .rows = n, .cols = n };

    f64 vData[] = { 2.5, 0.5, -0.5 };
    Matrix v = { .data = vData, .rows = n, .cols = 1 };

    f64 tData[n*n] = { 0 };
    Matrix t = { .data = tData, .rows = n, .cols = 1 };

    MCholeski( x, 'U' );

    MMul( t, x, v );

    MMulTriV( x, v, 'U' );


    ASSERT( MEqual( t, v, 1e-10) );

#undef n
}

void test_mmultriv_l( void )
{
#define n 3
    f64 xData[] = {
        4.0,  12.0, -16.0,
        12.0, 37.0, -43.0,
        -16.0, -43.0, 98.0
    };
    Matrix x = { .data = xData, .rows = n, .cols = n };

    f64 lData[n*n] = { 0 };
    Matrix l = { .data = lData, .rows = n, .cols = n };

    f64 vData[] = { 2.5, 0.5, -0.5 };
    Matrix v = { .data = vData, .rows = n, .cols = 1 };

    f64 tData[n*n] = { 0 };
    Matrix t = { .data = tData, .rows = n, .cols = 1 };

    MCholeski( x, 'L' );

    MMul( t, x, v );

    MMulTriV( x, v, 'L' );


    ASSERT( MEqual( t, v, 1e-10) );

#undef n
}



#endif

