// fwd.h
//
// A simple scalar-based implementation of forward mode automatic differentiation.
// This is for testing and applications where a small number of parameters
// can be expected.
//
// A C implementation using macros can be found here:
// https://github.com/Tuxonomics/ForwardModeAD/blob/master/src/forward_normal/fw_univariate.h


#include <cmath>

#if 0
    typedef double f64;
#endif


template <typename T>
struct Fwd {
    T val;
    T dot;


    Fwd( void ) {};

    Fwd( f64 val ) : val( val ) { dot = T(0.0); };

    Fwd( Fwd<f64> x ) : val( x.val ), dot( x.dot ) {};

    Fwd( Fwd<Fwd<f64>> x ) : val( x.val ), dot( x.dot ) {};

    Fwd( T val, T dot ) : val( val ), dot( dot ) {};


    explicit operator T& ()       { return val; }
    explicit operator T  () const { return val; }


    // unary operator overloads

    Fwd& operator=( f64 newVal );
    Fwd  operator-();

    Fwd  operator+() const
    {
        return *this;
    }


    // binary operator overloads

    Fwd operator+( Fwd rhs );
    Fwd operator+( f64 rhs );


    Fwd operator-( Fwd rhs );
    Fwd operator-( f64 rhs );


    Fwd operator*( Fwd rhs );
    Fwd operator*( f64 rhs );


    Fwd operator/( Fwd rhs );
    Fwd operator/( f64 rhs );


    Fwd& operator+=( Fwd rhs );
    Fwd& operator+=( f64 rhs );

    Fwd& operator-=( Fwd rhs );
    Fwd& operator-=( f64 rhs );

    Fwd& operator*=( Fwd rhs );
    Fwd& operator*=( f64 rhs );

    Fwd& operator/=( Fwd rhs );
    Fwd& operator/=( f64 rhs );
};


// ------------------------ //
// unary operator overloads //
// ------------------------ //

template <typename T>
Fwd<T>& Fwd<T>::operator=( f64 newVal )
{
    val = T(newVal);
    dot = T(0.0);
    return *this;
}

template <typename T>
Fwd<T> Fwd<T>::operator-()
{
    return Fwd( -val, -dot );
}


// ------------------------- //
// binary operator overloads //
// ------------------------- //

template <typename T>
Fwd<T> Fwd<T>::operator+( Fwd<T> rhs )
{
    return Fwd<T>( val + rhs.val, dot + rhs.dot );
}

template <typename T>
Fwd<T> Fwd<T>::operator+( f64 rhs )
{
    return Fwd<T>( val + rhs, dot );
}

template <typename T>
Fwd<T> operator+( f64 lhs, Fwd<T> rhs )
{
    return rhs + lhs;
}


template <typename T>
Fwd<T> Fwd<T>::operator-( Fwd<T> rhs )
{
    return Fwd<T>( val - rhs.val, dot - rhs.dot );
}

template <typename T>
Fwd<T> Fwd<T>::operator-( f64 rhs )
{
    return Fwd<T>( val - rhs, dot );
}

template <typename T>
Fwd<T> operator-( f64 lhs, Fwd<T> rhs )
{
    return Fwd<T>( lhs - rhs.val, -rhs.dot );
}


template <typename T>
Fwd<T> Fwd<T>::operator*( Fwd<T> rhs )
{
    return Fwd<T>( val * rhs.val, (val * rhs.dot) + (dot * rhs.val) );
}

template <typename T>
Fwd<T> Fwd<T>::operator*( f64 rhs )
{
    return Fwd<T>( val * rhs, dot * rhs );
}

template <typename T>
Fwd<T> operator*( f64 lhs, Fwd<T> rhs )
{
    return rhs * lhs;
}


template <typename T>
Fwd<T> Fwd<T>::operator/( Fwd<T> rhs )
{
    return Fwd<T>( val / rhs.val, ((dot*rhs.val) - (val*rhs.dot)) / (rhs.val*rhs.val));
}

template <typename T>
Fwd<T> Fwd<T>::operator/( f64 rhs )
{
    return Fwd<T>( val / rhs, dot / rhs );
}

template <typename T>
Fwd<T> operator/( f64 lhs, Fwd<T> rhs )
{
    return Fwd<T>( lhs / rhs.val, - (lhs*rhs.dot) / (rhs.val*rhs.val) );
}


template <typename T>
Fwd<T>& Fwd<T>::operator+=( Fwd<T> rhs )
{
    *this = *this + rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator+=( f64 rhs )
{
    *this = *this + rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator-=( Fwd<T> rhs )
{
    *this = *this - rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator-=( f64 rhs )
{
    *this = *this - rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator*=( Fwd<T> rhs )
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator*=( f64 rhs )
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator/=( Fwd<T> rhs )
{
    *this = *this / rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator/=( f64 rhs )
{
    *this = *this / rhs;
    return *this;
}


// ----------------------------- //
// comparison operator overloads //
// ----------------------------- //


template <typename T>
bool operator==( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val == rhs.val;
}

template <typename T>
bool operator==( Fwd<T> lhs, f64 rhs )
{
    return lhs.val == rhs;
}

template <typename T>
bool operator==( f64 lhs, Fwd<T> rhs )
{
    return lhs == rhs.val;
}


template <typename T>
bool operator!=( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val != rhs.val;
}

template <typename T>
bool operator!=( Fwd<T> lhs, f64 rhs )
{
    return lhs.val != rhs;
}

template <typename T>
bool operator!=( f64 lhs, Fwd<T> rhs )
{
    return lhs != rhs.val;
}


template <typename T>
bool operator<( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val < rhs.val;
}

template <typename T>
bool operator<( Fwd<T> lhs, f64 rhs )
{
    return lhs.val < rhs;
}

template <typename T>
bool operator<( f64 lhs, Fwd<T> rhs )
{
    return lhs < rhs.val;
}


template <typename T>
bool operator>( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val > rhs.val;
}

template <typename T>
bool operator>( Fwd<T> lhs, f64 rhs )
{
    return lhs.val > rhs;
}

template <typename T>
bool operator>( f64 lhs, Fwd<T> rhs )
{
    return lhs > rhs.val;
}


template <typename T>
bool operator<=( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val <= rhs.val;
}

template <typename T>
bool operator<=( Fwd<T> lhs, f64 rhs )
{
    return lhs.val <= rhs;
}

template <typename T>
bool operator<=( f64 lhs, Fwd<T> rhs )
{
    return lhs <= rhs.val;
}


template <typename T>
bool operator>=( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val >= rhs.val;
}

template <typename T>
bool operator>=( Fwd<T> lhs, f64 rhs )
{
    return lhs.val >= rhs;
}

template <typename T>
bool operator>=( f64 lhs, Fwd<T> rhs )
{
    return lhs >= rhs.val;
}



// ----------------------------- //
// elementary functions overload //
// ----------------------------- //


namespace std {

template <typename T>
Fwd<T> sqrt( Fwd<T> x )
{
    T tmp = std::sqrt( x.val );
    return Fwd<T>( tmp, (T(0.5)*x.dot) / tmp );
}

template <typename T>
Fwd<T> pow( Fwd<T> x, f64 a )
{
    T tmp = std::pow( x.val, a - 1.0 );
    return Fwd<T>( x.val * tmp, T(a) * tmp * x.dot );
}

template <typename T>
Fwd<T> sin( Fwd<T> x )
{
    return Fwd<T>( std::sin(x.val), std::cos(x.val)*x.dot );
}

template <typename T>
Fwd<T> cos( Fwd<T> x )
{
    return Fwd<T>( std::cos(x.val), -std::sin(x.val)*x.dot );
}

template <typename T>
Fwd<T> tan( Fwd<T> x )
{
    T tmp = std::cos(x.val);
    return Fwd<T>( std::tan(x.val), x.dot / (tmp*tmp) );
}

template <typename T>
Fwd<T> atan( Fwd<T> x )
{
    return Fwd<T>( std::atan(x.val), x.dot / (T(1.0) + (x.val*x.val)) );
}

template <typename T>
Fwd<T> exp( Fwd<T> x )
{
    T tmp = std::exp(x.val);
    return Fwd<T>( tmp, tmp * x.dot );
}

template <typename T>
Fwd<T> log( Fwd<T> x )
{
    T tmp = std::log(x.val);
    T dot = tmp != tmp ? tmp : (x.dot / x.val);
    return Fwd<T>( tmp, dot );
}

template <typename T>
Fwd<T> logabs( Fwd<T> x )
{
    return Fwd<T>( std::log(std::abs(x.val)), x.dot / x.val );
}

template <typename T>
Fwd<T> sinh( Fwd<T> x )
{
    return Fwd<T>( std::sinh(x.val), x.dot * std::cosh(x.val) );
}

template <typename T>
Fwd<T> cosh( Fwd<T> x )
{
    return Fwd<T>( std::cosh(x.val), x.dot * std::sinh(x.val) );
}

template <typename T>
Fwd<T> tanh( Fwd<T> x )
{
    T tmp = std::tanh(x.val);
    return Fwd<T>( tmp, x.dot * (T(1.0) - (tmp * tmp)) );
}

template <typename T>
Fwd<T> atanh( Fwd<T> x )
{
    return Fwd<T>( std::atanh(x.val), T(1.0) / (T(1.0) - (x.val*x.val)) );
}

} // end namespace std



// Gradient from arbitrary function / functor

template <typename F>
void FiniteDifferences( f64 *g, f64 *x, size_t n, f64 eps, F f )
{
    f64 f0 = f( x, n );

    for ( u32 i = 0; i < n; ++i ) {
        x[i] += eps;
        g[i] = ( f( x, n ) - f0 ) / eps;
        x[i] -= eps;
    }
}



#if TEST

using namespace std;

#if 0
bool Equal( f64 a, f64 b, f64 eps )
{
    if ( abs(a - b) > eps ) {
        return false;
    }
    else {
        return true;
    }
}
#endif


void test_fwd( void )
{
    Fwd<f64> x(1.0);

    f64 y = f64(x);

    ASSERT( x.val == y );

    Fwd<f64> z = x;

    ASSERT( x.val == z.val );
}

void test_fwd_add( void )
{
    Fwd<f64> x(1.0, 1.0);
    Fwd<f64> y(2.0);

    ASSERT( x < y && x <= y );

    Fwd<f64> z = x + y;

    ASSERT( z.val == x.val + y.val );
    ASSERT( z.dot == x.dot );

    Fwd<f64> zz = z + 4.0;

    ASSERT( zz.val == z.val + 4.0 );
    ASSERT( zz.dot == z.dot );

    Fwd<f64> zzz = 4.0 + z;

    ASSERT( zzz.val == z.val + 4.0 );
    ASSERT( zzz.dot == z.dot );
}

void test_fwd_sub( void )
{
    Fwd<f64> x(2.0, 1.0);
    Fwd<f64> y(1.0);

    Fwd<f64> z = x - y;

    ASSERT( z.val == x.val - y.val );
    ASSERT( z.dot == x.dot );

    Fwd<f64> zz = z - 4.0;

    ASSERT( zz.val == z.val - 4.0 );
    ASSERT( zz.dot == z.dot );

    Fwd<f64> zzz = 4.0 - z;

    ASSERT( zzz.val == 4.0 - z.val );
    ASSERT( zzz.dot == -z.dot );
}

void test_fwd_mul( void )
{
    Fwd<f64> x(2.0, 1.0);
    Fwd<f64> y(1.0);

    Fwd<f64> z = x * y;

    ASSERT( z.val == x.val * y.val );
    ASSERT( z.dot == y.val * x.dot );

    Fwd<f64> zz = z * 4.0;

    ASSERT( zz.val == z.val * 4.0 );
    ASSERT( zz.dot == z.dot * 4.0 );

    Fwd<f64> zzz = 4.0 * z;

    ASSERT( zzz.val == 4.0 * z.val );
    ASSERT( zzz.dot == 4.0 * z.dot );
}

void test_fwd_div( void )
{
    Fwd<f64> x(2.0, 1.0);
    Fwd<f64> y(3.0);

    ASSERT( y > x && y >= x );

    Fwd<f64> z = x / y;

    ASSERT( z.val == x.val / y.val );
    ASSERT( z.dot == x.dot / y.val );

    Fwd<f64> zz = z / 4.0;

    ASSERT( zz.val == z.val / 4.0 );
    ASSERT( zz.dot == z.dot / 4.0 );

    Fwd<f64> zzz = 4.0 / z;

    ASSERT( zzz.val == 4.0 / z.val );
    ASSERT( zzz.dot == - (4.0*z.dot) / (z.val*z.val) );
}

void test_incr_decr( void )
{
    Fwd<f64> x(2.0, 1.0);

    x += 3.0;

    ASSERT( x.val == 5.0 );
    ASSERT( x.dot == 1.0 );

    x -= 3.0;

    ASSERT( x.val == 2.0 );
    ASSERT( x.dot == 1.0 );
}

void test_sqrt( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = sqrt(x);

    ASSERT( z.val == sqrt(x.val) );
    ASSERT( z.dot == 1.0 / (2.0 * sqrt(x.val)) );
}

void test_pow( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = pow(x, 2.0);

    ASSERT( z.val == pow(x.val, 2.0) );
    ASSERT( z.dot == 2.0 * x.val );
}

void test_sin( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = sin(x);

    ASSERT( z.val == sin(x.val) );
    ASSERT( z.dot == cos(x.val) );
}

void test_cos( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = cos(x);

    ASSERT( z.val == cos(x.val) );
    ASSERT( z.dot == -sin(x.val) );
}

void test_tan( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z  = tan(x);
    Fwd<f64> zz = sin(x) / cos(x);

    ASSERT( Equal(z.val, zz.val, 1e-10) );
    ASSERT( Equal(z.dot, zz.dot, 1e-10) );
}

void test_atan( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = atan(x);

    ASSERT( Equal( z.val, atan(x.val), 1e-10) );
    ASSERT( Equal( z.dot, 1.0 / (1.0 + (x.val*x.val)), 1e-10) );
}

void test_exp( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = exp(x);

    ASSERT( Equal(z.val, exp(x.val), 1e-10) );
    ASSERT( Equal(z.dot, exp(x.val), 1e-10) );
}

void test_log( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = log(x);

    ASSERT( Equal(z.val, log(x.val), 1e-10) );
    ASSERT( Equal(z.dot, 1.0 / x.val, 1e-10) );
}

void test_logabs( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = logabs(-x);

    ASSERT( Equal(z.val, log(abs(x.val)), 1e-10) );
    ASSERT( Equal(z.dot, 1.0 / x.val, 1e-10) );
}

void test_sinh( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = sinh(x);

    ASSERT( Equal(z.val, sinh(x.val), 1e-10) );
    ASSERT( Equal(z.dot, cosh(x.val), 1e-10) );
}

void test_cosh( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = cosh(x);

    ASSERT( Equal(z.val, cosh(x.val), 1e-10) );
    ASSERT( Equal(z.dot, sinh(x.val), 1e-10) );
}

void test_tanh( void )
{
    Fwd<f64> x(2.0, 1.0);

    Fwd<f64> z = tanh(x);

    ASSERT( Equal(z.val, tanh(x.val), 1e-10) );
    ASSERT( Equal(z.dot, 1.0 - tanh(x.val)*tanh(x.val), 1e-10) );
}

void test_atanh( void )
{
    Fwd<f64> x(0.5, 1.0);

    Fwd<f64> z = atanh(x);

    ASSERT( Equal(z.val, atanh(x.val), 1e-10) );
    ASSERT( Equal(z.dot, 1.0 / (1.0 - x.val*x.val), 1e-10) );
}

#endif











