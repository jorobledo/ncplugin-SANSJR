#define HAS_Iq
#define FORM_VOL
#line 1 "./kernel_header.c"
#ifdef __OPENCL_VERSION__
# define USE_OPENCL
#elif defined(__CUDACC__)
# define USE_CUDA
#elif defined(_OPENMP)
# define USE_OPENMP
#endif

// Use SAS_DOUBLE to force the use of double even for float kernels
#define SAS_DOUBLE dou ## ble

// If opencl is not available, then we are compiling a C function
// Note: if using a C++ compiler, then define kernel as extern "C"
#ifdef USE_OPENCL

   #define USE_GPU
   #define pglobal global
   #define pconstant constant

   typedef int int32_t;

   #if defined(USE_SINCOS)
   #  define SINCOS(angle,svar,cvar) svar=sincos(angle,&cvar)
   #else
   #  define SINCOS(angle,svar,cvar) do {const double _t_=angle; svar=sin(_t_);cvar=cos(_t_);} while (0)
   #endif
   // Intel CPU on Mac gives strange values for erf(); on the verified
   // platforms (intel, nvidia, amd), the cephes erf() is significantly
   // faster than that available in the native OpenCL.
   #define NEED_ERF
   // OpenCL only has type generic math
   #define expf exp
   #ifndef NEED_ERF
   #  define erff erf
   #  define erfcf erfc
   #endif

#elif defined(USE_CUDA)

   #define USE_GPU
   #define local __shared__
   #define pglobal
   #define constant __constant__
   #define pconstant const
   #define kernel extern "C" __global__

   // OpenCL powr(a,b) = C99 pow(a,b), b >= 0
   // OpenCL pown(a,b) = C99 pow(a,b), b integer
   #define powr(a,b) pow(a,b)
   #define pown(a,b) pow(a,b)
   //typedef int int32_t;
   #if defined(USE_SINCOS)
   #  define SINCOS(angle,svar,cvar) sincos(angle,&svar,&cvar)
   #else
   #  define SINCOS(angle,svar,cvar) do {const double _t_=angle; svar=sin(_t_);cvar=cos(_t_);} while (0)
   #endif

#else // !USE_OPENCL && !USE_CUDA

   #define local
   #define pglobal
   #define constant const
   #define pconstant const

   #ifdef __cplusplus
      #include <cstdio>
      #include <cmath>
      using namespace std;
      #if defined(_MSC_VER)
         #include <limits>
         #include <float.h>
         #define kernel extern "C" __declspec( dllexport )
         inline double trunc(double x) { return x>=0?floor(x):-floor(-x); }
         inline double fmin(double x, double y) { return x>y ? y : x; }
         inline double fmax(double x, double y) { return x<y ? y : x; }
         #define isnan(x) _isnan(x)
         #define isinf(x) (!_finite(x))
         #define isfinite(x) _finite(x)
         #define NAN (std::numeric_limits<double>::quiet_NaN()) // non-signalling NaN
         #define INFINITY (std::numeric_limits<double>::infinity())
         #define NEED_ERF
         #define NEED_EXPM1
         #define NEED_TGAMMA
     #else
         #define kernel extern "C"
         #include <cstdint>
     #endif
     inline void SINCOS(double angle, double &svar, double &cvar) { svar=sin(angle); cvar=cos(angle); }
   #else // !__cplusplus
     #include <inttypes.h>  // C99 guarantees that int32_t types is here
     #include <stdio.h>
     #if defined(__TINYC__)
         typedef int int32_t;
         #include <math.h>
         // TODO: check isnan is correct
         inline double _isnan(double x) { return x != x; } // hope this doesn't optimize away!
         #undef isnan
         #define isnan(x) _isnan(x)
         // Defeat the double->float conversion since we don't have tgmath
         inline SAS_DOUBLE trunc(SAS_DOUBLE x) { return x>=0?floor(x):-floor(-x); }
         inline SAS_DOUBLE fmin(SAS_DOUBLE x, SAS_DOUBLE y) { return x>y ? y : x; }
         inline SAS_DOUBLE fmax(SAS_DOUBLE x, SAS_DOUBLE y) { return x<y ? y : x; }
         #define NEED_ERF
         #define NEED_EXPM1
         #define NEED_TGAMMA
         #define NEED_CBRT
         // expf missing from windows?
         #define expf exp
     #else
         #include <tgmath.h> // C99 type-generic math, so sin(float) => sinf
     #endif
     // MSVC doesn't support C99, so no need for dllexport on C99 branch
     #define kernel
     #define SINCOS(angle,svar,cvar) do {const double _t_=angle; svar=sin(_t_);cvar=cos(_t_);} while (0)
   #endif  // !__cplusplus
   // OpenCL powr(a,b) = C99 pow(a,b), b >= 0
   // OpenCL pown(a,b) = C99 pow(a,b), b integer
   #define powr(a,b) pow(a,b)
   #define pown(a,b) pow(a,b)

#endif // !USE_OPENCL

#if defined(NEED_CBRT)
   #define cbrt(_x) pow(_x, 0.33333333333333333333333)
#endif

#if defined(NEED_EXPM1)
   // TODO: precision is a half digit lower than numpy on mac in [1e-7, 0.5]
   // Run "explore/precision.py sas_expm1" to see this (may have to fiddle
   // the xrange for log to see the complete range).
   static SAS_DOUBLE expm1(SAS_DOUBLE x_in) {
      double x = (double)x_in;  // go back to float for single precision kernels
      // Adapted from the cephes math library.
      // Copyright 1984 - 1992 by Stephen L. Moshier
      if (x != x || x == 0.0) {
         return x; // NaN and +/- 0
      } else if (x < -0.5 || x > 0.5) {
         return exp(x) - 1.0;
      } else {
         const double xsq = x*x;
         const double p = (((
            +1.2617719307481059087798E-4)*xsq
            +3.0299440770744196129956E-2)*xsq
            +9.9999999999999999991025E-1);
         const double q = ((((
            +3.0019850513866445504159E-6)*xsq
            +2.5244834034968410419224E-3)*xsq
            +2.2726554820815502876593E-1)*xsq
            +2.0000000000000000000897E0);
         double r = x * p;
         r =  r / (q - r);
         return r+r;
       }
   }
#endif

// Standard mathematical constants:
//   M_E, M_LOG2E, M_LOG10E, M_LN2, M_LN10, M_PI, M_PI_2=pi/2, M_PI_4=pi/4,
//   M_1_PI=1/pi, M_2_PI=2/pi, M_2_SQRTPI=2/sqrt(pi), SQRT2, SQRT1_2=sqrt(1/2)
// OpenCL defines M_constant_F for float constants, and nothing if double
// is not enabled on the card, which is why these constants may be missing
#ifndef M_PI
#  define M_PI 3.141592653589793
#endif
#ifndef M_PI_2
#  define M_PI_2 1.570796326794897
#endif
#ifndef M_PI_4
#  define M_PI_4 0.7853981633974483
#endif
#ifndef M_E
#  define M_E 2.718281828459045091
#endif
#ifndef M_SQRT1_2
#  define M_SQRT1_2 0.70710678118654746
#endif

// Non-standard function library
// pi/180, used for converting between degrees and radians
// 4/3 pi for computing sphere volumes
// square and cube for computing squares and cubes
#ifndef M_PI_180
#  define M_PI_180 0.017453292519943295
#endif
#ifndef M_4PI_3
#  define M_4PI_3 4.18879020478639
#endif
double square(double x) { return x*x; }
double cube(double x) { return x*x*x; }
double sas_sinx_x(double x) { return x==0 ? 1.0 : sin(x)/x; }

// CRUFT: support old style models with orientation received qx, qy and angles

// To rotate from the canonical position to theta, phi, psi, first rotate by
// psi about the major axis, oriented along z, which is a rotation in the
// detector plane xy. Next rotate by theta about the y axis, aligning the major
// axis in the xz plane. Finally, rotate by phi in the detector plane xy.
// To compute the scattering, undo these rotations in reverse order:
//     rotate in xy by -phi, rotate in xz by -theta, rotate in xy by -psi
// The returned q is the length of the q vector and (xhat, yhat, zhat) is a unit
// vector in the q direction.
// To change between counterclockwise and clockwise rotation, change the
// sign of phi and psi.

#if 1
//think cos(theta) should be sin(theta) in new coords, RKH 11Jan2017
#define ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sn, cn) do { \
    SINCOS(phi*M_PI_180, sn, cn); \
    q = sqrt(qx*qx + qy*qy); \
    cn  = (q==0. ? 1.0 : (cn*qx + sn*qy)/q * sin(theta*M_PI_180));  \
    sn = sqrt(1 - cn*cn); \
    } while (0)
#else
// SasView 3.x definition of orientation
#define ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sn, cn) do { \
    SINCOS(theta*M_PI_180, sn, cn); \
    q = sqrt(qx*qx + qy*qy);\
    cn = (q==0. ? 1.0 : (cn*cos(phi*M_PI_180)*qx + sn*qy)/q); \
    sn = sqrt(1 - cn*cn); \
    } while (0)
#endif

#if 1
#define ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, xhat, yhat, zhat) do { \
    q = sqrt(qx*qx + qy*qy); \
    const double qxhat = qx/q; \
    const double qyhat = qy/q; \
    double sin_theta, cos_theta; \
    double sin_phi, cos_phi; \
    double sin_psi, cos_psi; \
    SINCOS(theta*M_PI_180, sin_theta, cos_theta); \
    SINCOS(phi*M_PI_180, sin_phi, cos_phi); \
    SINCOS(psi*M_PI_180, sin_psi, cos_psi); \
    xhat = qxhat*(-sin_phi*sin_psi + cos_theta*cos_phi*cos_psi) \
         + qyhat*( cos_phi*sin_psi + cos_theta*sin_phi*cos_psi); \
    yhat = qxhat*(-sin_phi*cos_psi - cos_theta*cos_phi*sin_psi) \
         + qyhat*( cos_phi*cos_psi - cos_theta*sin_phi*sin_psi); \
    zhat = qxhat*(-sin_theta*cos_phi) \
         + qyhat*(-sin_theta*sin_phi); \
    } while (0)
#else
// SasView 3.x definition of orientation
#define ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, cos_alpha, cos_mu, cos_nu) do { \
    q = sqrt(qx*qx + qy*qy); \
    const double qxhat = qx/q; \
    const double qyhat = qy/q; \
    double sin_theta, cos_theta; \
    double sin_phi, cos_phi; \
    double sin_psi, cos_psi; \
    SINCOS(theta*M_PI_180, sin_theta, cos_theta); \
    SINCOS(phi*M_PI_180, sin_phi, cos_phi); \
    SINCOS(psi*M_PI_180, sin_psi, cos_psi); \
    cos_alpha = cos_theta*cos_phi*qxhat + sin_theta*qyhat; \
    cos_mu = (-sin_theta*cos_psi*cos_phi - sin_psi*sin_phi)*qxhat + cos_theta*cos_psi*qyhat; \
    cos_nu = (-cos_phi*sin_psi*sin_theta + sin_phi*cos_psi)*qxhat + sin_psi*cos_theta*qyhat; \
    } while (0)
#endif

//# Beginning of rotational operation definitions

typedef struct {
          double R31, R32;
      } QACRotation;

typedef struct {
    double R11, R12;
    double R21, R22;
    double R31, R32;
} QABCRotation;

// Fill in the rotation matrix R from the view angles (theta, phi) and the
// jitter angles (dtheta, dphi).  This matrix can be applied to all of the
// (qx, qy) points in the image to produce R*[qx,qy]' = [qa,qc]'
static void
qac_rotation(
    QACRotation *rotation,
    double theta, double phi,
    double dtheta, double dphi)
{
    double sin_theta, cos_theta;
    double sin_phi, cos_phi;

    // reverse view matrix
    SINCOS(theta*M_PI_180, sin_theta, cos_theta);
    SINCOS(phi*M_PI_180, sin_phi, cos_phi);
    const double V11 = cos_phi*cos_theta;
    const double V12 = sin_phi*cos_theta;
    const double V21 = -sin_phi;
    const double V22 = cos_phi;
    const double V31 = sin_theta*cos_phi;
    const double V32 = sin_phi*sin_theta;

    // reverse jitter matrix
    SINCOS(dtheta*M_PI_180, sin_theta, cos_theta);
    SINCOS(dphi*M_PI_180, sin_phi, cos_phi);
    const double J31 = sin_theta;
    const double J32 = -sin_phi*cos_theta;
    const double J33 = cos_phi*cos_theta;

    // reverse matrix
    rotation->R31 = J31*V11 + J32*V21 + J33*V31;
    rotation->R32 = J31*V12 + J32*V22 + J33*V32;
}

// Apply the rotation matrix returned from qac_rotation to the point (qx,qy),
// returning R*[qx,qy]' = [qa,qc]'
static void
qac_apply(
    QACRotation *rotation,
    double qx, double qy,
    double *qab_out, double *qc_out)
{
    // Indirect calculation of qab, from qab^2 = |q|^2 - qc^2
    const double dqc = rotation->R31*qx + rotation->R32*qy;
    const double dqab_sq = -dqc*dqc + qx*qx + qy*qy;
    //*qab_out = sqrt(fabs(dqab_sq));
    *qab_out = dqab_sq > 0.0 ? sqrt(dqab_sq) : 0.0;
    *qc_out = dqc;
}

// Fill in the rotation matrix R from the view angles (theta, phi, psi) and the
// jitter angles (dtheta, dphi, dpsi).  This matrix can be applied to all of the
// (qx, qy) points in the image to produce R*[qx,qy]' = [qa,qb,qc]'
static void
qabc_rotation(
    QABCRotation *rotation,
    double theta, double phi, double psi,
    double dtheta, double dphi, double dpsi)
{
    double sin_theta, cos_theta;
    double sin_phi, cos_phi;
    double sin_psi, cos_psi;

    // reverse view matrix
    SINCOS(theta*M_PI_180, sin_theta, cos_theta);
    SINCOS(phi*M_PI_180, sin_phi, cos_phi);
    SINCOS(psi*M_PI_180, sin_psi, cos_psi);
    const double V11 = -sin_phi*sin_psi + cos_phi*cos_psi*cos_theta;
    const double V12 = sin_phi*cos_psi*cos_theta + sin_psi*cos_phi;
    const double V21 = -sin_phi*cos_psi - sin_psi*cos_phi*cos_theta;
    const double V22 = -sin_phi*sin_psi*cos_theta + cos_phi*cos_psi;
    const double V31 = sin_theta*cos_phi;
    const double V32 = sin_phi*sin_theta;

    // reverse jitter matrix
    SINCOS(dtheta*M_PI_180, sin_theta, cos_theta);
    SINCOS(dphi*M_PI_180, sin_phi, cos_phi);
    SINCOS(dpsi*M_PI_180, sin_psi, cos_psi);
    const double J11 = cos_psi*cos_theta;
    const double J12 = sin_phi*sin_theta*cos_psi + sin_psi*cos_phi;
    const double J13 = sin_phi*sin_psi - sin_theta*cos_phi*cos_psi;
    const double J21 = -sin_psi*cos_theta;
    const double J22 = -sin_phi*sin_psi*sin_theta + cos_phi*cos_psi;
    const double J23 = sin_phi*cos_psi + sin_psi*sin_theta*cos_phi;
    const double J31 = sin_theta;
    const double J32 = -sin_phi*cos_theta;
    const double J33 = cos_phi*cos_theta;

    // reverse matrix
    rotation->R11 = J11*V11 + J12*V21 + J13*V31;
    rotation->R12 = J11*V12 + J12*V22 + J13*V32;
    rotation->R21 = J21*V11 + J22*V21 + J23*V31;
    rotation->R22 = J21*V12 + J22*V22 + J23*V32;
    rotation->R31 = J31*V11 + J32*V21 + J33*V31;
    rotation->R32 = J31*V12 + J32*V22 + J33*V32;
}

// Apply the rotation matrix returned from qabc_rotation to the point (qx,qy),
// returning R*[qx,qy]' = [qa,qb,qc]'
static void
qabc_apply(
    QABCRotation *rotation,
    double qx, double qy,
    double *qa_out, double *qb_out, double *qc_out)
{
    *qa_out = rotation->R11*qx + rotation->R12*qy;
    *qb_out = rotation->R21*qx + rotation->R22*qy;
    *qc_out = rotation->R31*qx + rotation->R32*qy;
}

// ##### End of rotation operation definitions ######
#line 1 "./models/lib/polevl.c"
/*							polevl.c
 *							p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 * The function p1evl() assumes that C_N = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */


/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

static
double polevl( double x, pconstant double *coef, int N )
{

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }

    return ans;
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

static
double p1evl( double x, pconstant double *coef, int N )
{
    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return ans;
}

#line 1 "./models/lib/sas_J1.c"
/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

#if FLOAT_SIZE>4
//Cephes double pression function

constant double RPJ1[8] = {
    -8.99971225705559398224E8,
    4.52228297998194034323E11,
    -7.27494245221818276015E13,
    3.68295732863852883286E15,
    0.0,
    0.0,
    0.0,
    0.0 };

constant double RQJ1[8] = {
    6.20836478118054335476E2,
    2.56987256757748830383E5,
    8.35146791431949253037E7,
    2.21511595479792499675E10,
    4.74914122079991414898E12,
    7.84369607876235854894E14,
    8.95222336184627338078E16,
    5.32278620332680085395E18
    };

constant double PPJ1[8] = {
    7.62125616208173112003E-4,
    7.31397056940917570436E-2,
    1.12719608129684925192E0,
    5.11207951146807644818E0,
    8.42404590141772420927E0,
    5.21451598682361504063E0,
    1.00000000000000000254E0,
    0.0} ;


constant double PQJ1[8] = {
    5.71323128072548699714E-4,
    6.88455908754495404082E-2,
    1.10514232634061696926E0,
    5.07386386128601488557E0,
    8.39985554327604159757E0,
    5.20982848682361821619E0,
    9.99999999999999997461E-1,
    0.0 };

constant double QPJ1[8] = {
    5.10862594750176621635E-2,
    4.98213872951233449420E0,
    7.58238284132545283818E1,
    3.66779609360150777800E2,
    7.10856304998926107277E2,
    5.97489612400613639965E2,
    2.11688757100572135698E2,
    2.52070205858023719784E1 };

constant double QQJ1[8] = {
    7.42373277035675149943E1,
    1.05644886038262816351E3,
    4.98641058337653607651E3,
    9.56231892404756170795E3,
    7.99704160447350683650E3,
    2.82619278517639096600E3,
    3.36093607810698293419E2,
    0.0 };

static
double cephes_j1(double x)
{

    double w, z, p, q, abs_x, sign_x;

    const double Z1 = 1.46819706421238932572E1;
    const double Z2 = 4.92184563216946036703E1;

    // 2017-05-18 PAK - mathematica and mpmath use J1(-x) = -J1(x)
    if (x < 0) {
        abs_x = -x;
        sign_x = -1.0;
    } else {
        abs_x = x;
        sign_x = 1.0;
    }

    if( abs_x <= 5.0 ) {
        z = abs_x * abs_x;
        w = polevl( z, RPJ1, 3 ) / p1evl( z, RQJ1, 8 );
        w = w * abs_x * (z - Z1) * (z - Z2);
        return( sign_x * w );
    }

    w = 5.0/abs_x;
    z = w * w;
    p = polevl( z, PPJ1, 6)/polevl( z, PQJ1, 6 );
    q = polevl( z, QPJ1, 7)/p1evl( z, QQJ1, 7 );

    // 2017-05-19 PAK improve accuracy using trig identies
    // original:
    //    const double THPIO4 =  2.35619449019234492885;
    //    const double SQ2OPI = 0.79788456080286535588;
    //    double sin_xn, cos_xn;
    //    SINCOS(abs_x - THPIO4, sin_xn, cos_xn);
    //    p = p * cos_xn - w * q * sin_xn;
    //    return( sign_x * p * SQ2OPI / sqrt(abs_x) );
    // expanding p*cos(a - 3 pi/4) - wq sin(a - 3 pi/4)
    //    [ p(sin(a) - cos(a)) + wq(sin(a) + cos(a)) / sqrt(2)
    // note that sqrt(1/2) * sqrt(2/pi) = sqrt(1/pi)
    const double SQRT1_PI = 0.56418958354775628;
    double sin_x, cos_x;
    SINCOS(abs_x, sin_x, cos_x);
    p = p*(sin_x - cos_x) + w*q*(sin_x + cos_x);
    return( sign_x * p * SQRT1_PI / sqrt(abs_x) );
}

#else
//Single precission version of cephes
constant float JPJ1[8] = {
    -4.878788132172128E-009,
    6.009061827883699E-007,
    -4.541343896997497E-005,
    1.937383947804541E-003,
    -3.405537384615824E-002,
    0.0,
    0.0,
    0.0
    };

constant float MO1J1[8] = {
    6.913942741265801E-002,
    -2.284801500053359E-001,
    3.138238455499697E-001,
    -2.102302420403875E-001,
    5.435364690523026E-003,
    1.493389585089498E-001,
    4.976029650847191E-006,
    7.978845453073848E-001
    };

constant float PH1J1[8] = {
    -4.497014141919556E+001,
    5.073465654089319E+001,
    -2.485774108720340E+001,
    7.222973196770240E+000,
    -1.544842782180211E+000,
    3.503787691653334E-001,
    -1.637986776941202E-001,
    3.749989509080821E-001
    };

static
float cephes_j1f(float xx)
{

    float x, w, z, p, q, xn;

    const float Z1 = 1.46819706421238932572E1;


    // 2017-05-18 PAK - mathematica and mpmath use J1(-x) = -J1(x)
    x = xx;
    if( x < 0 )
        x = -xx;

    if( x <= 2.0 ) {
        z = x * x;
        p = (z-Z1) * x * polevl( z, JPJ1, 4 );
        return( xx < 0. ? -p : p );
    }

    q = 1.0/x;
    w = sqrt(q);

    p = w * polevl( q, MO1J1, 7);
    w = q*q;
    // 2017-05-19 PAK improve accuracy using trig identies
    // original:
    //    const float THPIO4F =  2.35619449019234492885;    /* 3*pi/4 */
    //    xn = q * polevl( w, PH1J1, 7) - THPIO4F;
    //    p = p * cos(xn + x);
    //    return( xx < 0. ? -p : p );
    // expanding cos(a + b - 3 pi/4) is
    //    [sin(a)sin(b) + sin(a)cos(b) + cos(a)sin(b)-cos(a)cos(b)] / sqrt(2)
    xn = q * polevl( w, PH1J1, 7);
    float cos_xn, sin_xn;
    float cos_x, sin_x;
    SINCOS(xn, sin_xn, cos_xn);  // about xn and 1
    SINCOS(x, sin_x, cos_x);
    p *= M_SQRT1_2*(sin_xn*(sin_x+cos_x) + cos_xn*(sin_x-cos_x));

    return( xx < 0. ? -p : p );
}
#endif

#if FLOAT_SIZE>4
#define sas_J1 cephes_j1
#else
#define sas_J1 cephes_j1f
#endif

//Finally J1c function that equals 2*J1(x)/x
static
double sas_2J1x_x(double x)
{
    return (x != 0.0 ) ? 2.0*sas_J1(x)/x : 1.0;
}

#line 1 "./models/lib/wrc_cyl.c"
/*
    Functions for WRC implementation of flexible cylinders. See
    W R Chen, P D Butler and L J Magid,
    Incorporating Intermicellar Interactions in the Fitting of
    SANS Data from Cationic Wormlike Micelles.
    Langmuir, 22(15) 2006 6539-6548
*/

static double
Rgsquare(double L, double b)
{
    const double x = L/b;
    // Use Horner's method to evaluate Pedersen eq 15:
    //     alpha^2 = [1.0 + (x/3.12)^2 + (x/8.67)^3] ^ (0.176/3)
    const double alphasq =
        pow(1.0 + x*x*(1.027284681130835e-01 + 1.534414548417740e-03*x),
            5.866666666666667e-02);
    return alphasq*L*b/6.0;
}

static double
Rgsquareshort(double L, double b)
{
    const double r = b/L;  // = 1/n_b in Pedersen ref.
    return Rgsquare(L, b) * (1.0 + r*(-1.5 + r*(1.5 + r*0.75*expm1(-2.0/r))));
}

static double
w_WR(double x)
{
    // Pedersen eq. 16:
    //    w = [1 + tanh((x-C4)/C5)]/2
    const double C4 = 1.523;
    const double C5 = 0.1477;
    return 0.5 + 0.5*tanh((x - C4)/C5);
}

static double
Sdebye(double qsq)
{
#if FLOAT_SIZE>4
#  define DEBYE_CUTOFF 0.25  // 1e-15 error
#else
#  define DEBYE_CUTOFF 1.0  // 4e-7 error
#endif

    if (qsq < DEBYE_CUTOFF) {
        const double x = qsq;
        // mathematica: PadeApproximant[2*Exp[-x^2] + x^2-1)/x^4, {x, 0, 8}]
        const double A1=1./15., A2=1./60, A3=0., A4=1./75600.;
        const double B1=2./5., B2=1./15., B3=1./180., B4=1./5040.;
        return ((((A4*x + A3)*x + A2)*x + A1)*x + 1.)
             / ((((B4*x + B3)*x + B2)*x + B1)*x + 1.);
    } else {
        return 2.*(expm1(-qsq) + qsq)/(qsq*qsq);
    }
}

static double
a_long(double q, double L, double b/*, double p1, double p2, double q0*/)
{
    const double p1 = 4.12;
    const double p2 = 4.42;
    const double q0 = 3.1;

    // Constants C1, ..., C5 derived from least squares fit (Pedersen, eq 13,16)
    const double C1 = 1.22;
    const double C2 = 0.4288;
    const double C3 = -1.651;
    const double C4 = 1.523;
    const double C5 = 0.1477;
    const double miu = 0.585;

    const double C = (L/b>10.0 ? 3.06*pow(L/b, -0.44) : 1.0);
    //printf("branch B-%d q=%g L=%g b=%g\n", C==1.0, q, L, b);
    const double r2 = Rgsquare(L,b);
    const double r = sqrt(r2);
    const double qr_b = q0*r/b;
    const double qr_b_sq = qr_b*qr_b;
    const double qr_b_4 = qr_b_sq*qr_b_sq;
    const double qr_b_miu = pow(qr_b, -1.0/miu);
    const double em1_qr_b_sq = expm1(-qr_b_sq);
    const double sech2 = 1.0/square(cosh((qr_b-C4)/C5));
    const double w = w_WR(qr_b);

    const double t1 = pow(q0, 1.0 + p1 + p2)/(b*(p1-p2));
    const double t2 = C/(15.0*L) * (
        + 14.0*b*b*em1_qr_b_sq/(q0*qr_b_sq)
        + 2.0*q0*r2*exp(-qr_b_sq)*(11.0 + 7.0/qr_b_sq));
    const double t11 = ((C3*qr_b_miu + C2)*qr_b_miu + C1)*qr_b_miu;
    const double t3 = r*sech2/(2.*C5)*t11;
    const double t4 = r*(em1_qr_b_sq + qr_b_sq)*sech2 / (C5*qr_b_4);
    const double t5 = -4.0*r*qr_b*em1_qr_b_sq/qr_b_4 * (1.0 - w);
    const double t10 = 2.0*(em1_qr_b_sq + qr_b_sq)/qr_b_4 * (1.0 - w); //=Sdebye*(1-w)
    const double t6 = 4.0*b/q0 * t10;
    const double t7 = r*((-3.0*C3*qr_b_miu -2.0*C2)*qr_b_miu -1.0*C1)*qr_b_miu/(miu*qr_b);
    const double t9 = C*b/L * (4.0 - exp(-qr_b_sq) * (11.0 + 7.0/qr_b_sq) + 7.0/qr_b_sq)/15.0;
    const double t12 = b*b*M_PI/(L*q0*q0) + t2 + t3 - t4 + t5 - t6 + t7*w;
    const double t13 = -b*M_PI/(L*q0) + t9 + t10 + t11*w;

    const double a1 = pow(q0,p1)*t13 - t1*pow(q0,-p2)*(t12 + b*p1/q0*t13);
    const double a2 = t1*pow(q0,-p1)*(t12 + b*p1/q0*t13);

    const double ans = a1*pow(q*b, -p1) + a2*pow(q*b, -p2) + M_PI/(q*L);
    return ans;
}

static double
_short(double r2, double exp_qr_b, double L, double b, double p1short,
        double p2short, double q0)
{
    const double qr2 = q0*q0 * r2;
    const double b3 = b*b*b;
    const double q0p = pow(q0, -4.0 + p1short);

    double yy = 1.0/(L*r2*r2) * b/exp_qr_b*q0p
        * (8.0*b3*L
           - 8.0*b3*exp_qr_b*L
           + 2.0*b3*exp_qr_b*L*p2short
           - 2.0*b*exp_qr_b*L*p2short*qr2
           + 4.0*b*exp_qr_b*L*qr2
           - 2.0*b3*L*p2short
           + 4.0*b*L*qr2
           - M_PI*exp_qr_b*qr2*q0*r2
           + M_PI*exp_qr_b*p2short*qr2*q0*r2);

    return yy;
}
static double
a_short(double qp, double L, double b
        /*double p1short, double p2short*/, double q0)
{
    const double p1short = 5.36;
    const double p2short = 5.62;

    const double r2 = Rgsquareshort(L,b);
    const double exp_qr_b = exp(r2*square(q0/b));
    const double pdiff = p1short - p2short;
    const double a1 = _short(r2,exp_qr_b,L,b,p1short,p2short,q0)/pdiff;
    const double a2= -_short(r2,exp_qr_b,L,b,p2short,p1short,q0)/pdiff;
    const double ans = a1*pow(qp*b, -p1short) + a2*pow(qp*b, -p2short) + M_PI/(qp*L);
    return ans;
}


static double
Sexv(double q, double L, double b)
{
    // Pedersen eq 13, corrected by Chen eq A.5, swapping w and 1-w
    const double C1=1.22;
    const double C2=0.4288;
    const double C3=-1.651;
    const double miu = 0.585;
    const double qr = q*sqrt(Rgsquare(L,b));
    const double qr_miu = pow(qr, -1.0/miu);
    const double w = w_WR(qr);
    const double t10 = Sdebye(qr*qr)*(1.0 - w);
    const double t11 = ((C3*qr_miu + C2)*qr_miu + C1)*qr_miu;

    return t10 + w*t11;
}

// Modified by Yun on Oct. 15,
static double
Sexv_new(double q, double L, double b)
{
    const double qr = q*sqrt(Rgsquare(L,b));
    const double qr2 = qr*qr;
    const double C = (L/b > 10.0) ? 3.06*pow(L/b, -0.44) : 1.0;
    const double t9 = C*b/L * (4.0 - exp(-qr2) * (11.0 + 7.0/qr2) + 7.0/qr2)/15.0;

    const double Sexv_orig = Sexv(q, L, b);

    // calculating the derivative to decide on the correction (cutoff) term?
    // Note: this is modified from WRs original code
    const double del=1.05;
    const double qdel = (Sexv(q*del,L,b) - Sexv_orig)/(q*(del - 1.0));

    if (qdel < 0) {
        //printf("branch A1-%d q=%g L=%g b=%g\n", C==1.0, q, L, b);
        return t9 + Sexv_orig;
    } else {
        //printf("branch A2-%d q=%g L=%g b=%g\n", C==1.0, q, L, b);
        const double w = w_WR(qr);
        const double t10 = Sdebye(qr*qr)*(1.0 - w);
        return t9 + t10;
    }
}


static double
Sk_WR(double q, double L, double b)
{
    const double Rg_short = sqrt(Rgsquareshort(L, b));
    double q0short = fmax(1.9/Rg_short, 3.0);
    double ans;


    if( L > 4*b ) { // L > 4*b : Longer Chains
        if (q*b <= 3.1) {
            ans = Sexv_new(q, L, b);
        } else { //q(i)*b > 3.1
            ans = a_long(q, L, b /*, p1, p2, q0*/);
        }
    } else { // L <= 4*b : Shorter Chains
        if (q*b <= q0short) { // q*b <= fmax(1.9/Rg_short, 3)
            //printf("branch C-%d q=%g L=%g b=%g\n", square(q*Rg_short)<DEBYE_CUTOFF, q, L, b);
            // Note that q0short is usually 3, but it will be greater than 3
            // small enough b, depending on the L/b ratio:
            //     L/b == 1 => b < 2.37
            //     L/b == 2 => b < 1.36
            //     L/b == 3 => b < 1.00
            //     L/b == 4 => b < 0.816
            // 2017-10-01 pkienzle: moved low q approximation into Sdebye()
            ans = Sdebye(square(q*Rg_short));
        } else {  // q*b > max(1.9/Rg_short, 3)
            //printf("branch D q=%g L=%g b=%g\n", q, L, b);
            ans = a_short(q, L, b /*, p1short, p2short*/, q0short);
        }
    }

    return ans;
}

#line 1 "./models/flexible_cylinder.c"
static double
form_volume(double length, double kuhn_length, double radius)
{
    return 1.0;
}

static double
Iq(double q,
   double length,
   double kuhn_length,
   double radius,
   double sld,
   double solvent_sld)
{
    const double contrast = sld - solvent_sld;
    const double cross_section = sas_2J1x_x(q*radius);
    const double volume = M_PI*radius*radius*length;
    const double flex = Sk_WR(q, length, kuhn_length);
    return 1.0e-4 * volume * square(contrast*cross_section) * flex;
}
