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

#line 1 "./models/lib/sas_J0.c"
/*							j0.c
 *
 *	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30       10000       4.4e-17     6.3e-18
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

#if FLOAT_SIZE>4
//Cephes double precission
double cephes_j0(double x);

 constant double PPJ0[8] = {
        7.96936729297347051624E-4,
        8.28352392107440799803E-2,
        1.23953371646414299388E0,
        5.44725003058768775090E0,
        8.74716500199817011941E0,
        5.30324038235394892183E0,
        9.99999999999999997821E-1,
        0.0
    };

 constant double PQJ0[8] = {
        9.24408810558863637013E-4,
        8.56288474354474431428E-2,
        1.25352743901058953537E0,
        5.47097740330417105182E0,
        8.76190883237069594232E0,
        5.30605288235394617618E0,
        1.00000000000000000218E0,
        0.0
    };

 constant double QPJ0[8] = {
        -1.13663838898469149931E-2,
        -1.28252718670509318512E0,
        -1.95539544257735972385E1,
        -9.32060152123768231369E1,
        -1.77681167980488050595E2,
        -1.47077505154951170175E2,
        -5.14105326766599330220E1,
        -6.05014350600728481186E0,
    };

 constant double QQJ0[8] = {
        /*  1.00000000000000000000E0,*/
        6.43178256118178023184E1,
        8.56430025976980587198E2,
        3.88240183605401609683E3,
        7.24046774195652478189E3,
        5.93072701187316984827E3,
        2.06209331660327847417E3,
        2.42005740240291393179E2,
    };

 constant double YPJ0[8] = {
        1.55924367855235737965E4,
        -1.46639295903971606143E7,
        5.43526477051876500413E9,
        -9.82136065717911466409E11,
        8.75906394395366999549E13,
        -3.46628303384729719441E15,
        4.42733268572569800351E16,
        -1.84950800436986690637E16,
 };


 constant double YQJ0[7] = {
        /* 1.00000000000000000000E0,*/
        1.04128353664259848412E3,
        6.26107330137134956842E5,
        2.68919633393814121987E8,
        8.64002487103935000337E10,
        2.02979612750105546709E13,
        3.17157752842975028269E15,
        2.50596256172653059228E17,
  };

 constant double RPJ0[8] = {
        -4.79443220978201773821E9,
        1.95617491946556577543E12,
        -2.49248344360967716204E14,
        9.70862251047306323952E15,
        0.0,
        0.0,
        0.0,
        0.0
  };

 constant double RQJ0[8] = {
        /* 1.00000000000000000000E0,*/
        4.99563147152651017219E2,
        1.73785401676374683123E5,
        4.84409658339962045305E7,
        1.11855537045356834862E10,
        2.11277520115489217587E12,
        3.10518229857422583814E14,
        3.18121955943204943306E16,
        1.71086294081043136091E18,
  };

double cephes_j0(double x)
{
    double w, z, p, q, xn;

    //const double TWOOPI = 6.36619772367581343075535E-1;
    const double SQ2OPI = 7.9788456080286535587989E-1;
    const double PIO4 = 7.85398163397448309616E-1;

    const double DR1 = 5.78318596294678452118E0;
    const double DR2 = 3.04712623436620863991E1;


    if( x < 0 )
        x = -x;

    if( x <= 5.0 ) {
        z = x * x;
        if( x < 1.0e-5 )
            return( 1.0 - z/4.0 );

        p = (z - DR1) * (z - DR2);
        p = p * polevl( z, RPJ0, 3)/p1evl( z, RQJ0, 8 );
        return( p );
    }

    w = 5.0/x;
    q = 25.0/(x*x);
    p = polevl( q, PPJ0, 6)/polevl( q, PQJ0, 6 );
    q = polevl( q, QPJ0, 7)/p1evl( q, QQJ0, 7 );
    xn = x - PIO4;

    double sn, cn;
    SINCOS(xn, sn, cn);
    p = p * cn - w * q * sn;

    return( p * SQ2OPI / sqrt(x) );
}
#else
//Cephes single precission
float cephes_j0f(float x);

 constant float MOJ0[8] = {
        -6.838999669318810E-002,
        1.864949361379502E-001,
        -2.145007480346739E-001,
        1.197549369473540E-001,
        -3.560281861530129E-003,
        -4.969382655296620E-002,
        -3.355424622293709E-006,
        7.978845717621440E-001
  };

 constant float PHJ0[8] = {
        3.242077816988247E+001,
        -3.630592630518434E+001,
        1.756221482109099E+001,
        -4.974978466280903E+000,
        1.001973420681837E+000,
        -1.939906941791308E-001,
        6.490598792654666E-002,
        -1.249992184872738E-001
  };

 constant float JPJ0[8] = {
        -6.068350350393235E-008,
        6.388945720783375E-006,
        -3.969646342510940E-004,
        1.332913422519003E-002,
        -1.729150680240724E-001,
        0.0,
        0.0,
        0.0
 };

float cephes_j0f(float x)
{
    float xx, w, z, p, q, xn;

    //const double YZ1 =  0.43221455686510834878;
    //const double YZ2 = 22.401876406482861405;
    //const double YZ3 = 64.130620282338755553;
    const float DR1 =  5.78318596294678452118;
    const float PIO4F = 0.7853981633974483096;

    if( x < 0 )
        xx = -x;
    else
        xx = x;

    // 2017-05-18 PAK - support negative x
    if( xx <= 2.0 ) {
        z = xx * xx;
        if( xx < 1.0e-3 )
            return( 1.0 - 0.25*z );

        p = (z-DR1) * polevl( z, JPJ0, 4);
        return( p );
    }

    q = 1.0/xx;
    w = sqrt(q);

    p = w * polevl( q, MOJ0, 7);
    w = q*q;
    xn = q * polevl( w, PHJ0, 7) - PIO4F;
    p = p * cos(xn + xx);
    return(p);
}
#endif

#if FLOAT_SIZE>4
#define sas_J0 cephes_j0
#else
#define sas_J0 cephes_j0f
#endif

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

#line 1 "./models/lib/sas_JN.c"
/*							jn.c
 *
 *	Bessel function of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, jn();
 *
 * y = jn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The ratio of jn(x) to j0(x) is computed by backward
 * recurrence.  First the ratio jn/jn-1 is found by a
 * continued fraction expansion.  Then the recurrence
 * relating successive orders is applied until j0 or j1 is
 * reached.
 *
 * If n = 0 or 1 the routine for j0 or j1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   range      # trials      peak         rms
 *    DEC       0, 30        5500       6.9e-17     9.3e-18
 *    IEEE      0, 30        5000       4.4e-16     7.9e-17
 *
 *
 * Not suitable for large n or x. Use jv() instead.
 *
 */

/*							jn.c
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

#if FLOAT_SIZE > 4

double cephes_jn( int n, double x );
double cephes_jn( int n, double x ) {

    // PAK: seems to be machine epsilon/2
    const double MACHEP = 1.11022302462515654042E-16;
    double pkm2, pkm1, pk, xk, r, ans;
    int k, sign;

    if( n < 0 ) {
        n = -n;
        if( (n & 1) == 0 )	/* -1**n */
            sign = 1;
        else
            sign = -1;
    } else {
        sign = 1;
    }

    if( x < 0.0 ) {
        if( n & 1 )
            sign = -sign;
    x = -x;
    }

    if( n == 0 )
        return( sign * cephes_j0(x) );
    if( n == 1 )
        return( sign * cephes_j1(x) );
    if( n == 2 )
        return( sign * (2.0 * cephes_j1(x) / x  -  cephes_j0(x)) );

    if( x < MACHEP )
        return( 0.0 );

    k = 53;
    pk = 2 * (n + k);
    ans = pk;
    xk = x * x;

    do {
        pk -= 2.0;
        ans = pk - (xk/ans);
    } while( --k > 0 );

    /* backward recurrence */

    pk = 1.0;

    ans = x/ans;
    pkm1 = 1.0/ans;

    k = n-1;
    r = 2 * k;

    do {
        pkm2 = (pkm1 * r  -  pk * x) / x;
        pk = pkm1;
        pkm1 = pkm2;
        r -= 2.0;
    } while( --k > 0 );

    if( fabs(pk) > fabs(pkm1) )
        ans = cephes_j1(x)/pk;
    else
        ans = cephes_j0(x)/pkm1;

    return( sign * ans );
}

#else

float cephes_jnf(int n, float x);
float cephes_jnf(int n, float x)
{
    // PAK: seems to be machine epsilon/2
    const double MACHEP = 5.9604645e-08;
    float pkm2, pkm1, pk, xk, r, ans;
    int k, sign;

    if( n < 0 ) {
        n = -n;
        if( (n & 1) == 0 ) /* -1**n */
            sign = 1;
        else
            sign = -1;
    } else {
        sign = 1;
    }

    if( x < 0.0 ) {
        if( n & 1 )
            sign = -sign;
        x = -x;
    }

    if( n == 0 )
        return( sign * cephes_j0f(x) );
    if( n == 1 )
        return( sign * cephes_j1f(x) );
    if( n == 2 )
        return( sign * (2.0 * cephes_j1f(x) / x  -  cephes_j0f(x)) );

    if( x < MACHEP )
        return( 0.0 );

    k = 24;
    pk = 2 * (n + k);
    ans = pk;
    xk = x * x;

    do {
        pk -= 2.0;
        ans = pk - (xk/ans);
    } while( --k > 0 );

    /* backward recurrence */

    pk = 1.0;

    const float xinv = 1.0/x;
    pkm1 = ans * xinv;
    k = n-1;
    r = (float )(2 * k);

    do {
        pkm2 = (pkm1 * r  -  pk * x) * xinv;
        pk = pkm1;
        pkm1 = pkm2;
        r -= 2.0;
    } while( --k > 0 );

    r = pk;
    if( r < 0 )
        r = -r;
    ans = pkm1;
    if( ans < 0 )
        ans = -ans;

    if( r > ans )  /* if( fabs(pk) > fabs(pkm1) ) */
        ans = sign * cephes_j1f(x)/pk;
    else
        ans = sign * cephes_j0f(x)/pkm1;
    return( ans );
}
#endif

#if FLOAT_SIZE>4
#define sas_JN cephes_jn
#else
#define sas_JN cephes_jnf
#endif

#line 1 "./models/lib/gauss76.c"
// Created by Andrew Jackson on 4/23/07

 #ifdef GAUSS_N
 # undef GAUSS_N
 # undef GAUSS_Z
 # undef GAUSS_W
 #endif
 #define GAUSS_N 76
 #define GAUSS_Z Gauss76Z
 #define GAUSS_W Gauss76Wt

// Gaussians
constant double Gauss76Wt[76] = {
	.00126779163408536,		//0
	.00294910295364247,
	.00462793522803742,
	.00629918049732845,
	.00795984747723973,
	.00960710541471375,
	.0112381685696677,
	.0128502838475101,
	.0144407317482767,
	.0160068299122486,
	.0175459372914742,		//10
	.0190554584671906,
	.020532847967908,
	.0219756145344162,
	.0233813253070112,
	.0247476099206597,
	.026072164497986,
	.0273527555318275,
	.028587223650054,
	.029773487255905,
	.0309095460374916,		//20
	.0319934843404216,
	.0330234743977917,
	.0339977794120564,
	.0349147564835508,
	.0357728593807139,
	.0365706411473296,
	.0373067565423816,
	.0379799643084053,
	.0385891292645067,
	.0391332242205184,		//30
	.0396113317090621,
	.0400226455325968,
	.040366472122844,
	.0406422317102947,
	.0408494593018285,
	.040987805464794,
	.0410570369162294,
	.0410570369162294,
	.040987805464794,
	.0408494593018285,		//40
	.0406422317102947,
	.040366472122844,
	.0400226455325968,
	.0396113317090621,
	.0391332242205184,
	.0385891292645067,
	.0379799643084053,
	.0373067565423816,
	.0365706411473296,
	.0357728593807139,		//50
	.0349147564835508,
	.0339977794120564,
	.0330234743977917,
	.0319934843404216,
	.0309095460374916,
	.029773487255905,
	.028587223650054,
	.0273527555318275,
	.026072164497986,
	.0247476099206597,		//60
	.0233813253070112,
	.0219756145344162,
	.020532847967908,
	.0190554584671906,
	.0175459372914742,
	.0160068299122486,
	.0144407317482767,
	.0128502838475101,
	.0112381685696677,
	.00960710541471375,		//70
	.00795984747723973,
	.00629918049732845,
	.00462793522803742,
	.00294910295364247,
	.00126779163408536		//75 (indexed from 0)
};

constant double Gauss76Z[76] = {
	-.999505948362153,		//0
	-.997397786355355,
	-.993608772723527,
	-.988144453359837,
	-.981013938975656,
	-.972229228520377,
	-.961805126758768,
	-.949759207710896,
	-.936111781934811,
	-.92088586125215,
	-.904107119545567,		//10
	-.885803849292083,
	-.866006913771982,
	-.844749694983342,
	-.822068037328975,
	-.7980001871612,
	-.77258672828181,
	-.74587051350361,
	-.717896592387704,
	-.688712135277641,
	-.658366353758143,		//20
	-.626910417672267,
	-.594397368836793,
	-.560882031601237,
	-.526420920401243,
	-.491072144462194,
	-.454895309813726,
	-.417951418780327,
	-.380302767117504,
	-.342012838966962,
	-.303146199807908,		//30
	-.263768387584994,
	-.223945802196474,
	-.183745593528914,
	-.143235548227268,
	-.102483975391227,
	-.0615595913906112,
	-.0205314039939986,
	.0205314039939986,
	.0615595913906112,
	.102483975391227,			//40
	.143235548227268,
	.183745593528914,
	.223945802196474,
	.263768387584994,
	.303146199807908,
	.342012838966962,
	.380302767117504,
	.417951418780327,
	.454895309813726,
	.491072144462194,		//50
	.526420920401243,
	.560882031601237,
	.594397368836793,
	.626910417672267,
	.658366353758143,
	.688712135277641,
	.717896592387704,
	.74587051350361,
	.77258672828181,
	.7980001871612,	//60
	.822068037328975,
	.844749694983342,
	.866006913771982,
	.885803849292083,
	.904107119545567,
	.92088586125215,
	.936111781934811,
	.949759207710896,
	.961805126758768,
	.972229228520377,		//70
	.981013938975656,
	.988144453359837,
	.993608772723527,
	.997397786355355,
	.999505948362153		//75
};

#line 1 "./models/pringle.c"
double form_volume(double radius, double thickness, double alpha, double beta);

double Iq(double q,
          double radius,
          double thickness,
          double alpha,
          double beta,
          double sld,
          double sld_solvent);


static
void _integrate_bessel(
    double radius,
    double alpha,
    double beta,
    double q_sin_psi,
    double q_cos_psi,
    double n,
    double *Sn,
    double *Cn)
{
    // translate gauss point z in [-1,1] to a point in [0, radius]
    const double zm = 0.5*radius;
    const double zb = 0.5*radius;

    // evaluate at Gauss points
    double sumS = 0.0;		// initialize integral
    double sumC = 0.0;		// initialize integral
    double r;
    for (int i=0; i < GAUSS_N; i++) {
        r = GAUSS_Z[i]*zm + zb;

        const double qrs = r*q_sin_psi;
        const double qrrc = r*r*q_cos_psi;

        double y = GAUSS_W[i] * r * sas_JN(n, beta*qrrc) * sas_JN(2*n, qrs);
        double S, C;
        SINCOS(alpha*qrrc, S, C);
        sumS += y*S;
        sumC += y*C;
    }

    *Sn = zm*sumS / (radius*radius);
    *Cn = zm*sumC / (radius*radius);
}

static
double _sum_bessel_orders(
    double radius,
    double alpha,
    double beta,
    double q_sin_psi,
    double q_cos_psi)
{
    //calculate sum term from n = -3 to 3
    //Note 1:
    //    S_n(-x) = (-1)^S_n(x)
    //    => S_n^2(-x) = S_n^2(x),
    //    => sum_-k^k Sk = S_0^2 + 2*sum_1^kSk^2
    //Note 2:
    //    better precision to sum terms from smaller to larger
    //    though it doesn't seem to make a difference in this case.
    double Sn, Cn, sum;
    sum = 0.0;
    for (int n=3; n>0; n--) {
      _integrate_bessel(radius, alpha, beta, q_sin_psi, q_cos_psi, n, &Sn, &Cn);
      sum += 2.0*(Sn*Sn + Cn*Cn);
    }
    _integrate_bessel(radius, alpha, beta, q_sin_psi, q_cos_psi, 0, &Sn, &Cn);
    sum += Sn*Sn+ Cn*Cn;
    return sum;
}

static
double _integrate_psi(
    double q,
    double radius,
    double thickness,
    double alpha,
    double beta)
{
    // translate gauss point z in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4;

    double sum = 0.0;
    for (int i = 0; i < GAUSS_N; i++) {
        double psi = GAUSS_Z[i]*zm + zb;
        double sin_psi, cos_psi;
        SINCOS(psi, sin_psi, cos_psi);
        double bessel_term = _sum_bessel_orders(radius, alpha, beta, q*sin_psi, q*cos_psi);
        double sinc_term = square(sas_sinx_x(q * thickness * cos_psi / 2.0));
        double pringle_kernel = 4.0 * sin_psi * bessel_term * sinc_term;
        sum += GAUSS_W[i] * pringle_kernel;
    }

    return zm * sum;
}

double form_volume(double radius, double thickness, double alpha, double beta)
{
    return M_PI*radius*radius*thickness;
}

static double
radius_from_excluded_volume(double radius, double thickness)
{
    return 0.5*cbrt(0.75*radius*(2.0*radius*thickness + (radius + thickness)*(M_PI*radius + thickness)));
}

static double
radius_effective(int mode, double radius, double thickness, double alpha, double beta)
{
    switch (mode) {
    default:
    case 1: // equivalent cylinder excluded volume
        return radius_from_excluded_volume(radius, thickness);
    case 2: // equivalent volume sphere
        return cbrt(M_PI*radius*radius*thickness/M_4PI_3);
    case 3: // radius
        return radius;
    }
}

double Iq(
    double q,
    double radius,
    double thickness,
    double alpha,
    double beta,
    double sld,
    double sld_solvent)
{
    double form = _integrate_psi(q, radius, thickness, alpha, beta);
    double contrast = sld - sld_solvent;
    double volume = M_PI*radius*radius*thickness;
    return 1.0e-4*form * square(contrast * volume);
}
