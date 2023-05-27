#define HAS_Iqabc
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
#line 1 "./models/lib/sas_3j1x_x.c"
/**
* Spherical Bessel function 3*j1(x)/x
*
* Used for low q to avoid cancellation error.
* Note that the values differ from sasview ~ 5e-12 rather than 5e-14, but
* in this case it is likely cancellation errors in the original expression
* using double precision that are the source.
*/
double sas_3j1x_x(double q);

// The choice of the number of terms in the series and the cutoff value for
// switching between series and direct calculation depends on the numeric
// precision.
//
// Point where direct calculation reaches machine precision:
//
//   single machine precision eps 3e-8 at qr=1.1 **
//   double machine precision eps 4e-16 at qr=1.1
//
// Point where Taylor series reaches machine precision (eps), where taylor
// series matches direct calculation (cross) and the error at that point:
//
//   prec   n eps  cross  error
//   single 3 0.28  0.4   6.2e-7
//   single 4 0.68  0.7   2.3e-7
//   single 5 1.18  1.2   7.5e-8
//   double 3 0.01  0.03  2.3e-13
//   double 4 0.06  0.1   3.1e-14
//   double 5 0.16  0.2   5.0e-15
//
// ** Note: relative error on single precision starts increase on the direct
// method at qr=1.1, rising from 3e-8 to 5e-5 by qr=1e3.  This should be
// safe for the sans range, with objects of 100 nm supported to a q of 0.1
// while maintaining 5 digits of precision.  For usans/sesans, the objects
// are larger but the q is smaller, so again it should be fine.
//
// See explore/sph_j1c.py for code to explore these ranges.

// Use 4th order series
#if FLOAT_SIZE>4
#define SPH_J1C_CUTOFF 0.1
#else
#define SPH_J1C_CUTOFF 0.7
#endif

double sas_3j1x_x(double q)
{
    // 2017-05-18 PAK - support negative q
    if (fabs(q) < SPH_J1C_CUTOFF) {
        const double q2 = q*q;
        return (1.0 + q2*(-3./30. + q2*(3./840. + q2*(-3./45360.))));// + q2*(3./3991680.)))));
    } else {
        double sin_q, cos_q;
        SINCOS(q, sin_q, cos_q);
        return 3.0*(sin_q/q - cos_q)/(q*q);
    }
}

#line 1 "./models/lib/gauss150.c"
// Created by Andrew Jackson on 4/23/07

 #ifdef GAUSS_N
 # undef GAUSS_N
 # undef GAUSS_Z
 # undef GAUSS_W
 #endif
 #define GAUSS_N 150
 #define GAUSS_Z Gauss150Z
 #define GAUSS_W Gauss150Wt


// Note: using array size 152 rather than 150 so that it is a multiple of 4.
// Some OpenCL devices prefer that vectors start and end on nice boundaries.
constant double Gauss150Z[152]={
  	-0.9998723404457334,
  	-0.9993274305065947,
  	-0.9983473449340834,
  	-0.9969322929775997,
  	-0.9950828645255290,
  	-0.9927998590434373,
  	-0.9900842691660192,
  	-0.9869372772712794,
  	-0.9833602541697529,
  	-0.9793547582425894,
  	-0.9749225346595943,
  	-0.9700655145738374,
  	-0.9647858142586956,
  	-0.9590857341746905,
  	-0.9529677579610971,
  	-0.9464345513503147,
  	-0.9394889610042837,
  	-0.9321340132728527,
  	-0.9243729128743136,
  	-0.9162090414984952,
  	-0.9076459563329236,
  	-0.8986873885126239,
  	-0.8893372414942055,
  	-0.8795995893549102,
  	-0.8694786750173527,
  	-0.8589789084007133,
  	-0.8481048644991847,
  	-0.8368612813885015,
  	-0.8252530581614230,
  	-0.8132852527930605,
  	-0.8009630799369827,
  	-0.7882919086530552,
  	-0.7752772600680049,
  	-0.7619248049697269,
  	-0.7482403613363824,
  	-0.7342298918013638,
  	-0.7198995010552305,
  	-0.7052554331857488,
  	-0.6903040689571928,
  	-0.6750519230300931,
  	-0.6595056411226444,
  	-0.6436719971150083,
  	-0.6275578900977726,
  	-0.6111703413658551,
  	-0.5945164913591590,
  	-0.5776035965513142,
  	-0.5604390262878617,
  	-0.5430302595752546,
  	-0.5253848818220803,
  	-0.5075105815339176,
  	-0.4894151469632753,
  	-0.4711064627160663,
  	-0.4525925063160997,
  	-0.4338813447290861,
  	-0.4149811308476706,
  	-0.3959000999390257,
  	-0.3766465660565522,
  	-0.3572289184172501,
  	-0.3376556177463400,
  	-0.3179351925907259,
  	-0.2980762356029071,
  	-0.2780873997969574,
  	-0.2579773947782034,
  	-0.2377549829482451,
  	-0.2174289756869712,
  	-0.1970082295132342,
  	-0.1765016422258567,
  	-0.1559181490266516,
  	-0.1352667186271445,
  	-0.1145563493406956,
  	-0.0937960651617229,
  	-0.0729949118337358,
  	-0.0521619529078925,
  	-0.0313062657937972,
  	-0.0104369378042598,
  	0.0104369378042598,
  	0.0313062657937972,
  	0.0521619529078925,
  	0.0729949118337358,
  	0.0937960651617229,
  	0.1145563493406956,
  	0.1352667186271445,
  	0.1559181490266516,
  	0.1765016422258567,
  	0.1970082295132342,
  	0.2174289756869712,
  	0.2377549829482451,
  	0.2579773947782034,
  	0.2780873997969574,
  	0.2980762356029071,
  	0.3179351925907259,
  	0.3376556177463400,
  	0.3572289184172501,
  	0.3766465660565522,
  	0.3959000999390257,
  	0.4149811308476706,
  	0.4338813447290861,
  	0.4525925063160997,
  	0.4711064627160663,
  	0.4894151469632753,
  	0.5075105815339176,
  	0.5253848818220803,
  	0.5430302595752546,
  	0.5604390262878617,
  	0.5776035965513142,
  	0.5945164913591590,
  	0.6111703413658551,
  	0.6275578900977726,
  	0.6436719971150083,
  	0.6595056411226444,
  	0.6750519230300931,
  	0.6903040689571928,
  	0.7052554331857488,
  	0.7198995010552305,
  	0.7342298918013638,
  	0.7482403613363824,
  	0.7619248049697269,
  	0.7752772600680049,
  	0.7882919086530552,
  	0.8009630799369827,
  	0.8132852527930605,
  	0.8252530581614230,
  	0.8368612813885015,
  	0.8481048644991847,
  	0.8589789084007133,
  	0.8694786750173527,
  	0.8795995893549102,
  	0.8893372414942055,
  	0.8986873885126239,
  	0.9076459563329236,
  	0.9162090414984952,
  	0.9243729128743136,
  	0.9321340132728527,
  	0.9394889610042837,
  	0.9464345513503147,
  	0.9529677579610971,
  	0.9590857341746905,
  	0.9647858142586956,
  	0.9700655145738374,
  	0.9749225346595943,
  	0.9793547582425894,
  	0.9833602541697529,
  	0.9869372772712794,
  	0.9900842691660192,
  	0.9927998590434373,
  	0.9950828645255290,
  	0.9969322929775997,
  	0.9983473449340834,
  	0.9993274305065947,
  	0.9998723404457334,
  	0., // zero padding is ignored
  	0.  // zero padding is ignored
};

constant double Gauss150Wt[152]={
  	0.0003276086705538,
  	0.0007624720924706,
  	0.0011976474864367,
  	0.0016323569986067,
  	0.0020663664924131,
  	0.0024994789888943,
  	0.0029315036836558,
  	0.0033622516236779,
  	0.0037915348363451,
  	0.0042191661429919,
  	0.0046449591497966,
  	0.0050687282939456,
  	0.0054902889094487,
  	0.0059094573005900,
  	0.0063260508184704,
  	0.0067398879387430,
  	0.0071507883396855,
  	0.0075585729801782,
  	0.0079630641773633,
  	0.0083640856838475,
  	0.0087614627643580,
  	0.0091550222717888,
  	0.0095445927225849,
  	0.0099300043714212,
  	0.0103110892851360,
  	0.0106876814158841,
  	0.0110596166734735,
  	0.0114267329968529,
  	0.0117888704247183,
  	0.0121458711652067,
  	0.0124975796646449,
  	0.0128438426753249,
  	0.0131845093222756,
  	0.0135194311690004,
  	0.0138484622795371,
  	0.0141714592928592,
  	0.0144882814685445,
  	0.0147987907597169,
  	0.0151028518701744,
  	0.0154003323133401,
  	0.0156911024699895,
  	0.0159750356447283,
  	0.0162520081211971,
  	0.0165218992159766,
  	0.0167845913311726,
  	0.0170399700056559,
  	0.0172879239649355,
  	0.0175283451696437,
  	0.0177611288626114,
  	0.0179861736145128,
  	0.0182033813680609,
  	0.0184126574807331,
  	0.0186139107660094,
  	0.0188070535331042,
  	0.0189920016251754,
  	0.0191686744559934,
  	0.0193369950450545,
  	0.0194968900511231,
  	0.0196482898041878,
  	0.0197911283358190,
  	0.0199253434079123,
  	0.0200508765398072,
  	0.0201676730337687,
  	0.0202756819988200,
  	0.0203748563729175,
  	0.0204651529434560,
  	0.0205465323660984,
  	0.0206189591819181,
  	0.0206824018328499,
  	0.0207368326754401,
  	0.0207822279928917,
  	0.0208185680053983,
  	0.0208458368787627,
  	0.0208640227312962,
  	0.0208731176389954,
  	0.0208731176389954,
  	0.0208640227312962,
  	0.0208458368787627,
  	0.0208185680053983,
  	0.0207822279928917,
  	0.0207368326754401,
  	0.0206824018328499,
  	0.0206189591819181,
  	0.0205465323660984,
  	0.0204651529434560,
  	0.0203748563729175,
  	0.0202756819988200,
  	0.0201676730337687,
  	0.0200508765398072,
  	0.0199253434079123,
  	0.0197911283358190,
  	0.0196482898041878,
  	0.0194968900511231,
  	0.0193369950450545,
  	0.0191686744559934,
  	0.0189920016251754,
  	0.0188070535331042,
  	0.0186139107660094,
  	0.0184126574807331,
  	0.0182033813680609,
  	0.0179861736145128,
  	0.0177611288626114,
  	0.0175283451696437,
  	0.0172879239649355,
  	0.0170399700056559,
  	0.0167845913311726,
  	0.0165218992159766,
  	0.0162520081211971,
  	0.0159750356447283,
  	0.0156911024699895,
  	0.0154003323133401,
  	0.0151028518701744,
  	0.0147987907597169,
  	0.0144882814685445,
  	0.0141714592928592,
  	0.0138484622795371,
  	0.0135194311690004,
  	0.0131845093222756,
  	0.0128438426753249,
  	0.0124975796646449,
  	0.0121458711652067,
  	0.0117888704247183,
  	0.0114267329968529,
  	0.0110596166734735,
  	0.0106876814158841,
  	0.0103110892851360,
  	0.0099300043714212,
  	0.0095445927225849,
  	0.0091550222717888,
  	0.0087614627643580,
  	0.0083640856838475,
  	0.0079630641773633,
  	0.0075585729801782,
  	0.0071507883396855,
  	0.0067398879387430,
  	0.0063260508184704,
  	0.0059094573005900,
  	0.0054902889094487,
  	0.0050687282939456,
  	0.0046449591497966,
  	0.0042191661429919,
  	0.0037915348363451,
  	0.0033622516236779,
  	0.0029315036836558,
  	0.0024994789888943,
  	0.0020663664924131,
  	0.0016323569986067,
  	0.0011976474864367,
  	0.0007624720924706,
  	0.0003276086705538,
  	0., // zero padding is ignored
  	0.  // zero padding is ignored
};

#line 1 "./models/lib/sphere_form.c"
double sphere_volume(double radius);
double sphere_form(double q, double radius, double sld, double solvent_sld);

double sphere_volume(double radius)
{
    return M_4PI_3*cube(radius);
}

double sphere_form(double q, double radius, double sld, double solvent_sld)
{
    const double fq = sphere_volume(radius) * sas_3j1x_x(q*radius);
    const double contrast = (sld - solvent_sld);
    return 1.0e-4*square(contrast * fq);
}


#line 1 "./models/bcc_paracrystal.c"
static double
bcc_Zq(double qa, double qb, double qc, double dnn, double d_factor)
{
    // Equations from Matsuoka 26-27-28, multiplied by |q|
    const double a1 = (-qa + qb + qc)/2.0;
    const double a2 = (+qa - qb + qc)/2.0;
    const double a3 = (+qa + qb - qc)/2.0;
    const double d_a = dnn/sqrt(0.75);

#if 1
    // Matsuoka 29-30-31
    //     Z_k numerator: 1 - exp(a)^2
    //     Z_k denominator: 1 - 2 cos(d a_k) exp(a) + exp(2a)
    // Rewriting numerator
    //         => -(exp(2a) - 1)
    //         => -expm1(2a)
    // Rewriting denominator
    //         => exp(a)^2 - 2 cos(d ak) exp(a) + 1)
    //         => (exp(a) - 2 cos(d ak)) * exp(a) + 1
    const double arg = -0.5*square(dnn*d_factor)*(a1*a1 + a2*a2 + a3*a3);
    const double exp_arg = exp(arg);
    const double Zq = -cube(expm1(2.0*arg))
        / ( ((exp_arg - 2.0*cos(d_a*a1))*exp_arg + 1.0)
          * ((exp_arg - 2.0*cos(d_a*a2))*exp_arg + 1.0)
          * ((exp_arg - 2.0*cos(d_a*a3))*exp_arg + 1.0));

#elif 0
    // ** Alternate form, which perhaps is more approachable
    //     Z_k numerator   => -[(exp(2a) - 1) / 2.exp(a)] 2.exp(a)
    //                     => -[sinh(a)] exp(a)
    //     Z_k denominator => [(exp(2a) + 1) / 2.exp(a) - cos(d a_k)] 2.exp(a)
    //                     => [cosh(a) - cos(d a_k)] 2.exp(a)
    //     => Z_k = -sinh(a) / [cosh(a) - cos(d a_k)]
    //            = sinh(-a) / [cosh(-a) - cos(d a_k)]
    //
    // One more step leads to the form in sasview 3.x for 2d models
    //            = tanh(-a) / [1 - cos(d a_k)/cosh(-a)]
    //
    const double arg = 0.5*square(dnn*d_factor)*(a1*a1 + a2*a2 + a3*a3);
    const double sinh_qd = sinh(arg);
    const double cosh_qd = cosh(arg);
    const double Zq = sinh_qd/(cosh_qd - cos(d_a*a1))
                    * sinh_qd/(cosh_qd - cos(d_a*a2))
                    * sinh_qd/(cosh_qd - cos(d_a*a3));
#else
    const double arg = 0.5*square(dnn*d_factor)*(a1*a1 + a2*a2 + a3*a3);
    const double tanh_qd = tanh(arg);
    const double cosh_qd = cosh(arg);
    const double Zq = tanh_qd/(1.0 - cos(d_a*a1)/cosh_qd)
                    * tanh_qd/(1.0 - cos(d_a*a2)/cosh_qd)
                    * tanh_qd/(1.0 - cos(d_a*a3)/cosh_qd);
#endif

    return Zq;
}


// occupied volume fraction calculated from lattice symmetry and sphere radius
static double
bcc_volume_fraction(double radius, double dnn)
{
    return 2.0*sphere_volume(sqrt(0.75)*radius/dnn);
    // note that sqrt(0.75) = root3/2 and sqrt(0.75)/dnn=1/d_a
    //Thus this is correct
}

static double
form_volume(double radius)
{
    return sphere_volume(radius);
}


static double Iq(double q, double dnn,
    double d_factor, double radius,
    double sld, double solvent_sld)
{
    // translate a point in [-1,1] to a point in [0, 2 pi]
    const double phi_m = M_PI;
    const double phi_b = M_PI;
    // translate a point in [-1,1] to a point in [0, pi]
    const double theta_m = M_PI_2;
    const double theta_b = M_PI_2;

    double outer_sum = 0.0;
    for(int i=0; i<GAUSS_N; i++) {
        double inner_sum = 0.0;
        const double theta = GAUSS_Z[i]*theta_m + theta_b;
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);
        const double qc = q*cos_theta;
        const double qab = q*sin_theta;
        for(int j=0;j<GAUSS_N;j++) {
            const double phi = GAUSS_Z[j]*phi_m + phi_b;
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);
            const double qa = qab*cos_phi;
            const double qb = qab*sin_phi;
            const double form = bcc_Zq(qa, qb, qc, dnn, d_factor);
            inner_sum += GAUSS_W[j] * form;
        }
        inner_sum *= phi_m;  // sum(f(x)dx) = sum(f(x)) dx
        outer_sum += GAUSS_W[i] * inner_sum * sin_theta;
    }
    outer_sum *= theta_m;
    const double Zq = outer_sum/(4.0*M_PI);
    const double Pq = sphere_form(q, radius, sld, solvent_sld);
    return bcc_volume_fraction(radius, dnn) * Pq * Zq;
    // note that until we can return non fitable values to the GUI this
    // can only be queried by a script. Otherwise we can drop the
    // bcc_volume_fraction as it is effectively included in "scale."
}


static double Iqabc(double qa, double qb, double qc,
    double dnn, double d_factor, double radius,
    double sld, double solvent_sld)
{
    const double q = sqrt(qa*qa + qb*qb + qc*qc);
    const double Zq = bcc_Zq(qa, qb, qc, dnn, d_factor);
    const double Pq = sphere_form(q, radius, sld, solvent_sld);
    return bcc_volume_fraction(radius, dnn) * Pq * Zq;
}
