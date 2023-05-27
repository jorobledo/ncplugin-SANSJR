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

#line 1 "./models/binary_hard_sphere.c"
double form_volume(void);

double Iq(double q,
    double lg_radius, double sm_radius,
    double lg_vol_frac, double sm_vol_frac,
    double lg_sld, double sm_sld, double solvent_sld
    );

void calculate_psfs(double qval,
    double r2, double nf2,
    double aa, double phi,
    double *s11, double *s22, double *s12
    );

double form_volume(void)
{
    return 1.0;
}

double Iq(double q,
    double lg_radius, double sm_radius,
    double lg_vol_frac, double sm_vol_frac,
    double lg_sld, double sm_sld, double solvent_sld)
{
    double r2,r1,nf2,phi,aa,rho2,rho1,rhos,inten;       //my local names
    double psf11,psf12,psf22;
    double phi1,phi2,phr,a3;
    double v1,v2,n1,n2,qr1,qr2,b1,b2,sc1,sc2;

    r2 = lg_radius;
    r1 = sm_radius;
    phi2 = lg_vol_frac;
    phi1 = sm_vol_frac;
    rho2 = lg_sld;
    rho1 = sm_sld;
    rhos = solvent_sld;


    phi = phi1 + phi2;
    aa = r1/r2;
    //calculate the number fraction of larger spheres (eqn 2 in reference)
    a3=aa*aa*aa;
    phr=phi2/phi;
    nf2 = phr*a3/(1.0-phr+phr*a3);
    // calculate the PSF's here
    calculate_psfs(q,r2,nf2,aa,phi,&psf11,&psf22,&psf12);

    // /* do form factor calculations  */

    v1 = M_4PI_3*r1*r1*r1;
    v2 = M_4PI_3*r2*r2*r2;

    n1 = phi1/v1;
    n2 = phi2/v2;

    qr1 = r1*q;
    qr2 = r2*q;

    sc1 = sas_3j1x_x(qr1);
    sc2 = sas_3j1x_x(qr2);
    b1 = r1*r1*r1*(rho1-rhos)*M_4PI_3*sc1;
    b2 = r2*r2*r2*(rho2-rhos)*M_4PI_3*sc2;
    inten = n1*b1*b1*psf11;
    inten += sqrt(n1*n2)*2.0*b1*b2*psf12;
    inten += n2*b2*b2*psf22;
    ///* convert I(1/A) to (1/cm)  */
    inten *= 1.0e8;
    ///*convert rho^2 in 10^-6A to A*/
    inten *= 1.0e-12;
    return(inten);
}


void calculate_psfs(double qval,
    double r2, double nf2,
    double aa, double phi,
    double *s11, double *s22, double *s12)
{
    //  variable qval,r2,nf2,aa,phi,&s11,&s22,&s12

    //   calculate constant terms
    double s2,v,a3,v1,v2,g11,g12,g22,wmv,wmv3,wmv4;
    double a1,a2i,a2,b1,b2,b12,gm1,gm12;
    double yy,ay,ay2,ay3,t1,t2,t3,f11,y2,y3,tt1,tt2,tt3;
    double c11,c22,c12,f12,f22,ttt1,ttt2,ttt3,ttt4,yl,y13;
    double t21,t22,t23,t31,t32,t33,t41,t42,yl3,wma3,y1;

    s2 = 2.0*r2;
//    s1 = aa*s2;  why is this never used?  check original paper?
    v = phi;
    a3 = aa*aa*aa;
    v1=((1.-nf2)*a3/(nf2+(1.-nf2)*a3))*v;
    v2=(nf2/(nf2+(1.-nf2)*a3))*v;
    g11=((1.+.5*v)+1.5*v2*(aa-1.))/(1.-v)/(1.-v);
    g22=((1.+.5*v)+1.5*v1*(1./aa-1.))/(1.-v)/(1.-v);
    g12=((1.+.5*v)+1.5*(1.-aa)*(v1-v2)/(1.+aa))/(1.-v)/(1.-v);
    wmv = 1/(1.-v);
    wmv3 = wmv*wmv*wmv;
    wmv4 = wmv*wmv3;
    a1=3.*wmv4*((v1+a3*v2)*(1.+v+v*v)-3.*v1*v2*(1.-aa)*(1.-aa)*(1.+v1+aa*(1.+v2))) + ((v1+a3*v2)*(1.+2.*v)+(1.+v+v*v)-3.*v1*v2*(1.-aa)*(1.-aa)-3.*v2*(1.-aa)*(1.-aa)*(1.+v1+aa*(1.+v2)))*wmv3;
    a2i=((v1+a3*v2)*(1.+v+v*v)-3.*v1*v2*(1.-aa)*(1.-aa)*(1.+v1+aa*(1.+v2)))*3*wmv4 + ((v1+a3*v2)*(1.+2.*v)+a3*(1.+v+v*v)-3.*v1*v2*(1.-aa)*(1.-aa)*aa-3.*v1*(1.-aa)*(1.-aa)*(1.+v1+aa*(1.+v2)))*wmv3;
    a2=a2i/a3;
    b1=-6.*(v1*g11*g11+.25*v2*(1.+aa)*(1.+aa)*aa*g12*g12);
    b2=-6.*(v2*g22*g22+.25*v1/a3*(1.+aa)*(1.+aa)*g12*g12);
    b12=-3.*aa*(1.+aa)*(v1*g11/aa/aa+v2*g22)*g12;
    gm1=(v1*a1+a3*v2*a2)*.5;
    gm12=2.*gm1*(1.-aa)/aa;
    //c
    //c   calculate the direct correlation functions and print results
    //c
    //  do 20 j=1,npts

    yy=qval*s2;
    //c   calculate direct correlation functions
    //c   ----c11
    ay=aa*yy;
    ay2 = ay*ay;
    ay3 = ay*ay*ay;
    t1=a1*(sin(ay)-ay*cos(ay));
    t2=b1*(2.*ay*sin(ay)-(ay2-2.)*cos(ay)-2.)/ay;
    t3=gm1*((4.*ay*ay2-24.*ay)*sin(ay)-(ay2*ay2-12.*ay2+24.)*cos(ay)+24.)/ay3;
    f11=24.*v1*(t1+t2+t3)/ay3;

    //c ------c22
    y2=yy*yy;
    y3=yy*y2;
    tt1=a2*(sin(yy)-yy*cos(yy));
    tt2=b2*(2.*yy*sin(yy)-(y2-2.)*cos(yy)-2.)/yy;
    tt3=gm1*((4.*y3-24.*yy)*sin(yy)-(y2*y2-12.*y2+24.)*cos(yy)+24.)/ay3;
    f22=24.*v2*(tt1+tt2+tt3)/y3;

    //c   -----c12
    yl=.5*yy*(1.-aa);
    yl3=yl*yl*yl;
    wma3 = (1.-aa)*(1.-aa)*(1.-aa);
    y1=aa*yy;
    y13 = y1*y1*y1;
    ttt1=3.*wma3*v*sqrt(nf2)*sqrt(1.-nf2)*a1*(sin(yl)-yl*cos(yl))/((nf2+(1.-nf2)*a3)*yl3);
    t21=b12*(2.*y1*cos(y1)+(y1*y1-2.)*sin(y1));
    t22=gm12*((3.*y1*y1-6.)*cos(y1)+(y1*y1*y1-6.*y1)*sin(y1)+6.)/y1;
    t23=gm1*((4.*y13-24.*y1)*cos(y1)+(y13*y1-12.*y1*y1+24.)*sin(y1))/(y1*y1);
    t31=b12*(2.*y1*sin(y1)-(y1*y1-2.)*cos(y1)-2.);
    t32=gm12*((3.*y1*y1-6.)*sin(y1)-(y1*y1*y1-6.*y1)*cos(y1))/y1;
    t33=gm1*((4.*y13-24.*y1)*sin(y1)-(y13*y1-12.*y1*y1+24.)*cos(y1)+24.)/(y1*y1);
    t41=cos(yl)*((sin(y1)-y1*cos(y1))/(y1*y1) + (1.-aa)/(2.*aa)*(1.-cos(y1))/y1);
    t42=sin(yl)*((cos(y1)+y1*sin(y1)-1.)/(y1*y1) + (1.-aa)/(2.*aa)*sin(y1)/y1);
    ttt2=sin(yl)*(t21+t22+t23)/(y13*y1);
    ttt3=cos(yl)*(t31+t32+t33)/(y13*y1);
    ttt4=a1*(t41+t42)/y1;
    f12=ttt1+24.*v*sqrt(nf2)*sqrt(1.-nf2)*a3*(ttt2+ttt3+ttt4)/(nf2+(1.-nf2)*a3);

    c11=f11;
    c22=f22;
    c12=f12;
    *s11=1./(1.+c11-(c12)*c12/(1.+c22));
    *s22=1./(1.+c22-(c12)*c12/(1.+c11));
    *s12=-c12/((1.+c11)*(1.+c22)-(c12)*(c12));

    return;
}
