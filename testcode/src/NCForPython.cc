#include "NCExtraTestUtils.hh"

// Expose crude (and horribly fragile) C functions which can be called from
// Python via the ctypes module. Updates to function names and signatures in
// this file should be carefully synchronised with code in the ../python/lib.py
// file.

extern "C" {

  void nctest_getmanyxsvalues( double sigma, double lambda_cutoff, double R, double length, unsigned array_size, const double* ekin_array, double* output_xs_array )
  {
    NCP::PhysicsModel pm(sigma,lambda_cutoff, R, length);
    for (unsigned i = 0; i < array_size; ++i)
      output_xs_array[i] = pm.calcDiffCrossSection(ekin_array[i]);
  }

  void nctest_samplemanyscatmu( double sigma, double lambda_cutoff, double R, double length, double ekin, unsigned nvalues, double* output_mu )
  {
    NCP::PhysicsModel pm(sigma,lambda_cutoff, R, length);
    auto rng = NC::getRNG();
    for (unsigned i = 0; i < nvalues; ++i)
      output_mu[i] = pm.sampleScatteringEvent(rng,ekin).mu;
  }

}
