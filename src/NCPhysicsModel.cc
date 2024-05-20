#include "NCPhysicsModel.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

// Include external C dependencies for sasview models
extern "C"
  {
    #include "../models/kernel_header.c"
    #include "../models/lib/polevl.c"
    #include "../models/lib/sas_J1.c"
    #include "../models/lib/gauss20.c"
    #include "../models/cylinder.c"
  }

bool NCP::PhysicsModel::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_<ourpluginname>
  //    <sigmavalue> <wavelength threshold value>
  //

  //Verify we have exactly one line and two words:
  if ( data.size() != 1 || data.at(0).size()!=4 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be two numbers on a single line");

  //Parse and validate values:
  double sigma, lambda_cutoff,R,length;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), sigma )
       || ! NC::safe_str2dbl( data.at(0).at(1), lambda_cutoff )
       || ! NC::safe_str2dbl( data.at(0).at(2), R)
       || ! NC::safe_str2dbl( data.at(0).at(3), length)
       || ! (sigma>0.0) || !(lambda_cutoff>=0.0) || !(R>=0.0) || (length>=0.0))
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" section (should be two positive floating point values)" );

  //Parsing done! Create and return our model:
  return PhysicsModel(sigma,lambda_cutoff,R,length);
}

NCP::PhysicsModel::PhysicsModel( double sigma, double lambda_cutoff, double R, double length)
  : m_sigma(sigma),
    m_cutoffekin(NC::wl2ekin(lambda_cutoff)),
    m_radius(R),
    m_length(length)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( m_sigma > 0.0 );
  nc_assert( m_cutoffekin > 0.0);
  nc_assert( m_radius > 0.0);
  nc_assert( m_length > 0.0);
}


double NCP::PhysicsModel::calcCrossSection( double Q ) const
{
  // Define scattering amplitude (F1) and it's form factor (F2); 
  double F1, F2;

  // To be user defined:
  double sld{4.0};
  double solvent_sld{1.0};
  double radius{20.0};
  double length{400.0};

  // Randomly oriented cylinder form factor calculator
  Fq(Q, &F1, &F2, sld, solvent_sld, radius, length);
  // Normalize form factor per particle volume
  F2 /= form_volume(radius, length);

  // 1/v * F(Q)
  return F2;
}

NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;

  if ( ! (neutron_ekin > m_cutoffekin) ) {
    //Special case: We are asked to sample a scattering event for a neutron
    //energy where we have zero cross section! Although in a real simulation we
    //would usually not expect this to happen, users with custom code might
    //still generate such calls. The only consistent thing to do when the cross
    //section is zero is to not change the neutron state parameters, which means:
    result.ekin_final = neutron_ekin;
    result.mu = 1.0;
    return result;
  }

  //Implement our actual model here. Of course it is trivial for the example
  //model. For a more realistic or complicated model, it might be that
  //additional helper classes or functions should be created and used, in order
  //to keep the code here manageable:

  result.ekin_final = neutron_ekin;//Elastic scattering

  double k{std::sqrt(NC::ekin2ksq(neutron_ekin))}; //modulus of wavevector

  double qmax{1.0}; // qmax for sampling, should be user defined maybe
  double qmin{1e-3}; // qmin for sampling, should be user defined also.
  double step{1e-2}; // steps to calculate Imax

  // Get de Maximum value of the XS to sample from it
  double Imax{std::numeric_limits<double>::lowest()};
  double maxX{qmin};
  for (double x = qmin; x<=qmax; x+=step){
    double val{NCP::PhysicsModel::calcCrossSection(x)};
    if (val > Imax){
      Imax = val;
      maxX = x;
    }
  }

  // Sample randomly following the XS distribution:
  // valid only for elastic scattering
  double q{}, xs{}, acceptance{};
  do {
    q = rng.generate()*qmax;
    acceptance = rng.generate()*Imax;
    xs = NCP::PhysicsModel::calcCrossSection(q);
  } while (acceptance>xs);
  result.mu = 1 - 0.5*NC::ncsquare(q/k);
  return result;
}


